# there were some issues in bbsAssistant package but it does not appear to be
# actively maintained. So install from my fork
# remotes::install_github("see24/bbsAssistant", ref = "fix-QualityCurrentID-filter")

library(bbsAssistant)
library(dplyr)
library(reproducible)
library(napops)

# Set parameters for script
start_yr <- 1991
end_yr <- 2023

quality_filter <- "RunType == 1"

juris <- "ontario"

# Bird Conservation Regions to filter to
bcrs <- c(7,8)

options(reproducible.cachePath = "analysis/data/cache_data")

# Step 1: #=====================================================================
# download BBS data, extract from files, filter to area and time range
# and based on quality, Get species letter codes, change format to one col per
# species and one row per stop/year/route
# Get the BBS data
# downloads the data the first time and will import it if already downloaded.
# using Cache from reproducible to cache long running functions
bbs_raw <- Cache(grab_bbs_data(sb_id = "625f151ed34e85fa62b7f926",
                               bbs_dir = "analysis/data/raw_data/bbs_2022"))

obs_ON_1991 <- Cache(munge_bbs_data(
  bbs_raw,
  states = region_codes %>% filter(State %in% juris) %>% pull(StateNum),
  zero.fill = FALSE,
  year.range = start_yr:end_yr,
  keep.stop.level.data = TRUE
))

# species napops
napop_sp <- napops::list_species()

# Use ECCC Avian Core database from
# https://ecollab.ncr.int.ec.gc.ca/theme/cws-scf/MBECAvianTaxonomicDatabase/Forms/AllItems.aspx
av_core <- readxl::read_xlsx("analysis/data/raw_data/ECCC_Avian_Taxonomic_Database/ECCC Avian Core 20230601.xlsx")

av_core <- av_core %>%
  select(Species_ID, Parent_ID, Full_Species, English_Name, Scientific_Name,
         Prev_Species_ID, BBS_Number) %>%
  # napops does not include subspecies so they should be grouped with parent
  mutate(Species_ID = coalesce(Parent_ID, Species_ID))

# Use av_core to connect AOU codes in bbs to 4 letter codes in napops
sp_lu_tbl <- bbs_raw$species_list %>%
  select(AOU, English_Common_Name, Scientific_Name, AOU4) %>%
  left_join(av_core, by = c(AOU = "BBS_Number")) %>%
  full_join(napop_sp, by = c(Species_ID = "Species")) %>%
  mutate(Scientific_Name = coalesce(Scientific_Name.x, Scientific_Name.y,
                                    Scientific_Name), .keep = "unused")

obs_ON_1991_2a <- obs_ON_1991 %>%
  filter(BCR %in% bcrs) %>%
  filter(!!rlang::parse_expr(quality_filter)) %>%
  # get 4 letter codes for sp names
  left_join(sp_lu_tbl %>% select(AOU, Species_ID), by = "AOU")

obs_ON_1991_2 <- obs_ON_1991_2a %>%
  # Some are missing so keep numeric codes for those
  mutate(sp_code = coalesce(Species_ID, as.character(AOU)), .keep = "unused") %>%
  select(-contains("Car"), -contains("Noise"), -RouteTotal, -TotalSpp,
         -contains("ObsFirstYear"), -Assistant, -RunType, -QualityCurrentID,
         -julian) %>%
  tidyr::pivot_longer(contains("Stop"), names_to = "Stop", values_to = "Abund") %>%
  filter(Abund != 0) %>%
  tidyr::pivot_wider(names_from = sp_code, values_from = Abund, values_fill = 0)

# flag species not matched to codes and check against napops species list as well
no_match_aou <- obs_ON_1991_2a %>% filter(is.na(Species_ID)) %>% pull(AOU) %>% unique()

filter(species_list, AOU %in% no_match_aou) %>% select(-Spanish_Common_Name) %>%
  {stopifnot(nrow(.) == 0)}

# Step 2: #=====================================================================
# Get locations of each stop on the routes from data supplied by CWS, not
# currently available online.

# Get names
# sf::st_layers("analysis/data/raw_data/BBS Routes_CurrentStops_Discontinued.gdb")

cur_stops <- sf::st_read("analysis/data/raw_data/BBS Routes_CurrentStops_Discontinued.gdb",
                         layer = "CURRENT_STOPS")
discon_stops <- sf::st_read("analysis/data/raw_data/BBS Routes_CurrentStops_Discontinued.gdb",
                         layer = "DISCONTINUED_STOPS")

# bbox is wrong in discon_stops
new_bb <- discon_stops %>% sf::st_drop_geometry() %>%
  summarise(across(c(POINT_X, POINT_Y),
                   lst(min = \(x){min(x, na.rm = TRUE)},
                       max = \(x){max(x, na.rm = TRUE)}))) %>%
  select(xmin = POINT_X_min, ymin = POINT_Y_min, xmax = POINT_X_max,
         ymax = POINT_Y_max) %>% unlist()

attr(new_bb, "class") <-  "bbox"

attr(sf::st_geometry(discon_stops), "bbox") <-  new_bb

ON_stops <- bind_rows(cur_stops, discon_stops %>% mutate(Year = as.integer(Year))) %>%
  sf::st_drop_geometry() %>%
  filter(Province == obs_ON_1991$StateNum %>% unique()) %>%
  # add country code to province route code for joining
  mutate(RTENO = Province*1000 + route +12400000) %>%
  # There are 2 rows that are identical
  distinct()

# check there is one record for each route, stop combo
test <- ON_stops %>% sf::st_drop_geometry() %>% group_by(RTENO, Stop) %>%
  mutate(n = n()) %>%
  filter(n > 1)

# keep the first occurrence of each route stop combo. Will keep current or first
# of records in discontinued
ON_stops <- ON_stops  %>% group_by(RTENO, Stop) %>%
  summarise(across(everything(), first))

obs_ON_1991_stops <- obs_ON_1991_2 %>%
  mutate(RTENO = as.numeric(RTENO), Stop = as.numeric(stringr::str_extract(Stop, "\\d\\d?"))) %>%
  left_join(ON_stops %>% select(RTENO, Stop, POINT_X, POINT_Y),
            by = c("RTENO", "Stop"))

# check for NA coords
missing <- obs_ON_1991_stops %>% filter(is.na(POINT_X)) %>% pull(Route) %>% unique()

# 6 discontinued routes don't have coords for each stop, just have start and end
# below I fill from the first stop for sunrise and take mean forest based on
# stops 1 and 50


# Step 3: #=====================================================================
# Covariates for calculating offsets
# Need time since sunrise (TSSR), ordinal day, and proportion forest, since BBS
# is always on roads don't need to do that

# Get time since sunrise for each stop on the route based on start and
# end times and number of stops

# get total time for route and assume each stop same length, use to get TSSR for each
obs_ON_1991_stops_2 <- obs_ON_1991_stops %>% ungroup() %>%
  mutate(timez = lutz::tz_lookup_coords(POINT_Y, POINT_X, method = "accurate")) %>%
  # fill in tz for ones that are NA coords
  arrange(RTENO, Year, Stop) %>%
  group_by(RTENO, Year) %>%
  tidyr::fill(timez) %>%
  group_by(timez) %>%
  mutate(StartTime = lubridate::ymd_hm(paste(Date, StartTime), tz = first(timez)),
         EndTime = lubridate::ymd_hm(paste(Date, EndTime), tz = first(timez)),
         totTime = EndTime - StartTime) %>%
  group_by(RTENO, Year) %>%
  mutate(nStops = n(), stopInc = 1:n(), incTime = totTime/nStops,
         stopTime1 = StartTime + (stopInc - 1) * incTime,
         stopTime2 = StartTime + (stopInc) * incTime) %>%
  ungroup() %>%
  mutate(sunrise = suncalc::getSunlightTimes(
    data = data.frame(date = Date, lon = POINT_X, lat = POINT_Y),
    keep = "sunrise")$sunrise) %>%
  group_by(RTENO, Year) %>%
  tidyr::fill(sunrise) %>%
  ungroup() %>%
  mutate(tssr = difftime(stopTime1, sunrise, units = "hours"),
         od = lubridate::yday(Date))

# Use tree cover, 1km^2 resolution data stored in borealbirds/qpad-offsets repo
# tree <- reproducible::prepInputs(url = "https://raw.githubusercontent.com/borealbirds/qpad-offsets/main/data/tree.tif",
#                                 destinationPath = "analysis/data/raw_data/",
#                                 writeTo = NULL)
#
# tree <- terra::subst(tree, from = c(254, 255), c(NA, 0))
#
# obs_tree <- terra::extract(tree, obs_ON_1991_stops_2 %>%
#                             filter(!is.na(POINT_X)) %>%
#                             sf::st_as_sf(coords = c("POINT_X", "POINT_Y"),
#                                          crs = 4326) %>%
#                             sf::st_transform(sf::st_crs(tree)) %>%
#                             terra::vect())

# Use NALCMS data same as napops paper. I am using 2020 data could in theory use
# version closest to actual year of data collection
url_nalcms_2020_Can <- "http://www.cec.org/files/atlas_layers/1_terrestrial_ecosystems/1_01_0_land_cover_2020_30m/can_land_cover_2020_30m_tif.zip"
nalcms_2020_Can <- reproducible::prepInputs(url = url_nalcms_2020_Can,
                                            destinationPath = "analysis/data/raw_data/",
                                            writeTo = NULL)

lcc_lu <- tibble::tribble(
  ~from, ~to, ~becomes,
  1,       2,        1,
  5,       6,        1)

# get boundary of bird data + ~2000m so window doesn't include edges
SA <- ON_stops %>% ungroup() %>%
  summarise(across(c(POINT_X, POINT_Y), lst(min, max), na.rm = TRUE)) %>%
  select(xmin = POINT_X_min, ymin = POINT_Y_min, xmax = POINT_X_max,
         ymax = POINT_Y_max) %>% unlist() %>%
  sf::st_bbox() %>%
  {. + 0.02} %>%
  sf::st_set_crs(sf::st_crs(cur_stops)) %>%
  sf::st_as_sfc() %>%
  sf::st_transform(sf::st_crs(nalcms_2020_Can))

# Crop and convert to prop forest in 5x5 window
for_cov <- Cache(nalcms_2020_Can %>%
  terra::crop(terra::vect(SA)) %>%
  terra::classify(rcl = lcc_lu, others = 0, right = NA))

for_cov150 <- Cache(terra::focal(for_cov, w = 5, fun = "mean"))

obs_for_cov <- obs_ON_1991_stops_2 %>% filter(!is.na(POINT_X)) %>%
  select(RTENO, Year, Stop, POINT_X, POINT_Y) %>%
  sf::st_as_sf(coords = c("POINT_X", "POINT_Y"),
               crs = sf::st_crs(cur_stops)) %>%
  sf::st_transform(sf::st_crs(nalcms_2020_Can)) %>%
  {terra::extract(for_cov150, ., bind = TRUE)} %>%
  sf::st_as_sf() %>%
  sf::st_drop_geometry() %>%
  rename(for_cov = focal_mean)

# join to forest cover and calculate average from first and last stop when stop
# coords are missing
obs_ON_1991_stops_3 <- obs_ON_1991_stops_2 %>%
  left_join(obs_for_cov, by = join_by(RTENO, Year, Stop)) %>%
  group_by(RTENO, Year) %>%
  mutate(mean_for_cov = mean(for_cov, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(for_cov = coalesce(for_cov, mean_for_cov))


# Step 4: #=====================================================================
# Calculate offsets using napops package

# Needs to be available in order for faster versions of functions to work
# preloading a table of species min, max median for each model type
sp_tbl <- list_species()

sum_covs <- function(df){
  df %>%
    summarise(across(-c(Species, Survey_Method), lst(min, max, median)))
}

all_cov_tbl <- Cache(sp_tbl %>% rowwise() %>%
                       mutate(sum_covs(covariates_removal(Species)),
                              sum_covs(covariates_distance(Species))))


# A is the known or estimated area of survey, p is availability given presence, q is detectability given avaibalility.
# offset=log(A) + log(p) + log(q)
obs_ON_1991_offsets_a <- obs_ON_1991_stops_3 %>%
  select(-c(
    Active, Latitude, Longitude, Stratum, BCR, RouteTypeID, RouteTypeDetailID,
    ObsN, Month, Day, StartTemp, EndTemp, TempScale, StartWind, EndWind,
    StartSky, EndSky, StartTime, EndTime, WindMean, TempMean, timez, totTime,
    nStops, stopInc, incTime, stopTime1, stopTime2, sunrise, mean_for_cov)) %>%
  mutate(tssr = as.numeric(tssr)) %>%
  tidyr::pivot_longer(-c(CountryNum, StateNum, Route, RPID, Year, RTENO, RouteName,
                  Date, Stop, POINT_X, POINT_Y, tssr, od, for_cov),
               names_to = "species", values_to = "count") %>%
  filter(count != 0)

obs_ON_1991_offsets <- obs_ON_1991_offsets_a %>% group_by(species) %>%
  mutate(
    on_road = TRUE,
    cr_est = cue_rate3(unique(species), model = "best", od, tssr,
                       pairwise = TRUE) %>% pull(CR_est),
    edr_est = edr2(unique(species), model = "best", road = on_road, forest = for_cov,
                   pairwise = TRUE) %>% pull(EDR_est),
    # doing these calculations manually since avail and percept funs were to slow
    avail_est = 1 - exp(-3 * cr_est),
    percept_est = (pi * edr_est^2 * (1 - exp(-400^2/edr_est^2)))/(pi * 400^2),
    # changing to 4 instead of actual distance of 400m makes offsets match QPAD better.
    area = pi*4^2,
    correction = area*avail_est*percept_est,
    offset = log(correction)
  )

# current version of avail and percept are very slow and throw errors if more
# than one species. Have submitted issues for speed up of avail (cue_rate)

# Many warnings are because napops only includs landbirds so many have NA offsets

test <- obs_ON_1991_offsets_a %>% ungroup() %>% slice(1:10) %>% group_by(species) %>%
  mutate(on_road = TRUE,
         edr_est = edr(unique(species), model = "best", road = on_road,
                       forest = for_cov, pairwise = TRUE) %>% pull(EDR_est) %>% .[1],
         avail_est = avail(unique(species), model = "best", od = od, tssr = tssr,
                           pairwise = TRUE, time = 3) %>% pull(p),
         percept_est = percept(unique(species), model = "best", road = on_road,
                              forest = for_cov, distance = 400,
                              pairwise = TRUE) %>% pull(q)) %>%
  ungroup() %>%
  mutate(my_avail = obs_ON_1991_offsets$avail_est[1:10],
         my_percept = obs_ON_1991_offsets$percept_est[1:10],
         my_edr = obs_ON_1991_offsets$edr_est[1:10]) %>%
  select(species, matches("avail"), matches("percept"), matches("edr"))


# TODO: investigate whether the offsets make sense. Ours have a much wider range
# than the ones from BAM did: (-8.6, 11.3) and (-4.3  2.3), respectively. Have
# investigated the formulas and checked for any difference with the native
# napops functions and found none. I also compared visually with the NA-POPS
# dashboard and the values for avail and percept seem similar.

bamEnv <- new.env()
load(file.path(in_dat_pth, "0_data/processed/BAMv6_RoFpackage_2022-01.RData"),
     envir = bamEnv)

sp_off <- get("off", envir = bamEnv)

# Try comparing mean offset by species
bam_off_plt <- colMeans(sp_off) %>% as_tibble(rownames = "species") %>%
  mutate(species = forcats::fct_reorder(species, value)) %>%
  ggplot(aes(species, value))+
  geom_col()+coord_flip()

na_pop_off_plt <- obs_ON_1991_offsets %>% select(RTENO, Year, Stop, species, offset) %>%
  tidyr::pivot_wider(names_from = species, values_from = offset) %>%
  select(-c(RTENO, Year, Stop)) %>%
  colMeans(na.rm = TRUE) %>% as_tibble(rownames = "species") %>%
  filter(!is.na(value), species != "CORE") %>%
  mutate(species = forcats::fct_reorder(species, value))%>%
  ggplot(aes(species, value))+
  geom_col()+
  coord_flip()

# CORE common redpole has an extremely low offset compared to the rest (~-8 vs
# ~6) other wise the napops offsets are much higher than the BAM ones but have a
# similar distribution and a few species I checked are on similar sides of it

off_comp_plt <- ggpubr::ggarrange(bam_off_plt+ggtitle("BAM"),
                                  na_pop_off_plt+ggtitle("NA-POPS"))
ggsave("analysis/figures/compare_offsets_BAM_NAPOP.pdf", height = 15, width = 7,
       units = "in")
