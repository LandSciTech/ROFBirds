# there were some issues in bbsAssistant pacakge but it does not appear to be
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

obs_ON_1991_2 <- obs_ON_1991 %>%
  filter(BCR %in% bcrs) %>%
  filter(!!rlang::parse_expr(quality_filter)) %>%
  # get 4 letter codes for sp names
  left_join(species_list %>% select(AOU, AOU4), by = "AOU") %>%
  # Some are missing so keep numeric codes for those
  mutate(sp_code = coalesce(AOU4, as.character(AOU)), .keep = "unused") %>%
  select(-contains("Car"), -contains("Noise"), -RouteTotal, -TotalSpp,
         -contains("ObsFirstYear"), -Assistant, -RunType, -QualityCurrentID,
         -julian) %>%
  tidyr::pivot_longer(contains("Stop"), names_to = "Stop", values_to = "Abund") %>%
  filter(Abund != 0) %>%
  tidyr::pivot_wider(names_from = sp_code, values_from = Abund, values_fill = 0)


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
  summarise(across(c(POINT_X, POINT_Y), lst(min, max), na.rm = TRUE)) %>%
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
# A is the known or estimated area of survey, p is availability given presence, q is detectability given avaibalility.
# offset=log(A) + log(p) + log(q)
obs_ON_1991_offsets <- obs_ON_1991_stops_3 %>%
  select(-c(
    Active, Latitude, Longitude, Stratum, BCR, RouteTypeID, RouteTypeDetailID,
    ObsN, Month, Day, StartTemp, EndTemp, TempScale, StartWind, EndWind,
    StartSky, EndSky, StartTime, EndTime, WindMean, TempMean, timez, totTime,
    nStops, stopInc, incTime, stopTime1, stopTime2, sunrise, mean_for_cov)) %>%
  mutate(tssr = as.numeric(tssr)) %>%
  tidyr::pivot_longer(-c(CountryNum, StateNum, Route, RPID, Year, RTENO, RouteName,
                  Date, Stop, POINT_X, POINT_Y, tssr, od, for_cov),
               names_to = "species", values_to = "count") %>%
  filter(count != 0) %>%
  slice(1:1000)

out <- obs_ON_1991_offsets %>% group_by(species) %>%
  mutate(on_road = TRUE,
    p = napops::avail(unique(species), model = "best", od, tssr, time = 3,
                           pairwise = TRUE),
         q = napops::percept(unique(species), model = "best", road = on_road,
                             forest = for_cov,
                             distance = 400, pairwise = TRUE))

# current version of avail and percept are very slow and throw errors if more
# than one species. one speed up issue submitted considering another one.


