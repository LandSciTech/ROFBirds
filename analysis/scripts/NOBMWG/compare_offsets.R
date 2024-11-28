library(tidyverse)
library(wildRtrax)
library(napops)

wt_auth()

projects <- wt_get_download_summary("PC")

dat_aru <- projects %>%
  # pick a project in boreal and with few tasks
  filter(project == "CWS-Ontario Atlas Digital Point Counts Boreal FMUs 2022") %>%
  pull(project_id) %>%
  wt_download_report(sensor_id = "ARU", reports = "main", weather_cols = FALSE)

dat_pc <- projects %>%
  # pick a project in boreal and with few tasks
  filter(project == "CWS-Ontario Atlas Digital Point Counts Boreal FMUs 2022") %>%
  pull(project_id) %>%
  wt_download_report(sensor_id = "PC", reports = "main", weather_cols = FALSE)

dat_pc <- dat_pc %>% wt_tidy_species(sensor = "PC", zerofill = TRUE)

dat_pc_clean <- dat_pc %>%
  wt_make_wide(sensor = "PC")

dat_aru_clean <- dat_aru %>% wt_tidy_species(sensor = "ARU", zerofill = TRUE) %>%
  wt_replace_tmtt() %>%
  wt_make_wide()

# in this case the results are the same whether you use the ARU data and clean
# as recommended or use the PC version except a few different or modified column
# names
# comp <- left_join(dat_aru_clean,
#           dat_pc_clean,
#           by = join_by(organization, project_id, location, location_id, location_buffer_m, longitude,
#                        latitude, recording_date_time == survey_date),
#           suffix = c("_aru", "_pc")) %>%
#   pivot_longer(cols = matches("[A-Z]", ignore.case = FALSE), names_to = c("species", "sensor"),
#                names_sep = "_", values_to = "count") %>%
#   pivot_wider(names_from = sensor, values_from = count) %>%
#   filter(aru != pc)

# QPAD offsets
dat_pc_off <- wt_qpad_offsets(dat_pc_clean, sensor = "PC")
# dat_aru_off <- wt_qpad_offsets(dat_aru_clean, sensor = "ARU")

# all(dat_aru_off == dat_pc_off)

# napops offsets
# fetch_data()

# Need to get covariates

# Use NALCMS data same as napops paper. I am using 2020 data could in theory use
# version closest to actual year of data collection
options(reproducible.cachePath = "analysis/data/cache_data")
url_nalcms_2020_Can <- "http://www.cec.org/files/atlas_layers/1_terrestrial_ecosystems/1_01_0_land_cover_2020_30m/can_land_cover_2020_30m_tif.zip"
nalcms_2020_Can <- reproducible::prepInputs(url = url_nalcms_2020_Can,
                                            destinationPath = "analysis/data/raw_data/",
                                            writeTo = NULL)

lcc_lu <- tibble::tribble(
  ~from, ~to, ~becomes,
  1,       2,        1,
  5,       6,        1)

# get boundary of bird data + ~2000m so window doesn't include edges
SA <- dat_pc %>% ungroup() %>%
  summarise(across(c(longitude, latitude), lst(min = \(x) min(x, na.rm = TRUE),
                                               max = \(x) max(x, na.rm = TRUE)))) %>%
  select(xmin = longitude_min, ymin = latitude_min, xmax = longitude_max,
         ymax = latitude_max) %>% unlist() %>%
  sf::st_bbox() %>%
  {. + 0.02} %>%
  sf::st_as_sfc() %>%
  sf::st_set_crs(4326) %>%
  sf::st_transform(sf::st_crs(nalcms_2020_Can))

# Crop and convert to prop forest in 5x5 window
for_cov <- reproducible::Cache(nalcms_2020_Can %>%
                   terra::crop(terra::vect(SA)) %>%
                   terra::classify(rcl = lcc_lu, others = 0, right = NA))
for_cov150 <- reproducible::Cache(terra::focal(for_cov, w = 5, fun = "mean"))

obs_for_cov <- dat_pc %>% filter(!is.na(latitude)) %>%
  select(project_id, location_id, survey_id, species_code, latitude, longitude) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326) %>%
  sf::st_transform(sf::st_crs(nalcms_2020_Can)) %>%
  {terra::extract(for_cov150, ., bind = TRUE)} %>%
  sf::st_as_sf() %>%
  sf::st_drop_geometry() %>%
  rename(for_cov = focal_mean)

dat_pc_cov <- dat_pc %>%
  filter(!is.na(latitude)) %>%
  left_join(obs_for_cov, by = join_by(project_id, location_id, survey_id, species_code)) %>%
  # time zone says UTC but is really unknown
  mutate(timez = lutz::tz_lookup_coords(latitude, longitude, method = "accurate")) %>%
  mutate(survey_date = lubridate::force_tzs(survey_date, timez, tzone_out = "UTC")) %>%
  mutate(sunrise = suncalc::getSunlightTimes(
  data = data.frame(date = lubridate::as_date(survey_date), lon = longitude, lat = latitude),
  keep = "sunrise")$sunrise) %>%
  mutate(tssr = difftime(survey_date, sunrise, units = "hours") %>% as.numeric(),
         od = lubridate::yday(survey_date))


# Needs to be available in order for faster versions of functions to work
# preloading a table of species min, max median for each model type
sp_tbl <- list_species()

sum_covs <- function(df){
  df %>%
    summarise(across(-c(Species, Survey_Method), lst(min, max, median)))
}

all_cov_tbl <- reproducible::Cache(sp_tbl %>% rowwise() %>%
                       mutate(sum_covs(covariates_removal(Species)),
                              sum_covs(covariates_distance(Species))))

source("R/napops_mod_functions.R")

dat_pc_off_napops <- dat_pc_cov %>% filter(!is.na(for_cov)) %>%
  mutate(on_road = FALSE) %>%
  group_by(species_code) %>%
  mutate(cr_est = cue_rate3(unique(species_code), model = "best", od, tssr,
                            pairwise = TRUE) %>% pull(CR_est),
         edr_est = edr2(unique(species_code), model = "best", road = on_road, forest = for_cov,
                        pairwise = TRUE) %>% pull(EDR_est))



nap_out <- dat_pc_cov %>%
  select(species_code, od, tssr, for_cov, detection_distance,
         survey_duration_method, survey_date) %>%
  slice(2) %>%
  mutate(
    on_road = FALSE,
    avail_off = avail(species_code, "best", od = od, tssr = tssr,
                      time = 5, pairwise = TRUE)$p,
    precep_off = percept(species_code, "best", forest = for_cov, road = on_road,
                         distance = 400, pairwise = TRUE)$q,
    cr_off = cue_rate(species_code, "best", od = od, tssr = tssr,
                      pairwise = TRUE)$CR_est[,1],
    edr_off = edr(species_code, "best", forest = for_cov, road = on_road,
                  pairwise = TRUE)$EDR_est[,1],
    off = log(pi*400^2)+log(avail_off)+log(precep_off))


wt_qpad_offsets(dat_pc_clean %>% slice(2), species = "REVI", sensor = "PC")
wt_qpad_offsets(dat_aru_clean %>% slice(2), species = "REVI", sensor = "ARU")
# Don't know what is going on the rasters for the covariates are not working

# when distance is INF then A is based on intercept from edr

# try with maxdis = 400
debug(wt_qpad_offsets)

qpad_off_x$MAXDUR <- 300
qpad_off_x$MAXDIS <- 400

wildRtrax:::.make_off("REVI", qpad_off_x)

nap_out %>% select(contains("off"))
qpad_off

# why is p 1? why is edr so different? Need to check the paper to compare the math notation.
