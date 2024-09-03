# ******************************************************************************
# ******************************************************************************
# Implementation of David Isle's Landscape Distribution Modeling ECCC workflow
# for Northern Ontario: Running models
# ******************************************************************************
# ******************************************************************************
# Includes Fitting models, generating predictions, mapping predictions and model
# validation through cross-validation
# ******************************************************************************
# Notable changes of method:
# ******************************************************************************

library(INLA)
library(inlabru)
library(tidyverse)
library(sf)
library(ggtext)


devtools::load_all(".")

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)

Study_Area_bound <- read_sf("analysis/data/derived_data/INLA_data/study_area.gpkg")

analysis_data <- readRDS("analysis/data/derived_data/INLA_data/analysis_data_package.rds")
# attach(analysis_data)

covars <- terra::rast("analysis/data/derived_data/INLA_data/covariate_stack.grd")

# Raster with target properties
target_raster <- covars$mat

# Atlas Squares # 10 km x 10 km grid
ONSquares <- st_make_grid(
  Study_Area_bound,
  cellsize = units::set_units(10*10,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE)%>%
  st_as_sf() %>%
  st_intersection(Study_Area_bound) %>%
  na.omit() %>%
  mutate(sq_id = 1:nrow(.)) %>%
  rename(geometry = x)

BCR_PROV <- st_read("analysis/data/raw_data/BCR_Terrestrial/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  group_by(BCR,PROVINCE_S)%>%
  summarize(geometry = st_union(geometry)) %>%
  st_transform(st_crs(Study_Area_bound))

ONGrid <- analysis_data$ONGrid %>%
  st_intersection(BCR_PROV) %>%
  dplyr::rename(geometry = x)

# Loop through species, fit models, generate maps
species_to_model <- analysis_data$species_to_model %>%
  arrange(desc(n_squares),desc(n_detections))

# adding a more southern dispersed species to compare
species_to_fit <- species_to_model %>%
  subset(Species_Code_BSC %in% c("SOSA","LEYE","GRYE", "DEJU"))

# This is a functionified version of what David had in 4_Analysis_INLA
for (sp_code in species_to_fit$Species_Code_BSC){

  sp_mod <- fit_inla(sp_code, analysis_data, st_crs(Study_Area_bound), Study_Area_bound)

  sp_pred <- predict_inla(analysis_data, sp_mod)

  map_inla_preds(sp_code, analysis_data, sp_pred, st_crs(Study_Area_bound), ONSquares)

}

# Evaluate model performance by splitting data into modeling and test sets
# Splits of interest include:
# # train on old data (eg atlas 2) test on far north
# # Spatial cross-validation to test overall performance David's
#   5_Model_Comparison_Crossval this seems to fit with all but one grid square
#   and test on last grid. Another option could be to make split based on
#   spatially stratified random sample so hopefully both samples have fairly
#   full representation.
