# ******************************************************************************
# ******************************************************************************
# Implementation of David Isle's Landscape Distribution Modeling ECCC workflow
# for Northern Ontario: Covariate download
# ******************************************************************************
# ******************************************************************************
# Includes downloading rasters and extracting data from GEE
# ******************************************************************************
# Notable changes of method:
# * Added GEE data extraction using BAMs variable list.
#   https://docs.google.com/spreadsheets/d/1XATuq8BOYC2KkbJObaturaf4NFD2ufvn/edit?usp=sharing&ouid=104837701987164094932&rtpof=true&sd=true
# ******************************************************************************
library(sf)
library(tidyverse)
library(terra)
library(exactextractr)
library(data.table) # used in gee_points_extraction
library(rgee)

devtools::load_all(".")
# Get data from previous scripts
Study_Area_bound <- read_sf("analysis/data/derived_data/INLA_data/study_area.gpkg")

analysis_data <- readRDS(file = "analysis/data/interim_data/analysis_data_use.rds")

all_surveys <- analysis_data$all_surveys
full_count_matrix <- analysis_data$full_count_matrix

rm(analysis_data)

# use reproducible package to Cache long running code
options(reproducible.cachePath = "analysis/data/cache_data")

# Prepare Covariates #===========================================================

## Download and crop #==========================================================

# Dataframe to store covariates
all_surveys_covariates <- all_surveys %>% dplyr::select(Obs_Index)

# 1 km x 1 km grid
# TODO check this transformation makes sense
ONGrid <- reproducible::Cache(st_make_grid(
  Study_Area_bound,
  cellsize = units::set_units(10*10,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE
) %>%
  st_as_sf() %>%
  st_intersection(Study_Area_bound) %>%
  na.omit() %>%
  # TODO: Not ideal, dropping several areas that are complex due to BCR intersection
  # Needed for conversion to sp
  st_set_geometry("geometry") %>%
  filter(st_geometry_type(geometry) == "POLYGON") %>%
  st_cast())

ONGrid$point_id <- 1:nrow(ONGrid)
ONGrid_centroid <- st_centroid(ONGrid)

# Add 20 km buffer so that covariates can extend slightly outside province if possible
SABoundary_buffer <- Study_Area_bound %>% st_buffer(20000)

# Download covariate layers
url_nalcms_2020_Can <- "http://www.cec.org/files/atlas_layers/1_terrestrial_ecosystems/1_01_0_land_cover_2020_30m/can_land_cover_2020_30m_tif.zip"
nalcms_2020_Can <- reproducible::Cache(reproducible::prepInputs(
  url = url_nalcms_2020_Can,
  destinationPath = "analysis/data/raw_data/",
  targetFile = "analysis/data/raw_data/CAN_NALCMS_landcover_2020_30m.tif",
  alsoExtract = NA,
  fun = "terra::rast",
  studyArea = SABoundary_buffer,
  # this doesn't quite work as expected, ends up in raw_data/NALCMS_NorthON
  writeTo = "analysis/data/derived_data/covariates/NALCMS_NorthON.tif",
  useCache = TRUE
))

terra::writeRaster(nalcms_2020_Can, "analysis/data/derived_data/covariates/NALCMS_NorthON.tif",
                   overwrite = TRUE)

url_adaptwest_norm <- "https://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6v73/normals/Normal_1961_1990_bioclim.zip"

clim_norm <- reproducible::prepInputs(
  url = url_adaptwest_norm,
  destinationPath = "analysis/data/raw_data/",
  targetFile = "Normal_1961_1990/Normal_1961_1990_bioclim/Normal_1961_1990_MAT.tif"
)

mat_SA <- terra::crop(clim_norm,
                      terra::vect(SABoundary_buffer %>% st_transform(st_crs(clim_norm))),
                      mask = TRUE, overwrite = TRUE,
                      filename = "analysis/data/derived_data/covariates/mat_NorthON.tif")

# Path to spatial covariates
covar_folder <- "analysis/data/derived_data/covariates/"

# TODO: get these variables?
# # Annual mean temperature
# AMT <- rast(paste0(covar_folder,"National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif")) %>%
#   crop(st_transform(ONBoundary_buffer,crs(.))) %>%
#   project(st_as_sf(ONBoundary_buffer), res = 250)
#
# # Land cover of Canada 2020 (not reprojected... faster to reproject grid)
# lcc2020 <- rast(paste0(covar_folder,"National/LandCoverCanada2020/landcover-2020-classification.tif")) %>%
#   crop(st_transform(ONBoundary_buffer,crs(.)))
#
# # Stand canopy closure
# SCC <- rast(paste0(covar_folder,"National/NationalForestInventory/NFI_MODIS250m_2011_kNN_Structure_Stand_CrownClosure_v1.tif")) %>%
#   crop(st_transform(ONBoundary_buffer,crs(.))) %>%
#   project(st_as_sf(ONBoundary_buffer), res = 250)


# Extract survey covariates within 1 km of location #===============

# I had created a stack but that requires same res which loses a lot of
# information... ignore for now

# # Create a stack of covariate rasters
# covars <- covar_folder %>% list.files(full.names = TRUE) %>%
#   set_names(list.files(covar_folder) %>% str_remove("_NorthON\\.tif")) %>%
#   map(rast)
#
# # TODO resampling to MAT is a bad idea... should keep at native resolution for extraction
# covars$NALCMS <- covars$NALCMS %>% resample(covars$mat, method = "near")
#
# covars <- rast(covars)
#
# terra::writeRaster(covars, "analysis/data/derived_data/INLA_data/covariate_stack.grd",
#                    overwrite = TRUE)

# Note Using terra very slow for fractional coverage

# Continuous covariates
all_surveys_1km <- all_surveys_covariates %>%
  st_buffer(1000) %>%
  mutate(#elevation_1km = exact_extract(elevation, ., 'mean'),
    AMT_1km = exact_extract(mat_SA, ., 'mean'))#,
#SCC_1km = exact_extract(SCC, ., 'mean'))


# Proportion of each land cover class
prop_LCC_1km <- exact_extract(nalcms_2020_Can, all_surveys_1km,"frac") %>%
  suppressWarnings()
names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km <- prop_LCC_1km %>% dplyr::select(sort(names(.)))

# TODO: this is not very safe could easily miss one in a new area, would be
# better as look up table

# Fix land cover class names
prop_LCC_1km <- prop_LCC_1km %>%
  mutate(
    Needleleaf_forest_1km = LCC_1_1km + LCC_2_1km,
    Mixed_forest_1km = LCC_5_1km + LCC_6_1km,
    Grass_shrub_1km = LCC_8_1km + LCC_10_1km,
    Barren_lichen_1km = LCC_13_1km + LCC_16_1km,
    Crop_1km = LCC_15_1km,
    Urban_1km = LCC_17_1km,
    Wetland_1km = LCC_14_1km,
    Water_1km = LCC_18_1km) %>%
  dplyr::select(Needleleaf_forest_1km:Water_1km)

all_surveys_1km <- bind_cols(all_surveys_1km, prop_LCC_1km) %>%
  st_drop_geometry()

# Join with survey dataset top get back to points
all_surveys_covariates <- left_join(all_surveys_covariates, all_surveys_1km,
                                    by = join_by(Obs_Index))

# Do GEE covariate extraction #=================================================
visit <- all_surveys %>%
  sf::st_drop_geometry() %>%
  # rename to fit WT format
  dplyr::mutate(
    project_id = Project_Name, location_id = survey_ID,
    latitude = Latitude, longitude = Longitude,
    year = lubridate::year(Date_Time), Obs_Index, .keep = "none"
  )

visit_yr <- visit %>%
  dplyr::select(project_id, location_id, latitude, longitude, year) %>%
  dplyr::mutate(year = as.integer(year)) %>%
  unique() %>%
  dplyr::mutate(id = paste(project_id, location_id, latitude, longitude, year, sep = "_")) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = 5072)

BAM_var_list <- read.csv(file.path("analysis/data/raw_data/", "Model_VariableList.csv"))

# To force download of all data set Complete to 0
# BAM_var_list$Complete <- 0
#
# write.csv(BAM_var_list, file.path("analysis/data/raw_data/", "Model_VariableList.csv"),
#           row.names = FALSE)

# initialize rgee
ee$Authenticate(auth_mode='notebook')
ee$Initialize(project='ee-sarahendicott-eccc')  # <-- EDIT THIS FOR YOUR PROJECT

# Optionally
# make a request to verify you are connected.
ee$String('Hello from the Earth Engine servers!')$getInfo()


# takes around 3 hours to run
tictoc::tic()
# Function uses GEE and the BAM variables list to extract BAM variables that are
# available from GEE including for year matched time series variables
gee_points_extract(meth_path = file.path("analysis/data/raw_data/", "Model_VariableList.csv"),
                   loc.yr = visit_yr,
                   out_save = "analysis/data/derived_data/INLA_data/survey_covariates")

tictoc::toc()

static_dat <- read.csv("analysis/data/derived_data/INLA_data/survey_covariates/data_covariates_GEE_static.csv")
match_dat <- read.csv("analysis/data/derived_data/INLA_data/survey_covariates/data_covariates_GEE_match.csv")

# Here following the national models code:
# Remove lat lon fields due to rounding errors that cause mismatches
# Format landcover classes as factors
# Zero out heights < 0.1 and NA the height cv values for those
# Remove 0 values for some layers that should be NAs
# Format landcover classes as factors
# Create eBird method and make method a factor
BAM_var_list <- read.csv(file.path("analysis/data/raw_data/", "Model_VariableList.csv"))

meth.use <- BAM_var_list %>%
  dplyr::filter(Complete==1)

# TODO: Steps for variables not yet extracted are commented out below. These
# should really be programmed from the variables spreadsheet so that they can be
# applied even when some variables are missing
visit.covs <- visit %>%
  mutate(id=paste(project_id, location_id, latitude, longitude, year, sep="_"),
         year=as.integer(year)) %>%
  left_join(static_dat %>% select(-c(latitude, longitude))) %>%
  left_join(match_dat %>% select(-c(latitude, longitude))) %>%
  dplyr::select(-id) %>%
  dplyr::select(all_of(colnames(visit)), any_of(meth.use$Label)) %>%
  mutate(#across(c(NLCD_1km, VLCE_1km),~ifelse (.x == 0, NA, .x)),
         across(c(#SCANFIheight_1km, SCANFIheight_5x5,
                  ETHheight_1km, ETHheight_5x5),
                  #LFheigth_1km, LFheigth_5x5),
                ~ifelse (.x < 0.1, 0, .x)),
         #SCANFIheightcv_1km = ifelse (SCANFIheight_1km < 0.1, NA, SCANFIheightcv_1km),
         #SCANFIheightcv_5x5= ifelse (SCANFIheight_5x5 < 0.1, NA, SCANFIheightcv_5x5),
         ETHheightcv_1km = ifelse (ETHheight_1km < 0.1, NA, ETHheightcv_1km),
         ETHheightcv_5x5 = ifelse (ETHheight_5x5 < 0.1, NA, ETHheightcv_5x5),
         #LFheigthcv_1km = ifelse (LFheigth_1km < 0.1, NA, LFheigthcv_1km),
         #LFheigthcv_5x5= ifelse (LFheigth_5x5 < 0.1, NA, LFheigthcv_5x5),
         #Test this
         # across(c(SCANFIheightcv_1km, SCANFIheightcv_5x5, ETHheightcv_1km, ETHheightcv_5x5,
         #         LFheigthcv_1km, LFheigthcv_5x5), ~ifelse (c(SCANFIheight_1km,
         #         SCANFIheight_5x5, ETHheight_1km, ETHheight_5x5, LFheigth_1km,
         #         LFheigth_5x5) < 0.1, NA, .x))
         across(c(MODISLCC_1km, MODISLCC_5x5),
                  #SCANFI_1km, NLCD_1km,
                  #ABoVE_1km, VLCE_1km),
                as.factor))

all_surveys_covariates <- all_surveys_covariates %>%
  left_join(visit.covs %>% select(-c(project_id, location_id, latitude, longitude, year)),
            by = "Obs_Index")

# Extract grid covariates within 1 km of location #=============================

# Continuous covariates
ONGrid_1km <- ONGrid_centroid %>%
  st_buffer(1000) %>%
  mutate(#elevation_1km = exact_extract(elevation, ., 'mean'),
    AMT_1km = exact_extract(mat_SA, ., 'mean'))#,
# SCC_1km = exact_extract(SCC, ., 'mean'))

# TODO: avoid this repetition
# Proportion of each land cover class
prop_LCC_1km <- exact_extract(nalcms_2020_Can, ONGrid_1km,"frac") %>% suppressWarnings()
names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km <- prop_LCC_1km %>% dplyr::select(sort(names(.)))

# Fix land cover class names
prop_LCC_1km <- prop_LCC_1km %>%
  mutate(
    Needleleaf_forest_1km = LCC_1_1km + LCC_2_1km,
    Mixed_forest_1km = LCC_5_1km + LCC_6_1km,
    Grass_shrub_1km = LCC_8_1km + LCC_10_1km,
    Barren_lichen_1km = LCC_13_1km + LCC_16_1km,
    Crop_1km = LCC_15_1km,
    Urban_1km = LCC_17_1km,
    Wetland_1km = LCC_14_1km,
    Water_1km = LCC_18_1km) %>%
  dplyr::select(Needleleaf_forest_1km:Water_1km)

ONGrid_1km <- bind_cols(ONGrid_1km, prop_LCC_1km)

# Re-initialize covariate list for grid extraction.
# To force download of all data set Complete to 0
BAM_var_list$Complete <- 0

write.csv(BAM_var_list, file.path("analysis/data/raw_data/", "Model_VariableList.csv"),
          row.names = FALSE)
tictoc::tic()
gee_points_extract(meth_path = file.path("analysis/data/raw_data/", "Model_VariableList.csv"),
                   loc.yr = ONGrid_centroid %>% mutate(year = 1, id = point_id),
                   year_get = "max",
                   out_save = "analysis/data/derived_data/INLA_data/grid_covariates")
tictoc::toc()

static_grid_dat <- read.csv("analysis/data/derived_data/INLA_data/grid_covariates/data_covariates_GEE_static.csv")
match_grid_dat <- read.csv("analysis/data/derived_data/INLA_data/grid_covariates/data_covariates_GEE_match.csv")

# same as above but dropped the commented out parts.
ONGrid_gee_covs <- ONGrid_centroid %>%
  left_join(static_grid_dat, by = join_by(point_id)) %>%
  left_join(match_grid_dat %>% select(-id, -year), by = join_by(point_id)) %>%
  dplyr::select(-id) %>%
  dplyr::select(all_of(colnames(ONGrid_centroid)), any_of(meth.use$Label)) %>%
  mutate(across(c(ETHheight_1km, ETHheight_5x5),
                ~ifelse (.x < 0.1, 0, .x)),
         ETHheightcv_1km = ifelse (ETHheight_1km < 0.1, NA, ETHheightcv_1km),
         ETHheightcv_5x5 = ifelse (ETHheight_5x5 < 0.1, NA, ETHheightcv_5x5),
         across(c(MODISLCC_1km, MODISLCC_5x5),
                as.factor))


# Join with ONGrid
ONGrid <- ONGrid %>% left_join(st_drop_geometry(ONGrid_1km), by = join_by(point_id)) %>%
  left_join(st_drop_geometry(ONGrid_gee_covs), by = join_by(point_id))

# Covariate PCA #===============================================================

# Conduct Principal Components Analysis on covariates to identify axes of major
# variation in habitat

# Remove surveys with no covariate information
surveys_to_remove <- st_drop_geometry(all_surveys_covariates)%>%
  select(-Obs_Index, -matches("MODISLCC")) %>%  # dropping Modis because can't have factor here
  select(-where(\(x) all(is.na(x))|all(x == 0)))

surveys_to_remove <- which(is.na(rowSums(surveys_to_remove))) # drops 109 surveys that have NA for any variable

if(length(surveys_to_remove) > 0){
  all_surveys_covariates <- all_surveys_covariates[-surveys_to_remove,]
  all_surveys <- all_surveys[-surveys_to_remove,]
  full_count_matrix <- full_count_matrix[-surveys_to_remove,]
}

covars_for_PCA <- all_surveys_covariates %>%
  st_drop_geometry() %>%
  select(-Obs_Index, -matches("MODISLCC")) %>%
  # dplyr::select(AMT_1km:Water_1km)
  select(-where(\(x){v <- var(x); is.na(v)|v == 0})) # drops variables with any NAs which may not be necessary

pca <- prcomp(covars_for_PCA, scale = TRUE)

# Interpretation of specific axes (e.g., axes 1 and 2)

summary(pca)   # Proportion variance explaind by axes
factoextra::fviz_eig(pca)  # Scree plot (first 5 axes explain 85% of variation in habitat between sites)
pca            # Variable loadings

factoextra::fviz_pca_var(pca,
             axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = viridis::viridis(10, direction = -1),
             repel = TRUE     # Avoid text overlapping
)

factoextra::fviz_pca_var(pca,
                         axes = c(1,3),
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = viridis::viridis(10, direction = -1),
                         repel = TRUE     # Avoid text overlapping
)

# Predict PCA values for each survey location and standardize (mean = 0, sd = 1)

all_surveys_PCA <- predict(pca, newdata = st_drop_geometry(all_surveys_covariates)[,names(pca$center)])
ONGrid_PCA <- predict(pca, newdata = st_drop_geometry(ONGrid)[,names(pca$center)])

# TODO: use across and scale
for (covar in colnames(all_surveys_PCA)){

  covar_mean <- mean(all_surveys_PCA[,covar],na.rm = TRUE)
  covar_sd <- sd(all_surveys_PCA[,covar],na.rm = TRUE)

  all_surveys_PCA[,covar] <- (as.data.frame(all_surveys_PCA)[,covar] - covar_mean)/covar_sd
  ONGrid_PCA[,covar] <- (as.data.frame(ONGrid_PCA)[,covar] - covar_mean)/covar_sd

}

all_surveys_covariates <- all_surveys_covariates %>%
  st_drop_geometry() %>%
  bind_cols(all_surveys_PCA)


all_surveys <- full_join(all_surveys,all_surveys_covariates)
ONGrid <- bind_cols(ONGrid,ONGrid_PCA)

# Save #========================================================================

analysis_data_package <- list(

  all_surveys = all_surveys, # Survey information (protocol, location, time, covariates)
  full_count_matrix = full_count_matrix, # counts of each species for each survey

  pca = pca,

  # Contains covariates on ON-wide grid
  ONGrid = ONGrid

)

saveRDS(analysis_data_package,"analysis/data/derived_data/INLA_data/analysis_data_covars.rds")
