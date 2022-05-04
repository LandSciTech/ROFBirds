# Prepare data for model prediction from raw data available in ROFSyncSim
library(tidyverse)
library(sf)

sourceData <- "C:/Users/endicotts/Documents/gitprojects/ROFSyncSim/ROFDemo_data"
in_dat_pth <- "analysis/data/derived_data"

rof <- read_sf(file.path(sourceData, "project_ranges.shp")) %>%
  filter(RANGE_NAME %in% c("Pagwachuan", "Nipigon","Missisa", "James Bay"))
# Use loadSpatialInputs from caribouMetrics to load the data and align crs and
# check rasters are aligned. This uses rasterizeLineDensity which is not
# necessary hear but is needed for caribou so convert from density to binary so
# can use the same data
inData <- caribouMetrics::loadSpatialInputs(
  rof,
  refRast = file.path(sourceData, "plc250.tif"),
  inputsList = list(roads = file.path(sourceData, "road_ORNMNRFROF2020.shp"),
                    esker = file.path(sourceData, "esker.shp")),
  convertToRast = c("esker", "roads")
)

inData$roads <- inData$roads > 0
names(inData$roads) <- "roads"
inData$esker <- inData$esker > 0
names(inData$esker) <- "esker"
plc_layers <- raster::layerize(inData$refRast)

# translate plc values to land cover class names
plc_classes <- read_csv(file.path(sourceData, "ROFBirdModels/plcClasses.csv"))

plc_classes <- plc_classes %>%
  mutate(Class = Class %>% str_replace_all("[^[:alpha:]]", "_") %>%
           str_replace_all("\\_+", "_"))

names(plc_layers) <- plc_classes %>%
  filter(Code %in% raster::unique(inData$refRast)) %>%
  pull(Class)

rastfw250 <- raster::focalWeight(plc_layers, 250, type = "Gauss")
rastfw750 <- raster::focalWeight(plc_layers, 750, type = "Gauss")

plc_layers250 <- map(1:raster::nlayers(plc_layers),
                     ~pfocal::pfocal(plc_layers[[.x]], rastfw250,
                                     transform_function = "MULTIPLY",
                                     reduce_function = "SUM",
                                     mean_divider = "KERNEL_COUNT"))

plc_layers250 <- map(plc_layers250, ~`names<-`(.x, paste0(names(.x), "_250")))

plc_layers750 <- map(1:raster::nlayers(plc_layers),
                     ~pfocal::pfocal(plc_layers[[.x]], rastfw750,
                                     transform_function = "MULTIPLY",
                                     reduce_function = "SUM",
                                     mean_divider = "KERNEL_COUNT"))

plc_layers750 <- map(plc_layers750, ~`names<-`(.x, paste0(names(.x), "_750")))

# make a raster stack of predictors to use for extracting
# BAM's models only use the 750 versions in the end.
pred_stk <- raster::stack(c(plc_layers750, inData$roads, inData$esker))

# Get point locations prepared by BAM
load(file.path(in_dat_pth, "0_data/processed/BAMv6_RoFpackage_2022-01.RData"))

bird_pts <- xx1

rm(xx1, xx2, BB, cn2)

bird_pts <- bird_pts %>% st_as_sf(coords = c("X", "Y")) %>% st_set_crs(4269) %>%
  select(PKEY_V4) %>% st_transform(st_crs(pred_stk))

ext_pt_data <- raster::extract(pred_stk, bird_pts, df = TRUE)

bird_pts <- bind_cols(bird_pts, ext_pt_data)

# get index of rows that are in bird_pts and are not NA. NA means the points did
# not overlap the raster because these rasters are only for the ROF
pkey_keep <- bird_pts %>% filter(!is.na(Ocean_750)) %>%
  pull(PKEY_V4)

inds <- which(rownames(off) %in% pkey_keep)

bird_pts <- select(bird_pts, -PKEY_V4, -ID) %>% st_drop_geometry()

cross_val_out <- map_dfr(SPP[1:2], run_cross_val,
                         samp_row_ind = inds,
                         samp_col_ind = 1, y = y, off =  off, pred_vars = bird_pts,
                         pred_var_nms = names(bird_pts),
                         save_dir = file.path(sourceData, "ROFBirdModels"))

