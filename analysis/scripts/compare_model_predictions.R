# compare and evaluate bird model predictions for the Hudson bay lowlands
library(dplyr)
library(terra)
library(sf)
library(purrr)

in_dat_pth <- "analysis/data/derived_data/models_to_compare"

# where caribou ranges are stored, only needed for small example
sourceData <- "C:/Users/endicotts/Documents/gitprojects/MissisaBooPaper/data/inputNV/caribouRanges"

# 4 letter species code and 6 letter species code for species to compare
sp4 <- "ALFL"
sp6 <- "aldfly"

# For example use Missisa as study area
study_area <- read_sf(file.path(sourceData, "Caribou_Range_Boundary.shp")) %>%
  filter(RANGE_NAME == "Missisa")

sa_bbox <- st_bbox(study_area) %>% st_as_sfc()

pats <- list(bateman = paste0(sp4, ".*suitability.tif$"),
             bam_rof = paste0(sp4, ".*mean.tif$"),
             ebird_st = paste0(sp6, "_abundance.*.tif"))

fls <- map(pats, \(x) list.files(in_dat_pth, x, recursive = TRUE))

rasts <- map(fls, \(x) rast(file.path(in_dat_pth, x)))

map(rasts, \(x) print(plot(x)))

# crop and then transform??? maybe just transform bbox and crop to that then transform raster.
# Or better to keep rasters in original format and just transform points used for comparison....?
