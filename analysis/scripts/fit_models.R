# Sarah's version of the bird modeling process

# Steps to the analysis
# 1 Get the environmental and bird data in a format that can be modelled
# 2 Run GBM once to determine which variables to use
# 3 Run GBM nBootstrap times with important variables
# 4 Make predictions for each bootstrapped model and take the mean

library(tidyverse)
library(fs)
library(sf)

in_dat_pth <- "analysis/data/derived_data"
sourceData <- "C:/Users/endicotts/Documents/gitprojects/ROFSyncSim/ROFDemo_data"

# For example use Missisa as study area
study_area <- read_sf(file.path(sourceData, "project_ranges.shp")) %>%
  filter(RANGE_NAME == "Missisa")

# PARAMETERS to set #====================================
# Number of bootstrap samples to use. A number between 2 and 250
nBootstrap <- 2

# Number of species to make models for. A number between 1 and 89
nSPP <- 2

# The relative influence at which to cut off the variables included in the final
# model. A number between 0 and 100. Make if higher for models with fewer
# variables
rel_inf_thold <- 5

# Load BAM inputs #=======================================

# unzip files because Google Drive broke them down into many zip files for download
#
# file_ls <- list.files(in_dat_pth, full.names = TRUE)
#
# walk(file_ls, unzip, overwrite = FALSE)

# load data that BAM prepared and rename objects to be literate
# data saved at line 196 Andy_scripts/01-2_create-brt-database.R
load(file.path(in_dat_pth, "0_data/processed/BAMv6_RoFpackage_2022-01.RData"))

# Rename objects to be more literate

# BB is 2916, 250 matrix each col seems to be numbers b/w 1 and 11053. It is 250
# different samples of points with the same number in hudson plains and boreal
# shield. The first is all points in hudson plains and a sample w/o replacement
# from boreal shield other columns are samples with replacement
pt_resamps_250 <- BB

# cn2 222 vector of variable names all are cols in xx2 but it has other cols too
# It is cols in xx2 that have SD above 0.0001 (84 removed) and correlation > 0.9
# (54 removed) when a pair are correlated the one in the row is kept a subset is
# then hard coded in below (This is from a reduced list of 110 variables that
# BAM came up with after discussions with the RoF team and our own internal
# discussions)
vars_sel <- cn2

# off 11053, 89 rownames PKEY point identifier also in xx1, colnames are sp names
# These are offsets.

# SPP 89 vector of species codes. Filtered to those with >20 detections

# xx1 11053, 25 meta data on points like PKEY, survey year, BCR,
#   aou YN, study area YN and dist (distance to study area boundary)
pt_meta <- xx1

# xx2 11053, 360 PKEY as rownames colnames are variables some not in cn2. Values
# of predictor variables extracted for each point
pt_vars <- xx2

# y sparse matrix PKEY as rownames, colnames are species, values are b/w 0 and 80
# Not clear where it comes from I think it gives species presence, possibly density

rm(xx1, xx2, BB, cn2)

vars_sel <- c("eskerpoint",
              "agriculture_G750.O", "bedrock_G750.O", "biomass2015.ntems",
              "bog_G750.O", "communities_G750.O", "coniftreed_G750.O", "decidtreed_G750.O",
              "disturbance_G750.O", "elev", "fen_G750.O", "G750LandCover_Veg_v1.grd",
              "G750LandCover_VegNonTreed_v1.grd", "G750LandCover_VegTreed_v1.grd",
              "G750Species_Abie_Bal_v1.grd", "G750Species_Acer_Neg_v1.grd",
              "G750Species_Acer_Pen_v1.grd", "G750Species_Acer_Rub_v1.grd",
              "G750Species_Acer_Sac_v1.grd", "G750Species_Acer_Sah_v1.grd",
              "G750Species_Acer_Spi_v1.grd", "G750Species_Acer_Spp_v1.grd",
              "G750Species_Alnu_Spp_v1.grd", "G750Species_Betu_All_v1.grd",
              "G750Species_Betu_Pap_v1.grd", "G750Species_Betu_Pop_v1.grd",
              "G750Species_Fagu_Gra_v1.grd", "G750Species_Frax_Ame_v1.grd",
              "G750Species_Frax_Nig_v1.grd", "G750Species_Frax_Pen_v1.grd",
              "G750Species_Genc_Spp_v1.grd", "G750Species_Genh_Spp_v1.grd",
              "G750Species_Lari_Lar_v1.grd", "G750Species_Pice_Abi_v1.grd",
              "G750Species_Pice_Gla_v1.grd", "G750Species_Pice_Mar_v1.grd",
              "G750Species_Pice_Rub_v1.grd", "G750Species_Pinu_Ban_v1.grd",
              "G750Species_Pinu_Con_v1.grd", "G750Species_Pinu_Res_v1.grd",
              "G750Species_Pinu_Str_v1.grd", "G750Species_Popu_Bal_v1.grd",
              "G750Species_Popu_Gra_v1.grd", "G750Species_Popu_Spp_v1.grd",
              "G750Species_Popu_Tre_v1.grd", "G750Species_Prun_Pen_v1.grd",
              "G750Species_Quer_Mac_v1.grd", "G750Species_Quer_Rub_v1.grd",
              "G750Species_Sali_Spp_v1.grd", "G750Species_Sorb_Ame_v1.grd",
              "G750Species_Thuj_Occ_v1.grd", "G750Species_Tili_Ame_v1.grd",
              "G750Species_Tsug_Can_v1.grd", "G750Species_Ulmu_Ame_v1.grd",
              "G750SpeciesGroups_Broadleaf_Spp_v1.grd", "G750SpeciesGroups_Needleleaf_Spp_v1.grd",
              "G750SpeciesGroups_Unknown_Spp_v1.grd", "G750Structure_Biomass_TotalDead_v1.grd",
              "G750Structure_Stand_Age_v1.grd", "G750Structure_Volume_Total_v1.grd",
              "heath_G750.O", "height2015.ntems", "LIDARheight", "marsh_G750.O",
              "mixedtreed_G750.O", "mudflat_G750.O", "openwater_G750.O", "road_yesno",
              "slope", "sparsetreed_G750.O", "swamp_G750.O", "TPI", "treecover",
              "turbidwater_G750.O", "volume2015.ntems", "ecozone")

pt_vars <- mutate(pt_vars, ecozone = ifelse(pt_meta$ecozone == "hudson_plain", 1, 0))

# Train models #==============================

# train a model with cross validation to get best ntrees and important variables
cross_val_out <- map_dfr(SPP[1:nSPP], run_cross_val,
                         samp_row_ind = pt_resamps_250[,1],
                         samp_col_ind = 1, y = y, off =  off, pred_vars = pt_vars,
                         pred_var_nms = vars_sel,
                         save_dir = file.path(in_dat_pth, "results"))

# save as a csv
cross_val_out %>% unnest(rel_inf) %>%
  write_csv(file.path(in_dat_pth, "results/cross_val_results.csv"))

# fit models for bootstrap samples using ntrees and important
# variables from cross val run
boot_out <- pmap(cross_val_out, run_boot,
                 rel_inf_thold = rel_inf_thold,
                 resamps = as.data.frame(pt_resamps_250[, 1:nBootstrap]),
                 y = y, off =  off, pred_vars = pt_vars,
                 save_dir = file.path(in_dat_pth, "results"))

# Make predictions #=============================

# Build boot_out from stored info
boot_out <- map(SPP[1:nSPP],
                ~list.files(file.path(in_dat_pth, "results"),
                            pattern = paste0(.x, "_brt_\\d.*_trees_samp.*rds"),
                            full.names = TRUE))

# Build cross_val_out from stored info
cross_val_out <- read_csv(file.path(in_dat_pth, "results/cross_val_results.csv"),
                           show_col_types = FALSE) %>%
  nest(rel_inf = c(var, rel.inf))


# get predictors actually used in models
vars_used <- cross_val_out %>% unnest(rel_inf) %>%
  filter(rel.inf > rel_inf_thold) %>%
  pull(var) %>% unique()

# make a stack of predictor variable rasters
fl <- list.files(file.path(in_dat_pth,"0_data/processed/prediction-rasters"),
                 full.names = TRUE)
names(fl) <- fl %>% path_file() %>%  path_ext_remove()

fl_use <- fl[which(names(fl) %in% vars_used)]

var_stack <- terra::rast(fl_use)

# store this stack for future use with models
terra::writeRaster(var_stack,
                    file.path(in_dat_pth, "results/rast_stack_all_imp_vars.grd"),
                    overwrite = TRUE)

# Crop to study area
var_stack <- terra::crop(var_stack,
                         study_area %>% st_transform(st_crs(var_stack)))

# make predictions for each bootstrap model version and then take the mean
# Note that predict.gbm says it does not include the offset could we make a
# raster of the offset? Would depend on time of day so probably not really.
pred_out <- map(boot_out, ~map(., ~terra::predict(var_stack, readRDS(.x),
             filename = str_replace(.x, "\\.rds", ".tif"),
             type = "response", overwrite = TRUE) %>%
               setNames(fs::path_file(.x) %>% fs::path_ext_remove())))

pred_out_mean <- map(
  pred_out,
  ~terra::app(terra::rast(.x), fun = mean,
                na.rm = TRUE,
                filename = file.path(in_dat_pth, "results",
                                     str_replace(names(.x[[1]]), "_samp_.",
                                                 "_boot_mean.tif")),
                overwrite = TRUE)
  )

pred_out_sd <- map(
  pred_out,
  ~terra::app(terra::rast(.x), fun = sd,
                na.rm = TRUE,
                filename = file.path(in_dat_pth, "results",
                                     str_replace(names(.x[[1]]), "_samp_.",
                                                 "_boot_sd.tif")),
                overwrite = TRUE)
)

pred_out_mean %>% setNames(cross_val_out$spp) %>% terra::rast() %>% terra::plot()
pred_out_sd %>% setNames(cross_val_out$spp) %>% terra::rast() %>% terra::plot()
