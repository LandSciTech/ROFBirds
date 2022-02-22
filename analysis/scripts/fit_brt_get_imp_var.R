# Sarah's version of the bird modelling process
library(tidyverse)
library(fs)

in_dat_pth <- "analysis/data/derived_data"
# unzip files because Google Drive broke them down into many zip files for download
#
# file_ls <- list.files(in_dat_pth, full.names = TRUE)
#
# walk(file_ls, unzip, overwrite = FALSE)

# Steps to the analysis
# 1 Get the environmental and bird data in a format that can be modelled
# 2 Run GBM once to determine which variables to use
# 3 Run GBM nBootstrap times with important variables
# 4 Make predictions for each bootstrapped model and take the mean

# load data that BAM prepared and rename objects to be literate
load(file.path(in_dat_pth, "0_data/processed/BAMv6_RoFpackage_2022-01.RData"))

# Rename objects to be more literate

# BB is 2916, 250 matrix each col seems to be numbers b/w 1 and 11053. It is 250
# different samples of points with the same number in hudson plains and boreal
# shield. The first is all points in hudson plains and a sample w/o replacement
# from boreal shield other columns are samples with replacement
pt_resamps_250 <- BB

# cn2 222 vector of variable names all are cols in xx2 but it has other cols too
# It is cols in xx2 that have SD above 0.0001 (84 removed) and correlation > 0.9
# (54 removed) when a pair are correlated the one in the row is kept
# a subset is then hard coded in below (there was no explanation of why)
vars_sel <- cn2

# off 11053, 89 rownames PKEY point identifier also in xx1, colnames are sp names
# These are offsets, don't know where they come from.

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
              "turbidwater_G750.O", "volume2015.ntems")

run_brt_xv <- function(spp, RATE=0.001) {
  i <- 1
  si <- BB[,i]
  if (sum(y[si, spp]) < 1)
    return(structure(sprintf("0 detections for %s", spp), class="try-error"))
  xi <- data.frame(
    count=as.numeric(y[si, spp]),
    offset=off[si, spp],
    ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
    xx2[si, cn2])
  out <- try(gbm.step(xi,
                      gbm.y = 1,
                      gbm.x = 3:ncol(xi),
                      offset = xi$offset,
                      family = "poisson",
                      tree.complexity = 3,
                      learning.rate = RATE,
                      bag.fraction = 0.5))
  if (!inherits(out, "try-error"))
    out$rof_settings <- list(RATE=RATE, spp=spp, i=i)
  out
}

res <- run_brt_xv("ALFL")

fl <- list.files(file.path(in_dat_pth,"0_data/processed/prediction-rasters"),
                 full.names = TRUE)
names(fl) <- fl %>% path_file() %>%  path_ext_remove()

fl_use <- fl[which(names(fl) %in% c(vars_sel, "ecozone"))]

var_stack <- raster::stack(fl_use)

# Note that predict.gbm says it does not include the offset
res_pred <- predict(var_stack, res,
                    filename = file.path(in_dat_pth, "results/test_pred.tif"))
