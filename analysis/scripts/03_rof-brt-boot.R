# R code for fitting bootstrapped BRT SDM models for birds in Far North Ontario



library(mefa4)
library(gbm)
library(dismo)
library(ggplot2)
library(segmented)

load("0_data/processed/BAMv6_RoFpackage_2022-01.Rdata")
cn2 <- c("eskerpoint",
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


# Function to run bootstrapped BRT models
run_brt1 <- function(i, spp) {
  RATE=0.001
  #  u <- rel_inf(res)
  #  vars <- c("count", "offset", as.character(u$var[u$rel.inf > 0]))
  #  ntree <- res$n.trees
  ntree <- 10000
  si <- BB[,i]
  if (sum(y[si, spp]) < 1)
    return(structure(sprintf("0 detections for %s", spp), class="try-error"))
  xi <- data.frame(
    count=as.numeric(y[si, spp]),
    offset=off[si, spp],
    ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
    xx2[si, cn2])# [,vars]
  out <- try(gbm::gbm(xi$count ~ . + offset(xi$offset),
                      data=xi[,-(1:2)],
                      n.trees = ntree,
                      interaction.depth = 3,
                      shrinkage = RATE,
                      bag.fraction = 0.5,
                      distribution = "poisson",
                      var.monotone = NULL,
                      keep.data = FALSE,
                      verbose = FALSE,
                      n.cores = 1))
  if (!inherits(out, "try-error"))
    out$rof_settings <- list(RATE=RATE, spp=spp, i=i)
  out
}


# Test on 2 spp with 2 bootstrapped runs, using the 
dir <- "2_pipeline/store/brt-boot-1"
if(!dir.exists(dir)){dir.create(dir)}
system.time({
  for (spp in SPP[1:2]) {
    cat("\n\n------------------------------", spp, "------------------------------\n\n")
    d1 <- paste0("2_pipeline/store/brt-boot-1/", spp)
    if(!dir.exists(d1)){dir.create(d1)}
    for(i in 1:2){
      res <- run_brt1(i, spp)
      save(res, file=paste0("2_pipeline/store/brt-boot-1/", spp, "/", i, ".RData"))
    }
  }
})






















