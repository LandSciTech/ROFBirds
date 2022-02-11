# R code for predicting density by land cover for birds in Far North Ontario

library(mefa4)
library(jsonlite)
library(sf)
library(rgeos)
library(raster)

#pl <- st_read("d:/bam/2021/rof/predictors-layers/ecoregions-clipped/BCR7and8Ontario.shp")

fl <- list.files("0_data/processed/prediction-rasters")
names(fl) <- gsub("\\.tif", "", fl)

#load("d:/bam/2021/rof/BAMv6_RoFpackage.RData")

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
compare_sets(names(fl), cn2)


## density by landcover

lc <- cn2[endsWith(cn2, ".O")]
lcc <- sapply(strsplit(lc, "_"), "[[", 1)

LC <- list()

for (k in seq_along(lc)[-c(1,2)]) {
  cat(lc[k], "\n")
  flush.console()
  
  r <- raster(paste0(
    "0_data/processed/prediction-rasters/",
    lc[k], ".tif"))
  LC[[lcc[k]]] <- values(r)
}
LC <- do.call(cbind, LC)
rs <- rowSums(LC, na.rm=TRUE)
summary(rs)
range(LC, na.rm=TRUE)
LC[is.na(LC)] <- 0
wm <- find_max(LC)
r <- raster(paste0("2_pipeline/store/brt-boot-pred-mosaic/ALFL.tif"))
wm$index[is.na(values(r))] <- NA

rlcDom <- r
values(rlcDom) <- wm$index
#writeRaster(rlcDom, file="d:/bam/2021/rof/dominant-landcover.tif")
levels(wm$index)

u0 <- mask(rlcDom, pl)
s <- !is.na(values(u0))
DD <- NULL
for (spp in SPP) {
  cat(spp, "\n")
  flush.console()
  
  ri <- raster(paste0("2_pipeline/store/brt-boot-pred-mosaic/", spp, ".tif"))
  q <- quantile(ri, 0.999)
  ri <- mask(ri, pl)
  values(ri)[!is.na(values(ri)) & values(ri) > q] <- q
  ag <- aggregate(data.frame(D=values(ri)[s]),
                  list(LCC=levels(wm$index)[values(u0)[s]]), summary)
  ag <- data.frame(Species=spp, ag)
  
  DD <- rbind(DD, ag)
}

write.csv(DD, row.names=FALSE, file="3_outputs/tables/SppDensityByOLCC_v2.csv")


# Make the denisty x landcover plots 
library(ggplot2)
spp.lc <- read.csv(file="3_outputs/tables/SppDensityByOLCC_v2.csv")

for(spp in SPP){
  dat <- spp.lc[which(spp.lc$Species == spp), ]
  dir <- "3_outputs/maps/brt-boot-pred-mosaic/"
  png(paste0(dir, "/", spp, "-dens-by-lcc.png"), width=1000, height=750)
  print(ggplot(data = dat, aes(x = as.factor(LCC), y = D.Mean)) + geom_bar(stat = "identity") + coord_flip() + 
    labs(title = spp, x = "Land cover\n", y = "\nDensity (males/ha)") + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid.minor = element_blank(), 
                       axis.text = element_text(size = rel(2)), axis.title = element_text(size = rel(3)), 
                       plot.title = element_text(size = 30), axis.line = element_line(colour = "black")))
  dev.off()
  
}










