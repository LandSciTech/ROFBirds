# R code for predicting bird density from bootstrapped BRT models, and making maps, for birds in Far North Ontario


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

r <- raster(paste0(
  "0_data/processed/prediction-rasters/",
  "elev.tif"))

re <- r
values(re)[!is.na(values(re))] <- 1
pl <- rgdal::readOGR("0_data/processed/shapefiles/BCR7and8Ontario.shp")
pl <- spTransform(pl, proj4string(re))
re <- mask(re, pl[1,])
values(re)[is.na(values(re))] <- 0
values(re)[is.na(values(r))] <- NA

dir <- "0_data/processed/predictor-layers"
if(!dir.exists(dir)){dir.create(dir)}

writeRaster(re, paste0(
  "0_data/processed/prediction-rasters/",
  "ecozone.tif"))

n <- length(values(r))
Chunks <- sort(sample(1:10, n, replace = TRUE))
table(Chunks)
dir1 <- "0_data/processed/predictor-layers/chunks2"
if(!dir.exists(dir1)){dir.create(dir1)}
save(Chunks, file="0_data/processed/predictor-layers/chunks2/chunks.RData")

ree <- raster("0_data/processed/prediction-rasters/esker_LCC_raster.tif")
ree <- projectRaster(ree, r, method="ngb")
table(values(ree), useNA="a")
values(ree)[!is.na(values(r)) & is.na(values(ree))] <- 0
writeRaster(ree, paste0(
  "0_data/processed/prediction-rasters/",
  "eskerpoint.tif"))

j <- 1
cn2x <- c("ecozone", cn2)
for (j in 1:10) {
  cat("Preparing chunk", j, "\n")
  flush.console()
  ss <- Chunks==j
  M <- matrix(0, sum(ss), length(cn2x))
  dimnames(M) <- list(which(ss), cn2x)
  M <- data.frame(M)
  
  for (k in cn2x) {
    cat(j, k, "\n")
    flush.console()
    
    r <- raster(paste0(
      "0_data/processed//prediction-rasters/",
      k, ".tif"))
    M[,k] <- values(r)[ss]
    
  }
  
  save(M, file=paste0("0_data/processed/predictor-layers/chunks2/variables-", j, ".RData"))
}


## Making chunked predictions for species. This loop creates a list of predictions for each chunk across all of the bootstrap
## replicates. For each chunk, an .RData file is saved that contains a list, where each list item is the predicted values for 
## that chunk at that bootstrap replicate. 


library(dismo)
library(gbm)
library(raster)

r <- raster(paste0(
  "0_data/processed/prediction-rasters/",
  "elev.tif"))
s <- !is.na(values(r))

load("0_data/processed/BAMv6_RoFpackage_2022-01.Rdata")

if (!dir.exists("2_pipeline/store/brt-boot-pred/"))
  dir.create("2_pipeline/store/brt-boot-pred/")
#spp <- "OVEN"
#j <- 1
CHUNKS <- 1:10
for (j in CHUNKS) {
  gc()
  load(paste0("0_data/processed/predictor-layers/chunks2/variables-", j, ".RData"))
  for (spp in SPP) {
    cat(j, spp, "\n")
    flush.console()
    dir <- paste0("2_pipeline/store/brt-boot-1/", spp)
    ls <- list.files(dir)
    p.list <- lapply(1:length(ls), function(x){
      load(file.path(dir, ls[x]))
      if (inherits(res, "gbm")) {
        if (!dir.exists(paste0("2_pipeline/store/brt-boot-pred/", spp)))
          dir.create(paste0("2_pipeline/store/brt-boot-pred/", spp))
        p <- predict.gbm(res, M, res$n.trees, type="response")
      }
      return(p)
    })
    save(p.list, file=paste0("2_pipeline/store/brt-boot-pred/", spp, "/", spp, "-chunk-", j, ".RData"))
  }
}


## Putting together the pieces into raster maps. This loops creates raster maps that are the mean predicted values across 
## the bootstrapped replicates for each spp
library(raster)
SPP <- list.dirs("2_pipeline/store/brt-boot-pred", full.names=FALSE)[-1]
load("0_data/processed/predictor-layers/chunks2/chunks.RData")
r <- raster(paste0(
  "0_data/processed/prediction-rasters/",
  "elev.tif"))
s <- !is.na(values(r))

for (spp in SPP) {
  cat(spp, "\n")
  flush.console()
  spppred <- list()
  for (j in 1:10) {
    f <- paste0("2_pipeline/store/brt-boot-pred/",
                spp, "/", spp, "-chunk-", j, ".RData")
    load(f)
    spppred[[j]] <- do.call(cbind, lapply(1:length(p.list), function(x) p.list[[x]]))
  }
  v <- rowMeans(do.call(rbind, spppred))
  ri <- r
  values(ri)[s] <- v[s]
  dir <- "2_pipeline/store/brt-boot-pred-mosaic"
  if(!dir.exists(dir)){dir.create(dir)}
  writeRaster(ri,
              paste0(dir, "/", spp, ".tif"),
              overwrite=TRUE)
}

## make png maps

library(sf)
library(rgeos)
library(raster)

# RoF boundary
r <- raster(paste0("2_pipeline/store/brt-boot-pred-mosaic/ALFL.tif"))
pl <- rgdal::readOGR(dsn = "0_data/processed", layer = "boundary")
pl <- spTransform(pl, proj4string(r))
#u <- mask(r, pl)

for (spp in SPP) {
  cat(spp, "\n")
  flush.console()
  
  ri <- raster(paste0("2_pipeline/store/brt-boot-pred-mosaic/", spp, ".tif"))
  q <- quantile(ri, 0.999)
  values(ri)[!is.na(values(ri)) & values(ri) > q] <- q
  u <- mask(ri, pl)
  N <- round(2 * sum(values(u), na.rm=TRUE) * 2.5^2 / 10^6, 3)
  
  dir <- "3_outputs/maps/brt-boot-pred-mosaic"
  if(!dir.exists(dir)){dir.create(dir)}
  png(paste0(dir, "/", spp, ".png"), width=500, height=500)
  op <- par(mfrow=c(1,1), mar=c(1,1,2,6))
  plot(ri, axes=FALSE, box=FALSE, col=hcl.colors(100, "Lajolla"),
       #        main=paste(spp, "Mean Density (males/ha)\nPopulation size =", N, "M inds."))
       main=paste(spp, "Mean Density (males/ha)"))
  plot(pl, add=TRUE, border=4)
  par(op)
  dev.off()
  
}
