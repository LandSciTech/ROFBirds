
library(tidyverse)
library(sf)
library(raster)
library(terra)


aou <- read.csv("0_data/raw/IBP-AOS-LIST21.csv")

acad <- read.csv("0_data/raw/ACAD Regional 2021.03.21-filtered.csv")


(which(acad$Common.Name %in% aou$COMMONNAME))

acad[-which(acad$Common.Name %in% aou$COMMONNAME), ]


dir <- "2_pipeline/store/brt-boot-pred-mosaic"
sp.list <- (list.files(dir, pattern = "tif") %>% str_split_fixed(pattern = "\\.", n = 2))[, 1]
sp.list <- sp.list[-grep("uncertainty", sp.list)]


sp.list[-which(sp.list %in% aou$COMMONNAME)]

aou_rof <- aou[which(aou$SPEC %in% sp.list), ]

acad_rof <- acad[which(acad$Common.Name %in% aou_rof$COMMONNAME), ]

acad_rof$spp_code <- aou$SPEC[match(acad_rof$Common.Name, aou$COMMONNAME)]


acad_rof_scores <- data.frame(acad_rof[, c("Common.Name", "spp_code", "Region", "RCS.b")])

ecozone <- raster("0_data/processed/prediction-rasters/ecozone.tif")
plot(ecozone)
re <- ecozone


pl <- st_read("0_data/processed/shapefiles/BCR7and8Ontario.shp")
pl <- st_transform(pl, proj4string(re))
re <- mask(re, pl)


pl$BCR <- as.numeric(pl$BCR)
bcr = rasterize(pl, re, "BCR")
writeRaster(bcr, paste0(
  "0_data/processed/",
  "bcr.tif"), overwrite = TRUE)
bcr.val <- values(bcr)

load("0_data/processed/predictor-layers/chunks2/chunks.RData")

sa <- st_read("0_data/processed/boundary.shp")
sa <- st_transform(sa, proj4string(re))

bcr <- mask(bcr, sa)



dir <- "2_pipeline/store/brt-boot"

pred7 <- list()
pred8 <- list()
for(i in 1:length(sp.list)){
  cat(paste0("--------------", sp.list[i], "--------------"))
  r <- raster(paste0("2_pipeline/store/brt-boot-pred-mosaic/", sp.list[i], ".tif")) 
  q <- quantile(r, 0.999)
  values(r)[!is.na(values(r)) & values(r) > q] <- q
  r <- mask(r, sa)
  values(r)[is.na(values(bcr))] <- NA
  v <- values(r)
  
  pred7[[i]] <- data.frame(spp = v[which(values(bcr) == 7)])
  colnames(pred7[[i]]) <- sp.list[i]
  pred8[[i]] <- data.frame(spp = v[which(values(bcr) == 8)])
  colnames(pred8[[i]]) <- sp.list[i]
}


pred7 <- do.call(cbind, pred7)
pred7[is.na(pred7)] <- 0

pred8 <- do.call(cbind, pred8)
pred8[is.na(pred8)] <- 0

pred_bcr <- rbind(data.frame(bcr = 7, pred7), data.frame(bcr = 8, pred8))

saveRDS(pred_bcr, file = "2_pipeline/store/pred_allspec_bcr.rds")


# Get the ACAD scores
acad7 <- acad_rof_scores %>% filter(Region == "BCR07") %>% dplyr::select(spp_code, RCS.b) %>% arrange(spp_code)
acad7 <- acad7[match(colnames(pred7), acad7$spp_code), ]
acad7$RCS.b[which(is.na(acad7$RCS.b))] <- 0
acad7$spp_code[which(is.na(acad7$spp_code))] <- colnames(pred7)[-which(colnames(pred7) %in% acad7$spp_code)]

score7 <- as.matrix(pred7) %*% acad7$RCS.b

acad8 <- acad_rof_scores %>% filter(Region == "BCR08") %>% dplyr::select(spp_code, RCS.b) %>% arrange(spp_code)
acad8 <- acad8[match(colnames(pred8), acad8$spp_code), ]
acad8$RCS.b[which(is.na(acad8$RCS.b))] <- 0

score8 <- as.matrix(pred8) %*% acad8$RCS.b
score8.rel <- score8/max(rbind(score7, score8), na.rm = T)

score7.rel <- score7/max(rbind(score7, score8), na.rm = T)


rv <- bcr
values(rv)[which(!is.na(values(rv)))] <- 0
values(rv)[which(values(bcr) == 7)] <- score7.rel
values(rv)[which(values(bcr) == 8)] <- score8.rel

rv <- mask(rv, sa)
plot(rv)

bcr_sa <- st_intersection(pl, sa)


png("3_outputs/maps/acad_scores.png", width=1000, height=1000)
op <- par(mfrow=c(1,1), mar=c(1,1,2,6))
plot(rv, axes=FALSE, box=FALSE)#, col=hcl.colors(100, "Lajolla"))
plot(bcr_sa$geometry, add=TRUE, border=4)
plot(sa$geometry, add=TRUE, border=6, colour = "black")
par(op)
dev.off()

saveRDS(data.frame(boreal = acad8, hp = acad7$RCS.b), file = "0_data/processed/acad_scores.Rds")
