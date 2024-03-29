# R code for predicting density by land cover for birds in Far North Ontario

library(mefa4)
library(jsonlite)
library(sf)
library(rgeos)
library(raster)

#pl <- st_read("d:/bam/2021/rof/predictors-layers/ecoregions-clipped/BCR7and8Ontario.shp")

fl <- list.files("0_data/processed/prediction-rasters")
names(fl) <- gsub("\\.tif", "", fl)

load("0_data/processed/BAMv6_RoFpackage_2022-01.Rdata")

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

SPP <- list.dirs("2_pipeline/store/brt-boot-pred", recursive = FALSE, full.names=FALSE)

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










