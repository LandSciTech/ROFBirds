#The purpose of this script is to clip a roads raster, classifying 
#whether each pixel is on/adjacent to a roadside. The raster will
#be clipped to Ontario.
library(raster)
library(dismo)
library(rpart)
library(maptools)
library(data.table)
library(rgdal)
library(dplyr)
library(ggplot2)
library(spatial)
library(rasterVis)
library(RColorBrewer)
#First, we need a tif file of the Ontario study area
#Use Ontario study area grid with 250-m resolution we created for Beaudoin
#layers
Ontariostudyarea<-raster("0_data/raw/spatial data/Canada shapefile/Ontariostudyareamask.tif")

## Road on/off 
road<- raster("0_data/raw/spatial data/RoadsOnOff Layer Ontario/roadonoff1.tif")
mr <- c(1, 2500000, 1,  NA, NA, 0)
rcroad <- matrix(mr, ncol=3, byrow=TRUE)
rrc <- reclassify(road,rcroad)
#   

croprrc <- crop(rrc, extent(Ontariostudyarea), snap="out")
croprrc.r<-resample(croprrc, Ontariostudyarea)
maskrrc <- mask(x=croprrc.r, mask=Ontariostudyarea)
plot(maskrrc)
writeRaster(maskrrc, filename="0_data/processed/RoadsOnOff Layer Ontario/roadonoff.Ontariostudyarea.tif", overwrite=TRUE)

