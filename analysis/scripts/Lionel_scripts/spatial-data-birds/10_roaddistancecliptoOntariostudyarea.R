#The purpose of this script is to clip a road distances raster, in which
#pixel values are the distances (in m) to the nearest pixel classified as a 
#roadside. The raster will be clipped to Ontario. The initial raster was constructed
#in ArcGIS from the raster stored in
#"0_data/processed/RoadsOnOff Layer Ontario/roadonoff.Ontariostudyarea.tif". The 
#intermediate product is what is stored in "0_data/raw/spatial data/Linear Disturbances".
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


#First, we need a tif file of the Ontario study area grid with 250-m resolution we created for Beaudoin
#layers
Ontariostudyarea<-raster("0_data/raw/spatial data/Canada shapefile/Ontariostudyareamask.tif")

## Tiff file with Euclidean distance to nearest road
roaddist<- raster("0_data/raw/spatial data/Linear Disturbances/NearestRoad.tif")


croprrc <- crop(roaddist, extent(Ontariostudyarea), snap="out")
croprrc.r<-resample(croprrc, Ontariostudyarea)
maskrrc <- mask(x=croprrc.r, mask=Ontariostudyarea)
plot(maskrrc)
writeRaster(maskrrc, filename="0_data/processed/Linear Disturbances/NearestRoad.CroppedOntario.tif", overwrite=TRUE)

