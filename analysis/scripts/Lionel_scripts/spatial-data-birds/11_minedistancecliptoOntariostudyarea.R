#The purpose of this script is to clip a mine distances raster, in which
#pixel values are the distances (in m) to the nearest mine. The raster will 
#be clipped to Ontario. The initial raster was constructed
#in ArcGIS from a geodatabase obtained from the Province of Ontario’s 
#Ministry of Energy, Northern Development and Mines (ENDM) and is stored in:
#"0_data/raw/spatial data/Mining Locations" along with the original geodatabase.
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
minedist<- raster("0_data/raw/spatial data/Mining Locations/NearestMine.tif")
lcc_crs <-"+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
minedist.p<-projectRaster(minedist, crs=lcc_crs)
# LCC projection wasn't quite the same as Ontario mask: now it will be  

croprrc <- crop(minedist.p, extent(Ontariostudyarea), snap="out")
croprrc.r<-resample(croprrc, Ontariostudyarea)
maskrrc <- mask(x=croprrc.r, mask=Ontariostudyarea)
plot(maskrrc)
writeRaster(maskrrc, filename="0_data/processed/Mining Locations/NearestMine.CroppedOntario.tif", overwrite=TRUE)

