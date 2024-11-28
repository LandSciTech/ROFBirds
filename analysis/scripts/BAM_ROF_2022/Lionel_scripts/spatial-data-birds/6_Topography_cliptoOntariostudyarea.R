#The purpose of this script is to clip topographic variables from
#a continental scale to Ontario and resample clipped layers from 1-km to 250-m resolution.
#Data were originally obtained from a website described in:
#Amatulli et al. 2018. A suite of global, cross-scale topographic 
#variables for environmental and biodiversity modeling. Scientific 
#data 5, no. 1: 1-15.

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

#Terrain Ruggedness Index (TRI) is the mean of the absolute 
#differences in elevation between a focal cell and its 8 surrounding cells.
#It quantifies the total elevation change across the 3×3 cells. 
#Flat areas have a value of zero whereas mountain areas with steep ridges 
#have positive values, which can be greater than 2000m in the Himalaya area. 
TRI <- raster("0_data/raw/spatial data/Topography/triLCC.tif")

cropTRI <- crop(TRI, extent(Ontariostudyarea), snap="out")
cropTRI.r<-resample(cropTRI, Ontariostudyarea)
maskTRI <- mask(x=cropTRI.r, mask=Ontariostudyarea)
writeRaster(maskTRI, filename="0_data/processed/Topography/triLCC.Ontariostudyarea.tif", overwrite=TRUE)

#Topographic Position Index (TPI) is the difference between 
#the elevation of a focal cell and the mean of its 8 surrounding 
#cells. Positive and negative values correspond to ridges and 
#valleys, respectively, while zero values correspond generally 
#to flat areas (with the exception of a special case where a 
#focal cell with a value 5 can have surrounding cells with 
#values of 4×1 and 4×9, resulting in a TPI of 0). 
TPI <- raster("0_data/raw/spatial data/Topography/tpiLCC.tif")

cropTPI <- crop(TPI, extent(Ontariostudyarea), snap="out")
cropTPI.r<-resample(cropTPI, Ontariostudyarea)
maskTPI <- mask(x=cropTPI.r, mask=Ontariostudyarea)
writeRaster(maskTPI, filename="0_data/processed/Topography/tpiLCC.Ontariostudyarea.tif", overwrite=TRUE)

#Slope and aspect are decomposed into 3-dimensional vector 
#components (in the x, y, and z directions) using standard 
#trigonometric operators, and by calculating the resultant 
#vector magnitude within a user-specified moving window size 
slope <- raster("0_data/raw/spatial data/Topography/slopeLCC.tif")

cropslope <- crop(slope, extent(Ontariostudyarea), snap="out")
cropslope.r<-resample(cropslope, Ontariostudyarea)
maskslope <- mask(x=cropslope.r, mask=Ontariostudyarea)
writeRaster(maskslope, filename="0_data/processed/Topography/slopeLCC.Ontariostudyarea.tif", overwrite=TRUE)

#Topographic Wetness Index (TWI) or Compound Topographic Index (CTI)
#https://www.researchgate.net/post/Which-is-the-correct-way-to-calculate-the-topographic-wetness-index-with-ArcGis-Dextop
#This predictor was not available from Amatulli et al. (2018)
#and was not used in Ring of Fire bird models. The Boreal Avian Modelling
#Project has generated national and regional CTI rasters however,
#in case there is interest in using these layers in the future.

#Roughness is expressed as the largest inter-cell difference 
#of a focal cell and its 8 surrounding cells. A raster is available
#for clipping to Ontario but was not used in the Ring of Fire bird models.

#Land form categories
lf <- raster("0_data/raw/spatial data/Topography/lf_lcc1.tif")
croplandform <- crop(lf, extent(Ontariostudyarea), snap="out")
croplandform.r<-resample(croplandform, Ontariostudyarea)
masklandform <- mask(x=croplandform.r, mask=Ontariostudyarea)
writeRaster(masklandform, filename="0_data/processed/Topography/landformLCC.Ontariostudyarea.tif", overwrite=TRUE)

