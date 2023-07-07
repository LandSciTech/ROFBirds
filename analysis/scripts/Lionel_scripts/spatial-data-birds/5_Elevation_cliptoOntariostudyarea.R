#The purpose of this script is to clip elevation data from a
#continental 1-km resolution layer to just Ontario's extent
#and resample Ontario elevation to 250-m resolution. The original raster
#was generated using Climate NA software by Wang et al. (2016)
#(https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156720) 
#and the raster is stored at: 
#https://adaptwest.databasin.org/pages/adaptwest-climatena/.
#Note: Original data source is described as having a Lambert
#Aximuthal Equal Area projection, but projection of the downloaded
#raster appears to already have been in Lambert Conformal Conic.
library(raster)
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

elev <- raster("0_data/raw/spatial data/North America ASCII DEM Elevation/ClimateNA_DEM.asc")

cropelev <- crop(elev, extent(Ontariostudyarea), snap="out")
#crop raster to just Ontario to reduce memory requirements.
cropelev.r<-resample(cropelev, Ontariostudyarea)
#resamples elevation from 1-km resolution to 250-m resolution
maskelev <- mask(x=cropelev.r, mask=Ontariostudyarea)
writeRaster(maskelev, filename="0_data/processed/Elevation/elevOntariostudyarea.tif", overwrite=TRUE)
