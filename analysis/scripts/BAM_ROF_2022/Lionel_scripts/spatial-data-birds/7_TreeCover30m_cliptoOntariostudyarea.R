#The purpose of this script is to extract canopy cover from a 
#global layer to Ontario and resample it to 250-m resolution.
#Tree cover (density, Global 30m Landsat Tree Canopy Version 4)
#https://landsat.gsfc.nasa.gov/article/global-30m-landsat-tree-canopy-version-4-released
#The most recent canopy cover estimates are for 2019.

#Note: the original data was stored in an ArcGIS geodatabase.
#At some point between 2019 and now, the geodatabase was opened 
#in ArcGIS and reprojected to Lambert Conformal Conic. Both
#the original geodatabase and the intermediate product made in
#ArcGIS are provided here.
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
Ontariostudyarea<-raster("0_data/raw/Canada shapefile/Ontariostudyareamask.tif")

tree <- raster("0_data/raw/spatial data/Global LANDSAT 30TM Tree Cover Layer/treecoverlcc1.tif")

croptree <- crop(tree, extent(Ontariostudyarea), snap="out")
croptree.r<-resample(croptree, Ontariostudyarea)
masktree <- mask(x=croptree.r, mask=Ontariostudyarea)
plot(masktree)
#Note: tree cover values >250 should be reset to 0
values(masktree) = ifelse(values(masktree)>250,0,values(masktree))
plot(masktree)
writeRaster(masktree, filename="0_data/processed/Global LANDSAT 30TM Tree Cover Layer/treecoverLCC.Ontariostudyarea.tif", overwrite=TRUE)

