#The purpose of this script is to clip all Ontario raster
#layers used in data extraction to point counts to a
#smaller area encompassing the Ring of Fire region
#for prediction purposes. Outputs are stored in a folder
#called "0_data/processed/prediction rasters/" where the 
#contents can be easily assembled into a prediction raster
#stack or used individually as needed to obtain predicted 
#densities of birds from models. 

# Load packages and data
library(colorspace)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(gridExtra)
library(maptools)
library(raster)
library(reshape2)
library(rgdal)
library(RColorBrewer)
#library(rgeos)
library(RODBC)
library(sf)
library(sp)
library(tidyr)

lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
# lat_long <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
#GIS data extraction
#Import BCR 7 & 8 and make sure it is in LCC projection
BCR.7.8<- readOGR("0_data/raw/spatial data/BCR_Terrestrial_master/BCR7and8Ontario.shp")
BCR.7.8.LCC<-spTransform(BCR.7.8, CRS=lcc_crs)

##Elevation (1 km)
elev <- raster('0_data/processed/Elevation/elevOntariostudyarea.tif')
#already projected in lcc_crs
elev.crop<-crop(elev, BCR.7.8.LCC)
writeRaster(elev.crop, filename="0_data/processed/prediction rasters/elev.tif")


## Tree cover (density, Global 30m Landsat Tree Canopy Version 4)
#https://landsat.gsfc.nasa.gov/article/global-30m-landsat-tree-canopy-version-4-released
treecover<-raster("0_data/processed/Global LANDSAT 30TM Tree Cover Layer/treecoverLCC.Ontariostudyarea.tif")
plot(treecover)#the colours look like they are backwards: areas I'd
treecover.crop<-crop(treecover, BCR.7.8.LCC)
writeRaster(treecover.crop, filename="0_data/processed/prediction rasters/treecover.tif")

## Forest height (this layer was not processed prior to clipping to BCR 7/8)
#link: https://landscape.jpl.nasa.gov/
LIDARheight<-raster("0_data/raw/spatial data/ForestHeightLidar/Simard_Pinto_3DGlobalVeg_L3C.tif")
#not projected in lcc_crs
crs(LIDARheight)<-"+proj=longlat +datum=WGS84 +no_defs"
LIDARheight.LCC<-projectRaster(from=LIDARheight, to=elev, crs=lcc_crs) 
LIDARheight.LCC.crop<-crop(LIDARheight.LCC, BCR.7.8.LCC)
writeRaster(LIDARheight.LCC.crop, filename="0_data/processed/prediction rasters/LIDARheight.tif")

## Road on/off 
road<- raster("0_data/raw/spatial data/RoadsOnOff Layer Ontario/roadonoff.Ontariostudyarea.tif")
mr <- c(1, 2500000, 1,  NA, NA, 0)
rcroad <- matrix(mr, ncol=3, byrow=TRUE)
rrc <- reclassify(road,rcroad)#projected in lcc_crs
rrc.crop<-crop(rrc, BCR.7.8.LCC)
writeRaster(rrc.crop, filename="0_data/processed/prediction rasters/road_yesno.tif")


## Landform and topoedaphic 
#setwd("E:/CWS Wood Thrush Contract/CHID-GWmodel/")
#Terrain Ruggedness Index (TRI) is the mean of the absolute differences in elevation between a focal cell and its 8 surrounding cells. It quantifies the total elevation change across the 3×3 cells. Flat areas have a value of zero whereas mountain areas with steep ridges have positive values, which can be greater than 2000m in the Himalaya area. 
#Topographic Position Index (TPI) is the difference between the elevation of a focal cell and the mean of its 8 surrounding cells. Positive and negative values correspond to ridges and valleys, respectively, while zero values correspond generally to flat areas (with the exception of a special case where a focal cell with a value 5 can have surrounding cells with values of 4×1 and 4×9, resulting in a TPI of 0). 
#Roughness is expressed as the largest inter-cell difference of a focal cell and its 8 surrounding cells.
#Slope and aspect are decomposed into 3-dimensional vector components (in the x, y, and z directions) using standard trigonometric operators, and by calculating the resultant vector magnitude within a user-specified moving window size 

TPI <- raster("0_data/processed/Topography/tpiLCC.Ontariostudyarea.tif")
crs.TPI<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
TPI.crop<-crop(TPI, BCR.7.8.LCC)
writeRaster(TPI.crop, filename="0_data/processed/prediction rasters/TPI.tif")

TRI <- raster("0_data/processed/Topography/triLCC.Ontariostudyarea.tif")
crs.TRI<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
TRI.crop<-crop(TRI, BCR.7.8.LCC)
writeRaster(TRI.crop, filename="0_data/processed/prediction rasters/TRI.tif")

slope <- raster("0_data/processed/Topography/slopeLCC.Ontariostudyarea.tif")
crs.slope<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
slope.crop<-crop(slope, BCR.7.8.LCC)
writeRaster(slope.crop, filename="0_data/processed/prediction rasters/slope.tif")


#Ontario Land Cover Layers - could be used for points missing FNLC data
olchab <- list.files("0_data/processed/Ontario Land Cover/",pattern="OLC.tif$")
for (i in 1:length(olchab)) { 
  st.fnlchab <- raster(paste0("0_data/processed/Ontario Land Cover/", olchab[i]))
  #crs(st.fnlchab)<-"+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +datum=NAD83 +units=m +no_defs"
  olc.local.LCC<-projectRaster(from=st.fnlchab, to=elev, crs=lcc_crs) 
  olc.local.LCC.crop<-crop(olc.local.LCC, BCR.7.8.LCC)
  names(olc.local.LCC.crop) <- gsub("OLC","local.O",names(st.fnlchab))
  writeRaster(olc.local.LCC.crop, filename=paste0("0_data/processed/prediction rasters/",names(olc.local.LCC.crop),".tif"),overwrite=TRUE)
}

for (i in 1:length(olchab)) { 
  st.fnlchab <- raster(paste0("0_data/processed/Ontario Land Cover/", olchab[i]))
  #crs(st.fnlchab)<-"+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +datum=NAD83 +units=m +no_defs"
  olc.local.LCC<-projectRaster(from=st.fnlchab, to=elev, crs=lcc_crs) 
  ## sigma = 250m
  fw250<-focalWeight(x=olc.local.LCC,d=250,type="Gauss")
  olc.local.LCC_Gauss250<-focal(olc.local.LCC,w=fw250,na.rm=TRUE)
  olc.local.LCC_Gauss250.crop<-crop(olc.local.LCC_Gauss250, BCR.7.8.LCC)
  names(olc.local.LCC_Gauss250.crop) <- gsub("OLC","_G250.O",names(st.fnlchab))
  writeRaster(olc.local.LCC_Gauss250.crop, filename=paste0("0_data/processed/prediction rasters/",names(olc.local.LCC_Gauss250.crop),".tif"),overwrite=TRUE)
}

for (i in 1:length(olchab)) { 
  st.fnlchab <- raster(paste0("0_data/processed/Ontario Land Cover/", olchab[i]))
  #crs(st.fnlchab)<-"+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +datum=NAD83 +units=m +no_defs"
  olc.local.LCC<-projectRaster(from=st.fnlchab, to=elev, crs=lcc_crs) 
  ## sigma = 750m
  fw750<-focalWeight(x=olc.local.LCC,d=750,type="Gauss")
  olc.local.LCC_Gauss750<-focal(olc.local.LCC,w=fw750,na.rm=TRUE)
  olc.local.LCC_Gauss750.crop<-crop(olc.local.LCC_Gauss750, BCR.7.8.LCC)
  names(olc.local.LCC_Gauss750.crop) <- gsub("OLC","_G750.O",names(st.fnlchab))
  writeRaster(olc.local.LCC_Gauss750.crop, filename=paste0("0_data/processed/prediction rasters/",names(olc.local.LCC_Gauss750.crop),".tif"),overwrite=TRUE)
}


#Landscape-scale rasters for NALC 2005 layers
NALC2005 <- list.files("0_data/processed/North American Land Cover 2005 MODIS/",pattern="Ontarioclip.tif$")
for (i in 1:length(NALC2005)) { 
  na2005 <- raster(paste0("0_data/processed/North American Land Cover 2005 MODIS/",NALC2005[i]))
  names(na2005) <- gsub("2005LCC.Ontarioclip",".nalc2005local",names(na2005))
  na2005.crop<-crop(na2005, BCR.7.8.LCC)
  writeRaster(na2005.crop, filename=paste0("0_data/processed/prediction rasters/",names(na2005.crop),".tif"),overwrite=TRUE)
}

# obtain weighted sums of neighourhood cells using Gaussian filter with sigma=250, and 750m for Beaudoin and CTI layers, save outputs as rasters
NALC2005 <- list.files("0_data/processed/prediction rasters/",pattern="local.tif$")
for (i in 1:length(NALC2005)) {
  na2005 <- raster(paste0("0_data/processed/prediction rasters/",NALC2005[i]))
  ## sigma = 250m
  fw250<-focalWeight(x=na2005,d=250,type="Gauss")
  na2005_Gauss250<-focal(na2005,w=fw250,na.rm=TRUE)
  names(na2005_Gauss250)<-gsub("local","_G250",names(na2005))
  na2005_Gauss250.crop<-crop(na2005_Gauss250, BCR.7.8.LCC)
  writeRaster(na2005_Gauss250.crop, filename=paste0("0_data/processed/prediction rasters/",names(na2005_Gauss250.crop),".tif"),overwrite=TRUE)
}

for (i in 1:length(NALC2005)) {
  na2005 <- raster(paste0("0_data/processed/prediction rasters/",NALC2005[i]))
  ## sigma = 750m
  fw750<-focalWeight(x=na2005,d=750,type="Gauss")
  na2005_Gauss750<-focal(na2005,w=fw750,na.rm=TRUE)
  names(na2005_Gauss750)<-gsub("local","_G750",names(na2005))
  na2005_Gauss750.crop<-crop(na2005_Gauss750, BCR.7.8.LCC)
  writeRaster(na2005_Gauss750.crop, filename=paste0("0_data/processed/prediction rasters/",names(na2005_Gauss750.crop),".tif"),overwrite=TRUE)
}


#Beaudoin 2011 layers
b2011 <- list.files("0_data/processed/Beaudoin/2011/2_Beaudoin 2011 TIFFs Ontario/",pattern=".tif$")
for (i in 1:length(b2011)) { 
  bs2011 <- raster(paste0("0_data/processed/Beaudoin/2011/2_Beaudoin 2011 TIFFs Ontario/",b2011[i]))
  names(bs2011) <- gsub("Ontario2011_250m_","local",names(bs2011))
  bs2011.crop<-crop(bs2011, BCR.7.8.LCC)
  writeRaster(bs2011.crop, filename=paste0("0_data/processed/prediction rasters/",names(bs2011.crop),".tif"),overwrite=TRUE)
}

for (i in 1:length(b2011)) {
  bs2011 <- raster(paste0("0_data/processed/Beaudoin/2011/2_Beaudoin 2011 TIFFs Ontario/",b2011[i]))
  ## sigma = 250m
  fw250<-focalWeight(x=bs2011,d=250,type="Gauss")
  bs2011_Gauss250<-focal(bs2011,w=fw250,na.rm=TRUE)
  names(bs2011_Gauss250)<-gsub("Ontario2011_250m_","G250",names(bs2011))
  bs2011_Gauss250.crop<-crop(bs2011_Gauss250, BCR.7.8.LCC)
  writeRaster(bs2011_Gauss250.crop, filename=paste0("0_data/processed/prediction rasters/",names(bs2011_Gauss250.crop),".tif"),overwrite=TRUE)
}


for (i in 1:length(b2011)) {
  bs2011 <- raster(paste0("0_data/processed/Beaudoin/2011/2_Beaudoin 2011 TIFFs Ontario/",b2011[i]))
  ## sigma = 750m
  fw750<-focalWeight(x=bs2011,d=750,type="Gauss")
  bs2011_Gauss750<-focal(bs2011,w=fw750,na.rm=TRUE)
  names(bs2011_Gauss750)<-gsub("Ontario2011_250m_","G750",names(bs2011))
  bs2011_Gauss750.crop<-crop(bs2011_Gauss750, BCR.7.8.LCC)
  writeRaster(bs2011_Gauss750.crop, filename=paste0("0_data/processed/prediction rasters/",names(bs2011_Gauss750.crop),".tif"),overwrite=TRUE)
}

#Recent Derived 30-m Landsat national layers clipped to Ontario
#and resampled to 250 m
#These might be used to update some of the variables from the Beaudoin 
#layer

firemask <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/firemaskOntario.tif")
firemask.LCC<-projectRaster(from=firemask, to=elev, CRS=lcc_crs)
firemask.crop<-crop(firemask.LCC, BCR.7.8.LCC)
writeRaster(firemask.crop, filename="0_data/processed/prediction rasters/firemask.ntems.tif")

fireyear <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/fireyearOntario.tif")
fireyear.LCC<-projectRaster(from=fireyear, to=elev, CRS=lcc_crs)
fireyear.crop<-crop(fireyear.LCC, BCR.7.8.LCC)
writeRaster(fireyear.crop, filename="0_data/processed/prediction rasters/fireyear.ntems.tif")

fireseverity <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/fireseverityOntario.tif")
fireseverity.LCC<-projectRaster(from=fireseverity, to=elev, CRS=lcc_crs)
fireseverity.crop<-crop(fireseverity.LCC, BCR.7.8.LCC)
writeRaster(fireseverity.crop, filename="0_data/processed/prediction rasters/fireseverity.ntems.tif")

harvestmask <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/harvestmask250Ontario.tif")
harvestmask.LCC<-projectRaster(from=harvestmask, to=elev, CRS=lcc_crs)
harvestmask.crop<-crop(harvestmask.LCC, BCR.7.8.LCC)
writeRaster(harvestmask.crop, filename="0_data/processed/prediction rasters/harvestmask.ntems.tif")

harvestyear <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/harvestyear250Ontario.tif")
harvestyear.LCC<-projectRaster(from=harvestyear, to=elev, CRS=lcc_crs)
harvestyear.crop<-crop(harvestyear.LCC, BCR.7.8.LCC)
writeRaster(harvestyear.crop, filename="0_data/processed/prediction rasters/harvestyear.ntems.tif")

biomass2015 <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/forestbiomass2015Ontario.tif")
biomass2015.LCC<-projectRaster(from=biomass2015, to=elev, CRS=lcc_crs)
biomass2015.crop<-crop(biomass2015.LCC, BCR.7.8.LCC)
writeRaster(biomass2015.crop, filename="0_data/processed/prediction rasters/biomass2015.ntems.tif")

volume2015 <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/forestvolume2015Ontario.tif")
volume2015.LCC<-projectRaster(from=volume2015, to=elev, CRS=lcc_crs)
volume2015.crop<-crop(volume2015.LCC, BCR.7.8.LCC)
writeRaster(volume2015.crop, filename="0_data/processed/prediction rasters/volume2015.ntems.tif")

height2015 <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/meanforestheight2015Ontario.tif")
height2015.LCC<-projectRaster(from=height2015, to=elev, CRS=lcc_crs)
height2015.crop<-crop(height2015.LCC, BCR.7.8.LCC)
writeRaster(height2015.crop, filename="0_data/processed/prediction rasters/height2015.ntems.tif")

#Exceedance Layers
exceedance.05 <- raster("0_data/processed/Acidity/Exceedance05.Ontario.tif")
exceedance.05.crop<-crop(exceedance.05, BCR.7.8.LCC)
writeRaster(exceedance.05.crop, filename="0_data/processed/prediction rasters/exceedance.05.tif")

exceedance.20 <- raster("0_data/processed/Acidity/Exceedance20.Ontario.tif")
exceedance.20.crop<-crop(exceedance.20, BCR.7.8.LCC)
writeRaster(exceedance.20.crop, filename="0_data/processed/prediction rasters/exceedance.20.tif")

#Use one or more of these shapefiles to get point count
#distances to/occurrence in or out of study area, to 
#weight point count sampling probability

load("0_data/processed/RoF_SpaDESspatialpointcountdata_Feb18.RData")
getdistances<-SScombo[,c("X","Y","PKEY_V4","SS_V4","survey_year")]
write.csv(getdistances, file="0_data/processed/RoF_SpaDESspatialpointcountdatalocations.csv")
load("0_data/raw/gis_for_lionel.RData")
writeOGR(rof_circle, dsn="0_data/processed", 
         layer="RingOfFireCircle", 
         driver="ESRI Shapefile")
writeOGR(caribou_ranges, dsn="0_data/processed", 
         layer="CaribouRanges", 
         driver="ESRI Shapefile")
writeOGR(far_north_boundary, dsn="0_data/processed", 
         layer="FarNorthBoundary", 
         driver="ESRI Shapefile")

#Use Near tool in ArcGIS to get distances from point counts to boundary of study area
load("0_data/processed/RoF_SpaDESspatialpointcountdata_Feb18.RData")
distancefromstudyarea<-read.csv("0_data/raw/spatial data/PointDistanceFromRoFStudyAreaBoundary.csv")
str(distancefromstudyarea)
dists<-distancefromstudyarea[,c("PKEY_V4","NEAR_DIST")]
SScombo.final.distadded<-merge(SScombo.Beaudoinadded,dists,by=c("PKEY_V4"))
save(SScombo, 
     SScombo2, 
     SScomboA, 
     SScomboB, 
     SScombo.Beaudoinadded, 
     SScombo.final.distadded,
     file="0_data/processed/RoF_SpaDESspatialpointcountdata_Feb27.RData")

