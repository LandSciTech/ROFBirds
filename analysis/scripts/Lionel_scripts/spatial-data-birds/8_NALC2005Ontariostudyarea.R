#The purpose of this script is to take MODIS satellite data for 
#land cover classes at 250-m resolution and clip it from continental
#scale to Ontario.

# These land cover data are from 2005.
# For more information on MODIS, see http://modis.gsfc.nasa.gov/
 
#The layers will be further clipped to just Bird Conservation Regions
#7 and 8 in a later script. While this is less efficient than directly
#clipping the national layers to just BCRs 7 and 8 directly, this two-
#step clipping process is how the R project was originally done.

#The clipped land cover rasters are then used to
#create 3 raster stacks. One stack ("local") containing the original values 
#of each variable in each 250-m cell. The other two stacks contain mean
#values of each variable at larger spatial scales or "landscape" scales,
#approximately within 250 m and 750 m of each cell. Essentially buffers of
#those distances are generated around each 250-m cell in a moving window
#analysis and the mean value of each variable within the cells in each 
#buffer is calculated. However, instead of using a "hard" buffer boundary,
#I used a "soft" buffer in which cells closer to each 250-m cell have 
#more weight or influence on the final estimated mean landscape-scale value
#of each variable, while more distant cells are downweighted. The buffer
#distances of 250 m and 750 m define where the down-weighting starts becoming
#significant. Because weight declines with distance as a Gaussian function,
#this "soft" buffer process is called Gaussian filtering. You will see
#examples of Gaussian filtering in other scripts in this R project.



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
Ontariostudyarea<-raster("0_data/raw/RingOfFireSpaDESSpatialData/Canada shapefile/Ontariostudyareamask.tif")

##MODIS-based landcover (250-m)
nalc2005<-raster("0_data/raw/spatial data/North American Land Cover 2005 MODIS/NA_LandCover_2005_LCC.img")
   
cropnalc2005 <- crop(nalc2005, extent(Ontariostudyarea), snap="out")
cropnalc2005.r<-resample(cropnalc2005, Ontariostudyarea)
masknalc2005 <- mask(x=cropnalc2005.r, mask=Ontariostudyarea)
plot(masknalc2005)
writeRaster(masknalc2005, filename="0_data/processed/North American Land Cover 2005 MODIS/NALC2005LCC.Ontariostudyarea.tif", overwrite=TRUE)
#create rasters for individual cover types

masknalc2005<-raster("0_data/processed/North American Land Cover 2005 MODIS/NALC2005LCC.Ontariostudyarea.tif")
#   The following table describes the display of land cover classification in the 
# .img file:
#   >Value  Class                                     RGB values
# >----------------------------------------------------------------------------
# >1      Temperate or sub-polar needleleaf forest  0     0.24  0
# >2      Sub-polar taiga needleleaf forest         0.58  0.61  0.44
# >3      Tropical or sub-tropical broadleaf 
# >         evergreen forest                        0     0.39  0
# >4      Tropical or sub-tropical broadleaf 
# >         deciduous forest                        0.12  0.67  0.02
# >5      Temperate or sub-polar broadleaf 
# >         deciduous forest                        0.08  0.55  0.24
# >6      Mixed forest                              0.36  0.46  0.17
# >7      Tropical or sub-tropical shrubland        0.7   0.62  0.18
# >8      Temperate or sub-polar shrubland          0.7   0.54  0.2
# >9      Tropical or sub-tropical grassland        0.91  0.86  0.37
# >10     Temperate or sub-polar grassland          0.88  0.81  0.54
# >11     Sub-polar or polar shrubland-lichen-moss  0.61  0.46  0.33
# >12     Sub-polar or polar grassland-lichen-moss  0.73  0.83  0.56
# >13     Sub-polar or polar barren-lichen-moss     0.25  0.54  0.45
# >14     Wetland                                   0.42  0.64  0.54
# >15     Cropland                                  0.9   0.68  0.4
# >16     Barren lands                              0.66  0.67  0.68
# >17     Urban                                     0.86  0.13  0.15
# >18     Water                                     0.3   0.44  0.64
# >19     Snow and Ice                              1     0.98  1
# >
temp.needleleaf<-masknalc2005
values(temp.needleleaf) = ifelse(values(temp.needleleaf)==1,1,0)
plot(temp.needleleaf)#, xlim=c(1000000,1200000),ylim=c(6000000,6200000))
writeRaster(temp.needleleaf, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/temp.needleleaf2005LCC.Ontarioclip.tif", overwrite=TRUE)

taiga<-masknalc2005
values(taiga) = ifelse(values(taiga)==2,1,0)
plot(taiga)
writeRaster(taiga, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/taiga2005LCC.Ontarioclip.tif", overwrite=TRUE)

temp.broadleaf<-masknalc2005
values(temp.broadleaf) = ifelse(values(temp.broadleaf)==5,1,0)
plot(temp.broadleaf)#
writeRaster(temp.broadleaf, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/temp.broadleaf2005LCC.Ontarioclip.tif", overwrite=TRUE)

temp.mixed<-masknalc2005
values(temp.mixed) = ifelse(values(temp.mixed)==6,1,0)
plot(temp.mixed)#
writeRaster(temp.mixed, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/temp.mixed2005LCC.Ontarioclip.tif", overwrite=TRUE)

temp.shrub<-masknalc2005
values(temp.shrub) = ifelse(values(temp.shrub)==8,1,0)
plot(temp.shrub)#
writeRaster(temp.shrub, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/temp.shrub2005LCC.Ontarioclip.tif", overwrite=TRUE)

temp.grass<-masknalc2005
values(temp.grass) = ifelse(values(temp.grass)==10,1,0)
plot(temp.grass)#
writeRaster(temp.grass, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/temp.grass2005LCC.Ontarioclip.tif", overwrite=TRUE)

polar.shrub<-masknalc2005
values(polar.shrub) = ifelse(values(polar.shrub)==11,1,0)
plot(polar.shrub)#
writeRaster(polar.shrub, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/polar.shrub2005LCC.Ontarioclip.tif", overwrite=TRUE)

polar.grass<-masknalc2005
values(polar.grass) = ifelse(values(polar.grass)==12,1,0)
plot(polar.grass)#
writeRaster(polar.grass, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/polar.grass2005LCC.Ontarioclip.tif", overwrite=TRUE)

polar.barren<-masknalc2005
values(polar.barren) = ifelse(values(polar.barren)==13,1,0)
plot(polar.barren)#
writeRaster(polar.barren, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/polar.barren2005LCC.Ontarioclip.tif", overwrite=TRUE)

wetland<-masknalc2005
values(wetland) = ifelse(values(wetland)==14,1,0)
plot(wetland)#
writeRaster(wetland, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/wetland2005LCC.Ontarioclip.tif", overwrite=TRUE)

cropland<-masknalc2005
values(cropland) = ifelse(values(cropland)==15,1,0)
plot(cropland)#
writeRaster(cropland, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/cropland2005LCC.Ontarioclip.tif", overwrite=TRUE)

barren<-masknalc2005
values(barren) = ifelse(values(barren)==16,1,0)
plot(barren)#
writeRaster(barren, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/barren2005LCC.Ontarioclip.tif", overwrite=TRUE)

urban<-masknalc2005
values(urban) = ifelse(values(urban)==17,1,0)
plot(urban)#, xlim=c(1000000,1200000),ylim=c(6000000,6200000))
writeRaster(urban, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/urban2005LCC.Ontarioclip.tif", overwrite=TRUE)

water<-masknalc2005
values(water) = ifelse(values(water)==18,1,0)
plot(water)#
writeRaster(water, filename="0_data/processed/North American Land Cover 2005 MODIS/indiv layers/water2005LCC.Ontarioclip.tif", overwrite=TRUE)

#Landscape-scale rasters for NALC 2005 layers
NALC2005 <- list.files("0_data/processed/North American Land Cover 2005 MODIS/indiv layers/",pattern="Ontarioclip.tif$")
na2005 <- stack(raster(paste0("0_data/processed/North American Land Cover 2005 MODIS/indiv layers/",NALC2005[1])))
for (i in 2:length(NALC2005)) { na2005 <- addLayer(na2005, raster(paste0("0_data/processed/North American Land Cover 2005 MODIS/indiv layers/",NALC2005[i])))}
names(na2005) <- gsub("2005LCC.Ontarioclip",".nalc2005local",names(na2005))
writeRaster(na2005, filename="0_data/processed/North American Land Cover 2005 MODIS/NALC 2005 local and landscape GRD/NALC2005.local.grd", format="raster",overwrite=TRUE)


# obtain weighted sums of neighourhood cells using Gaussian filter with sigma=250, and 750m for Beaudoin and CTI layers, save outputs as rasters
na2005<-stack("0_data/processed/North American Land Cover 2005 MODIS/NALC 2005 local and landscape GRD/NALC2005.local.grd")

## sigma = 250m
fw250<-focalWeight(x=na2005,d=250,type="Gauss")
na2005_Gauss250<-stack(focal(na2005[[1]],w=fw250,na.rm=TRUE))
names(na2005_Gauss250)<-gsub("local","_G250",names(na2005)[[1]])
for(i in 2:nlayers(na2005)){
  na2005_Gauss250<-addLayer(na2005_Gauss250,focal(na2005[[i]],w=fw250,na.rm=TRUE))
  names(na2005_Gauss250)[i]<-gsub("local","_G250",names(na2005)[[i]])
}
na2005_Gauss250<-brick(na2005_Gauss250)
writeRaster(na2005_Gauss250, filename="0_data/processed/North American Land Cover 2005 MODIS/NALC 2005 local and landscape GRD/NALC2005.G250.grd", format="raster",overwrite=TRUE)


# ## sigma = 750m
fw750<-focalWeight(x=na2005,d=750,type="Gauss")
na2005_Gauss750<-stack(focal(na2005[[1]],w=fw750,na.rm=TRUE))
names(na2005_Gauss750)<-gsub("local","_G750",names(na2005)[[1]])
for(i in 2:nlayers(na2005)){
  na2005_Gauss750<-addLayer(na2005_Gauss750,focal(na2005[[i]],w=fw750,na.rm=TRUE))
  names(na2005_Gauss750)[i]<-gsub("local","_G750",names(na2005)[[i]])
}
na2005_Gauss750<-brick(na2005_Gauss750)
writeRaster(na2005_Gauss750, filename="0_data/processed/North American Land Cover 2005 MODIS/NALC 2005 local and landscape GRD/NALC2005.G750.grd", format="raster",overwrite=TRUE)
