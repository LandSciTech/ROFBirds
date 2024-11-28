#Beaudoin layers (Beaudoin et al. 2014): 
#These are national raster layers of forest stand age and structure
#(species cover, biomass, canopy cover, volume) developed from 250-m
#MODIS satellite imagery combined with Permanent Sample Plot forestry
#data. Estimates of forest age and structure are at 250-m resolution
#and are available for the years 2001 and 2011. The 2001 estimates 
#will be extracted to point counts from earlier years while the 2011
#estimates will be extracted to point counts from later years in a
#later R script in this R project.

#The purpose of this script is to clip the national layers to Ontario.
#The layers will be further clipped to just Bird Conservation Regions
#7 and 8 in a later script. While this is less efficient than directly
#clipping the national layers to just BCRs 7 and 8 directly, this two-
#step clipping process is how the R project was originally done.

#The clipped Ontario forest age and structure rasters are then used to
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

#The Beaudoin layers are used for predictor values in the Ring of Fire 
#bird models. Values from the Beaudoin layers have also been used to fill in gaps in the 
#vegetation data in SpaDES, which by default would be CASFRI or maybe other VRI data. 
memory.limit(size=56000)
#First run this command to increase the amount of memory that can be
#allotted to R processing. You will need a lot of memory to create, operate on
#and save raster stacks. The number can be increased to allot more memory
#but will be limited by what is available on the computer where this
#R project is being used. Obviously a computer with more memory is better.
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

#First, we need a tif file of the Ontario study area for Ring of Fire models
library(sf)
lcc_crs <-"+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0
+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
prov <- st_read("0_data/raw/spatial data/Canada shapefile/gpr_000a11a_e/gpr_000a11a_e.shp")
str(prov)
ONTarea<-prov[prov$PRENAME=="Ontario",]
ONTsf <- st_transform(ONTarea, lcc_crs)
ggplot() + 
  geom_sf(data = ONTsf, size = 0.5, color = "black", fill = "white") + 
  ggtitle("Ontario study area") + 
  coord_sf()

extent(ONTsf)
#e. Now use the raster package to convert to a raster and set its CRS
#r.0E = raster()#default
#extent(r.0E) = extent(bcsf)
#projection(r.0E) = CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0
#+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
#res(r.0E)<-250#set resolution to 250 m, same as Beaudoin
#rp <- rasterize(bcsf, r.0E, 'PRNAME')

b2011 <- list.files("0_data/raw/spatial data/Beaudoin/2011/",pattern="tif$")
bs2011 <- stack(raster(paste0("0_data/raw/spatial data/Beaudoin/2011/",b2011[1])))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(paste0("0_data/raw/spatial data/Beaudoin/2011/",b2011[i])))}
names(bs2011) <- gsub("NFI_MODIS250m_2011_kNN_","",names(bs2011))
#dropping the first part of raster layer name to make name shorter

crop2011 <- crop(bs2011, extent(ONTsf), snap="out")
#Once national layers are read into R, they are clipped with crop  
#and mask functions to a shapefile of Ontario. Clipping to a smaller
#area will reduce memory requirements for later operations.
rp<-rasterize(ONTsf, crop2011 )
writeRaster(rp, filename="0_data/raw/spatial data/Canada shapefile/Ontariostudyareamask.tif", overwrite=TRUE)

rp<-raster("0_data/raw/spatial data/Canada shapefile/Ontariostudyareamask.tif")
mask2011 <- mask(x=crop2011, mask=rp)
writeRaster(mask2011, filename="0_data/processed/Beaudoin/2011/1_Beaudoin 2011 GRD Ontario/Ontario2011_250m.grd", overwrite=TRUE, format="raster", bylayer=TRUE, suffix=names(mask2011))
#writes each file to a separate GRD file

#create TIFF files from GRD files
grd2011 <- list.files("0_data/processed/Beaudoin/2011/1_Beaudoin 2011 GRD Ontario/",pattern="grd$")
for (i in grd2011){
  ras<-stack(paste0("0_data/processed/Beaudoin/2011/1_Beaudoin 2011 GRD Ontario/",i))
  writeRaster(ras, filename=paste0("0_data/processed/Beaudoin/2011/2_Beaudoin 2011 TIFFs Ontario/",i,".tif"))
}

ontario2011.1<- stack("0_data/processed/Beaudoin/2011/2_Beaudoin 2011 TIFFs Ontario/Ontario2011_250m_LandCover_NonVeg_v1.grd")
plot(ontario2011.1)#[1] "LandCover_NonVeg_v1"
rm(crop2011, mask2011, bs2011)
gc()#It may be necessary to restart R to crop and mask the next set of Beaudoin layers for 2001
#due to memory requirements

rp<-raster("0_data/raw/spatial data/Canada shapefile/Ontariostudyareamask.tif")
b2001 <- list.files("0_data/raw/spatial data/Beaudoin/2001/",pattern="tif$")
bs2001 <- stack(raster(paste0("0_data/raw/spatial data/Beaudoin/2001/",b2001[1])))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(paste0("0_data/raw/spatial data/Beaudoin/2001/",b2001[i])))}
names(bs2001) <- gsub("NFI_MODIS250m_2001_kNN_","",names(bs2001))

crop2001 <- crop(bs2001, rp, snap="out")#rp already created; use instead of extent(bcsf)
mask2001 <- mask(x=crop2001, mask=rp)
writeRaster(mask2001, filename="0_data/processed/Beaudoin/2001/1_Beaudoin 2001 GRD Ontario/Ontario2001_250m.tif", overwrite=TRUE, format="raster", bylayer=TRUE, suffix=names(mask2001))
#writes each file to a separate file

#create TIFF files
grd2001 <- list.files("0_data/processed/Beaudoin/2001/1_Beaudoin 2001 GRD Ontario/",pattern="grd$")
for (i in grd2001){
  ras<-stack(paste0("0_data/processed/Beaudoin/2001/1_Beaudoin 2001 GRD Ontario/",i))
  writeRaster(ras, filename=paste0("0_data/processed/Beaudoin/2001/2_Beaudoin 2001 TIFFs Ontario/",i,".tif"))
}


#Landscape-scale rasters for Beaudoin 2001 layers
b2001 <- list.files("0_data/processed/Beaudoin/2001/2_Beaudoin 2001 TIFFs Ontario/",pattern=".tif$")
bs2001 <- stack(raster(paste0("0_data/processed/Beaudoin/2001/2_Beaudoin 2001 TIFFs Ontario/",b2001[1])))
for (i in 2:length(b2001)) { bs2001 <- addLayer(bs2001, raster(paste0("0_data/processed/Beaudoin/2001/2_Beaudoin 2001 TIFFs Ontario/",b2001[i])))}
names(bs2001) <- gsub("Ontario2001_250m_","B2001_local",names(bs2001))
writeRaster(bs2001, filename="0_data/processed/Beaudoin/2001/Beaudoin 2001 local and landscape GRD/Beaudoin2001.local.grd", format="raster",overwrite=TRUE)

# obtain weighted sums of neighourhood cells using Gaussian filter with sigma=250, and 750m for Beaudoin and CTI layers, save outputs as rasters
bs2001<-stack("0_data/processed/Beaudoin/2001/Beaudoin 2001 local and landscape GRD/Beaudoin2001.local.grd")

## sigma = 250m
fw250<-focalWeight(x=bs2001,d=250,type="Gauss")
bs2001_Gauss250<-stack(focal(bs2001[[1]],w=fw250,na.rm=TRUE))
names(bs2001_Gauss250)<-gsub("local","_G250",names(bs2001)[[1]])
for(i in 2:nlayers(bs2001)){
  bs2001_Gauss250<-addLayer(bs2001_Gauss250,focal(bs2001[[i]],w=fw250,na.rm=TRUE))
  names(bs2001_Gauss250)[i]<-gsub("local","_G250",names(bs2001)[[i]])
}
bs2001_Gauss250<-brick(bs2001_Gauss250)
writeRaster(bs2001_Gauss250, filename="0_data/processed/Beaudoin/2001/Beaudoin 2001 local and landscape GRD/Beaudoin2001.G250.grd", format="raster",overwrite=TRUE)


# ## sigma = 750m
fw750<-focalWeight(x=bs2001,d=750,type="Gauss")
bs2001_Gauss750<-stack(focal(bs2001[[1]],w=fw750,na.rm=TRUE))
names(bs2001_Gauss750)<-gsub("local","_G750",names(bs2001)[[1]])
for(i in 2:nlayers(bs2001)){
  bs2001_Gauss750<-addLayer(bs2001_Gauss750,focal(bs2001[[i]],w=fw750,na.rm=TRUE))
  names(bs2001_Gauss750)[i]<-gsub("local","_G750",names(bs2001)[[i]])
}
bs2001_Gauss750<-brick(bs2001_Gauss750)
writeRaster(bs2001_Gauss750, filename="0_data/processed/Beaudoin/2001/Beaudoin 2001 local and landscape GRD/Beaudoin2001.G750.grd", format="raster",overwrite=TRUE)

#Landscape-scale rasters for Beaudoin 2011 layers
b2011 <- list.files("0_data/processed/Beaudoin/2011/2_Beaudoin 2011 TIFFs Ontario/",pattern=".tif$")
bs2011 <- stack(raster(paste0("0_data/processed/Beaudoin/2011/2_Beaudoin 2011 TIFFs Ontario/",b2011[1])))
for (i in 2:length(b2011)) { bs2011 <- addLayer(bs2011, raster(paste0("0_data/processed/Beaudoin/2011/2_Beaudoin 2011 TIFFs Ontario/",b2011[i])))}
names(bs2011) <- gsub("Ontario2011_250m_","B2011_local",names(bs2011))
writeRaster(bs2011, filename="0_data/processed/Beaudoin/2011/Beaudoin 2011 local and landscape GRD/Beaudoin2011.local.grd", format="raster",overwrite=TRUE)

# obtain weighted sums of neighourhood cells using Gaussian filter with sigma=250, and 750m for Beaudoin and CTI layers, save outputs as rasters
bs2011<-stack("0_data/processed/Beaudoin/2011/Beaudoin 2011 local and landscape GRD/Beaudoin2011.local.grd")

## sigma = 250m
fw250<-focalWeight(x=bs2011,d=250,type="Gauss")
bs2011_Gauss250<-stack(focal(bs2011[[1]],w=fw250,na.rm=TRUE))
names(bs2011_Gauss250)<-gsub("local","_G250",names(bs2011)[[1]])
for(i in 2:nlayers(bs2011)){
  bs2011_Gauss250<-addLayer(bs2011_Gauss250,focal(bs2011[[i]],w=fw250,na.rm=TRUE))
  names(bs2011_Gauss250)[i]<-gsub("local","_G250",names(bs2011)[[i]])
}
bs2011_Gauss250<-brick(bs2011_Gauss250)
writeRaster(bs2011_Gauss250, filename="0_data/processed/Beaudoin/2011/Beaudoin 2011 local and landscape GRD/Beaudoin2011.G250.grd", format="raster",overwrite=TRUE)


# ## sigma = 750m
fw750<-focalWeight(x=bs2011,d=750,type="Gauss")
bs2011_Gauss750<-stack(focal(bs2011[[1]],w=fw750,na.rm=TRUE))
names(bs2011_Gauss750)<-gsub("local","_G750",names(bs2011)[[1]])
for(i in 2:nlayers(bs2011)){
  bs2011_Gauss750<-addLayer(bs2011_Gauss750,focal(bs2011[[i]],w=fw750,na.rm=TRUE))
  names(bs2011_Gauss750)[i]<-gsub("local","_G750",names(bs2011)[[i]])
}
bs2011_Gauss750<-brick(bs2011_Gauss750)
writeRaster(bs2011_Gauss750, filename="0_data/processed/Beaudoin/2011/Beaudoin 2011 local and landscape GRD/Beaudoin2011.G750.grd", format="raster",overwrite=TRUE)
