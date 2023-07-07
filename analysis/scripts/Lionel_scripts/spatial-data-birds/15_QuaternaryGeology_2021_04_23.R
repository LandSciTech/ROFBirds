#The purpose of this script is to classify point counts according 
#to if they are located on or close to an esker (=1) or not (=0), and
#if point counts are located in uplands (=1) or not (=0).

#The original data source is a Quaternary Geology geodatabase available
#from: 
#http://www.geologyontario.mndm.gov.on.ca/mndmaccess/mndm_dir.asp?type=pub&id=EDS014-REV

#I extracted eskers as a shapefile in ArcGIS ("eskers.shp") and used
#it to create a 250-m esker raster in ArcGIS ("eskers_LCC_raster")

memory.limit(size=56000)
library(data.table)
library(dismo)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(raster)
library(rpart)
library(maptools)
library(rgdal)
theme_set(theme_bw())
library(sf)
lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
latlong <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

load("0_data/processed/RoF_SpaDESspatialpointcountdata_Feb27.RData")
ls()
SScombo2<-SpatialPointsDataFrame(coords=SScombo.final.distadded[,c("X","Y")],data=SScombo.Beaudoinadded, proj4string = latlong)

eskerLCC<-raster("0_data/raw/spatial data/Quaternary Geology/esker_LCC_raster/esker_LCC_raster.tif")
#250-m resolution
eskerLCC[is.na(eskerLCC)]<-0
plot(eskerLCC)

SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
SScombo.final.distadded <- cbind(SScombo.final.distadded,"eskerLCC250"=raster::extract(eskerLCC, SScombo2@coords))
nrow(SScombo.final.distadded[SScombo.final.distadded$eskerLCC==0,])
nrow(SScombo.final.distadded[SScombo.final.distadded$eskerLCC==1,])
#rastereskerpoints<-SScombo.final.distadded[SScombo.final.distadded$eskerLCC==1,]
#BorealEskerStudypoints<-SScombo.final.distadded[substr(SScombo.final.distadded$SS_V4,1,9)=="NONBORESK",]
#write.csv(rastereskerpoints, file="rastereskerpoints.csv")
#write.csv(BorealEskerStudypoints, file="BorealEskerStudypoints.csv")

#So if we use this raster I created from esker polygons selected from a shapefile,
#there are 403 points located in pixels assigned to eskers, but none of
#these points are from the Northern Boreal Esker study.
#There are 88 other points from the Northern Boreal Esker study that I am
#assuming are eskers, but none of them were classified as within an esker
#pixel from the esker raster I made.

eskers.shp<-st_read("0_data/raw/spatial data/Quaternary Geology/eskers.shp")
eskers.sf <- st_transform(eskers.shp, latlong)
eskers.sf <- st_transform(eskers.sf, lcc_crs)

prov <- st_read("0_data/raw/spatial data/Canada shapefile/gpr_000a11a_e/gpr_000a11a_e.shp")
str(prov)
Ontario<-prov[prov$PRENAME=="Ontario",]
Ontariosf <- st_transform(Ontario, lcc_crs)

ggplot() + 
  geom_sf(data = Ontariosf, size = 0.5, color = "black", fill = "white") + 
  geom_sf(data = eskers.sf, size = 0.5, color = "blue", fill = "blue") + 
  ggtitle("Eskers") + 
  coord_sf()

#I assessed whether using a finer-scale esker raster would result in 
#more esker point counts, including the Northern Boreal Esker project.
#I speculated that my original 250-m raster was too coarse and that
#many eskers might have been small enough that they weren't identified
#as eskers in the 250-m pixels.

#I used the esker shapefile to create a 30-m raster, which I then
#used to create a proportion-esker-within-150 m raster at 30-m resolution.
#I then averaged and resampled this raster at coarser resolution, thinking
#that a larger number of pixels with point counts would be identified as
#being on eskers.

eskers.shp<-readOGR("0_data/raw/spatial data/Quaternary Geology/eskers.shp")

eskers.shp$esker<-ifelse(eskers.shp$FEATURE=="esker or area of eskers; direction of flow know or assumed",1,0)
eskers.shp<-spTransform(eskers.shp, CRS=lcc_crs)
#create sf object from data

#create a fine-scale esker raster
library(fasterize)
#use Ontario raster as a template
Ontariostudyareamask<-raster("0_data/raw/spatial data/Canada shapefile/Ontariostudyareamask.tif")
#250 m resolution
#change resolution to 30 m
Ontario30<-Ontariostudyareamask
res(Ontario30)<-30
plot(Ontario30)

#30-m esker raster
r.eskers30<-rasterize(
  eskers.shp,
  Ontario30,
  field = "esker",
  fun="max"#,
  #background = NA_real_,
  #by = NULL
)
plot(r.eskers30)
writeRaster(r.eskers30, filename="0_data/processed/Quaternary Geology/eskers30.tif", overwrite=TRUE)
#saves binary esker raster at 30-m resolution within study area

#150-m Gaussian filter 
## sigma = 150m
r.eskers30[is.na(r.eskers30[])] <- 0
fw150<-focalWeight(x=r.eskers30,d=150,type="Gauss")
r.eskers30.G150<-focal(r.eskers30,w=fw150,na.rm=TRUE)
writeRaster(r.eskers30.G150, filename="0_data/processed/Quaternary Geology/eskers30_G150.tif", overwrite=TRUE)
plot(r.eskers30.G150)

r.eskers30.G150<-raster("0_data/processed/Quaternary Geology/eskers30_G150.tif", overwrite=TRUE)
ras.aggregate <- aggregate(r.eskers30.G150, fact=5)
res(ras.aggregate)

#aggregate at 250 m resolution
ras<-ras.aggregate
#change resolution from 150x150 to 250x250
resampleFactor <- 250/150
inCols <- ncol(ras)
inRows <- nrow(ras)
resampledRaster <- raster(ncol=(inCols / resampleFactor), nrow=(inRows / resampleFactor))
extent(resampledRaster) <- extent(ras)
res(resampledRaster)<-250#set to exactly 250 m before resampling occurs
resampledRaster <- resample(ras,resampledRaster,datatype="INT1U",method='bilinear')
writeRaster(resampledRaster, filename="0_data/processed/Quaternary Geology/prop.esker_250m.tif", overwrite=TRUE)
#resampled raster projection reset to latlong

resampledRaster<-raster("0_data/processed/Quaternary Geology/prop.esker_250m.tif")
projection(resampledRaster)<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
resampledRaster[is.na(resampledRaster)]<-0

SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
SScombo.final.distadded <- cbind(SScombo.final.distadded,"prop.eskerLCC250"=raster::extract(resampledRaster, SScombo2@coords))
nrow(SScombo.final.distadded[SScombo.final.distadded$prop.eskerLCC250==0,])#28430
nrow(SScombo.final.distadded[SScombo.final.distadded$prop.eskerLCC250>0,])#0
points.with.some.esker<-SScombo.final.distadded[SScombo.final.distadded$prop.eskerLCC250>0,]
write.csv(points.with.some.esker, file="points.with.some.esker.csv")
#creating a 30-m raster, then getting proportion of eskers within "150 m",
#then aggregating to 150 m and resampling to 250 m resulted in 2173 point counts 
#classified as esker points.
#But none of the points from the Northern Boreal Esker were assigned as eskers.
#There are 1446 BBS points, 486 Ontario Breeding Bird Atlas points, 
#113 WAP points, and the rest are RLMBP.
SScombo.final.distadded$prob.esker.missed<-ifelse(substr(SScombo.final.distadded$SS_V4,1,9)=="NONBORESK",1,0)
SScombo.final.distadded$eskerpoint<-ifelse((SScombo.final.distadded$prob.esker.missed+
                                              SScombo.final.distadded$prop.eskerLCC250+
                                              SScombo.final.distadded$eskerLCC250)>0,1,0)

#a point will be classified as an esker point if any of the three variables in the above sum is greater than zero.  
nrow(SScombo.final.distadded[SScombo.final.distadded$eskerpoint==1,])
#2238 point count surveys classified as on eskers, including the 88 reclassified NONBORESK point counts
SScombo.eskerdata.added<-SScombo.final.distadded
save(SScombo.eskerdata.added, file="0_data/processed/RoF_SpaDESspatialpointcountdata_April17.RData")
#counts intersecting with a pixel containing at least a bit of esker. Why?

load("0_data/processed/RoF_SpaDESspatialpointcountdata_April17.RData")
ls()
SScombo2<-SpatialPointsDataFrame(coords=SScombo.eskerdata.added[,c("X","Y")],data=SScombo.eskerdata.added, proj4string = latlong)
SScombo2<-spTransform(SScombo2, CRS=lcc_crs)

prov <- st_read("0_data/raw/RingOfFireSpaDESSpatialData/Canada shapefile/gpr_000a11a_e/gpr_000a11a_e.shp")
str(prov)
Ontario<-prov[prov$PRENAME=="Ontario",]
Ontariosf <- st_transform(Ontario, lcc_crs)

#surficial geology shapefile
geology<-st_read("0_data/raw/spatial data/Quaternary Geology/EDS014-REV/GIS_DATA/Quaternary/geology_ll.shp")
str(geology)
geologysf <- st_transform(geology, latlong)
geologysf <- st_transform(geology, lcc_crs)
upland<-geologysf[geologysf$NUMBER==20|geologysf$NUMBER==22|geologysf$NUMBER==23|geologysf$NUMBER==25|geologysf$NUMBER==27,]
upland$upland<-1

ggplot() + 
  geom_sf(data = Ontariosf, size = 0.5, color = "black", fill = "white") + 
  geom_sf(data = upland, size = 0.5, color = "blue", fill = "blue") + 
  ggtitle("Upland Deposits") + 
  coord_sf()

#create an upland raster
library(fasterize)
#use Ontario study area raster as a template
Ontariostudyareamask<-raster("0_data/raw/RingOfFireSpaDESSpatialData/Canada shapefile/Ontariostudyareamask.tif")
#250 m resolution

#250-m upland raster
r.uplands<-rasterize(
  upland,
  Ontariostudyareamask,
  field = "upland",
  fun="max"#,
  #background = NA_real_,
  #by = NULL
)
plot(r.uplands)
writeRaster(r.uplands, filename="0_data/processed/Quaternary Geology/uplands.tif", overwrite=TRUE)
#saves binary uplands raster at 250-m resolution within study area

r.uplands[is.na(r.uplands)]<-0
SScombo.eskerdata.added <- cbind(SScombo.eskerdata.added,"upland"=raster::extract(r.uplands, SScombo2@coords))
nrow(SScombo.eskerdata.added[is.na(SScombo.eskerdata.added$upland),])#21111
SScombo.eskerdata.added[SScombo.eskerdata.added$upland==1,]#7325
save(SScombo.eskerdata.added, file="0_data/processed/RoF_SpaDESspatialpointcountdata_April23.RData")

#Erin wanted to know how many upland and esker point counts and non-esker point counts were in the study region
#and outside the Ring of Fire study region but in the broader study area
