#The purpose of this script is to take the spatial data clipped
#to Ontario and extract it to point counts for modeling. 

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


#Extracting from Access Database in R

#Ignore next section commented out: being kept for now

#I set up connection to Access database "BAM-V6-USE"
#The path below is specific to my laptop
#Path: E:\RoF-SpaDES-birds\0_data\raw\point counts\BAM-V6-Use
#Instructions for setting up: https://www.r-bloggers.com/2013/01/getting-access-data-into-r/
con <- odbcConnect("BAM-V6-USE")

#Look up any queries
queries<- sqlQuery(con, "BAM-V6-USE")

#Look up tables
tbls <- sqlTables(con)
tbls$TABLE_NAME
#I see a "location" table

usethisquery <- sqlFetch(con, "BAM-V6-USE")
str(usethisquery)#2393759 obs. of  23 variables
save(usethisquery, file="0_data/raw/point counts/BAM_v6.RData")
#The query that has been saved consists of all survey observations
#in the BAM database in long format. It does not include offsets,
#which would be estimated separately. This script is only for extracting
#spatial data to point counts so we need to summarize to stations.

load("0_data/raw/point counts/BAM_v6.RData")
ls()
str(usethisquery)
theme_set(theme_bw())

#Create a reduced file containing survey-year locations
PKEY.all<-unique(usethisquery[,c("dataset_code","dataset_name","location_name_V6","latitude","longitude","PKEY_V6","survey_year")])
str(PKEY.all)#247873 obs. of  7 variables
PKEY.all<-PKEY.all[!is.na(PKEY.all$latitude),]
nrow(PKEY.all)#247692 obs.
write.csv(PKEY.all, file="0_data/raw/point counts/AllPKEY.csv")

#indicates that X and Y coordinates for the
#national model data are currently projected as lat-long
#lat-long coordinates are used for for generating QPAD offsets
#(that will be for another script)
#projection of points needs to be originally defined to match 
#the coordinates in the data (so, lat-long) but can be changed
#to other projections later on.
lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
lat_long <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
sf <- sf::st_as_sf(PKEY.all, coords = c("longitude","latitude"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

basemap.laea<-st_read("0_data/raw/spatial data/NA_PoliticalDivisions/data/bound_p/boundary_p_v2.shp") ## Basemap for point plots
#this shapefile is in a Lambert Azimuthal Equal Area Projection
basemap.laea <- st_set_crs(basemap.laea, "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
basemap.new <- st_transform(basemap.laea, lat_long)
str(basemap.new)
levels(as.factor(basemap.new$NAME))
studyarea<-basemap.new[basemap.new$NAME=="Ontario",]

ggplot() +
  geom_sf(data = studyarea, size = 0.25, color = "blue", fill = "cyan1") +
  ggtitle("Study Area for Ring of Fire") +
  coord_sf()

#Now that the study area is in same the projection as the BAM points,
#I can use this shapefile to filter out just points within the study area
#but it needs to be a single polygon. Alternatively, I may want to extract
#data from the polygons in the study area to the points, then keep the
#points that get data (occur within the polygon)

## intersect polygons with points, keeping the information from both
pli = st_intersection(sf, studyarea)
nrow(pli)#81117 point count locations
str(pli)

ggplot() +
  geom_sf(data = studyarea, size = 0.25, color = "blue", fill = "cyan1") +
  geom_sf(data = pts.farnorth, size = 1, color = "darkred", fill = "red") +
  ggtitle("Study Area for Ring of Fire") +
  coord_sf()

## transform into a 'data.frame' by removing the geometry
#st_geometry(pli) = NULL
#head(pli)
st_write(pli, "0_data/raw/point counts/OntarioBAMXYandVisitsNoBBS.csv", layer_options = "GEOMETRY=AS_XY")
#coordinates saved will be lat/long


#BAM data for Ontario
bam<-read.csv("0_data/raw/point counts/OntarioBAMpointcountXYandVisitsNoBBS.csv")
str(bam)
bamXY<-unique(bam[,c("dataset_code","location_name_V6","X","Y","PKEY_V6","survey_year")])
#note that "longitude" and "latitude" original names of coords.
#were changed to "X" and "Y" when writing "pli" to a CSV
str(bamXY)#54813 rows

#Projects in Ontario area. Again, BBS data will have to be added.
levels(as.factor(bamXY$dataset_code))

#Breeding Bird Survey points for Ontario
#This is a file in long format with information about BBS survey-years at 
#individual stops. This file was created in another R project, which was
#used to add lat/long, summarize species abundance, generate QPAD offsets 
#for each species at each stop in each year visited to raw data downloaded
#from the BBS website. That R project will also be made available with
#the R project for extracting spatial data to Ring of Fire point counts.

#Offsets and abundances were generated separately for BBS data from BAM 
#data because BBS data are updated yearly on the web but not regularly
#incorporated into the BAM database.

load("0_data/raw/point counts/BBS.91.19.ON.stopsFeb2021.RData")
#contains a data frame "m5" from which I select specific columns
#and create the same variable names as in the BAM site locations file I made.
#That way, the BBS and BAM data can be added together in the same file 
#and data extraction to point counts is done to them together.
m5$dataset_code<-"BBSON"
m5$location_name_V6<-m5$SS
m5$X<-m5$lon
m5$Y<-m5$lat
m5$PKEY_V6<-paste0(m5$SS,":",m5$Year)
m5$survey_year<-m5$Year
bbs1<-unique(m5[,c("dataset_code","location_name_V6","X","Y","PKEY_V6","survey_year")])

str(bbs1)#tibble
bbs1<-data.frame(bbs1)
str(bbs1)#124332 obs. of  6 variables

Compiled<-rbind(bamXY, bbs1)
nrow(Compiled)#179145
write.csv(Compiled, file="0_data/raw/point counts/AllBAM_BBS_XY_OntarioVisits.csv")


#create spatial data frame from regular data frame SScombo to extract GIS data
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
latlong <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


basemap.laea<-st_read("0_data/raw/spatial data/NA_PoliticalDivisions/data/bound_p/boundary_p_v2.shp") ## Basemap for point plots
basemap.laea <- st_set_crs(basemap.laea, "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
basemap.new <- st_transform(basemap.laea, latlong)
str(basemap.new)
levels(as.factor(basemap.new$NAME))
studyarea<-basemap.new[basemap.new$NAME=="Ontario",]

ggplot() +
  geom_sf(data = studyarea, size = 0.25, color = "blue", fill = "cyan1") +
  ggtitle("Study Area for Ring of Fire") +
  coord_sf()

#Now that the study area is in same the projection as the BAM points,
#I can use this shapefile to filter out just points within the study area
#but it needs to be a single polygon. Alternatively, I may want to extract
#data from the polygons in the study area to the points, then keep the
#points that get data (occur within the polygon)

sf <- sf::st_as_sf(Compiled, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

## intersect polygons with points, keeping the information from both
pli = st_intersection(sf, studyarea)
nrow(pli)#175278 point count visits in Ontario
str(pli)

ggplot() +
  geom_sf(data = studyarea, size = 0.25, color = "blue", fill = "cyan1") +
  geom_sf(data = pli, size = 1, color = "darkred", fill = "red") +
  ggtitle("Study Area for Ring of Fire") +
  coord_sf()

## transform into a 'data.frame' by removing the geometry
#st_geometry(pli) = NULL
#head(pli)
st_write(pli, "0_data/raw/point counts/OntarioBAMBBSXYandVisits.csv", layer_options = "GEOMETRY=AS_XY")
#coordinates saved will be lat/long

OntarioPoints<-read.csv("0_data/raw/point counts/OntarioBAMBBSXYandVisits.csv", header=TRUE)
sf <- sf::st_as_sf(OntarioPoints, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
#GIS data extraction
#First, add Bird Conservation Region ID and
#filter to just points in BCR 7 and 8
BCRs<- st_read("0_data/raw/spatial data/BCR_Terrestrial_master/BCR_Terrestrial_master.shp")
BCRs <- st_set_crs(BCRs, "+proj=longlat +datum=NAD83 +no_defs")
BCRs.new <- st_transform(BCRs, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

BCR.7.8<-BCRs.new[BCRs.new$BCR==7|BCRs.new$BCR==8,]
pli = st_intersection(sf, BCR.7.8)
nrow(pli)#28430

st_write(pli, "0_data/raw/point counts/OntarioBAMBBSXYandVisits_BCR7and8.csv", layer_options = "GEOMETRY=AS_XY")

SScombo<-read.csv("0_data/raw/point counts/OntarioBAMBBSXYandVisits_BCR7and8.csv", header=TRUE)
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
latlong <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

SScombo2<-SpatialPointsDataFrame(coords=SScombo[,c("X","Y")],data=SScombo,proj4string = latlong)

##Elevation (1 km)
mydem <- raster('0_data/processed/Elevation/elevOntariostudyarea.tif')
SScombo2<-spTransform(SScombo2,proj4string(mydem))

SScombo <- cbind(SScombo,"elev"=raster::extract(mydem,SScombo2@coords))
nrow(SScombo[is.na(SScombo$elev),])#2
save(SScombo, SScombo2, file="0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb18.RData")

load("0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb18.RData")

## Tree cover (density, Global 30m Landsat Tree Canopy Version 4)
#https://landsat.gsfc.nasa.gov/article/global-30m-landsat-tree-canopy-version-4-released
treecover<-raster("0_data/processed/Global LANDSAT 30TM Tree Cover Layer/treecoverLCC.Ontariostudyarea.tif")
plot(treecover)#the colours look like they are backwards: areas I'd
#expect to have less forest cover have been coloured with higher values
SScombo2<-spTransform(SScombo2,proj4string(treecover))
SScombo <- cbind(SScombo,"treecover"=raster::extract(treecover,SScombo2@coords))
SScombo$treecover[which(SScombo$treecover>250)]<-0 # converting values >250 (which represent no forest data) into 0s

## Forest height (not pre-processed to Ontario prior to data extraction)
#link: https://landscape.jpl.nasa.gov/ 
height<-raster("0_data/raw/spatial data/ForestHeightLidar/Simard_Pinto_3DGlobalVeg_L3C.tif")
SScombo2<-spTransform(SScombo2,proj4string(height))
SScombo <- cbind(SScombo,"LIDARheight"=raster::extract(height,SScombo2@coords))

## Road on/off 
road<- raster("0_data/processed/RoadsOnOff Layer Ontario/roadonoff.Ontariostudyarea.tif")
mr <- c(1, 2500000, 1,  NA, NA, 0)
rcroad <- matrix(mr, ncol=3, byrow=TRUE)
rrc <- reclassify(road,rcroad)
SScombo2<-spTransform(SScombo2,proj4string(road))
SScombo <- cbind(SScombo,"road_yesno"=raster::extract(rrc,SScombo2@coords))  
#At this point stop and check the file for missing predictor values
nrow(SScombo)#28430
nrow(SScombo[is.na(SScombo$nalc),])#0
nrow(SScombo[is.na(SScombo$treecover),])#1
nrow(SScombo[is.na(SScombo$height),])#0
nrow(SScombo[is.na(SScombo$road),])#0

## Nearest road
roaddist<-raster("0_data/processed/Linear Disturbances/NearestRoad.CroppedOntario.tif")
SScombo2<-spTransform(SScombo2,proj4string(roaddist))
SScombo <- cbind(SScombo,"road_dist_m"=raster::extract(roaddist,SScombo2@coords))
###Road distance: I may need to subtract 250 (or 125: centre of cell) from all distances to get zero distances
#for roadside points. Need to know how to treat NA values too.
SScombo$road_dist_m.adj<-SScombo$road_dist_m-125
#consider setting BBS road distances to 0
nrow(SScombo[is.na(SScombo$road_dist_m),])#8333
###CAUTION!####


## Nearest mine
minedist<-raster("0_data/processed/Mining Locations/NearestMine.CroppedOntario.tif")
SScombo2<-spTransform(SScombo2,proj4string(minedist))
SScombo <- cbind(SScombo,"mine_dist_m"=raster::extract(minedist,SScombo2@coords))
nrow(SScombo[is.na(SScombo$mine_dist_m),])#60

## Topography
TPI <- raster("0_data/processed/Topography/tpiLCC.Ontariostudyarea.tif")
crs.TPI<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
SScombo2<-spTransform(SScombo2,CRS=crs.TPI)
SScombo <- cbind(SScombo,"TPI"=raster::extract(TPI,SScombo2@coords))

TRI <- raster("0_data/processed/Topography/triLCC.Ontariostudyarea.tif")
crs.TRI<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
SScombo2<-spTransform(SScombo2,CRS=crs.TRI)
SScombo <- cbind(SScombo,"TRI"=raster::extract(TRI,SScombo2@coords))

slope <- raster("0_data/processed/Topography/slopeLCC.Ontariostudyarea.tif")
crs.slope<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
SScombo2<-spTransform(SScombo2,CRS=crs.slope)
SScombo <- cbind(SScombo,"slope"=raster::extract(slope,SScombo2@coords))
save(SScombo, SScombo2, file="0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb18.RData")
#rm(height, clim2010, climate.tifs, climate2010, climtif, LCC, LCC2, TRI, TPI, treecover, slope)

#Far North Land Cover Layers
fnlc.local<-stack("0_data/processed/Far North Land Cover/FNLC local and landscape/FNLC.lcc.local.gri")
#fnlc.local's projection:
#+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(fnlc.local)){
  SScombo <- cbind(SScombo, raster::extract(fnlc.local[[i]],as.matrix(cbind(SScombo2$X,SScombo2$Y)))) #
  names(SScombo)[ncol(SScombo)] <- names(fnlc.local)[i]
}
nrow(SScombo[!is.na(SScombo$coniftreedlocal),])#20505
nrow(SScombo[!is.na(SScombo$mixedtreedlocal),])#20505
nrow(SScombo[!is.na(SScombo$decidswamplocal),])#20505
SScombo[!is.na(SScombo$coniftreedlocal),]

FNLC.lcc.G250<-stack("0_data/processed/Far North Land Cover/FNLC local and landscape/FNLC.lcc.G250.gri")
#fnlc.lcc.G250's projection:
#+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(FNLC.lcc.G250)){
  SScombo <- cbind(SScombo, raster::extract(FNLC.lcc.G250[[i]],as.matrix(cbind(SScombo2$X,SScombo2$Y)))) 
  names(SScombo)[ncol(SScombo)] <- names(FNLC.lcc.G250)[i]
}

nrow(SScombo[!is.na(SScombo$coniftreed_G250),])#20611
nrow(SScombo[!is.na(SScombo$coniferswamp_G250),])#20611
nrow(SScombo[!is.na(SScombo$decidswamp_G250),])#20611

FNLC.lcc.G750<-stack("0_data/processed/Far North Land Cover/FNLC local and landscape/FNLC.lcc.G750.gri")
#fnlc.lcc.G750's projection:
#+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(FNLC.lcc.G750)){
  SScombo <- cbind(SScombo, raster::extract(FNLC.lcc.G750[[i]],as.matrix(cbind(SScombo2$X,SScombo2$Y)))) # 
  names(SScombo)[ncol(SScombo)] <- names(FNLC.lcc.G750)[i]
}

nrow(SScombo[!is.na(SScombo$coniftreed_G750),])#21387
nrow(SScombo[!is.na(SScombo$coniferswamp_G750),])#21387
nrow(SScombo[!is.na(SScombo$decidswamp_G750),])#21387
#load("0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb18.RData")

#most point count locations in Ontario do not occur in the Far North
#so most point counts receive a value of NA


#Ontario Land Cover Layers - could be used for points missing FNLC data
olc.local<-stack("0_data/processed/Ontario Land Cover/OLC local and landscape/OLC.lcc.local.gri")
#olc.local's projection:
#+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +datum=NAD83 +units=m +no_defs
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS="+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(olc.local)){
  SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded, raster::extract(olc.local[[i]],as.matrix(cbind(SScombo.Beaudoinadded2$X,SScombo.Beaudoinadded2$Y)))) #
  names(SScombo.Beaudoinadded)[ncol(SScombo.Beaudoinadded)] <- names(olc.local)[i]
}

OLC.lcc.G250<-stack("0_data/processed/Ontario Land Cover/OLC local and landscape/OLC.lcc.G250.gri")
#olc.lcc.G250's projection:
#+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +datum=NAD83 +units=m +no_defs
SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(OLC.lcc.G250)){
  SScombo <- cbind(SScombo, raster::extract(OLC.lcc.G250[[i]],as.matrix(cbind(SScombo2$X,SScombo2$Y)))) 
  names(SScombo)[ncol(SScombo)] <- names(OLC.lcc.G250)[i]
}


OLC.lcc.G750<-stack("0_data/processed/Ontario Land Cover/OLC local and landscape/OLC.lcc.G750.gri")
#olc.lcc.G750's projection:
#+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +datum=NAD83 +units=m +no_defs
SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=0 +lon_0=-85 +lat_1=44.5 +lat_2=53.5 +x_0=930000 +y_0=6430000 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(OLC.lcc.G750)){
  SScombo <- cbind(SScombo, raster::extract(OLC.lcc.G750[[i]],as.matrix(cbind(SScombo2$X,SScombo2$Y)))) # 
  names(SScombo)[ncol(SScombo)] <- names(OLC.lcc.G750)[i]
}

save(SScombo, SScombo2, file="0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb18.RData")

#North American Land Cover Layers 2005
nalc.local<-stack("0_data/processed/North American Land Cover 2005 MODIS/NALC 2005 local and landscape GRD/NALC2005.local.gri")
#nalc.local's projection:
#+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs
SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

for(i in 1:nlayers(nalc.local)){
  SScombo <- cbind(SScombo, raster::extract(nalc.local[[i]],as.matrix(cbind(SScombo2$X,SScombo2$Y)))) #
  names(SScombo)[ncol(SScombo)] <- names(nalc.local)[i]
}
nrow(SScombo[!is.na(SScombo$cropland.nalc2005local),])#31309
nrow(SScombo[!is.na(SScombo$urban.nalc2005local),])#31309
nrow(SScombo[!is.na(SScombo$temp.mixed.nalc2005local),])#31309
SScombo[!is.na(SScombo$coniftreedlocal),]

nalc.2005.G250<-stack("0_data/processed/North American Land Cover 2005 MODIS/NALC 2005 local and landscape GRD/NALC2005.G250.gri")
#nalc.2005.G250's projection:
#+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs
SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

for(i in 1:nlayers(nalc.2005.G250)){
  SScombo <- cbind(SScombo, raster::extract(nalc.2005.G250[[i]],as.matrix(cbind(SScombo2$X,SScombo2$Y)))) 
  names(SScombo)[ncol(SScombo)] <- names(nalc.2005.G250)[i]
}

nrow(SScombo[!is.na(SScombo$temp.broadleaf.nalc2005_G250),])#31312

nalc.2005.G750<-stack("0_data/processed/North American Land Cover 2005 MODIS/NALC 2005 local and landscape GRD/NALC2005.G750.gri")
#nalc.2005.G750's projection:
#+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs
SScombo2<-spTransform(SScombo2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

for(i in 1:nlayers(nalc.2005.G750)){
  SScombo <- cbind(SScombo, raster::extract(nalc.2005.G750[[i]],as.matrix(cbind(SScombo2$X,SScombo2$Y)))) 
  names(SScombo)[ncol(SScombo)] <- names(nalc.2005.G750)[i]
}

nrow(SScombo[!is.na(SScombo$temp.broadleaf.nalc2005_G750),])#31308

#At this point, SScombo should be split into two data frames
SScomboA<-SScombo[SScombo$survey_year<2007,]#extract Beaudoin 2001 variables to these points
SScomboA2<-SpatialPointsDataFrame(coords=SScomboA[,c("X","Y")],data=SScomboA,proj4string = latlong)

SScomboB<-SScombo[SScombo$survey_year>2006,]#extract Beaudoin 2011 variables to these points
SScomboB2<-SpatialPointsDataFrame(coords=SScomboB[,c("X","Y")],data=SScomboB,proj4string = latlong)

#Beaudoin 2001 layers
b2001.local<-stack("0_data/processed/Beaudoin/2001/Beaudoin 2001 local and landscape GRD/Beaudoin2001.local.gri")
#b2001.local's projection:
#+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
SScomboA2<-spTransform(SScomboA2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(b2001.local)){
  SScomboA <- cbind(SScomboA, raster::extract(b2001.local[[i]],as.matrix(cbind(SScomboA2$X,SScomboA2$Y)))) #
  names(SScomboA)[ncol(SScomboA)] <- gsub("B2001_","",names(b2001.local)[i])
}
nrow(SScomboA[!is.na(SScomboA$localStructure_Volume_Total_v1.grd),])#12642


b2001.G250<-stack("0_data/processed/Beaudoin/2001/Beaudoin 2001 local and landscape GRD/Beaudoin2001.G250.gri")
#b2001.G250's projection:
#+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
SScomboA2<-spTransform(SScomboA2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(b2001.G250)){
  SScomboA <- cbind(SScomboA, raster::extract(b2001.G250[[i]],as.matrix(cbind(SScomboA2$X,SScomboA2$Y)))) 
  names(SScomboA)[ncol(SScomboA)] <- gsub("B2001__","",names(b2001.G250)[i])
}

nrow(SScomboA[!is.na(SScomboA$G250Structure_Volume_Total_v1.grd),])#12643

b2001.G750<-stack("0_data/processed/Beaudoin/2001/Beaudoin 2001 local and landscape GRD/Beaudoin2001.G750.gri")
#b2001.2005.G750's projection:
#+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
SScomboA2<-spTransform(SScomboA2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(b2001.G750)){
  SScomboA <- cbind(SScomboA, raster::extract(b2001.G750[[i]],as.matrix(cbind(SScomboA2$X,SScomboA2$Y)))) 
  names(SScomboA)[ncol(SScomboA)] <- gsub("B2001__","",names(b2001.G750)[i])
}

nrow(SScomboA[!is.na(SScomboA$G750Structure_Volume_Total_v1.grd),])#12643
save(SScombo, SScombo2, SScomboA, file="0_data/processed/RoF_SpaDESspatialpointcountdata_Feb18.RData")

#Beaudoin 2011 layers
b2011.local<-stack("0_data/processed/Beaudoin/2011/Beaudoin 2011 local and landscape GRD/Beaudoin2011.local.gri")
#b2011.local's projection:
#+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
SScomboB2<-spTransform(SScomboB2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(b2011.local)){
  SScomboB <- cbind(SScomboB, raster::extract(b2011.local[[i]],as.matrix(cbind(SScomboB2$X,SScomboB2$Y)))) #
  names(SScomboB)[ncol(SScomboB)] <- gsub("B2011_","",names(b2011.local)[i])
}
nrow(SScomboB[!is.na(SScomboB$localStructure_Volume_Total_v1.grd),])#15787


b2011.G250<-stack("0_data/processed/Beaudoin/2011/Beaudoin 2011 local and landscape GRD/Beaudoin2011.G250.gri")
#b2011.G250's projection:
#+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
SScomboB2<-spTransform(SScomboB2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(b2011.G250)){
  SScomboB <- cbind(SScomboB, raster::extract(b2011.G250[[i]],as.matrix(cbind(SScomboB2$X,SScomboB2$Y)))) 
  names(SScomboB)[ncol(SScomboB)] <- gsub("B2011__","",names(b2011.G250)[i])
}

nrow(SScomboB[!is.na(SScomboB$G250Structure_Volume_Total_v1.grd),])#15787

b2011.G750<-stack("0_data/processed/Beaudoin/2011/Beaudoin 2011 local and landscape GRD/Beaudoin2011.G750.gri")
#b2011.2005.G750's projection:
#+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
SScomboB2<-spTransform(SScomboB2,CRS="+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(i in 1:nlayers(b2011.G750)){
  SScomboB <- cbind(SScomboB, raster::extract(b2011.G750[[i]],as.matrix(cbind(SScomboB2$X,SScomboB2$Y)))) 
  names(SScomboB)[ncol(SScomboB)] <- gsub("B2011__","",names(b2011.G750)[i])
}

nrow(SScomboB[!is.na(SScomboB$G750Structure_Volume_Total_v1.grd),])#15787
save(SScombo, SScombo2, SScomboA, SScomboB, file="0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb18.RData")

SScombo.Beaudoinadded<-bind_rows(SScomboA,SScomboB)
#create SScombo.Beaudoinadded spatial points data frame
SScombo.Beaudoinadded2<-SpatialPointsDataFrame(coords=SScombo.Beaudoinadded[,c("X","Y")],data=SScombo.Beaudoinadded, proj4string = latlong)

#Recent Derived 30-m Landsat national layers clipped to Ontario
#and resampled to 250 m
#These might be used to update some of the variables from the Beaudoin 
#layer

firemask <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/firemaskOntario.tif")
crs.firemask<-"+proj=lcc +lat_0=49 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.firemask)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"firemask.ntems"=raster::extract(firemask,SScombo.Beaudoinadded2@coords))
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$firemask.ntems),])
#1212 points within an area that was burned before or after the point count
SScombo.Beaudoinadded$firemask.ntems[is.na(SScombo.Beaudoinadded$firemask.ntems)]<-0
nrow(SScombo.Beaudoinadded[is.na(SScombo.Beaudoinadded$firemask.ntems),])#0

fireyear <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/fireyearOntario.tif")
crs.fireyear<-"+proj=lcc +lat_0=49 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.fireyear)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"fireyear.ntems"=raster::extract(fireyear,SScombo.Beaudoinadded2@coords))
#fire year uses 1900 as reference year
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$fireyear.ntems),])
#1212 points within an area that was burned before or after the point count
SScombo.Beaudoinadded$fireyear.ntems[is.na(SScombo.Beaudoinadded$fireyear.ntems)]<-0#should correspond to year 1900
nrow(SScombo.Beaudoinadded[is.na(SScombo.Beaudoinadded$fireyear.ntems),])#0

fireseverity <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/fireseverityOntario.tif")
crs.fireseverity<-"+proj=lcc +lat_0=49 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.fireseverity)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"fireseverity.ntems"=raster::extract(fireseverity,SScombo.Beaudoinadded2@coords))
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$fireseverity.ntems),])
#1212 points within an area that was burned before or after the point count
SScombo.Beaudoinadded$fireseverity.ntems[is.na(SScombo.Beaudoinadded$fireseverity.ntems)]<-0
nrow(SScombo.Beaudoinadded[is.na(SScombo.Beaudoinadded$fireseverity.ntems),])#0

harvestmask <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/harvestmask250Ontario.tif")
crs.harvestmask<-"+proj=lcc +lat_0=49 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.harvestmask)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"harvestmask.ntems"=raster::extract(harvestmask,SScombo.Beaudoinadded2@coords))
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$harvestmask.ntems),])
#9159 points within an area that was cut before or after the point count
SScombo.Beaudoinadded$harvestmask.ntems[is.na(SScombo.Beaudoinadded$harvestmask.ntems)]<-0
nrow(SScombo.Beaudoinadded[is.na(SScombo.Beaudoinadded$harvestmask.ntems),])#0

harvestyear <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/harvestyear250Ontario.tif")
crs.harvestyear<-"+proj=lcc +lat_0=49 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.harvestyear)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"harvestyear.ntems"=raster::extract(harvestyear,SScombo.Beaudoinadded2@coords))
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$harvestyear.ntems),])
#8788 points within an area that was cut before or after the point count
#harvest year uses 1900 as reference year
SScombo.Beaudoinadded$harvestyear.ntems[is.na(SScombo.Beaudoinadded$harvestyear.ntems)]<-0
nrow(SScombo.Beaudoinadded[is.na(SScombo.Beaudoinadded$harvestyear.ntems),])#0
#maybe reclassify harvest year as > some threshold number

biomass2015 <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/forestbiomass2015Ontario.tif")
crs.biomass2015<-"+proj=lcc +lat_0=49 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.biomass2015)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"biomass2015.ntems"=raster::extract(biomass2015,SScombo.Beaudoinadded2@coords))
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$biomass2015.ntems),])
#27526 points with forest biomass; treat NA values as 0 forest biomass?
SScombo.Beaudoinadded$biomass2015.ntems[is.na(SScombo.Beaudoinadded$biomass2015.ntems)]<-0#because maybe there IS no forest biomass here
nrow(SScombo.Beaudoinadded[is.na(SScombo.Beaudoinadded$biomass2015.ntems),])#0

volume2015 <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/forestvolume2015Ontario.tif")
crs.volume2015<-"+proj=lcc +lat_0=49 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.volume2015)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"volume2015.ntems"=raster::extract(volume2015,SScombo.Beaudoinadded2@coords))
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$volume2015.ntems),])
#27464 points with forest volume; treat NA values as 0 forest volume?
SScombo.Beaudoinadded$volume2015.ntems[is.na(SScombo.Beaudoinadded$volume2015.ntems)]<-0
nrow(SScombo.Beaudoinadded[is.na(SScombo.Beaudoinadded$volume2015.ntems),])#0

height2015 <- raster("0_data/raw/spatial data/NTEMS Derived Landsat Layers for 2015 Plus Harvest Fire Wetlands/meanforestheight2015Ontario.tif")
crs.height2015<-"+proj=lcc +lat_0=49 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.height2015)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"height2015.ntems"=raster::extract(height2015,SScombo.Beaudoinadded2@coords))
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$height2015.ntems),])
#27464 points with forest height; treat NA values as 0 forest height?
SScombo.Beaudoinadded$height2015.ntems[is.na(SScombo.Beaudoinadded$height2015.ntems)]<-0
nrow(SScombo.Beaudoinadded[is.na(SScombo.Beaudoinadded$height2015.ntems),])#0

#Exceedance Layers
exceedance.05 <- raster("0_data/processed/Acidity/Exceedance05.Ontario.tif")
crs.exceedance.05<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.exceedance.05)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"exceedance.05"=raster::extract(exceedance.05,SScombo.Beaudoinadded2@coords))
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$exceedance.05),])
#28150 points with exceedance values; leave NA as is

exceedance.20 <- raster("0_data/processed/Acidity/Exceedance20.Ontario.tif")
crs.exceedance.20<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
SScombo.Beaudoinadded2<-spTransform(SScombo.Beaudoinadded2,CRS=crs.exceedance.20)
SScombo.Beaudoinadded <- cbind(SScombo.Beaudoinadded,"exceedance.20"=raster::extract(exceedance.20,SScombo.Beaudoinadded2@coords))
nrow(SScombo.Beaudoinadded[!is.na(SScombo.Beaudoinadded$exceedance.20),])
#28150 points with exceedance values; leave NA as is
save(SScombo, SScombo2, SScomboA, SScomboB, SScombo.Beaudoinadded, file="0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb18.RData")

#Use one or more of these shapefiles to get point count
#distances to/occurrence in or out of study area, to 
#weight point count sampling probability

load("0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb18.RData")
getdistances<-SScombo[,c("X","Y","PKEY_V4","SS_V4","survey_year")]
write.csv(getdistances, file="0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdatalocations.csv")

load("0_data/raw/spatial data/Ring of Fire Boundaries/gis_for_project.RData")
writeOGR(rof_circle, dsn="0_data/processed/Ring of Fire Boundaries", 
         layer="RingOfFireCircle", 
         driver="ESRI Shapefile")
writeOGR(caribou_ranges, dsn="0_data/processed/Ring of Fire Boundaries", 
         layer="CaribouRanges", 
         driver="ESRI Shapefile")
writeOGR(far_north_boundary, dsn="0_data/processed/Ring of Fire Boundaries", 
         layer="FarNorthBoundary", 
         driver="ESRI Shapefile")

#Use Near tool in ArcGIS to get distances from point counts to boundary of study area and stored 
#the result in "PointDistanceFromRoFStudyAreaBoundary.csv", which was
#copied to "0_data/processed/point counts with data".

load("0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb18.RData")
distancefromstudyarea<-read.csv("0_data/processed/point counts with data/PointDistanceFromRoFStudyAreaBoundary.csv")
str(distancefromstudyarea)
dists<-distancefromstudyarea[,c("PKEY_V4","NEAR_DIST")]
SScombo.final.distadded<-merge(SScombo.Beaudoinadded,dists,by=c("PKEY_V4"))
save(SScombo, 
     SScombo2, 
     SScomboA, 
     SScomboB, 
     SScombo.Beaudoinadded, 
     SScombo.final.distadded,
     file="0_data/processed/point counts with data/RoF_SpaDESspatialpointcountdata_Feb27.RData")



# #reset working directory to R project folder
# lf <- raster("0_data/raw/spatial data/Topography/landformLCC.Ontariostudyarea.tif")
# #use crs of lf
# crs.lf<-"+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
# SScombo2<-spTransform(SScombo2,CRS=crs.lf)
# SScombo <- cbind(SScombo,"landform"=raster::extract(lf,SScombo2@coords))
# ###Land Form being assigned floating values for some reason; is R averaging land forms?
# 
# ###CAUTION!####
# #Land form SEEMS to have been assigned below, though there are NA values.
# 
# SScombo$landform[which(SScombo$landform>0&SScombo$landform<1)]<-0
# SScombo$landform<-as.integer(SScombo$landform)
# lf_classes<-data.frame(value=0:9,label=factor(c("water","valley","hilltop.in.valley","headwaters", "ridges.and.peaks","plain","local.ridge.in.plain","local.valley.in.plain","gentle.slopes","steep.slopes")))
# lfclass<-lf_classes$label[match(SScombo$landform,lf_classes$value)]
# SScombo$landform<-lfclass


