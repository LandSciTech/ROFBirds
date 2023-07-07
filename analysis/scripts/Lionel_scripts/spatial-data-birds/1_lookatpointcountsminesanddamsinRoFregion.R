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
load("0_data/raw/spatial data/Ring of Fire Boundaries/gis_for_project.RData")
ls()
#[1] "aou"                       "caribou_ranges"            "far_north_boundary"       
#[4] "far_north_boundary_claims" "on"                        "rof_circle" 

final.rof.aoi <- st_read("0_data/raw/spatial data/Ring of Fire Boundaries/ROF_RA_def.shp")

#all spatial polygon data frames in which crs = 
#"crs         : +proj=longlat +datum=NAD83 +no_defs" 
#From GNM-National scripts on Github
lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
onsf <- sf::st_as_sf(on)
sf1 <- st_set_crs(onsf, "+proj=longlat +datum=NAD83 +no_defs")
sf1 <- st_transform(sf1, lcc_crs)

caribousf <- sf::st_as_sf(caribou_ranges)
sf2 <- st_set_crs(caribousf, "+proj=longlat +datum=NAD83 +no_defs")
sf2 <- st_transform(sf2, lcc_crs)

final.rof.sf <- sf::st_as_sf(final.rof.aoi)#Ontario MNR Lambert
final.rof.sf2 <- st_set_crs(final.rof.sf, "+proj=lcc +lat_1=44.5 +lat_2=53.5 +lat_0=0 +lon_0=-85 +x_0=930000 +y_0=6430000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
final.rof.sf2 <- st_transform(final.rof.sf2, lcc_crs)

ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf2, size = 0.5, color = "darkgreen", fill = "green") + 
  geom_sf(data = final.rof.sf2, size = 0.5, color = "darkred", fill = NA) + 
  ggtitle("Ontario Caribou Ranges & Final AOI") + 
  coord_sf()

aousf <- sf::st_as_sf(aou)
sf3 <- st_set_crs(aousf, "+proj=longlat +datum=NAD83 +no_defs")
sf3 <- st_transform(sf3, lcc_crs)

ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf3, size = 0.5, color = "darkred", fill = "red") + 
  ggtitle("Ontario Area of Undertaking for Ontario") + 
  coord_sf()

claimssf <- sf::st_as_sf(far_north_boundary_claims)
sf4 <- st_set_crs(claimssf, "+proj=longlat +datum=NAD83 +no_defs")
sf4 <- st_transform(sf4, lcc_crs)

ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf4, size = 0.5, color = "darkred", fill = "red") + 
  ggtitle("Ontario Far North Boundary Claims") + 
  coord_sf()

boundarysf <- sf::st_as_sf(far_north_boundary)
sf5 <- st_set_crs(boundarysf, "+proj=longlat +datum=NAD83 +no_defs")
sf5 <- st_transform(sf5, lcc_crs)

ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf5, size = 0.5, color = "darkred", fill = "red") + 
  ggtitle("Ontario Far North Boundary") + 
  coord_sf()

ringsf <- sf::st_as_sf(rof_circle)
sf6 <- st_set_crs(ringsf, "+proj=longlat +datum=NAD83 +no_defs")
sf6 <- st_transform(sf6, lcc_crs)

ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf6, size = 0.5, color = "darkred", fill = NA) + 
  ggtitle("Ontario Ring of Fire") + 
  coord_sf()

#Now load BBS-ON data and extract to each region
load("0_data/raw/point counts/BBS.91.19.ON.stopsFeb2021.RData")
ls()
#m5
names(m5)
nrow(m5)#7420205
bbs<-unique(m5[,c("SS","Year","lat","lon")])
nrow(bbs)#124332
bbs.sf<-sf::st_as_sf(bbs, coords=c("lon","lat"))
p1 <- st_set_crs(bbs.sf, "+proj=longlat +datum=NAD83 +no_defs")
p1 <- st_transform(p1, lcc_crs)

#filter BBS points to different shapefiles and plot them
bbs.farnorth = 
  p1 %>% 
  filter(st_contains(sf5, ., sparse = FALSE))
nrow(bbs.farnorth)#2546
P.bbs.farnorth<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf5, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = bbs.farnorth, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Ontario Far North - BBS") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(bbs.farnorth), " surveys"))
ggsave(P.bbs.farnorth, file="0_data/raw/point counts/BBSpoints_OntarioFarNorth.png", units="in", width=5, height=5)

#need to dissolve caribou range polygons before using them as filter
sf2$area<-st_area(sf2)
sf2D <-
  sf2 %>%
  summarise(area = sum(area))

bbs.caribou = 
  p1 %>% 
  filter(st_contains(sf2D, ., sparse = FALSE))
nrow(bbs.caribou)#2266
P.bbs.caribou<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf2D, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = bbs.caribou, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Ontario Caribou Ranges - BBS") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(bbs.caribou), " surveys"))
ggsave(P.bbs.caribou, file="0_data/raw/point counts/BBSpoints_OntarioCaribouRanges.png", units="in", width=5, height=5)

bbs.aou = 
  p1 %>% 
  filter(st_contains(sf3, ., sparse = FALSE))
nrow(bbs.aou)#67221
P.bbs.aou<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf3, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = bbs.aou, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Ontario AOU - BBS") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(bbs.aou), " surveys"))
ggsave(P.bbs.aou, file="0_data/raw/point counts/BBSpoints_OntarioAOU.png", units="in", width=5, height=5)


bbs.final.rof.aoi = 
  p1 %>% 
  filter(st_contains(final.rof.sf2, ., sparse = FALSE))
nrow(bbs.final.rof.aoi)#2262
P.bbs.final.rof.aoi<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = final.rof.sf2, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = bbs.final.rof.aoi, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Final Ring of Fire Area of Interest - BBS") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(bbs.final.rof.aoi), " surveys"))
ggsave(P.bbs.final.rof.aoi, file="0_data/raw/point counts/BBSpoints_FinalRingOfFireAOI.png", units="in", width=5, height=5)


#now look at Forest Bird Monitoring Program
F1<-read.csv("0_data/raw/point counts/FBMP.Before2000.csv",header=TRUE)
F2<-read.csv("0_data/raw/point counts/FBMP.In2000.DistTimeInt.csv",header=TRUE)
F3<-read.csv("0_data/raw/point counts/FBMP.In2000.UnlimDist.csv",header=TRUE)
F4<-read.csv("0_data/raw/point counts/FBMP.After2000.csv",header=TRUE)
FBMP<-bind_rows(F1,F2,F3,F4)
names(FBMP)
nrow(FBMP)#520410
FBMPU<-unique(FBMP[,c("SiteNumb_Station","Latitude","Longitude","DateTime")])
nrow(FBMPU)#27411
FBMPU<-FBMPU[!is.na(FBMPU$Latitude),]
nrow(FBMPU)#22415

FBMP.sf<-sf::st_as_sf(FBMPU, coords=c("Longitude","Latitude"))
p2 <- st_set_crs(FBMP.sf, "+proj=longlat +datum=NAD83 +no_defs")
p2 <- st_transform(p2, lcc_crs)


#filter FBMP points to different shapefiles and plot them
fbmp.farnorth = 
  p2 %>% 
  filter(st_contains(sf5, ., sparse = FALSE))
nrow(fbmp.farnorth)#0
P.fbmp.farnorth<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf5, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = fbmp.farnorth, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Ontario Far North - FBMP") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(fbmp.farnorth), " surveys"))
ggsave(P.fbmp.farnorth, file="0_data/raw/point counts/FBMPpoints_OntarioFarNorth.png", units="in", width=5, height=5)

fbmp.caribou = 
  p2 %>% 
  filter(st_contains(sf2D, ., sparse = FALSE))
nrow(fbmp.caribou)#0
P.fbmp.caribou<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf2D, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = fbmp.caribou, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Ontario Caribou Ranges - FBMP") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(fbmp.caribou), " surveys"))
ggsave(P.fbmp.caribou, file="0_data/raw/point counts/FBMPpoints_OntarioCaribouRanges.png", units="in", width=5, height=5)

fbmp.aou = 
  p2 %>% 
  filter(st_contains(sf3, ., sparse = FALSE))
nrow(fbmp.aou)#7819
P.fbmp.aou<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf3, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = fbmp.aou, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Ontario AOU - FBMP") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(fbmp.aou), " surveys"))
ggsave(P.fbmp.aou, file="0_data/raw/point counts/FBMPpoints_OntarioAOU.png", units="in", width=5, height=5)


#BAM Database V.5
#PKEY1<-read.csv("National_PKEY_V5_dec3_batch1_before2004.csv",header=TRUE)
#PKEY2<-read.csv("National_PKEY_V5_dec3_batch2_after2004.csv",header=TRUE)
#BAM1<-rbind(PKEY1,PKEY2)#read.csv("National_PKEY_V5.csv",header=TRUE)
#SS1<-read.csv("National_XY_V5_dec3_SSeastof96_exceptON.FBMP.MNBBA.MNNFB.csv",header=TRUE)
#SS2<-read.csv("National_XY_V5_dec3_SSeastof96_ONonly.csv",header=TRUE)
#BAM2<-rbind(SS1,SS2)#read.csv("National_XY_V5_ALL.csv",header=TRUE)
#BAM<-merge(BAM1, BAM2, by=c("SS"))
#nrow(BAM)
#245789 - smaller because it doesn't include BBS data
#includes Ontario Breeding Bird Atlas data, which is most or
#all of data in Far North and Caribou Ranges

#BAM data for Ontario (V.6 database)
BAM<-read.csv("0_data/raw/point counts/OntarioBAMpointcountXYandVisitsNoBBS.csv")
str(BAM)
BAM<-BAM[!is.na(BAM$X),]
nrow(BAM)#54813

BAM.sf<-sf::st_as_sf(BAM, coords=c("X","Y"))
p3 <- st_set_crs(BAM.sf, "+proj=longlat +datum=NAD83 +no_defs")
p3 <- st_transform(p3, lcc_crs)

#filter BAM points to different shapefiles and plot them
bam.farnorth = 
  p3 %>% 
  filter(st_contains(sf5, ., sparse = FALSE))
nrow(bam.farnorth)#3128
P.bam.farnorth<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf5, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = bam.farnorth, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Ontario Far North - BAM V.6") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(bam.farnorth), " surveys"))
ggsave(P.bam.farnorth, file="0_data/raw/point counts/BAMpoints_OntarioFarNorth.png", units="in", width=5, height=5)

bam.caribou = 
  p3 %>% 
  filter(st_contains(sf2D, ., sparse = FALSE))
nrow(bam.caribou)#2326
P.bam.caribou<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf2D, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = bam.caribou, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Ontario Caribou Ranges - BAM V.6") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(bam.caribou), " surveys"))
ggsave(P.bam.caribou, file="0_data/raw/point counts/BAMpoints_OntarioCaribouRanges.png", units="in", width=5, height=5)


bam.aou = 
  p3 %>% 
  filter(st_contains(sf3, ., sparse = FALSE))
nrow(bam.aou)#31707
P.bam.aou<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = sf3, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = bam.aou, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Ontario AOU - BAM V.6") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(bam.aou), " surveys"))
ggsave(P.bam.aou, file="0_data/raw/point counts/BAMpoints_OntarioAOU.png", units="in", width=5, height=5)

bam.final.rof.aoi = 
  p3 %>% 
  filter(st_contains(final.rof.sf2, ., sparse = FALSE))
nrow(bam.final.rof.aoi)#2532
P.bam.aou<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = final.rof.sf2, size = 0.25, color = "black", fill = "white") + 
  geom_sf(data = bam.final.rof.aoi, size = 1, color = "darkred", fill = "red") + 
  ggtitle("Final Ring of Fire Area of Interest - BAM V.6") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nrow(bam.final.rof.aoi), " surveys"))
ggsave(P.bam.aou, file="0_data/raw/point counts/BAMpoints_FinalRingOfFireAOI.png", units="in", width=5, height=5)


#Create New Far North-Caribou Ranges Study Area
#Union and Dissolve sf5 and sf2d
sf7<-st_union(sf5, sf2D)


bam.FarNorthCaribou = 
  p3 %>% 
  filter(st_contains(sf7, ., sparse = FALSE))
nrow(bam.FarNorthCaribou)#4709

bbs.FarNorthCaribou = 
  p1 %>% 
  filter(st_contains(sf7, ., sparse = FALSE))
nrow(bbs.FarNorthCaribou)#4317

P.FarNorthCaribou<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "darkgreen", fill = "green") + 
  geom_sf(data = sf7, size = 0.5, color = "black", fill = "white") + 
  geom_sf(data = sf2D, size = 0.5, color = "black", fill = "yellow") + 
  geom_sf(data = sf6, size = 0.75, color = "darkred", fill = NA) + 
  geom_sf(data = bam.FarNorthCaribou, size = 1, color = "red", fill = "red") + 
  geom_sf(data = bbs.FarNorthCaribou, size = 1, color = "blue", fill = "blue") + 
  ggtitle("Ontario Far North + Caribou Range") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nlevels(as.factor(bam.FarNorthCaribou$SS))+nlevels(as.factor(bbs.FarNorthCaribou$SS)), " sites"))+
  annotate("text", x = 1250000, y = 7400000, label = paste0(nlevels(as.factor(bam.caribou$SS))+nlevels(as.factor(bbs.caribou$SS)), " within caribou ranges"))+
  annotate("text", x = 1250000, y = 7300000, label = paste0((nlevels(as.factor(bam.FarNorthCaribou$SS))+nlevels(as.factor(bbs.FarNorthCaribou$SS)))-(nlevels(as.factor(bam.caribou$SS))+nlevels(as.factor(bbs.caribou$SS))), " outside caribou ranges"))+
  annotate("text", x = 1250000, y = 7200000, label = paste0(nlevels(as.factor(bam.FarNorthCaribou$SS)), " BAM (red)"))+
  annotate("text", x = 1250000, y = 7100000, label = paste0(nlevels(as.factor(bbs.FarNorthCaribou$SS)), " BBS (blue)"))
ggsave(P.FarNorthCaribou, file="0_data/raw/point counts/NumberOfSites_OntarioFarNorthCaribou.png", units="in", width=6, height=5)

P.FinalROF.AOI<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "darkgreen", fill = "green") + 
  geom_sf(data = final.rof.sf2, size = 0.5, color = "black", fill = "white") + 
  #geom_sf(data = sf2D, size = 0.5, color = "black", fill = "yellow") + 
  geom_sf(data = sf6, size = 0.75, color = "darkred", fill = NA) + 
  geom_sf(data = bam.final.rof.aoi, size = 1, color = "red", fill = "red") + 
  geom_sf(data = bbs.final.rof.aoi, size = 1, color = "blue", fill = "blue") + 
  ggtitle("Final Ring of Fire Area of Interest") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nlevels(as.factor(bam.final.rof.aoi$SS))+nlevels(as.factor(bbs.final.rof.aoi$SS)), " sites"))+
  annotate("text", x = 1250000, y = 7400000, label = paste0(nlevels(as.factor(bam.caribou$SS))+nlevels(as.factor(bbs.caribou$SS)), " within caribou ranges"))+
  annotate("text", x = 1250000, y = 7300000, label = paste0((nlevels(as.factor(bam.final.rof.aoi$SS))+nlevels(as.factor(bbs.final.rof.aoi$SS)))-(nlevels(as.factor(bam.caribou$SS))+nlevels(as.factor(bbs.caribou$SS))), " outside caribou ranges"))+
  annotate("text", x = 1250000, y = 7200000, label = paste0(nlevels(as.factor(bam.final.rof.aoi$SS)), " BAM (red)"))+
  annotate("text", x = 1250000, y = 7100000, label = paste0(nlevels(as.factor(bbs.final.rof.aoi$SS)), " BBS (blue)"))
ggsave(P.FinalROF.AOI, file="0_data/raw/point counts/NumberOfSites_FinalRingOfFireAOI.png", units="in", width=6, height=5)


#Now get drill hole locations
DR1<-read.csv("0_data/raw/spatial data/Mining Locations/Drill_holes1980andearlier.csv",header=TRUE)
DR2<-read.csv("0_data/raw/spatial data/Mining Locations/Drill_holes1981_1995.csv",header=TRUE)
DR3<-read.csv("0_data/raw/spatial data/Mining Locations/Drill_holes1996andlater.csv",header=TRUE)
drill.holes<-rbind(DR1,DR2,DR3)
names(drill.holes)
levels(as.factor(drill.holes$LL.Datum))
drill.holes<-drill.holes[!is.na(drill.holes$Lng.Dec),]
drill.holes<-drill.holes[!drill.holes$Lng.Dec=="",]
drill.holes<-drill.holes[!is.na(drill.holes$Lat.Dec),]
drill.holes<-drill.holes[!drill.holes$Lat.Dec=="",]
drill.holes.sf<-sf::st_as_sf(drill.holes, coords=c("Lng.Dec","Lat.Dec"))
dh1 <- st_set_crs(drill.holes.sf, "+proj=longlat +datum=NAD83 +no_defs")
dh1 <- st_transform(dh1, lcc_crs)

dh1.On = 
  dh1 %>% 
  filter(st_contains(sf1, ., sparse = FALSE))
P.DrillHoles<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = dh1.On, size = 0.5, color = "darkred", fill = NA) + 
  ggtitle("Ontario Drill Holes") + 
  coord_sf()
ggsave(P.DrillHoles, file="0_data/raw/point counts/P.DrillHoles.png", units="in", width=6, height=5)


#Now get abandoned mine locations
AM<-read.csv("0_data/raw/spatial data/Mining Locations/AMIS_2018_12/EXCEL_spreadsheets/ENDM_AMIS_SITES_NOV2018.csv",header=TRUE)
names(AM)
AM<-AM[!is.na(AM$SITES.LONGITUDE.DD),]
AM<-AM[!AM$SITES.LONGITUDE.DD=="",]
AM<-AM[!is.na(AM$SITES.LATITUDE.DD),]
AM<-AM[!AM$SITES.LATITUDE.DD=="",]
AM.sf<-sf::st_as_sf(AM, coords=c("SITES.LONGITUDE.DD","SITES.LATITUDE.DD"))
am1 <- st_set_crs(AM.sf, "+proj=longlat +datum=NAD83 +no_defs")
am1 <- st_transform(am1, lcc_crs)

am1.On = 
  am1 %>% 
  filter(st_contains(sf1, ., sparse = FALSE))
P.AbandonedMines<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  geom_sf(data = am1.On, size = 0.5, color = "darkred", fill = NA) + 
  ggtitle("Ontario Abandoned Mines") + 
  coord_sf()
ggsave(P.AbandonedMines, file="0_data/raw/point counts/P.AbandonedMines.png", units="in", width=6, height=5)

dh1.sf7 = 
  dh1 %>% 
  filter(st_contains(sf7, ., sparse = FALSE))

am1.sf7 = 
  am1 %>% 
  filter(st_contains(sf7, ., sparse = FALSE))

P.FarNorthCaribou.WithMines<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "darkgreen", fill = "green") + 
  geom_sf(data = sf7, size = 0.5, color = "black", fill = "white") + 
  geom_sf(data = sf2D, size = 0.5, color = "black", fill = "yellow") + 
  geom_sf(data = dh1.sf7, size = 0.5, color = "violet", fill = "violet") + 
  geom_sf(data = am1.sf7, size = 0.5, color = "violet", fill = "violet") + 
  geom_sf(data = sf6, size = 0.75, color = "darkred", fill = NA) + 
  geom_sf(data = bam.FarNorthCaribou, size = 1, color = "red", fill = "red") + 
  geom_sf(data = bbs.FarNorthCaribou, size = 1, color = "blue", fill = "blue") + 
  ggtitle("Ontario Far North + Caribou Range") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nlevels(as.factor(bam.FarNorthCaribou$SS))+nlevels(as.factor(bbs.FarNorthCaribou$SS)), " sites"))+
  annotate("text", x = 1250000, y = 7400000, label = paste0(nlevels(as.factor(bam.caribou$SS))+nlevels(as.factor(bbs.caribou$SS)), " within caribou ranges"))+
  annotate("text", x = 1250000, y = 7300000, label = paste0((nlevels(as.factor(bam.FarNorthCaribou$SS))+nlevels(as.factor(bbs.FarNorthCaribou$SS)))-(nlevels(as.factor(bam.caribou$SS))+nlevels(as.factor(bbs.caribou$SS))), " outside caribou ranges"))+
  annotate("text", x = 1250000, y = 7200000, label = paste0(nlevels(as.factor(bam.FarNorthCaribou$SS)), " BAM (red)"))+
  annotate("text", x = 1250000, y = 7100000, label = paste0(nlevels(as.factor(bbs.FarNorthCaribou$SS)), " BBS (blue)"))+
  annotate("text", x = 1350000, y = 7000000, label = paste0(nrow(dh1.sf7)+nrow(am1.sf7), " mines (violet)"))
ggsave(P.FarNorthCaribou.WithMines, file="0_data/raw/point counts/NumberOfSites_OntarioFarNorthCaribou_WithMines.png", units="in", width=6, height=5)

#Now Hydro
hydro <- st_read("0_data/raw/spatial data/Ontario_Hydro_Network__OHN__-_Hydrographic_Poly-shp/Ontario_Hydro_Network__OHN__-_Hydrographic_Poly.shp")
str(hydro)
levels(as.factor(hydro$HYDROGRAPH))
hydro$geometry#in lat-long WGS84
DAMS<-hydro[hydro$HYDROGRAPH=="Dam",]
DAMSsf <- st_set_crs(DAMS, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
damsf <- st_transform(DAMSsf, lcc_crs)
dam.buf <- st_buffer(damsf, dist = 30000)

ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "black", fill = "white") + 
  geom_sf(data = sf7, size = 0.5, color = "black", fill = "white")+ 
  geom_sf(data = dam.buf, color = "yellow", fill = "yellow")    

#Get number of point counts within 30 km of dams in study area
#Dams did not end up being included in the 2020/2021 Ring of Fire bird models
P.FarNorthCaribou.WithDams<-ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "darkgreen", fill = "green") + 
  geom_sf(data = sf7, size = 0.5, color = "black", fill = "white") + 
  geom_sf(data = sf2D, size = 0.5, color = "black", fill = "yellow") + 
  geom_sf(data = dam.buf, size = 0.5, color = "violet", fill = "violet") + 
  geom_sf(data = sf6, size = 0.75, color = "darkred", fill = NA) + 
  geom_sf(data = bam.FarNorthCaribou, size = 1, color = "red", fill = "red") + 
  geom_sf(data = bbs.FarNorthCaribou, size = 1, color = "blue", fill = "blue") + 
  ggtitle("Ontario Far North + Caribou Range") + 
  coord_sf()+
  annotate("text", x = 1250000, y = 7500000, label = paste0(nlevels(as.factor(bam.FarNorthCaribou$SS))+nlevels(as.factor(bbs.FarNorthCaribou$SS)), " point count sites in study area"))+
  annotate("text", x = 1250000, y = 7400000, label = paste0(nlevels(as.factor(bam.caribou$SS))+nlevels(as.factor(bbs.caribou$SS)), " within caribou ranges"))+
  annotate("text", x = 1250000, y = 7300000, label = paste0((nlevels(as.factor(bam.FarNorthCaribou$SS))+nlevels(as.factor(bbs.FarNorthCaribou$SS)))-(nlevels(as.factor(bam.caribou$SS))+nlevels(as.factor(bbs.caribou$SS))), " outside caribou ranges"))+
  annotate("text", x = 1250000, y = 7200000, label = paste0(nlevels(as.factor(bam.FarNorthCaribou$SS)), " BAM (red)"))+
  annotate("text", x = 1250000, y = 7100000, label = paste0(nlevels(as.factor(bbs.FarNorthCaribou$SS)), " BBS (blue)"))+
  annotate("text", x = 1350000, y = 7000000, label = paste0(nrow(damsf), " dam locations (violet)"))
ggsave(P.FarNorthCaribou.WithDams, file="0_data/raw/point counts/NumberOfSites_OntarioFarNorthCaribou_WithDams.png", units="in", width=6, height=5)

#Get number of point counts within 500, 2500, 5000 m of mines and dams in
#Far North/Caribou area
dh1.sf7.500 <- st_buffer(dh1.sf7, dist = 500)
dh1.sf7.500 %>%
  st_set_geometry(NULL) %>%
  glimpse()
dh1.sf7.500$area<-st_area(dh1.sf7.500)
dh1.sf7.500D<-dh1.sf7.500 %>%
  summarise(area=sum(area))

bbs.FarNorthCaribou.Drillholes.500 = 
  bbs.FarNorthCaribou %>% 
  filter(st_contains(dh1.sf7.500D, ., sparse = FALSE))
nlevels(as.factor(bbs.FarNorthCaribou.Drillholes.500$SS))#13 sites (75 surveys)
bam.FarNorthCaribou.Drillholes.500 = 
  bam.FarNorthCaribou %>% 
  filter(st_contains(dh1.sf7.500D, ., sparse = FALSE))
nlevels(as.factor(bam.FarNorthCaribou.Drillholes.500$SS))#102 sites (122 surveys)

dh1.sf7.2500 <- st_buffer(dh1.sf7, dist = 2500)
dh1.sf7.2500 %>%
  st_set_geometry(NULL) %>%
  glimpse()
dh1.sf7.2500$area<-st_area(dh1.sf7.2500)
dh1.sf7.2500D<-dh1.sf7.2500 %>%
  summarise(area=sum(area))

bbs.FarNorthCaribou.Drillholes.2500 = 
  bbs.FarNorthCaribou %>% 
  filter(st_contains(dh1.sf7.2500D, ., sparse = FALSE))
nlevels(as.factor(bbs.FarNorthCaribou.Drillholes.2500$SS))#128 sites (768 surveys)
bam.FarNorthCaribou.Drillholes.2500 = 
  bam.FarNorthCaribou %>% 
  filter(st_contains(dh1.sf7.2500D, ., sparse = FALSE))
nlevels(as.factor(bam.FarNorthCaribou.Drillholes.2500$SS))#461 sites (572 surveys)

dh1.sf7.5000 <- st_buffer(dh1.sf7, dist = 5000)
dh1.sf7.5000 %>%
  st_set_geometry(NULL) %>%
  glimpse()
dh1.sf7.5000$area<-st_area(dh1.sf7.5000)
dh1.sf7.5000D<-dh1.sf7.5000 %>%
  summarise(area=sum(area))

bbs.FarNorthCaribou.Drillholes.5000 = 
  bbs.FarNorthCaribou %>% 
  filter(st_contains(dh1.sf7.5000D, ., sparse = FALSE))
nlevels(as.factor(bbs.FarNorthCaribou.Drillholes.5000$SS))#235 sites (1491 surveys)
bam.FarNorthCaribou.Drillholes.5000 = 
  bam.FarNorthCaribou %>% 
  filter(st_contains(dh1.sf7.5000D, ., sparse = FALSE))
nlevels(as.factor(bam.FarNorthCaribou.Drillholes.5000$SS))#731 sites (905 surveys)


am1.sf7.500 <- st_buffer(am1.sf7, dist = 500)
am1.sf7.500 %>%
  st_set_geometry(NULL) %>%
  glimpse()
am1.sf7.500$area<-st_area(am1.sf7.500)
am1.sf7.500D<-am1.sf7.500 %>%
  summarise(area=sum(area))

bbs.FarNorthCaribou.AbandMines.500 = 
  bbs.FarNorthCaribou %>% 
  filter(st_contains(am1.sf7.500D, ., sparse = FALSE))
nlevels(as.factor(bbs.FarNorthCaribou.AbandMines.500$SS))#1 site (8 surveys)
bam.FarNorthCaribou.AbandMines.500 = 
  bam.FarNorthCaribou %>% 
  filter(st_contains(am1.sf7.500D, ., sparse = FALSE))
nlevels(as.factor(bam.FarNorthCaribou.AbandMines.500$SS))#16 sites (16 surveys)

am1.sf7.2500 <- st_buffer(am1.sf7, dist = 2500)
am1.sf7.2500 %>%
  st_set_geometry(NULL) %>%
  glimpse()
am1.sf7.2500$area<-st_area(am1.sf7.2500)
am1.sf7.2500D<-am1.sf7.2500 %>%
  summarise(area=sum(area))

bbs.FarNorthCaribou.AbandMines.2500 = 
  bbs.FarNorthCaribou %>% 
  filter(st_contains(am1.sf7.2500D, ., sparse = FALSE))
nlevels(as.factor(bbs.FarNorthCaribou.AbandMines.2500$SS))#17 sites (70 surveys)
bam.FarNorthCaribou.AbandMines.2500 = 
  bam.FarNorthCaribou %>% 
  filter(st_contains(am1.sf7.2500D, ., sparse = FALSE))
nlevels(as.factor(bam.FarNorthCaribou.AbandMines.2500$SS))#129 sites (151 surveys)

am1.sf7.5000 <- st_buffer(am1.sf7, dist = 5000)
am1.sf7.5000 %>%
  st_set_geometry(NULL) %>%
  glimpse()
am1.sf7.5000$area<-st_area(am1.sf7.5000)
am1.sf7.5000D<-am1.sf7.5000 %>%
  summarise(area=sum(area))

bbs.FarNorthCaribou.AbandMines.5000 = 
  bbs.FarNorthCaribou %>% 
  filter(st_contains(am1.sf7.5000D, ., sparse = FALSE))
nlevels(as.factor(bbs.FarNorthCaribou.AbandMines.5000$SS))#46 sites (177 surveys)
bam.FarNorthCaribou.AbandMines.5000 = 
  bam.FarNorthCaribou %>% 
  filter(st_contains(am1.sf7.5000D, ., sparse = FALSE))
nlevels(as.factor(bam.FarNorthCaribou.AbandMines.5000$SS))#254 sites (318 surveys)

damsf.sf7 = 
  damsf %>% 
  filter(st_contains(sf7, ., sparse = FALSE))

damsf.sf7.500 <- st_buffer(damsf.sf7, dist = 500)
damsf.sf7.500 %>%
  st_set_geometry(NULL) %>%
  glimpse()
damsf.sf7.500$area<-st_area(damsf.sf7.500)
damsf.sf7.500D<-damsf.sf7.500 %>%
  summarise(area=sum(area))

bbs.FarNorthCaribou.Dams.500 = 
  bbs.FarNorthCaribou %>% 
  filter(st_contains(damsf.sf7.500D, ., sparse = FALSE))
nlevels(as.factor(bbs.FarNorthCaribou.Dams.500$SS))#1 site (11 surveys)
bam.FarNorthCaribou.Dams.500 = 
  bam.FarNorthCaribou %>% 
  filter(st_contains(damsf.sf7.500D, ., sparse = FALSE))
nlevels(as.factor(bam.FarNorthCaribou.Dams.500$SS))#0

damsf.sf7.2500 <- st_buffer(damsf.sf7, dist = 2500)
damsf.sf7.2500 %>%
  st_set_geometry(NULL) %>%
  glimpse()
damsf.sf7.2500$area<-st_area(damsf.sf7.2500)
damsf.sf7.2500D<-damsf.sf7.2500 %>%
  summarise(area=sum(area))

bbs.FarNorthCaribou.Dams.2500 = 
  bbs.FarNorthCaribou %>% 
  filter(st_contains(damsf.sf7.2500D, ., sparse = FALSE))
nlevels(as.factor(bbs.FarNorthCaribou.Dams.2500$SS))#6 sites (64 surveys)
bam.FarNorthCaribou.Dams.2500 = 
  bam.FarNorthCaribou %>% 
  filter(st_contains(damsf.sf7.2500D, ., sparse = FALSE))
nlevels(as.factor(bam.FarNorthCaribou.Dams.2500$SS))#0

damsf.sf7.5000 <- st_buffer(damsf.sf7, dist = 5000)
damsf.sf7.5000 %>%
  st_set_geometry(NULL) %>%
  glimpse()
damsf.sf7.5000$area<-st_area(damsf.sf7.5000)
damsf.sf7.5000D<-damsf.sf7.5000 %>%
  summarise(area=sum(area))

bbs.FarNorthCaribou.Dams.5000 = 
  bbs.FarNorthCaribou %>% 
  filter(st_contains(damsf.sf7.5000D, ., sparse = FALSE))
nlevels(as.factor(bbs.FarNorthCaribou.Dams.5000$SS))#12 sites (127 surveys)
bam.FarNorthCaribou.Dams.5000 = 
  bam.FarNorthCaribou %>% 
  filter(st_contains(damsf.sf7.5000D, ., sparse = FALSE))
nlevels(as.factor(bam.FarNorthCaribou.Dams.5000$SS))#12 sites (12 surveys)
