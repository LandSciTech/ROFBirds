#Once bird abundance has been summarized for individual BBS stops,
#if the BBS data are analyzed with other point count data in a single
#analysis (really, the only reason to ever analyze individual
#BBS stops as replicates), detection offset variables should be used
#in models to account for differences in survey methods and conditions
#that could affect detection probability of species.

#The QPAD detection offset method is described in:
#Sólymos, Péter, et al. "Calibrating indices of avian 
#density from non‐standardized survey data: making the
#most of a messy situation." Methods in Ecology and 
#Evolution 4.11 (2013): 1047-1058.

#There is an R package (QPAD) that is used to calculate
#detection offsets. It is possible to use QPAD to calculate
#offsets based on a specific set of point counts; however,
#it is usually easier to use the model coefficients stored
#in the most recent version of the QPAD package, from the best
#detection model for each species using the point count data
#that the QPAD package was last run on. In the latter case,
#functions stored within the "0_data/recurring/" folder (developed by 
#P. Solymos) in this R project can be used to obtain detection
#offsets with a very small amount of code.

## offsets

if (!requireNamespace("QPAD")) {
  if (!requireNamespace("remotes"))
    install.packages("remotes")
  remotes::install_github("psolymos/QPAD")
}
if (!requireNamespace("sp"))
  install.packages("sp")
if (!requireNamespace("maptools"))
  install.packages("maptools")
if (!requireNamespace("raster"))
  install.packages("raster")
if (!requireNamespace("intrval"))
  install.packages("intrval")

library(QPAD)
library(maptools)
library(intrval)
library(raster)

load_BAM_QPAD(version = 3)
if (getBAMversion() != "3")
  stop("This script requires BAM version 3")

source("0_data/recurring-master/recurring-master/offset/functions.R")
ls()
#[1] "make_off" "make_x"
rlcc <- raster("0_data/recurring-master/recurring-master/offset/data/lcc.tif")
rtree <- raster("0_data/recurring-master/recurring-master/offset/data/tree.tif")
rtz <- raster("0_data/recurring-master/recurring-master/offset/data/utcoffset.tif")
rd1 <- raster("0_data/recurring-master/recurring-master/offset/data/seedgrow.tif")
crs <- proj4string(rtree)

#Example: one species, one location
spp <- "OVEN"
## https://en.wikipedia.org/wiki/ISO_8601
dt <- "2019-06-07" # ISO 8601 in YYYY-MM-DD (0-padded)
tm <- "05:20" # ISO 8601 in hh:mm (24 hr clock, 0-padded)
lon <- -113.4938 # longitude WGS84 (EPSG: 4326)
lat <- 53.5461 # latitude WGS84 (EPSG: 4326)
dur <- 10 # mins
dis <- 100 # meters


x <- make_x(dt, tm, lon, lat, dur, dis)
o <- make_off(spp, x)
x
o

prov<-c("ON")
#Remember, if you are filtering data to additional/other provinces (coordinates unavailable for
#U.S. states), "prov<-c("ON")" will need to be edited.

for (i in prov){
  X<-read.csv(paste0("0_data/processed/8_BBS.",i,".readygetoffsets.csv"))
  str(X)
  print(nrow(X))
  X<-X[!is.na(X$tm),]
  X<-X[!is.na(X$dt),]
  X<-X[!is.na(X$lon),]
  X<-X[!is.na(X$lat),]
  X<-X[!is.na(X$dur),]
  X<-X[!is.na(X$dis),]
  nchar(X$tm)
  nrow(X[is.na(X$lat),])
  X$tm<-ifelse(nchar(X$tm)==4,paste0("0",X$tm),X$tm)#e.g. "5:47" to "05:47"
  X<-X[!X$lat<40,]#remove points that are too far south
  X<-X[!X$tm=="00:00",]
  x <- make_x(dt=X$dt, tm=X$tm,
              lon=X$lon, lat=X$lat,
              dur=X$dur, dis=X$dis, key=X$SS)
  
  s1 <- intersect(colnames(X), getBAMspecieslist())
  
  o1 <- matrix(0, nrow(x), length(s1))
  rownames(o1) <- rownames(X)
  colnames(o1) <- s1
  
  for (spp in s1) {
    cat(spp, "\n")
    flush.console()
    o1[,spp] <- make_off(spp, x)$offset
  }
  str(o1)
  sum(is.na(o1))
  
  o1.pkeyadded<-cbind(PKEY_V6=X$PKEY_V6, location_name_V6=X$location_name_V6, o1)
  write.csv(o1.pkeyadded, file=paste0("0_data/processed/9_BBS_",i,"_offsetspervisit.csv"))
  save(o1.pkeyadded, file=paste0("0_data/processed/9_BBS_",i,"_offsetspervisit.RData"))
  
}

