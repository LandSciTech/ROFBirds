#The purpose of this script is to generate rasters of potential 
#acid rain impacts on forests in the Ring of Fire Region, in case
#those impacts have any effect on bird abundance in the Ring of Fire
#bird models. Data for acid rain impacts come from GIS layers of  
#critical loads of acid deposition for soils, for either 95% or 80% 
#protection of the most sensitive tree species. The more positive
#the values within these "exceedance" layers, the stronger the predicted
#impacts of acid rain on trees. There is a soil exceedance layer for
#5 % of trees being impacted, along with a soil exceedance layer for
#20 % of trees being impacted. Higher values within the 5% layer mean
#the top 5 % of tree species are more likely to be affected. Higher values
#within the 20 % layer mean that the top 20 % of tree species are more
#likely to be affected. Check out the "README" file in "0_data/raw/spatial
#data/Acidity" for more information on the source of these exceedance layers.
#There may be updates to these layers by the time the bird models are
#ready to be updated and there are links to other GIS layers that potentially
#measure acid rain impacts.
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
lcc_crs <-"+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0
+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

prov <- st_read("0_data/raw/spatial data/Canada shapefile/gpr_000a11a_e/gpr_000a11a_e.shp")
str(prov)
Ontario<-prov[prov$PRENAME=="Ontario",]
Ontariosf <- st_transform(Ontario, lcc_crs)
ggplot() + 
  geom_sf(data = Ontariosf, size = 0.5, color = "black", fill = "white") + 
  ggtitle("Ontario study area") + 
  coord_sf()

#need to make sure any provincial boundaries 
#are dissolved before using 
#study area as a shape-file filter for
#clipping other shape-files
Ontariosf$area<-st_area(Ontariosf)
Ontariosf.D <-
  Ontariosf %>%
  summarise(area = sum(area))

#Exceedance measurements for critical loads sufficient to kill/damage 5 % of an area's trees
exceed.05 <- st_read("0_data/raw/spatial data/Acidity/soilCL05_exc_2010/soilCL05_exc_2010.shp")

#all spatial polygon data frames in which crs = 
#"crs         : +proj=longlat +datum=NAD83 +no_defs" 
#From GNM-National scripts on Github
lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
exceed.05.sf <- sf::st_as_sf(exceed.05)
sf1 <- st_set_crs(exceed.05.sf, "+proj=longlat +datum=NAD83 +no_defs")
sf1 <- st_transform(sf1, lcc_crs)

ggplot() + 
  geom_sf(data = sf1, size = 0.5, color = "blue", fill = "cyan1") + 
  ggtitle("Exceedance Map - CL 5 %") + 
  coord_sf()
#This layer is mostly Canada wide except for parts of Prairie
#Provinces and southeast Ontario

#filter exceedance shape-files to just Ontario study area
exceed.05.Ontario = 
  sf1 %>% 
  filter(st_contains(Ontariosf.D, ., sparse = FALSE))
ggplot() + 
  geom_sf(data = exceed.05.Ontario, size = 0.25, color = "blue", fill = "cyan1") + 
  geom_sf(data = Ontariosf, size = 0.25, color = "darkred", fill = NA) + 
  ggtitle("Ontario study area - Exceedance CL 5 %") + 
  coord_sf()

#now look at data: Exceedance
exceed.05.Ontario$Exceedance6655<-exceed.05.Ontario$Exceedance+6656
#Exceedance values made positive by adding 6656 to original values
#because colours can't be assigned to negative numbers
ggplot() + 
  geom_sf(data = exceed.05.Ontario, size = 0, color = exceed.05.Ontario$Exceedance6655, fill = exceed.05.Ontario$Exceedance6655) + 
  geom_sf(data = Ontariosf, size = 0.25, color = "darkred", fill = NA) + 
  ggtitle("Ontario study area - Exceedance CL 5 %") + 
  coord_sf()

#Now create rasters for each variable of interest in the shape-file
library(fasterize)
#use Ontario study area raster as a template
Ontariostudyareamask<-raster("0_data/raw/spatial data/Canada shapefile/Ontariostudyareamask.tif")

r.exceed.05.Ontario.exceedance<-fasterize(
  exceed.05.Ontario,
  Ontariostudyareamask,
  field = "Exceedance",
  fun="max",
  background = NA_real_,
  by = NULL
)
plot(r.exceed.05.Ontario.exceedance)
writeRaster(r.exceed.05.Ontario.exceedance, filename="0_data/processed/Acidity/Exceedance05.Ontario.tif", overwrite=TRUE)
#saves exceedance values for critical loads of nitrogen and sulfur,
#above which the top 5 % of sensitive tree species are negatively affected, as a raster layer.
#Note, since the raw data are for a national layer, the top 5 % refers to
#the most sensitive tree species nation-wide, not within the Ring of Fire 
#region.

r.exceed.05.Ontario.ndep<-fasterize(
  exceed.05.Ontario,
  Ontariostudyareamask,
  field = "ndep",
  fun="max",
  background = NA_real_,
  by = NULL
)
plot(r.exceed.05.Ontario.ndep)
writeRaster(r.exceed.05.Ontario.ndep, filename="0_data/processed/Acidity/DepositedN.Ontario.tif", overwrite=TRUE)
#saves total deposited nitrogen as a raster layer

r.exceed.05.Ontario.sdep<-fasterize(
  exceed.05.Ontario,
  Ontariostudyareamask,
  field = "sdep",
  fun="max",
  background = NA_real_,
  by = NULL
)
plot(r.exceed.05.Ontario.sdep)
writeRaster(r.exceed.05.Ontario.sdep, filename="0_data/processed/Acidity/DepositedS.Ontario.tif", overwrite=TRUE)
#saves total deposited sulfur as a raster layer

r.exceed.05.Ontario.CL05N<-fasterize(
  exceed.05.Ontario,
  Ontariostudyareamask,
  field = "CLMAXN_5",
  fun="max",
  background = NA_real_,
  by = NULL
)
plot(r.exceed.05.Ontario.CL05N)
writeRaster(r.exceed.05.Ontario.CL05N, filename="0_data/processed/Acidity/CLMAXN_5percentTrees.Ontario.tif", overwrite=TRUE)
#saves max critical load of nitrogen before 5 % of trees are
#negatively affected, as raster layer

r.exceed.05.Ontario.CL05S<-fasterize(
  exceed.05.Ontario,
  Ontariostudyareamask,
  field = "CLMAXS_5",
  fun="max",
  background = NA_real_,
  by = NULL
)
plot(r.exceed.05.Ontario.CL05S)
writeRaster(r.exceed.05.Ontario.CL05S, filename="0_data/processed/Acidity/CLMAXS_5percentTrees.Ontario.tif", overwrite=TRUE)
#saves max critical load of sulfur before 5 % of trees are
#negatively affected, as raster layer


#Now do the same steps for the exceedance shape file based
#on critical loads before 20 % of an area's trees are negatively
#affected
#Exceedance measurements for critical loads sufficient to kill/damage 5 % of an area's trees
exceed.20 <- st_read("0_data/raw/spatial data/Acidity/soilCL20_exc_2010/soilCL20_exc_2010.shp")

#all spatial polygon data frames in which crs = 
#"crs         : +proj=longlat +datum=NAD83 +no_defs" 
#From GNM-National scripts on Github
lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
exceed.20.sf <- sf::st_as_sf(exceed.20)
sf2 <- st_set_crs(exceed.20.sf, "+proj=longlat +datum=NAD83 +no_defs")
sf2 <- st_transform(sf2, lcc_crs)

ggplot() + 
  geom_sf(data = sf2, size = 0.5, color = "blue", fill = "cyan1") + 
  ggtitle("Exceedance Map - CL 20 %") + 
  coord_sf()
#This layer is mostly Canada wide except for parts of Prairie
#Provinces and southeast Ontario

#filter exceedance shape-files to just Ontario
exceed.20.Ontario = 
  sf2 %>% 
  filter(st_contains(Ontariosf.D, ., sparse = FALSE))
ggplot() + 
  geom_sf(data = exceed.20.Ontario, size = 0.25, color = "blue", fill = "cyan1") + 
  geom_sf(data = Ontariosf, size = 0.25, color = "darkred", fill = NA) + 
  ggtitle("Ontario study area - Exceedance CL 20 %") + 
  coord_sf()

#now look at data: Exceedance
exceed.20.Ontario$Exceedance11796<-exceed.20.Ontario$Exceedance+11797
#Exceedance values made positive by adding 6656 to original values
#because colours can't be assigned to negative numbers
ggplot() + 
  geom_sf(data = exceed.20.Ontario, size = 0, color = exceed.20.Ontario$Exceedance11796, fill = exceed.20.Ontario$Exceedance11796) + 
  geom_sf(data = Ontariosf, size = 0.25, color = "darkred", fill = NA) + 
  ggtitle("Ontario study area - Exceedance CL 20 %") + 
  coord_sf()

r.exceed.20.Ontario.exceedance<-fasterize(
  exceed.20.Ontario,
  Ontariostudyareamask,
  field = "Exceedance",
  fun="max",
  background = NA_real_,
  by = NULL
)
plot(r.exceed.20.Ontario.exceedance)
writeRaster(r.exceed.20.Ontario.exceedance, filename="0_data/processed/Acidity/Exceedance20.Ontario.tif", overwrite=TRUE)
#saves exceedance values for critical loads of nitrogen and sulfur,
#above which 20 % of trees are negatively affected, as a raster layer

r.exceed.20.Ontario.CL20N<-fasterize(
  exceed.20.Ontario,
  Ontariostudyareamask,
  field = "CLMAXN_20",
  fun="max",
  background = NA_real_,
  by = NULL
)
plot(r.exceed.20.Ontario.CL20N)
writeRaster(r.exceed.20.Ontario.CL20N, filename="0_data/processed/Acidity/CLMAXN_20percentTrees.Ontario.tif", overwrite=TRUE)
#saves max critical load of nitrogen before 20 % of trees are
#negatively affected, as raster layer

r.exceed.20.Ontario.CL20S<-fasterize(
  exceed.20.Ontario,
  Ontariostudyareamask,
  field = "CLMAXS_20",
  fun="max",
  background = NA_real_,
  by = NULL
)
plot(r.exceed.20.Ontario.CL20S)
writeRaster(r.exceed.20.Ontario.CL20S, filename="0_data/processed/Acidity/CLMAXS_20percentTrees.Ontario.tif", overwrite=TRUE)
#saves max critical load of sulfur before 20 % of trees are
#negatively affected, as raster layer

