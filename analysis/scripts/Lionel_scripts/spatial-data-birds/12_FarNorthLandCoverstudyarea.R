#The purpose of this script is to obtain land cover data from the Far North Land 
#Cover product from 
#https://geohub.lio.gov.on.ca/documents/01783626bc1a49a6ba5e18db0620890a/about, and
#obtain local-scale and landscape-scale values of
#the proportion of land occupied by different Far North cover types. The initial raster 
#was constructed in ArcGIS from a geodatabase obtained from the link above and had to be
#constructed from three separate raster layers that were in different projections. 
#Those three layers (FNLC_utm15, FNLC_utm16, FNLC_utm17) and the initial merged raster 
#(FNLCmergeLCC) are stored in:
#"0_data/raw/spatial data/Far North Land Cover" along with the original geodatabase
#"FarNorthLandCover.zip".

#Once the three raster layers were merged, they were aggregated and resampled 
#in R from 30-m to 250-m resolution.

#The reason for using these data, although they don't cover all
#of the BCR 7 and BCR 8 modelling area, is because of the fine-
#scale distinction among different wetland types. For point counts
#outside of the Far North Land Cover layer boundaries other
#land cover data are used to assign values.
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


#First, we need a tif file of the Ontario study area
#Use Ontario study area grid with 250-m resolution we created for Beaudoin
#layers
Ontariostudyarea<-raster("0_data/raw/spatial data/Canada shapefile/Ontariostudyareamask.tif")

##Far North Land Cover layer (merged 3 separate reprojected rasters resampled to 250 m from 30 m)
FNLCmergeLCC<-raster("0_data/raw/spatial data/Far North Land Cover/FNLCmergeLCC.tif")
plot(FNLCmergeLCC)

#create rasters for individual cover types
openwater<-FNLCmergeLCC
values(openwater) = ifelse(values(openwater)==1,1,0)
plot(openwater)
writeRaster(openwater, filename="0_data/processed/Far North Land Cover/openwaterFNLC.tif", overwrite=TRUE)

turbidwater<-FNLCmergeLCC
values(turbidwater) = ifelse(values(turbidwater)==2,1,0)
plot(turbidwater)
writeRaster(turbidwater, filename="0_data/processed/Far North Land Cover/turbidwaterFNLC.tif", overwrite=TRUE)

mudflat<-FNLCmergeLCC
values(mudflat) = ifelse(values(mudflat)==3,1,0)
plot(mudflat)
writeRaster(mudflat, filename="0_data/processed/Far North Land Cover/mudflatFNLC.tif", overwrite=TRUE)

intertidal<-FNLCmergeLCC
values(intertidal) = ifelse(values(intertidal)==4,1,0)
plot(intertidal)
writeRaster(intertidal, filename="0_data/processed/Far North Land Cover/intertidalmarshFNLC.tif", overwrite=TRUE)

supertidal<-FNLCmergeLCC
values(supertidal) = ifelse(values(supertidal)==5,1,0)
plot(supertidal)
writeRaster(supertidal, filename="0_data/processed/Far North Land Cover/supertidalmarshFNLC.tif", overwrite=TRUE)

freshwater<-FNLCmergeLCC
values(freshwater) = ifelse(values(freshwater)==6,1,0)
plot(freshwater)
writeRaster(freshwater, filename="0_data/processed/Far North Land Cover/freshwatermarshFNLC.tif", overwrite=TRUE)

heath<-FNLCmergeLCC
values(heath) = ifelse(values(heath)==7,1,0)
plot(heath)
writeRaster(heath, filename="0_data/processed/Far North Land Cover/heathFNLC.tif", overwrite=TRUE)

thicketswamp<-FNLCmergeLCC
values(thicketswamp) = ifelse(values(thicketswamp)==8,1,0)
plot(thicketswamp)
writeRaster(thicketswamp, filename="0_data/processed/Far North Land Cover/thicketswampFNLC.tif", overwrite=TRUE)

coniferswamp<-FNLCmergeLCC
values(coniferswamp) = ifelse(values(coniferswamp)==9,1,0)
plot(coniferswamp)
writeRaster(coniferswamp, filename="0_data/processed/Far North Land Cover/coniferswampFNLC.tif", overwrite=TRUE)

decidswamp<-FNLCmergeLCC
values(decidswamp) = ifelse(values(decidswamp)==10,1,0)
plot(decidswamp)
writeRaster(decidswamp, filename="0_data/processed/Far North Land Cover/decidswampFNLC.tif", overwrite=TRUE)

openfen<-FNLCmergeLCC
values(openfen) = ifelse(values(openfen)==11,1,0)
plot(openfen)
writeRaster(openfen, filename="0_data/processed/Far North Land Cover/openfenFNLC.tif", overwrite=TRUE)

treedfen<-FNLCmergeLCC
values(treedfen) = ifelse(values(treedfen)==12,1,0)
plot(treedfen)
writeRaster(treedfen, filename="0_data/processed/Far North Land Cover/treedfenFNLC.tif", overwrite=TRUE)

openbog<-FNLCmergeLCC
values(openbog) = ifelse(values(openbog)==13,1,0)
plot(openbog)
writeRaster(openbog, filename="0_data/processed/Far North Land Cover/openbogFNLC.tif", overwrite=TRUE)

treedbog<-FNLCmergeLCC
values(treedbog) = ifelse(values(treedbog)==14,1,0)
plot(treedbog)
writeRaster(treedbog, filename="0_data/processed/Far North Land Cover/treedbogFNLC.tif", overwrite=TRUE)

sparsetreed<-FNLCmergeLCC
values(sparsetreed) = ifelse(values(sparsetreed)==15,1,0)
plot(sparsetreed)
writeRaster(sparsetreed, filename="0_data/processed/Far North Land Cover/sparsetreedFNLC.tif", overwrite=TRUE)

decidtreed<-FNLCmergeLCC
values(decidtreed) = ifelse(values(decidtreed)==16,1,0)
plot(decidtreed)
writeRaster(decidtreed, filename="0_data/processed/Far North Land Cover/decidtreedFNLC.tif", overwrite=TRUE)

mixedtreed<-FNLCmergeLCC
values(mixedtreed) = ifelse(values(mixedtreed)==17,1,0)
plot(mixedtreed)
writeRaster(mixedtreed, filename="0_data/processed/Far North Land Cover/mixedtreedFNLC.tif", overwrite=TRUE)

coniftreed<-FNLCmergeLCC
values(coniftreed) = ifelse(values(coniftreed)==18,1,0)
plot(coniftreed)
writeRaster(coniftreed, filename="0_data/processed/Far North Land Cover/coniftreedFNLC.tif", overwrite=TRUE)

dist.nonwoody<-FNLCmergeLCC
values(dist.nonwoody) = ifelse(values(dist.nonwoody)==19,1,0)
plot(dist.nonwoody)
writeRaster(dist.nonwoody, filename="0_data/processed/Far North Land Cover/dist.nonwoodyFNLC.tif", overwrite=TRUE)

dist.treed<-FNLCmergeLCC
values(dist.treed) = ifelse(values(dist.treed)==20,1,0)
plot(dist.treed)
writeRaster(dist.treed, filename="0_data/processed/Far North Land Cover/dist.treedFNLC.tif", overwrite=TRUE)

sand.grav.mines<-FNLCmergeLCC
values(sand.grav.mines) = ifelse(values(sand.grav.mines)==21,1,0)
plot(sand.grav.mines)
writeRaster(sand.grav.mines, filename="0_data/processed/Far North Land Cover/sand.grav.minesFNLC.tif", overwrite=TRUE)

bedrock<-FNLCmergeLCC
values(bedrock) = ifelse(values(bedrock)==22,1,0)
plot(bedrock)
writeRaster(bedrock, filename="0_data/processed/Far North Land Cover/bedrockFNLC.tif", overwrite=TRUE)

communities<-FNLCmergeLCC
values(communities) = ifelse(values(communities)==23,1,0)
plot(communities)
writeRaster(communities, filename="0_data/processed/Far North Land Cover/communitiesFNLC.tif", overwrite=TRUE)

agriculture<-FNLCmergeLCC
values(agriculture) = ifelse(values(agriculture)==24,1,0)
plot(agriculture)
writeRaster(agriculture, filename="0_data/processed/Far North Land Cover/agricultureFNLC.tif", overwrite=TRUE)

fnlchab <- list.files("0_data/processed/Far North Land Cover/",pattern="FNLC.tif$")
st.fnlchab <- stack(raster(paste0("0_data/processed/Far North Land Cover/",fnlchab[1])))
for (i in 2:length(fnlchab)) { st.fnlchab <- addLayer(st.fnlchab, raster(paste0("0_data/processed/Far North Land Cover/",fnlchab[i])))}
names(st.fnlchab) <- gsub("FNLC","local",names(st.fnlchab))
writeRaster(st.fnlchab, filename="0_data/processed/Far North Land Cover/FNLC local and landscape/FNLC.lcc.local.grd", format="raster",overwrite=TRUE)

# obtain weighted sums of neighourhood cells using Gaussian filter with sigma=250, and 750m for Beaudoin and CTI layers, save outputs as rasters
st.fnlchabB<-stack("0_data/processed/Far North Land Cover/FNLC local and landscape/FNLC.lcc.local.grd")

## sigma = 250m
fw250<-focalWeight(x=st.fnlchabB,d=250,type="Gauss")
st.fnlchab_Gauss250<-stack(focal(st.fnlchabB[[1]],w=fw250,na.rm=TRUE))
names(st.fnlchab_Gauss250)<-gsub("local","_G250",names(st.fnlchabB)[[1]])
for(i in 2:nlayers(st.fnlchabB)){
  st.fnlchab_Gauss250<-addLayer(st.fnlchab_Gauss250,focal(st.fnlchabB[[i]],w=fw250,na.rm=TRUE))
  names(st.fnlchab_Gauss250)[i]<-gsub("local","_G250",names(st.fnlchabB)[[i]])
}
st.fnlchab_Gauss250<-brick(st.fnlchab_Gauss250)
writeRaster(st.fnlchab_Gauss250, filename="0_data/processed/Far North Land Cover/FNLC local and landscape/FNLC.lcc.G250.grd", format="raster",overwrite=TRUE)


# ## sigma = 750m
fw750<-focalWeight(x=st.fnlchabB,d=750,type="Gauss")
st.fnlchab_Gauss750<-stack(focal(st.fnlchabB[[1]],w=fw750,na.rm=TRUE))
names(st.fnlchab_Gauss750)<-gsub("local","_G750",names(st.fnlchabB)[[1]])
for(i in 2:nlayers(st.fnlchabB)){
  st.fnlchab_Gauss750<-addLayer(st.fnlchab_Gauss750,focal(st.fnlchabB[[i]],w=fw750,na.rm=TRUE))
  names(st.fnlchab_Gauss750)[i]<-gsub("local","_G750",names(st.fnlchabB)[[i]])
}
st.fnlchab_Gauss750<-brick(st.fnlchab_Gauss750)
writeRaster(st.fnlchab_Gauss750, filename="0_data/processed/Far North Land Cover/FNLC local and landscape/FNLC.lcc.G750.grd", format="raster",overwrite=TRUE)
