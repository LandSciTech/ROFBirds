#The purpose of this script is to obtain land cover data from the Ontario Land 
#Cover Complilation product (v.0.2) from 
#https://geohub.lio.gov.on.ca/documents/lio::ontario-land-cover-compilation-v-2-0/about, and
#obtain local-scale and landscape-scale values of
#the proportion of land occupied by different land cover types. The initial raster 
#was obtained in ArcGIS from a geodatabase obtained from the link above. 
#The initial layer was aggregated and resampled from 30-m to 250-m 
#resolution. Initial and intermediate layers are stored in:
#"0_data/raw/spatial data/Ontario Land Cover" along with the original geodatabase
#"Ontario_Land_Cover_Compilation_Version2.gdb".

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

##Ontario Land Cover layer (merged 3 separate raster products (incl.
#Far North layer, and SOLRIS layer) resampled to 250 m from 30 m)
#The purpose of using the Ontario Land Cover layer is those point counts
#that occur outside the Far North Layer and therefore lack data if we
#rely on the Far North data alone. The Ontario layer will be processed in 
#the same way, then if point counts are lacking FNLC values, values will
#be taken from the Ontario Land Cover.

OLC<-raster("0_data/raw/spatial data/Ontario Land Cover Compilation - FGDB/OLCC_V2_250m.tif")
plot(OLC)

#create rasters for individual cover types
openwater<-OLC
values(openwater) = ifelse(values(openwater)==1,1,0)
plot(openwater)
writeRaster(openwater, filename="0_data/processed/Ontario Land Cover/openwaterOLC.tif", overwrite=TRUE)

turbidwater<-OLC
values(turbidwater) = ifelse(values(turbidwater)==2,1,0)
plot(turbidwater)
writeRaster(turbidwater, filename="0_data/processed/Ontario Land Cover/turbidwaterOLC.tif", overwrite=TRUE)

shoreline<-OLC
values(shoreline) = ifelse(values(shoreline)==3,1,0)
plot(shoreline)
writeRaster(shoreline, filename="0_data/processed/Ontario Land Cover/shorelineOLC.tif", overwrite=TRUE)

mudflat<-OLC
values(mudflat) = ifelse(values(mudflat)==4,1,0)
plot(mudflat)
writeRaster(mudflat, filename="0_data/processed/Ontario Land Cover/mudflatOLC.tif", overwrite=TRUE)

marsh<-OLC
values(marsh) = ifelse(values(marsh)==5,1,0)
plot(marsh)
writeRaster(marsh, filename="0_data/processed/Ontario Land Cover/marshOLC.tif", overwrite=TRUE)

swamp<-OLC
values(swamp) = ifelse(values(swamp)==6,1,0)
plot(swamp)
writeRaster(swamp, filename="0_data/processed/Ontario Land Cover/swampOLC.tif", overwrite=TRUE)

fen<-OLC
values(fen) = ifelse(values(fen)==7,1,0)
plot(fen)
writeRaster(fen, filename="0_data/processed/Ontario Land Cover/fenOLC.tif", overwrite=TRUE)

bog<-OLC
values(bog) = ifelse(values(bog)==8,1,0)
plot(bog)
writeRaster(bog, filename="0_data/processed/Ontario Land Cover/bogOLC.tif", overwrite=TRUE)

heath<-OLC
values(heath) = ifelse(values(heath)==10,1,0)
plot(heath)
writeRaster(heath, filename="0_data/processed/Ontario Land Cover/heathOLC.tif", overwrite=TRUE)

sparsetreed<-OLC
values(sparsetreed) = ifelse(values(sparsetreed)==11,1,0)
plot(sparsetreed)
writeRaster(sparsetreed, filename="0_data/processed/Ontario Land Cover/sparsetreedOLC.tif", overwrite=TRUE)

treedupland<-OLC
values(treedupland) = ifelse(values(treedupland)==12,1,0)
plot(treedupland)
writeRaster(treedupland, filename="0_data/processed/Ontario Land Cover/treeduplandOLC.tif", overwrite=TRUE)

decidtreed<-OLC
values(decidtreed) = ifelse(values(decidtreed)==13,1,0)
plot(decidtreed)
writeRaster(decidtreed, filename="0_data/processed/Ontario Land Cover/decidtreedOLC.tif", overwrite=TRUE)

mixedtreed<-OLC
values(mixedtreed) = ifelse(values(mixedtreed)==14,1,0)
plot(mixedtreed)
writeRaster(mixedtreed, filename="0_data/processed/Ontario Land Cover/mixedtreedOLC.tif", overwrite=TRUE)

coniftreed<-OLC
values(coniftreed) = ifelse(values(coniftreed)==15,1,0)
plot(coniftreed)
writeRaster(coniftreed, filename="0_data/processed/Ontario Land Cover/coniftreedOLC.tif", overwrite=TRUE)

disturbance<-OLC
values(disturbance) = ifelse(values(disturbance)==18,1,0)
plot(disturbance)
writeRaster(disturbance, filename="0_data/processed/Ontario Land Cover/disturbanceOLC.tif", overwrite=TRUE)

sand.grav.mines<-OLC
values(sand.grav.mines) = ifelse(values(sand.grav.mines)==25,1,0)
plot(sand.grav.mines)
writeRaster(sand.grav.mines, filename="0_data/processed/Ontario Land Cover/sand.grav.minesOLC.tif", overwrite=TRUE)

bedrock<-OLC
values(bedrock) = ifelse(values(bedrock)==26,1,0)
plot(bedrock)
writeRaster(bedrock, filename="0_data/processed/Ontario Land Cover/bedrockOLC.tif", overwrite=TRUE)

communities<-OLC
values(communities) = ifelse(values(communities)==27,1,0)
plot(communities)
writeRaster(communities, filename="0_data/processed/Ontario Land Cover/communitiesOLC.tif", overwrite=TRUE)

agriculture<-OLC
values(agriculture) = ifelse(values(agriculture)==28,1,0)
plot(agriculture)
writeRaster(agriculture, filename="0_data/processed/Ontario Land Cover/agricultureOLC.tif", overwrite=TRUE)

fnlchab <- list.files("0_data/processed/Ontario Land Cover/",pattern="OLC.tif$")
st.fnlchab <- stack(raster(paste0("0_data/processed/Ontario Land Cover/",fnlchab[1])))
for (i in 2:length(fnlchab)) { st.fnlchab <- addLayer(st.fnlchab, raster(paste0("0_data/processed/Ontario Land Cover/",fnlchab[i])))}
names(st.fnlchab) <- gsub("OLC","local.O",names(st.fnlchab))
writeRaster(st.fnlchab, filename="0_data/processed/Ontario Land Cover/OLC local and landscape/OLC.lcc.local.grd", format="raster",overwrite=TRUE)
#added O on end of variable names to distinguish from FNLC variable names

# obtain weighted sums of neighourhood cells using Gaussian filter with sigma=250, and 750m, save outputs as rasters
st.fnlchabB<-stack("0_data/processed/Ontario Land Cover/OLC local and landscape/OLC.lcc.local.grd")

## sigma = 250m
fw250<-focalWeight(x=st.fnlchabB,d=250,type="Gauss")
st.fnlchab_Gauss250<-stack(focal(st.fnlchabB[[1]],w=fw250,na.rm=TRUE))
names(st.fnlchab_Gauss250)<-gsub("local","_G250",names(st.fnlchabB)[[1]])
for(i in 2:nlayers(st.fnlchabB)){
  st.fnlchab_Gauss250<-addLayer(st.fnlchab_Gauss250,focal(st.fnlchabB[[i]],w=fw250,na.rm=TRUE))
  names(st.fnlchab_Gauss250)[i]<-gsub("local","_G250",names(st.fnlchabB)[[i]])
}
st.fnlchab_Gauss250<-brick(st.fnlchab_Gauss250)
writeRaster(st.fnlchab_Gauss250, filename="0_data/processed/Ontario Land Cover/OLC local and landscape/OLC.lcc.G250.grd", format="raster",overwrite=TRUE)


# ## sigma = 750m
fw750<-focalWeight(x=st.fnlchabB,d=750,type="Gauss")
st.fnlchab_Gauss750<-stack(focal(st.fnlchabB[[1]],w=fw750,na.rm=TRUE))
names(st.fnlchab_Gauss750)<-gsub("local","_G750",names(st.fnlchabB)[[1]])
for(i in 2:nlayers(st.fnlchabB)){
  st.fnlchab_Gauss750<-addLayer(st.fnlchab_Gauss750,focal(st.fnlchabB[[i]],w=fw750,na.rm=TRUE))
  names(st.fnlchab_Gauss750)[i]<-gsub("local","_G750",names(st.fnlchabB)[[i]])
}
st.fnlchab_Gauss750<-brick(st.fnlchab_Gauss750)
writeRaster(st.fnlchab_Gauss750, filename="0_data/processed/Ontario Land Cover/OLC local and landscape/OLC.lcc.G750.grd", format="raster",overwrite=TRUE)
