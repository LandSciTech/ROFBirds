# R code for creating the bird-habitat database for running BRT SDM models for birds in Far North Ontario  


library(mefa4)
library(jsonlite)
library(sf)
library(rgeos)
library(raster)
library(spdep)
library(rgdal)


# Get the list of species for teh ring of fire
tab <- fromJSON("https://borealbirds.github.io/ring-of-fire/species/index.json")
rownames(tab) <- tab$id


# Attached is the file with the locations of ontario bam points north of bcr-12 (SS_V4)
# with designations for whether they are in the study area (study_areaYN),
# which ecozone they are in (ecozone), and whether or not they are in the area of the
# undertaking (aouYN).

# I assume you have the updated habitat data from Lionel,
# where the Ontario land cover goes from columns 101 to 154.
# use SScombo.final.distadded

## Load the list of BAM points in the RoF study area
load("0_data/processed/rof_bam_pts_ss.Rdata")
xx1 <- rof_bam_pts_ss
rm(rof_bam_pts_ss)
load("0_data/raw/point_counts_with_data/RoF_SpaDESspatialpointcountdata_April23.RData")
xx2 <- SScombo.eskerdata.added
rm(SScombo.eskerdata.added)
gc()

with(xx2, table(lc=eskerLCC250, pt=eskerpoint))


# Load the BAM V6 data
load("0_data/raw/BAMv6_ONBBS.RData")

SPP <- intersect(intersect(tab$id, colnames(y1)), colnames(y2))

# RoF study area boundary
poly <- st_read("0_data/processed/boundary.shp")
# BAM xy is NAD83, EPSG:4269
xy <- st_as_sf(xx1, coords = c("X", "Y"), crs = 4269)
xy <- st_transform(xy, st_crs(poly))

# Plot the locations and study area boundary
plot(xy$geometry, pch=".", col=ifelse(xx1$ecozone=="hudson_plain", 2, 4))
plot(poly$geometry, add=TRUE, border="grey", col="#00000022")

# Calculate distance from study area boundary
di <- st_distance(poly, xy) # units m because of crs
di <- di[1,] / 1000 # unit in km
xx1$dist <- as.numeric(di)

# Calculate distance from Hudson Plains BCR
ecozones <- st_transform(st_as_sf(readOGR(dsn = "0_data/raw/shapefiles", layer = "ecozones")), crs = st_crs(xy))
hp <- ecozones[which(ecozones$ZONE_NAME == "Hudson Plain"), ]
dh <- st_distance(hp, xy)
dh <- dh[1, ]/1000
xx1$dist_hp <- as.numeric(dh)

nn <- dnearneigh(st_coordinates(xy), d1 = 0, d2 = 2500, row.names = xy$SS_V4)
nn1 <- sapply(1:length(nn), function(x) length(nn[[x]]))
xx1$num_neigh <- nn1

rownames(xx1) <- xx1$PKEY_V4
rownames(xx2) <- xx2$PKEY_V4
PKEY <- rownames(xx1)

# Create the offsets file
off <- rbind(o1[,SPP], o2[,SPP])[PKEY,]
y <- rbind(y1[,SPP], y2[,SPP])[PKEY,]
xx2 <- xx2[PKEY,22:ncol(xx2)]

# Make sure treecover, elev, and offsets aren't NA
keep <- !is.na(xx2$treecover) & !is.na(xx2$elev) & rowSums(is.na(off)) == 0
xx1 <- xx1[keep,]
xx2 <- xx2[keep,]
y <- y[keep,]
off <- off[keep,]
xx2 <- xx2[,colSums(is.na(xx2))==0]


## high correlations & constant variables
get_cn <- function(z, rmax=0.9) {
  SD <- apply(z, 2, sd, na.rm=TRUE)
  COR <- cor(z[,SD > 0.0001], use="pairwise.complete.obs")
  cr <- mefa:::stack.dist(as.dist(COR), dim.names = TRUE)
  cr <- cr[!is.na(cr$dist),]
  cr <- cr[order(abs(cr$dist), decreasing = TRUE),]
  cr[,1] <- as.character(cr[,1])
  cr[,2] <- as.character(cr[,2])
  while(any(abs(cr$dist) > rmax)) {
    i <- cr[1,2]
    j <- cr[,1] == i | cr[,2] == i
    cr <- cr[!j,]
  }
  union(as.character(cr[,1]), as.character(cr[,2]))
}

# Initial variable list reduced to 110 through discussion with the RoF team, created as cn1
cn1 <- c("eskerpoint",
         "elev", "treecover", "LIDARheight", "road_yesno", "TPI", "TRI",
         "slope", "agriculture_G750.O", "bedrock_G750.O", "bog_G750.O",
         "communities_G750.O", "coniftreed_G750.O", "decidtreed_G750.O",
         "disturbance_G750.O", "fen_G750.O", "heath_G750.O", "marsh_G750.O",
         "mixedtreed_G750.O", "mudflat_G750.O", "openwater_G750.O", "shoreline_G750.O",
         "sparsetreed_G750.O", "swamp_G750.O", "treedupland_G750.O", "turbidwater_G750.O",
         "G750LandCover_NonVeg_v1.grd", "G750LandCover_Veg_v1.grd", "G750LandCover_VegNonTreed_v1.grd",
         "G750LandCover_VegTreed_v1.grd", "G750Species_Abie_Bal_v1.grd",
         "G750Species_Abie_Spp_v1.grd", "G750Species_Acer_Neg_v1.grd",
         "G750Species_Acer_Pen_v1.grd", "G750Species_Acer_Rub_v1.grd",
         "G750Species_Acer_Sac_v1.grd", "G750Species_Acer_Sah_v1.grd",
         "G750Species_Acer_Spi_v1.grd", "G750Species_Acer_Spp_v1.grd",
         "G750Species_Alnu_Rub_v1.grd", "G750Species_Alnu_Spp_v1.grd",
         "G750Species_Betu_All_v1.grd", "G750Species_Betu_Pap_v1.grd",
         "G750Species_Betu_Pop_v1.grd", "G750Species_Betu_Spp_v1.grd",
         "G750Species_Carp_Car_v1.grd", "G750Species_Cary_Cor_v1.grd",
         "G750Species_Fagu_Gra_v1.grd", "G750Species_Frax_Ame_v1.grd",
         "G750Species_Frax_Nig_v1.grd", "G750Species_Frax_Pen_v1.grd",
         "G750Species_Genc_Spp_v1.grd", "G750Species_Genh_Spp_v1.grd",
         "G750Species_Jugl_Cin_v1.grd", "G750Species_Jugl_Nig_v1.grd",
         "G750Species_Juni_Vir_v1.grd", "G750Species_Lari_Lar_v1.grd",
         "G750Species_Lari_Spp_v1.grd", "G750Species_Malu_Spp_v1.grd",
         "G750Species_Ostr_Vir_v1.grd", "G750Species_Pice_Abi_v1.grd",
         "G750Species_Pice_Eng_v1.grd", "G750Species_Pice_Gla_v1.grd",
         "G750Species_Pice_Mar_v1.grd", "G750Species_Pice_Rub_v1.grd",
         "G750Species_Pice_Spp_v1.grd", "G750Species_Pinu_Alb_v1.grd",
         "G750Species_Pinu_Ban_v1.grd", "G750Species_Pinu_Con_v1.grd",
         "G750Species_Pinu_Mon_v1.grd", "G750Species_Pinu_Pon_v1.grd",
         "G750Species_Pinu_Res_v1.grd", "G750Species_Pinu_Spp_v1.grd",
         "G750Species_Pinu_Str_v1.grd", "G750Species_Pinu_Syl_v1.grd",
         "G750Species_Popu_Bal_v1.grd", "G750Species_Popu_Gra_v1.grd",
         "G750Species_Popu_Spp_v1.grd", "G750Species_Popu_Tre_v1.grd",
         "G750Species_Popu_Tri_v1.grd", "G750Species_Prun_Pen_v1.grd",
         "G750Species_Prun_Ser_v1.grd", "G750Species_Quer_Alb_v1.grd",
         "G750Species_Quer_Mac_v1.grd", "G750Species_Quer_Rub_v1.grd",
         "G750Species_Quer_Spp_v1.grd", "G750Species_Sali_Spp_v1.grd",
         "G750Species_Sorb_Ame_v1.grd", "G750Species_Thuj_Occ_v1.grd",
         "G750Species_Tili_Ame_v1.grd", "G750Species_Tsug_Can_v1.grd",
         "G750Species_Tsug_Spp_v1.grd", "G750Species_Ulmu_Ame_v1.grd",
         "G750SpeciesGroups_Broadleaf_Spp_v1.grd", "G750SpeciesGroups_Needleleaf_Spp_v1.grd",
         "G750SpeciesGroups_Unknown_Spp_v1.grd", "G750Structure_Biomass_Branch_v1.grd",
         "G750Structure_Biomass_Foliage_v1.grd", "G750Structure_Biomass_StemBark_v1.grd",
         "G750Structure_Biomass_StemWood_v1.grd", "G750Structure_Biomass_TotalDead_v1.grd",
         "G750Structure_Biomass_TotalLiveAboveGround_v1.grd", "G750Structure_Stand_Age_v1.grd",
         "G750Structure_Stand_CrownClosure_v1.grd", "G750Structure_Stand_Height_v1.grd",
         "G750Structure_Volume_Merch_v1.grd", "G750Structure_Volume_Total_v1.grd",
         "biomass2015.ntems", "volume2015.ntems", "height2015.ntems")


# Run the get_cn function to remove corellated and constant variables from the cn1 list
cn2 <- sort(get_cn(xx2[,cn1]))


# Filter spp with >= 20 detections
SPP <- colnames(y)[colSums(y) >= 20]
y <- y[,SPP]
off <- off[,SPP]


# resample to have the same n in Boreal Shield than in Hudson Plain
n0 <- sum(xx1$ecozone == "hudson_plain")
i1 <- which(xx1$ecozone == "hudson_plain")
i2 <- which(xx1$ecozone != "hudson_plain")

# Weight the selection probability of i2 data by distance from the Hudson Plains and local point density 
i2.pts <- xx1[i2, ]
i2.dist <-  1 - (i2.pts$dist_hp - min(i2.pts$dist_hp))/(max(i2.pts$dist_hp) - min(i2.pts$dist_hp))
i2.dens <- 1 - (i2.pts$num_neigh - min(i2.pts$num_neigh))/(max(i2.pts$num_neigh) - min(i2.pts$num_neigh))

ww <- i2.dist*i2.dens


resample_fun <- function(i) {
  if (i == 1) {
    j1 <- i1
    j2 <- sample(i2, n0, replace=FALSE, prob = ww)
  } else {
    j1 <- sample(i1, n0, replace=FALSE)
    j2 <- sample(i2, n0, replace=FALSE, prob = ww)
    
  }
  c(j1, j2)
}

# Resample 250 times
B <- 250
set.seed(1)
BB <- sapply(1:B, resample_fun)

save(y, off, xx1, xx2, cn2, SPP, BB, file="0_data/processed/BAMv6_RoFpackage_2022-04.RData")



# Create table of detections x ecozone 
spp_det_table <- do.call(rbind, lapply(1:length(SPP), function(x){
  spp <- SPP[x] 
  det <- ifelse(as.numeric(y[, spp]) > 0, 1, 0)
  b <- sum(det[which(xx1$ecozone!="hudson_plain")])
  h <- sum(det[which(xx1$ecozone=="hudson_plain")])
  data.frame(spp = spp, boreal = b, hudson_plain = h, total = sum(det))
}))

saveRDS(spp_det_table, file = "3_outputs/tables/spp_det_table.Rds")





