# R code for creating the bird database for modeling spp. density in far North Ontario 


library(data.table)
library(Matrix)
library(googledrive)
library(sqldf)
library(rgdal)
library(sf)
library(ggplot2)
library(tidyverse)


load("0_data/raw/BAMv6_ONBBS.RData")


load("0_data/raw/point_counts_with_data/RoF_SpaDESspatialpointcountdata_Feb18.RData")

geo <- "0_data/raw/shapefiles"

l <- ogrListLayers(geo)

# Load maps, convert to sf objects, and transform to the NAD83 / Ontario MNR Lambert CRS
t <- readOGR(dsn = geo, layer = "ROF_RA_def")
b <- readOGR(dsn = geo, layer = "bcr_canada_lcc")

ra <- st_as_sf(readOGR(geo, "ROF_RA_def"))
cr <- st_transform(st_as_sf(readOGR(geo, "rof_caribou_ranges")), crs = st_crs(ra))
pts <- st_cast(st_transform(st_as_sf(SScombo2), crs = st_crs(ra)), "POINT")
on <- st_transform(st_as_sf(readOGR(dsn = geo, layer = "ontario")), crs = st_crs(ra))
aou <- st_transform(st_as_sf(readOGR(dsn = geo, layer = "aou_ontario")), crs = st_crs(ra))
ecozones <- st_transform(st_as_sf(readOGR(dsn = geo, layer = "ecozones")), crs = st_crs(ra))
zones <- ecozones[which(ecozones$ZONE_NAME == "Boreal Shield" | ecozones$ZONE_NAME == "Hudson Plain"), ]
zones_on <- st_intersection(zones, on, sparse = F)
bcr_on <- st_intersection(st_transform(st_as_sf(b), crs = st_crs(ra)), on, sparse = F)
bcr_rof <- bcr_on[which(bcr_on$BCR %in% c(7, 8)), ]

# Remove duplicated points (those with >1 survey) 
pts_ss <- pts[-which(duplicated(pts$SS_V4)), ]

pts_shield <- pts_ss[which(st_intersects(pts_ss, zones[which(zones$ZONE_NAME == "Boreal Shield"), ], sparse = F)), ]
pts_plain <- pts_ss[which(st_intersects(pts_ss, zones[which(zones$ZONE_NAME == "Hudson Plain"), ], sparse = F)), ]
pts_aou <- pts_ss[which(st_intersects(pts_ss, aou, sparse = FALSE)), ]

pts_ra <- pts_ss[st_intersects(pts_ss, ra, sparse = FALSE), ]
pts_ra_aou <- pts_ra[which(st_intersects(pts_ra, aou, sparse = FALSE)), ]
pts_ra_north <- pts_ra[-which(st_intersects(pts_ra, aou, sparse = FALSE)), ]
pts_ra_shield <- pts_ra[which(st_intersects(pts_ra, zones[which(zones$ZONE_NAME == "Boreal Shield"), ], sparse = F)), ]
pts_ra_plain <- pts_ra[which(st_intersects(pts_ra, zones[which(zones$ZONE_NAME == "Hudson Plain"), ], sparse = F)), ]

save(ra, cr, pts, on, aou, zones_on, pts_ss, pts_shield, pts_plain, pts_aou, pts_ra, pts_ra_aou, pts_ra_north, pts_ra_shield, pts_ra_plain, SScombo.Beaudoinadded, file = file.path("0_data/processed/points_spatial_data.Rdata"))

st_write(ra, "0_data/processed/boundary.shp")
st_write(bcr_rof, "0_data/processed/shapefiles/BCR7and8Ontario.shp", append = F)

rof_bam_pts_ss <- st_drop_geometry(pts_ss)
rof_bam_pts_ss$study_areaYN <- ifelse(st_intersects(pts_ss, ra, sparse = FALSE), 1, 0)
rof_bam_pts_ss$ecozone <- ifelse(st_intersects(pts_ss, zones[which(zones$ZONE_NAME == "Boreal Shield"), ], sparse = F), "boreal_sheild", "hudson_plain")
rof_bam_pts_ss$aouYN <- ifelse(st_intersects(pts_ss, aou, sparse = FALSE), 1, 0)

save(rof_bam_pts_ss, file = file.path("0_data/processed/rof_bam_pts_ss.Rdata"))



png(filename = "3_outputs/maps/bcr_map.png",
    width = 2000, height = 2000, units = "px", pointsize = 0.5,
    bg = "white", res = 300)
ggplot(data = bcr_on[which(bcr_on$BCR %in% c(7, 8, 12)), ]) + geom_sf(aes(fill = BCR)) + geom_sf(data = pts_ss) + 
  geom_sf(data = ra, fill = NA, size = 1.5, colour = "red") + theme_bw() + theme(panel.grid = element_blank())
dev.off()








