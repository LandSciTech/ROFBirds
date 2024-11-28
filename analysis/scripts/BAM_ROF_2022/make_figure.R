# figure for report
library(sf)
library(tmap)
library(dplyr)

bcrs <- read_sf("analysis/data/derived_data/0_data/raw/shapefiles/bcr_canada_lcc.shp")

mfz <- read_sf("analysis/data/derived_data/0_data/raw/shapefiles/aou_ontario.shp")

ecozones <- read_sf("analysis/data/derived_data/0_data/raw/shapefiles/ecozones.shp")

study_area <- read_sf("analysis/data/derived_data/0_data/raw/shapefiles/ROF_RA_def.shp")

on_prov <- read_sf("analysis/data/derived_data/0_data/raw/shapefiles/ontario.shp") %>%
  st_transform(st_crs(bcrs))

load("analysis/data/derived_data/0_data/processed/BAMv6_RoFpackage_2022-01.RData")

pts <- xx1 %>% st_as_sf(coords = c("X", "Y")) %>% st_set_crs(4269) %>%
  select(PKEY_V4) %>% st_transform(st_crs(bcrs))

bcr_7812 <- bcrs %>% filter(BCR %in% c(7, 8, 12)) %>% st_intersection(on_prov) %>%
 filter(!(WATER == 1 & BCR == 8)) %>%
  sfheaders::sf_remove_holes() %>%
  st_cast("MULTIPOLYGON") %>%
  st_cast("POLYGON", do_split = TRUE) %>%
  mutate(area = st_area(.)) %>%
  filter(area > units::set_units(100000000, "m2"))

mfz <- sfheaders::sf_remove_holes(mfz)

hudP_borSh <- ecozones %>%
  st_transform(st_crs(bcrs)) %>%
  filter(ZONE_NAME %in% c("Hudson Plain", "Boreal Shield")) %>%
  st_intersection(on_prov)

fig <- tm_shape(hudP_borSh)+
  tm_fill(col = "ZONE_NAME", title = "Ecozone", alpha = 0.6,
          palette = c("#f4a582", "#92c5de"))+
  tm_shape(bcr_7812, is.master = TRUE) +
  tm_borders(col = "black", lwd = 1.5)+
  tm_shape(mfz)+
  tm_borders(col = "#0571b0", lwd = 1.5)+
  tm_shape(study_area)+
  tm_borders(col = "black", lwd = 2)+
  tm_shape(pts)+
  tm_dots(col = "#ca0020")+
  tm_shape(filter(bcr_7812, area > units::set_units(1000000000, "m2"))) +
  tm_text("BCR", size= 1, fontface = "bold")+
  tm_scale_bar(position = c("left","bottom"), text.size=1)+
  tm_compass(position = c("left", "top"), size=3)+
  tm_layout(frame=F, legend.position = c(0.75, 0.4),
            legend.title.size = 0.8, legend.height = -0.55)+
  # tm_add_legend(type = c("line"),
  #               labels = c("Managed Forest\nBoundary"),
  #               col = c("#0571b0"), lwd = 1.5)+
  tm_add_legend(type = c("line"),
                col = c("black", "black", "#0571b0"), lwd = c(2, 1, 1),
                labels = c("Study Area", "Bird Conservation\nRegion",
                           "Managed Forest\nBoundary"))+
  tm_add_legend(type = "symbol",
                labels = "     Count locations",
                size = 0.5,
                col = "#ca0020",
                shape = 19)

tmap_save(fig, "analysis/figures/birdPts_BCR_Ecozone_StudyArea.png", dpi = 600)

