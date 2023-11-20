# compare and evaluate bird model predictions for the Hudson bay lowlands
library(dplyr)
library(terra)
library(sf)
library(purrr)
library(tidyr)
library(ggplot2)
library(stringr)

devtools::load_all(".")

in_dat_pth1 <- "analysis/data/derived_data/models_to_compare"

# where caribou ranges are stored, only needed for small example
sourceData <- "C:/Users/endicotts/Documents/gitprojects/MissisaBooPaper/data/inputNV/caribouRanges"

code6 <- list.files(file.path(in_dat_pth1, "Bateman_2020_sdms"),
                    pattern = "_breeding.*_suitability.tif",
                    recursive = TRUE) %>%
  str_extract("(bore.*boreal_forests_...._)(......\\d?)(_breed.*tif)", group = 2)

# spelling difference for burrowing owl
# setdiff(code6, ebirdst_runs$species_code)
code6 <- str_replace(code6, "borowl", "burowl") %>% str_subset("grgowl", negate = TRUE)



# 4 letter species code and 6 letter species code for species to compare
sp4 <- "BOWA" # will need a look up with this and sp6 to include bam models
sp6 <- "bohwax"

# number of points to compare across models
samp_size <- 1000

# For example use Missisa as study area
study_area <- read_sf(file.path("analysis/data/derived_data/0_data/raw/shapefiles", "ecozones.shp")) %>%
  filter(ZONE_NAME == "Hudson Plain")

sa_bbox <- st_bbox(study_area) %>% st_as_sfc()

# Use grid to sample a bunch of points to compare, how to standardize the
# predictions across different models with different outputs??
samp_pts1 <- terra::spatSample(terra::vect(study_area), size = samp_size,
                              method = "regular")

do_sp_compare <- function(sp6, sp4 = "XXXX", samp_pts, in_dat_pth){

  cat(paste0("## ", sp6))

  pats <- list(bateman = paste0(sp6, ".*_breeding.*suitability.tif$"),
               bam_rof = paste0(sp4, ".*mean.tif$"),
               ebird_st = paste0(sp6, "_abundance.*.tif"))

  fls <- map(pats, \(x) list.files(in_dat_pth, x, recursive = TRUE)) %>%
    # drop empty
    compact()

  rasts <- map(fls, \(x) rast(file.path(in_dat_pth, x)))

  # for some species the ebird file has multiple layers, keep breeding or resident
  rasts$ebird_st <- rasts$ebird_st[[which(names(rasts$ebird_st) %in% c("resident", "breeding"))]]

  crs_use <- terra::crs(rasts$bateman)
  rasts_crop <- map(rasts, \(x) crop_transform(x, bbx = sa_bbox, crs_out = crs_use))

  oldpar <- par(mfrow = c(2,2))
  iwalk(rasts_crop, \(x, nm) print(plot(x, main = nm)))
  par(oldpar)


  mod_samps <- map(rasts, \(x) extract_transform(x, samp_pts)[,2]) %>% bind_cols() %>%
    # Bateman models were multiplied by 10000 to store as integer
    mutate(ID = 1:n(), bateman = bateman/10000) %>%
    pivot_longer(-ID, names_to = "model", values_to = "prediction") %>%
    group_by(model) %>%
    mutate(rel_pred = scale(prediction, center = FALSE)[,1],
           rank_pred = rank(prediction, na.last = "keep")) %>%
    group_by(ID) %>%
    mutate(rank_sd = sd(rank_pred, na.rm = T))

  print(ggplot(mod_samps, aes(ID, prediction, col = model))+
    geom_point())

  print(ggplot(mod_samps, aes(ID, rel_pred, col = model))+
    geom_point())

  print(ggplot(mod_samps, aes(ID, rank_pred, col = model))+
    geom_point())

  print(ggplot(mod_samps, aes(ID, rank_sd, col = model))+
    geom_point())

}

walk(code6, \(x) do_sp_compare(x, samp_pts = samp_pts1, in_dat_pth = in_dat_pth1))
