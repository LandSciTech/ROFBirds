#' ---
#' title: "Compare and evaluate bird model predictions for the Hudson bay lowlands"
#' author: "Sarah Endicott"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    toc_float: true
#' ---

#'
#' #### Notes on differences between the models:
#' The Bateman 2020 models are SDMs that predict habitat suitability, the BAM
#' national models predict density in number of signing males per ha and the
#' ebird models predict "relative abundance" defined as the the count by an
#' expert eBirder on a 1 hour, 2 kilometer traveling checklist at the optimal
#' time of day. The eBird models also don't make predictions in areas where
#' there was < 0.05% spatial coverage in a 3km grid cell.
#'
#' Maps show the prediction rasters for each species. Graphs compare the values
#' from each model at 1000 points using four different measures. Prediction is
#' the raw value, scaled_pred is the prediction divided by the root-mean-square,
#' rank_pred is the rank of each value relative to other sampled points after
#' removing points that are NA for any model, and rank_sd gives the standard
#' deviation of ranks across models for each point. Points were projected on the
#' fly to avoid re-projecting the rasters
#'
#'

#+ include=FALSE
knitr::opts_chunk$set(echo = FALSE)

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

# codes from Bateman file names for boreal forest birds (there are more available)
# include both 6 letter codes and 4 letter codes so use to create look up
code6 <- list.files(file.path(in_dat_pth1, "Bateman_2020_sdms"),
                    pattern = "_breeding.*_suitability.tif",
                    recursive = TRUE) %>%
  str_extract("(.*_[A-Z]{4}_)(......?\\d?)(_breed.*tif)", group = 2)

code4 <- list.files(file.path(in_dat_pth1, "Bateman_2020_sdms"),
                    pattern = "_breeding.*_suitability.tif",
                    recursive = TRUE) %>%
  str_extract("(.*_)([A-Z]{4})(_.*)(_breed.*tif)", group = 2)

sp_codes <- tibble(code4, code6)

# currently BAM nat only has species from A-E, have emailed to check
sp_BAM_nat <- list.files(file.path(in_dat_pth1, "BAM_v4_National"),
                         pattern = "pred.*CAN-Mean.tif",
                         recursive = TRUE) %>%
  str_extract("(pred-)(....)(-CAN-Mean.tif)", group = 2) %>%
  tibble::as_tibble_col(column_name = "bam_nat")

all_sp_codes <- full_join(sp_codes, sp_BAM_nat, by = c(code4 = "bam_nat"),
                          keep = TRUE) %>%
  left_join(ebirdst::ebirdst_runs %>% select(species_code),
            by = c(code6 = "species_code"), keep = TRUE)

# Only keep species in both Bateman boreal/eastern forests and generalists,
# ebird status and trends, and BAM national
sp_code_use <- all_sp_codes %>% filter(if_all(everything(), ~!is.na(.x)))

# number of points to compare across models
samp_size <- 1000

# For example use Missisa as study area
study_area <- read_sf(file.path("analysis/data/derived_data/0_data/raw/shapefiles",
                                "ecozones.shp")) %>%
  filter(ZONE_NAME == "Hudson Plain")

# sa_bbox <- st_bbox(study_area) %>% st_as_sfc()

# Use grid to sample a bunch of points to compare, how to standardize the
# predictions across different models with different outputs??
samp_pts1 <- terra::spatSample(terra::vect(study_area), size = samp_size,
                              method = "regular")

# test one species in all models
# do_sp_compare("aldfly", "ALFL", samp_pts = samp_pts1, in_dat_pth = in_dat_pth1)

#+ results='asis', fig.height=6, fig.width=7, warning=FALSE

out <- map2(sp_code_use$code6, sp_code_use$code4,
            \(x, y) do_sp_compare(x, y, samp_pts = samp_pts1, in_dat_pth = in_dat_pth1)) %>%
  bind_rows(.id = "species")


#' ## Table
out %>% group_by(species) %>%
  summarise(mean_rank_sd = mean(rank_sd, na.rm = TRUE) %>% round(2)) %>%
  DT::datatable()

write.csv(out, "analysis/data/derived_data/model_compare.csv")

# rmarkdown::render("analysis/scripts/compare_model_predictions.R",
#                   knit_root_dir = here::here(), envir = new.env())

# problem loading ebird in row 38
