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
#' Maps show the prediction rasters for each species. The first set of graphs
#' compare the values from each model at 1000 points using four different
#' measures. Prediction is the raw value, scaled_pred is the prediction divided
#' by the root-mean-square, rank_pred is the rank of each value relative to
#' other sampled points after removing points that are NA for any model, and
#' rank_sd gives the standard deviation of ranks across models for each point.
#' Points were projected on the fly to avoid re-projecting the rasters.
#'
#' The second set of graphs compare mean observed abundance across multiple
#' recordings from the same location to the model predictions. I know that this
#' doesn't fully make sense because the predictions are in different units and
#' were trained on differently prepared data but I figured we would still expect
#' them to be correlated or somewhat related as a starting point.
#'
#' The first table gives the mean of the rank_sds at each sample point for each
#' species. The second table gives model performance measures for the validation
#' data. corr, intercept, slope, r.squared, and nobs for a linear model of mean
#' observed abundance ~ prediction. auc and rmse use the binary observation data
#' and the model predictions. For the same reasons mentioned above this is just
#' a starting point.
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

terraOptions(progress=0)

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

# There was a problem with chispa raster so I downloaded the one from the
# website manually

# currently BAM nat only has species from A-E, have emailed to check
sp_BAM_nat <- list.files(file.path(in_dat_pth1, "BAM_v4_National"),
                         pattern = "pred.*CAN-Mean.tif",
                         recursive = TRUE) %>%
  str_extract("(pred-)(....)(-CAN-Mean.tif)", group = 2) %>%
  tibble::as_tibble_col(column_name = "bam_nat")

sp_BAM_rof <- list.files(file.path(in_dat_pth1, "BAM_ROF"),
                         pattern = "^.....tif",
                         recursive = TRUE) %>%
  str_extract("(....)(.tif)", group = 1) %>%
  tibble::as_tibble_col(column_name = "bam_rof")

all_sp_codes <- full_join(sp_codes, sp_BAM_nat, by = c(code4 = "bam_nat"),
                          keep = TRUE) %>%
  full_join(sp_BAM_rof, by = c(code4 = "bam_rof"), keep = TRUE) %>%
  left_join(ebirdst::ebirdst_runs %>% select(species_code),
            by = c(code6 = "species_code"), keep = TRUE) %>%
  arrange(code6)

bat_grps <- read.csv(file.path(in_dat_pth1, "Bateman_2020_sdms/spp_names_codes_group.csv"))

mis_codes <- all_sp_codes %>% filter(is.na(code6), !is.na(bam_rof)) %>% pull(bam_rof)

bat_grps %>% filter(Code %in% mis_codes) %>% count(Group)

# Only keep species in both Bateman boreal/eastern forests and generalists,
# ebird status and trends, and BAM national
sp_code_use <- all_sp_codes %>% filter(if_all(everything(), ~!is.na(.x)))

# number of points to compare across models
samp_size <- 1000

# For example use Missisa as study area
study_area <- read_sf("analysis/data/derived_data/0_data/processed/shapefiles/BCR7and8Ontario.shp")

study_area <- summarise(study_area, name = "study_area")

# sa_bbox <- st_bbox(study_area) %>% st_as_sfc()

# Use grid to sample a bunch of points to compare, how to standardize the
# predictions across different models with different outputs??
samp_pts1 <- terra::spatSample(terra::vect(study_area), size = samp_size,
                              method = "regular")

# validation data
val_dat1 <- read.csv(file.path(in_dat_pth1, "cwsON_ARUs_validation.csv"))
val_dat1 <- st_as_sf(val_dat1, coords = c("longitude", "latitude"), crs = 4326)

# test one species in all models
# do_sp_compare("aldfly", "ALFL", samp_pts = samp_pts1, val_dat = val_dat1,
#               in_dat_pth = in_dat_pth1)

#+ results='asis', fig.height=6, fig.width=7, warning=FALSE
do_sp_compare_pos <- possibly(do_sp_compare, otherwise = data.frame(ID = NA), quiet = FALSE)

out <- map2(sp_code_use$code6[1:2], sp_code_use$code4[1:2],
            \(x, y) do_sp_compare_pos(x, y, samp_pts = samp_pts1, val_dat = val_dat1,
                                      in_dat_pth = in_dat_pth1))

#' ## Tables
#+ tab.cap="Variability of ranks of sample points between models measured as the mean standard deviation of rank across the sampled cells."
mod_samps <- out %>% map("mod_samps") %>% bind_rows() %>%  group_by(sp4, sp6)
mod_samps %>%
  summarise(mean_rank_sd = mean(rank_sd, na.rm = TRUE) %>% round(2), .groups = "drop") %>%
  DT::datatable()

#+ tab.cap="Comparison of model performance on validation data set. Note that the model predictions are not all in the same units but we might still expect them to be related. Not sure what the best way to address this is."
mod_perf <- out %>% map("mod_perf") %>% bind_rows()
mod_perf %>%
  select(model, sp4, sp6, corr, intercept, slope, r.squared, nobs, auc, rmse) %>%
  mutate(across(where(is.numeric), \(x)round(x, 3))) %>%
  DT::datatable()

write.csv(mod_samps, "analysis/data/derived_data/mod_samps.csv", row.names = FALSE)
write.csv(mod_perf, "analysis/data/derived_data/mod_perf.csv", row.names = FALSE)
# rmarkdown::render("analysis/scripts/compare_model_predictions.R",
#                   knit_root_dir = here::here(), envir = new.env())

# TODO figure out why third graph isn't printing except for last species. Doesn't matter which graph it is either.
