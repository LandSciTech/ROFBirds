# get ebird maps from status and trends data using ebirdst package

library(ebirdst)
library(stringr)
library(purrr)
library(dplyr)
# set_ebirdst_access_key("q74koe64s91l")

in_dat_pth <- "analysis/data/derived_data/models_to_compare"

# get all bird data that is included in the Audubon data we have

code6 <- list.files(file.path(in_dat_pth, "Bateman_2020_sdms"),
           pattern = "_breeding.*_suitability.tif",
           recursive = TRUE) %>%
  str_extract("(bore.*boreal_forests_...._)(......\\d?)(_breed.*tif)", group = 2)

# spelling difference for burrowing owl
# setdiff(code6, ebirdst_runs$species_code)
code6 <- str_replace(code6, "borowl", "burowl") %>% str_subset("grgowl", negate = TRUE)

purrr::map(code6, \(x){ebirdst_download(x, pattern = "abundance_.*_mean_hr",
                                        path = file.path(in_dat_pth, "eBird_SandT"))})

