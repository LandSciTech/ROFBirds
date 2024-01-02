# get ebird maps from status and trends data using ebirdst package

library(ebirdst)
library(stringr)
library(purrr)
library(dplyr)
library(googledrive)

# set_ebirdst_access_key("q74koe64s91l")

in_dat_pth <- "analysis/data/derived_data/models_to_compare"

# get all bird data that is included in the Audubon data we have

code6 <- list.files(file.path(in_dat_pth, "Bateman_2020_sdms"),
           pattern = "_breeding.*_suitability.tif",
           recursive = TRUE) %>%
  str_extract("(.*_[A-Z]{4}_)(......?\\d?)(_breed.*tif)", group = 2)

# Only keep species in both Bateman boreal forests and ebird status and trends
code6 <- intersect(code6, ebirdst::ebirdst_runs$species_code)




purrr::map(code6, \(x){ebirdst_download(x, pattern = "abundance_.*_mean_hr",
                                        path = file.path(in_dat_pth, "eBird_SandT"))})



# data from google drive for BAM national models didn't have permissions to
# download programmatically so did manually from here
# https://drive.google.com/drive/folders/1exWa6vfhGo1DNUL4ei2baDz77as7jYzY?usp=sharing
