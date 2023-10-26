# get ebird maps from status and trends data using ebirdst package

library(ebirdst)
set_ebirdst_access_key("q74koe64s91l")

in_dat_pth <- "analysis/data/derived_data/models_to_compare"

path <- ebirdst_download(species = "alder flycatcher", pattern = "abundance_full-year_mean_hr",
                         path = file.path(in_dat_pth, "eBird_SandT"))
