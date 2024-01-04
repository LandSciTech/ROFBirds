library(wildRtrax)
library(tidyverse)
library(reproducible)

options(reproducible.cachePath = "analysis/data/cache_data",
        reproducible.cacheSpeed = "fast",
        reproducible.useTerra = TRUE)

# Note that you need to use 'WT_USERNAME' and 'WT_PASSWORD'
Sys.setenv(WT_USERNAME = 'guest', WT_PASSWORD = 'Apple123')

wt_auth()

# Download all of the published ARU data from CWS-ONT
aru_projects <- wt_get_download_summary(sensor_id = "ARU") %>%
  as_tibble()

# CWS-ONT seem to be the only ARUs but there are point counts from other
# projects: BAM, ON-ATLAS (BAM), Gov-ON, NRCAN, and NL-UofG But none of these
# are published. Most is included in what we already have from BAM and or the
# BBS

# Currently only downloading one because it takes a long time.
# cwsON_arus <- aru_projects %>%
#   filter(status == "Published - Public",
#          grepl('^CWS-Ontario Boreal Shield', project)) %>%
#   dplyr::mutate(data = purrr::map(
#     .x = project_id,
#     .f = ~wt_download_report(project_id = .x, sensor_id = "ARU",
#                              weather_cols = T, reports = "main"))) %>%
#   select(-project_id, -organization) %>%
#   unnest(data) %>%
#   unnest(cols = c(organization_id, project, sensor, tasks, status))

# write.csv(cwsON_arus, file = "analysis/data/raw_data/CWS_ARUs/CWSOntarioBorealShieldLowlandsTransition2022.csv", row.names = FALSE)

cwsON_arus <- read.csv("analysis/data/raw_data/CWS_ARUs/CWSOntarioBorealShieldLowlandsTransition2022.csv")

spp_tbl <- wt_get_species()

cwsON_arus_sp <- left_join(cwsON_arus, spp_tbl %>% select(species_code, species_class),
                           by = "species_code") %>%
  wt_tidy_species(remove = c("mammal", "amphibian", "abiotic", "insect", "unknown", "NONE"),
                  zerofill = TRUE)

# What should we use for validation? Could it be the presence absence and just
# assume that greater predicted abundance is also greater probability of
# occurrence? Or modify abundance/density predictions to be presence absence in
# a different way...?  Or should we use abundance for models that try to predict
# it? For BAM would then want to do QPAD steps but for ebird they have a
# different relative abundance. Ideally would follow all the same data prep
# steps as the training data... Or could look at correlation between obs
# abundance and the prediction?


cwsON_arus_wide <- cwsON_arus_sp %>%  wt_replace_tmtt(calc = "round") %>%
  wt_make_wide(sound = "all")

cwsON_arus_qpad <- wt_qpad_offsets(cwsON_arus_wide, together = F)

write.csv(cwsON_arus_qpad, "analysis/data/derived_data/models_to_compare/cwsON_ARUs_validation.csv", row.names = FALSE)
