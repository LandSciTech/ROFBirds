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

cwsON_arus <- Cache(aru_projects %>%
                      filter(status == "Published - Public",
                             grepl('^CWS-Ontario', project)) %>%
                      dplyr::mutate(data = purrr::map(
                        .x = project_id,
                        .f = ~wt_download_report(project_id = .x, sensor_id = "ARU",
                                                 weather_cols = T, reports = "main"))) %>%
                      select(-project_id, -organization) %>%
                      unnest(data))

spp_tbl <- wt_get_species()

cwsON_arus_sp <- left_join(cwsON_arus, spp_tbl %>% select(species_code, species_class),
                           by = "species_code") %>%
  wt_tidy_species(remove = c("mammal", "amphibian", "abiotic", "insect", "unknown", "NONE"),
                  zerofill = F)



