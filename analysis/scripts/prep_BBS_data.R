# there were some issues in bbsAssistant pacakge but it does not appear to be
# actively maintained. So install from my fork
# remotes::install_github("see24/bbsAssistant", ref = "fix-QualityCurrentID-filter")

library(bbsAssistant)
library(dplyr)
library(reproducible)

# Set parameters for script
start_yr <- 1991
end_yr <- 2023

quality_filter <- "RunType == 1"

juris <- "ontario"

options(reproducible.cachePath = "analysis/data/cache_data")

# Not going to replicate creating offsets at this point since we don't need to
# but I will explain my understanding of how they are calculated.

# Step 1: #=====================================================================
# download BBS data, extract from files, filter to area and time range
# and based on quality, Get species letter codes, change format to one col per
# species and one row per stop/year/route
# Get the BBS data
# downloads the data the first time and will import it if already downloaded.
# using Cache from reproducible to cache long running functions
bbs_raw <- Cache(grab_bbs_data(sb_id = "625f151ed34e85fa62b7f926",
                               bbs_dir = "analysis/data/raw_data/bbs_2022"))

obs_ON_1991 <- Cache(munge_bbs_data(
  bbs_raw,
  states = region_codes %>% filter(State %in% juris) %>% pull(StateNum),
  zero.fill = FALSE,
  year.range = start_yr:end_yr,
  keep.stop.level.data = TRUE
))

obs_ON_1991_2 <- obs_ON_1991 %>%
  filter(!!rlang::parse_expr(quality_filter)) %>%
  # get 4 letter codes for sp names
  left_join(species_list %>% select(AOU, AOU4), by = "AOU") %>%
  # Some are missing so keep numeric codes for those
  mutate(sp_code = coalesce(AOU4, as.character(AOU)), .keep = "unused") %>%
  select(-contains("Car"), -contains("Noise")) %>%
  tidyr::pivot_longer(Stop1:Stop50, names_to = "Stop", values_to = "Abund") %>%
  filter(Abund != 0) %>%
  tidyr::pivot_wider(names_from = sp_code, values_from = Abund, values_fill = 0)


# Step 2: #=====================================================================
# Get locations of each stop on the routes from data supplied by CWS, not
# currently available online.

# Step 3: #=====================================================================
# Get time since sunrise for each stop on the route based on start and
# end times and number of stops

# Step 4: #=====================================================================
# Calculate offsets. Requires data and functions from this github
# repository: https://github.com/borealbirds/qpad-offsets/tree/main. This has
# been updated by BAM from the one they used originally. This uses species specific
# coefficients stored in the QPAD project to determine the offset using raster
# datasets provided in the repository.




