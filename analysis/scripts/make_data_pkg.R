# Make a data set that could be shared with others for modelling birds

# This just includes the observations and offsets and the point locations.
# Other meta data about the survey may be available but was not included by BAM.

# Combine pt meta data and observations and save as csv also save offsets

in_dat_pth <- "analysis/data/derived_data"

# load data that BAM prepared and rename objects to be literate
# data saved at line 196 Andy_scripts/01-2_create-brt-database.R
load(file.path(in_dat_pth, "0_data/processed/BAMv6_RoFpackage_2022-01.RData"),
     temp_env <- new.env())

BAM_package <- as.list(temp_env)

rm(temp_env)

pt_meta <- BAM_package$xx1

sp_dens_obs <- BAM_package$y

# change from sparse matrix to data frame
sp_dens_obs <- as.matrix(sp_dens_obs) %>% as_tibble(rownames = "PKEY")

sp_dens_obs_pts <- inner_join(pt_meta %>%
                                select(X, Y, PKEY_V4, survey_year, STATEABB, BCR,
                                       BCRNAME, ecozone, aou = aouYN) %>%
                                rename_with(tolower),
                              sp_dens_obs, by = c(pkey_v4 = "PKEY"))

write.csv(sp_dens_obs_pts, "analysis/data/derived_data/ONBBS_BAM_bird_density_BCR_7-8_1991-2019.csv",
          row.names = FALSE)

sp_off <- BAM_package$off

as_tibble(sp_off, rownames = "pkey_v4") %>%
  write.csv("analysis/data/derived_data/ONBBS_BAM_bird_offsets_BCR_7-8_1991-2019.csv",
            row.names = FALSE)
