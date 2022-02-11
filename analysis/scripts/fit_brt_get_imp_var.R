# Sarah's version of the modeling scripts from BAM

library(tidyverse)

in_dat_pth <- "analysis/data/derived_data"

# Step to analysis
# Get bird data and environmental variables into format for modelling
# Fit one BRT per species to determine important variables
# Fit n boostrap BRTs per species
# Make predictions using each bootstrap model and summarise

# load the prepared data from BAM scripts
load(file.path(in_dat_pth, "0_data/processed/BAMv6_RoFpackage_2022-01.RData"))

# rename the object loaded so they are literate
# first need to figure out what they are!!
