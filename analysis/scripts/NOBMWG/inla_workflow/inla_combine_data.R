# ******************************************************************************
# ******************************************************************************
# Implementation of David Isle's Landscape Distribution Modeling ECCC workflow
# for Northern Ontario: Combining data
# ******************************************************************************
# ******************************************************************************
# Includes determining offsets and getting eBird ranges
# ******************************************************************************
# Notable changes of method:
# *
# ******************************************************************************



library(tidyverse)
library(naturecounts)
library(napops)
library(ebirdst)
library(sf)
devtools::load_all(".")

# use reproducible package to Cache long running code
options(reproducible.cachePath = "analysis/data/cache_data")

Study_Area_bound <- read_sf("analysis/data/derived_data/INLA_data/study_area.gpkg")

analysis_data <- readRDS("analysis/data/derived_data/INLA_data/analysis_data_covars.rds")

all_surveys <- analysis_data$all_surveys
full_count_matrix <- analysis_data$full_count_matrix

BSC_species <- naturecounts::search_species_code() %>%
  # drops duplicate species_ids keeping the first.
  distinct(species_id, .keep_all = TRUE) %>%
  rename(BSC_spcd = BSCDATA, species_scientific_name = scientific_name) %>%
  dplyr::select(species_id, BSC_spcd,species_scientific_name,english_name,french_name) %>%
  unique()

# Identify species to run analysis for #========================================

# Process species names / labels (from Birds Canada)

ON_spcd <- search_species_code() %>% rename(spcd = BSCDATA, CommonName = english_name)
countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

# Process Species Names so they fit
ON_spcd$Label <- NA
ON_spcd$CommonName[ON_spcd$CommonName=="Rock Pigeon (Feral Pigeon)"] <- "Rock Pigeon"

for (i in 1:nrow(ON_spcd)) {
  Name <- ON_spcd$CommonName[i]
  if(nchar(Name) > 13){
    if(countSpaces(Name)>0){
      ON_spcd$Label[i] <- gsub(" "," \n",Name)
    }

  }
  else {
    ON_spcd$Label[i] <- ON_spcd$CommonName[i]
  }

}

# ---
# Only fit models for species detected in at least 50 atlas squares
# ---

n_detections <- full_count_matrix
rownames(n_detections) <- all_surveys$sq_id
n_detections <- n_detections %>%
  reshape2::melt() %>%
  dplyr::rename(sq_id = Var1, Species_Code_BSC = Var2, detected = value) %>%
  subset(detected>0) %>%
  group_by(Species_Code_BSC) %>%
  summarize(n_squares = length(unique(sq_id)),
            n_detections = sum(detected>0))

species_to_model <- left_join(n_detections,BSC_species,
                              by = c("Species_Code_BSC" = "BSC_spcd")) %>%

  # Remove erroneous observations
  subset(english_name %!in% c("passerine sp.",
                              "duck sp.",
                              "new world sparrow sp.",
                              "new world warbler sp.",
                              "blackbird sp.",
                              "gull sp.",
                              "Greater/Lesser Yellowlegs",
                              "Scolopacidae sp.",
                              "vireo sp.",
                              "Catharus sp.",
                              "Worthen's Sparrow",
                              "woodpecker sp.",
                              "Alder/Willow Flycatcher (Traill's Flycatcher)",
                              "Aythya sp.",
                              "Bat sp.",
                              "bird sp.",
                              "chickadee sp.",
                              "finch sp.",
                              "Golden-winged/Blue-winged Warbler",
                              "goose sp.",
                              "hawk sp.",
                              "new world flycatcher sp.",
                              "Northern Flicker (Yellow-shafted)",
                              "owl sp.",
                              "Philadelphia/Red-eyed Vireo",
                              "swallow sp.",
                              "swan sp.",
                              "tern sp.",
                              "Ula-ai-hawane",
                              "Yellow-rumped Warbler (Myrtle)"))

dim(species_to_model)
species_to_model$english_name %>% sort()

# ******************************************************************
# PART 4: IDENTIFY SPECIES WITH DETECTABILITY OFFSETS AVAILABLE # ------

napops_species <- list_species() %>% rename(Species_Code_NAPOPS = Species,
                                            Common_Name_NAPOPS = Common_Name,
                                            Scientific_Name_NAPOPS = Scientific_Name)

species_to_model <- left_join(species_to_model,napops_species[,c("Species_Code_NAPOPS","Common_Name_NAPOPS","Removal","Distance")],
                              by = c("Species_Code_BSC" = "Species_Code_NAPOPS"))


species_to_model$offset_exists <- FALSE
species_to_model$EDR <- NA
species_to_model$cue_rate <- NA
species_to_model$log_offset_5min <- 0


# Extract QPAD offsets if available
for (i in 1:nrow(species_to_model)){

  sp = species_to_model$Species_Code_BSC[i]
  print(sp)
  offset_exists <- FALSE
  sp_napops <- subset(napops_species,Species_Code_NAPOPS == sp)

  if (nrow(sp_napops)>0){
    if (sp_napops$Removal == 1 & sp_napops$Distance == 1){

      species_to_model$offset_exists[i] <- TRUE
      species_to_model$cue_rate[i] <- cue_rate(species = sp,od = 153, tssr = 0, model = 1)[3] %>% as.numeric()
      species_to_model$EDR[i] <- edr(species = sp,road = FALSE, forest = 0.5,model = 1)[3] %>% as.numeric()

      # Calculate A and p, which jointly determine offset
      A_metres <- c(pi*species_to_model$EDR[i]^2)
      p <- 1-exp(-5*species_to_model$cue_rate[i])

      species_to_model$log_offset_5min[i] <- log(A_metres * p)

    }
  }
}


# Get eBird Ranges #============================================================
# Download, trim, and save ebird range limits for species

# Need to set access key (only once) before downloading ranges
# ebirdst::set_ebirdst_access_key()

# Ensure path is correctly set
# usethis::edit_r_environ()

# This should read:
# EBIRDST_KEY='ntqm1ha68fov'
# EBIRDST_DATA_DIR='C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Landbird-Distribution-Modeling-ECCC/Data/Spatial/eBird/'

# # Spatial layers for ON
# ON_BCR <- st_read("../../../Data/Spatial/National/BCR/BCR_Terrestrial_master.shp")  %>%
#   subset(PROVINCE_S == "ONTARIO") %>%
#   st_make_valid() %>%
#   dplyr::select(BCR, PROVINCE_S)
#
# ONBoundary <- ON_BCR %>% st_union()

# List to contain sf objects for species ranges
species_ranges <- list()

# wasn't acutally filtered above...
species_to_model <- species_to_model %>% filter(n_squares > 50)

for (sp_code in species_to_model$Species_Code_BSC){
  print(sp_code)
  if (file.exists("analysis/data/derived_data/eBird_ranges_ON.RDS")){
    species_ranges <- readRDS("analysis/data/derived_data/eBird_ranges_ON.RDS")
    message("loading range from file")
  }

  # Check is species already in the list
  if (sp_code %in% names(species_ranges)) next

  # Download and trim ebird range for this species

  species_name = ON_spcd$CommonName[which(ON_spcd$spcd == sp_code)]
  species_label = ON_spcd$Label[which(ON_spcd$spcd == sp_code)]

  # Check if species is available
  check <- get_species(species_name)
  if (length(check)==0){
    print(paste0(sp_code, " not available in eBird"))
    next
  }

  if (length(check)>0 & is.na(check)){
    print(paste0(sp_code, " not available in eBird"))
    next
  }

  ebirdst_download_status(species_name, download_ranges = TRUE, download_abundance = FALSE,
                          dry_run = FALSE, pattern = "_27km_")

  path <- get_species_path(species_name)

  range <- load_ranges(species = species_name, resolution = "27km")

  range <- range %>% subset(season %in% c("resident","breeding")) %>%
    st_transform(.,crs = st_crs(Study_Area_bound)) %>%
    st_union() %>%
    st_crop(Study_Area_bound)

  species_ranges[[sp_code]] <- range
  saveRDS(species_ranges,file="analysis/data/derived_data/eBird_ranges_ON.RDS")

  # TODO: Finish download of more species
} # close species loop


# ******************************************************************
# ******************************************************************
# Save
# ******************************************************************
# ******************************************************************

analysis_data_package <- list(

  all_surveys = all_surveys, # Survey information (protocol, location, time, etc)
  full_count_matrix = full_count_matrix, # counts of each species for each survey

  pca = analysis_data$pca,

  # Contains covariates on ON-wide grid
  ONGrid = analysis_data$ONGrid,

  # Species to include in analysis (and number of detections in point counts)
  species_to_model = species_to_model,

  # Species codes and common names
  ON_spcd = ON_spcd,

  # eBird range limits
  species_ranges = species_ranges

)

saveRDS(analysis_data_package,"analysis/data/derived_data/INLA_data/analysis_data_package.rds")

write.csv(analysis_data_package$species_to_model,
          file = "analysis/data/derived_data/INLA_data/species_to_model.csv",
          row.names = FALSE)
