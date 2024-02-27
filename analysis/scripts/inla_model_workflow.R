# run settings #=============================
test_mode <- TRUE
if(test_mode){
  warning("Running in test mode with small data set.")
}

# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------
my_packs = c('tidyverse',
             'openxlsx',
             'RColorBrewer',
             'viridis',
             'ggrepel',
             'scales',
             'wildRtrax',
             'lubridate',
             'sf')

lapply(my_packs, library, character.only = TRUE)
rm(my_packs)

# currently using main but might need development branch
#remotes::install_github("ABbiodiversity/wildRtrax@development")

# load custom functions
devtools::load_all(".")

wt_auth() # Need your wildtrax username

options(reproducible.cachePath = "analysis/data/cache_data")

# Species list #================================================================
# Reconcile species lists between Birds Canada and NatureCounts
# ------------------------------------------------

BSC_species <- naturecounts::search_species_code() %>%
  distinct(species_id, .keep_all = TRUE) %>%
  rename(BSC_spcd = BSCDATA, species_scientific_name = scientific_name) %>%
  dplyr::select(BSC_spcd,species_scientific_name,english_name,french_name) %>%
  unique() %>%
  mutate(index = 1:nrow(.))

WT_species <- wildRtrax::wt_get_species() %>%
  # assuming we only want birds
  filter(species_class == "AVES") %>%
  rename(WT_spcd = species_code) %>%
  dplyr::select(WT_spcd,species_scientific_name) %>%
  unique()

all_species <- full_join(BSC_species,WT_species) %>% relocate(BSC_spcd,WT_spcd) %>% select(-index) %>% unique()

all_species <- subset(all_species, species_scientific_name %!in% c("NULL NULL"," "))

rm(WT_species)

# wildTrax Download #===========================================================
# Download WildTrax datasets within a study area boundary

# -----------------------------------------------
# Polygon delineating study area boundary
# -----------------------------------------------

if(!file.exists("analysis/data/raw_dataBCR_Terrestrial/BCR_Terrestrial_master.shp")){
  download.file(url = "https://birdscanada.org/download/gislab/bcr_terrestrial_shape.zip?_ga=2.123276619.1787785417.1706559261-142409804.1705510123",
                destfile = "analysis/data/raw_data/bcr_terrestrial_shape.zip", mode = "wb")
  unzip("analysis/data/raw_data/bcr_terrestrial_shape.zip",
        exdir = "analysis/data/raw_data")
}

Study_Area <- st_read("analysis/data/raw_data/BCR_Terrestrial/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  dplyr::select(BCR, PROVINCE_S)

# Only keep BCR 7 and 8 for N ON project
Study_Area <- Study_Area %>% filter(BCR %in% c(7,8))

# ---------------------------------------------------------
# Identify projects that are part of the Study Area
# ---------------------------------------------------------

# Use wildRtrax conversion from ARU to PC
# ARU_projects <- wt_get_download_summary(sensor_id = 'ARU') %>% subset(sensor == "ARU") %>% arrange(project)
PC_projects <- wt_get_download_summary(sensor_id = 'PC') %>% arrange(project)

# Subset projects for testing
if(test_mode){
  PC_projects <- PC_projects %>% filter(status == "Published - Public", tasks <1000) %>%
    filter(str_detect(project, "Ontario")) %>%
    slice(1:10)
}

# ---------------------------------------------------------
# Download projects with 'PC' data (includes projects converted from ARU)
# ---------------------------------------------------------

PC_fulldat <- reproducible::Cache(map(PC_projects$project_id, ~wt_dl_loc_year(.x, sens_id = "PC")) %>%
  bind_rows())

write.csv(PC_fulldat, file = "analysis/data/interim_data/WildTrax_PC_locations.csv", row.names = FALSE)

# ---------------------------------------------------------
# Subset to data within Study Area Boundary
# ---------------------------------------------------------

PC_fulldat <- read.csv(file = "analysis/data/interim_data/WildTrax_PC_locations.csv")

# Remove sites outside study area
PC_sf <- PC_fulldat %>%
  na.omit() %>%
  st_as_sf(coords=c("longitude","latitude"),crs=4326, remove = FALSE) %>%
  st_transform(crs = st_crs(Study_Area))%>%
  st_intersection(Study_Area)

# Summarize number of surveys from each project
PC_summary <- PC_sf %>%
  as.data.frame() %>%
  group_by(project_id) %>%
  summarize(n_surveys = n()) %>%
  subset(n_surveys > 5) %>%
  left_join(PC_projects[,c("project_id","organization","project")],
            by = join_by(project_id)) %>%
  arrange(desc(n_surveys)) %>%
  select(n_surveys,organization,project)

# Save summary of projects within the study area boundary
write.csv(PC_summary,file = "analysis/data/interim_data/WildTrax_PC_summary.csv",row.names = FALSE)

PC_summary  <- read.csv(file = "analysis/data/interim_data/WildTrax_PC_summary.csv")

# ---------------------------------------------------------
# Download species records and recording info from each Point Count survey
# ---------------------------------------------------------

# Note that all location information is stored in point count dataframe
PC_counts <- data.frame()
for (i in 1:nrow(PC_summary)){
  project_name = PC_summary$project[i]
  PID = subset(PC_projects, project == project_name)$project_id
  counts <- wt_download_report(project_id = PID, sensor_id = "PC",
                               reports = c("main"), weather_cols = FALSE) %>%
    mutate(individual_count = as.character(individual_count))
  PC_counts <- bind_rows(PC_counts,counts)
}

# ---------------------------------------------------------
# Convert WildTrax species codes to Birds Canada species codes when they disagree
# ---------------------------------------------------------

species <- subset(all_species,WT_spcd %in% unique(PC_counts$species_code))

species_to_fix <- subset(species,WT_spcd != BSC_spcd)
for (i in 1:nrow(species_to_fix)){
  PC_counts$species_code[which(PC_counts$species_code == species_to_fix$WT_spcd[i])] <- species_to_fix$BSC_spcd[i]
}

# Remove species with no BSC code
cat("Species removed:")
subset(PC_counts, species_code %!in% all_species$BSC_spcd)$species_common_name %>%
  table() %>% sort(decreasing = TRUE) %>% as.data.frame()

PC_counts <- subset(PC_counts, species_code %in% all_species$BSC_spcd)

PC_counts <- PC_counts %>%
  # make it work with PC by adding cols name same as ARU
  mutate(recording_date_time = survey_date,
         observer_id = observer,
         individual_count = ifelse(individual_count == "TMTC", "TMTT", individual_count)) %>%
  wt_replace_tmtt() %>%
  # remove extra columns
  select(-observer_id, -recording_date_time)



# ---------------------------------------------------------
# Dataframe to track information associated with each survey
# ---------------------------------------------------------

# 1 row per survey
PC_surveys <- PC_counts %>%
  select(-detection_distance,-detection_time,-species_code,-species_common_name,-species_scientific_name,-individual_count,-detection_heard,-detection_seen,-detection_comments) %>%
  unique()

# ---------------------------------------------------------
# Date of each survey
# ---------------------------------------------------------

# Three that have no time are made NA by this so I am not converting
# PC_surveys$survey_date <- lubridate::ymd_hms(PC_surveys$survey_date)

# ---------------------------------------------------------
# Append total duration of each survey
# ---------------------------------------------------------

table(PC_surveys$survey_duration_method)
method_definitions <- rbind(c(survey_duration_method = "0-1-2-3-4-5-6-7-8-9-10min", Survey_Duration_Minutes = 10),
                            c(survey_duration_method = "0-10min", Survey_Duration_Minutes = 10),
                            c(survey_duration_method = "0-3-10min", Survey_Duration_Minutes = 10),
                            c(survey_duration_method = "0-3-5-10min", Survey_Duration_Minutes = 10),
                            c(survey_duration_method = "0-3-5min", Survey_Duration_Minutes = 5),
                            c(survey_duration_method = "0-3min", Survey_Duration_Minutes = 3),
                            c(survey_duration_method = "0-5-10min", Survey_Duration_Minutes = 5),
                            c(survey_duration_method = "0-9min", Survey_Duration_Minutes = 9),
                            c(survey_duration_method = "0-5min", Survey_Duration_Minutes = 5),
                            c(survey_duration_method = "0-1min", Survey_Duration_Minutes = 1)) %>%
  as.data.frame()

if(any(PC_surveys$survey_duration_method %!in% method_definitions$survey_duration_method)){
  missing_method <- PC_surveys %>%
    filter(survey_duration_method %!in% method_definitions$survey_duration_method) %>%
    distinct(survey_duration_method)
  stop("survey duration method missing for ",
       paste0(missing_method$survey_duration_method, collapse = ", "))
}

PC_surveys <- left_join(PC_surveys, method_definitions,
                        by = join_by(survey_duration_method))


table(PC_surveys$Survey_Duration_Minutes, useNA = "always")


# ---------------------------------------------------------
# Append maximum distance of each survey
# ---------------------------------------------------------

table(PC_surveys$survey_distance_method)
method_definitions <- rbind(c(survey_distance_method = "0m-50m-100m", Max_Distance_Metres = 100),
                            c(survey_distance_method = "0m-50m-100m-150m-INF", Max_Distance_Metres = Inf),
                            c(survey_distance_method = "0m-50m-100m-INF", Max_Distance_Metres = Inf),
                            c(survey_distance_method = "0m-INF", Max_Distance_Metres = Inf)) %>%
  as.data.frame()

if(any(PC_surveys$survey_distance_method %!in% method_definitions$survey_distance_method)){
  missing_method <- PC_surveys %>%
    filter(survey_distance_method %!in% method_definitions$survey_distance_method) %>%
    distinct(survey_distance_method)
  stop("survey distance method missing for ",
       paste0(missing_method$survey_distance_method, collapse = ", "))
}


PC_surveys <- left_join(PC_surveys, method_definitions,
                        join_by(survey_distance_method))

# ---------------------------------------------------------
# Remove records outside the study area boundary
# ---------------------------------------------------------

PC_surveys <- PC_surveys %>%
  subset(!is.na(latitude) & !is.na(longitude)) %>%
  st_as_sf(coords=c("longitude","latitude"),crs=4326, remove = FALSE) %>%
  st_transform(crs = st_crs(Study_Area)) %>%
  st_set_agr("constant") %>%
  st_intersection(st_set_agr(Study_Area, "constant"))

# ---------------------------------------------------------
# Create a matrix of counts for each species
# ---------------------------------------------------------

PC_counts <- subset(PC_counts, survey_id %in% PC_surveys$survey_id)
PC_counts$survey_id <- factor(PC_counts$survey_id, levels = PC_surveys$survey_id)

# Total number of individuals detected per survey (sum of individual counts)
PC_counts <- PC_counts %>%
  group_by(species_code,survey_id) %>%
  summarize(total_count = sum(as.numeric(individual_count))) %>%
  as.data.frame() %>%
  pivot_wider(names_from = species_code,
              values_from = total_count,
              values_fill = 0,
              id_expand = TRUE,
              names_expand = TRUE) %>%
  arrange(survey_id)

# Same ordering as PC_surveys
mean(PC_counts$survey_id == PC_surveys$survey_id) # should be 1

# TODO stopped here

# *********************************************************
# *********************************************************
# COMBINE ARU AND POINT COUNTS INTO SINGLE DATAFRAME
# RENAME AND SELECT RELEVANT COLUMNS
# *********************************************************
# *********************************************************

# Combine SPT and SPM into single dataframe
ARU_recordings_combined <- bind_rows(ARU_recordings_SPT,ARU_recordings_SPM) %>%

  rename(Project_Name = project,
         survey_ID = recording_id,
         Date_Time = recording_date_time,
         Latitude = latitude,
         Longitude = longitude) %>%

  mutate(Data_Source = "WildTrax",
         Survey_Type = paste0("ARU_",Transcription_Method),
         Survey_Duration_Minutes = Duration_Seconds/60,
         Max_Distance_Metres = Inf) %>%

  dplyr::select(Data_Source,
                Project_Name,
                Survey_Type,
                survey_ID,
                Latitude,
                Longitude,
                Date_Time,
                Survey_Duration_Minutes,
                Max_Distance_Metres)

# Point Counts
PC_surveys <- PC_surveys %>%

  rename(Project_Name = project,
         survey_ID = survey_id,
         Date_Time = survey_date,
         Latitude = latitude,
         Longitude = longitude) %>%

  mutate(Data_Source = "WildTrax",
         Survey_Type = "Point_Count",
         Survey_Duration_Minutes = as.numeric(Survey_Duration_Minutes),
         Max_Distance_Metres = as.numeric(Max_Distance_Metres)) %>%

  dplyr::select(Data_Source,
                Project_Name,
                Survey_Type,
                survey_ID,
                Latitude,
                Longitude,
                Date_Time,
                Survey_Duration_Minutes,
                Max_Distance_Metres)

# Combined
WT_surveyinfo <- bind_rows(ARU_recordings_combined, PC_surveys)


# *********************************************************
# *********************************************************
# Count matrix (combined)
# *********************************************************
# *********************************************************

# All species in WT database
WT_species <- wildRtrax::wt_get_species() %>% subset(species_class == "AVES")
WT_species_codes <- WT_species$species_code

# Fill in matrix
WT_matrix <- matrix(0,nrow=nrow(WT_surveyinfo), ncol = length(WT_species$species_code),
                    dimnames = list(NULL,WT_species$species_code))

for (spp in WT_species_codes){

  if (spp %in% colnames(ARU_counts_SPT)) WT_matrix[which(WT_surveyinfo$Survey_Type == "ARU_SPT"),which(colnames(WT_matrix) == spp)] <- as.data.frame(ARU_counts_SPT)[,spp]
  if (spp %in% colnames(ARU_counts_SPM)) WT_matrix[which(WT_surveyinfo$Survey_Type == "ARU_SPM"),which(colnames(WT_matrix) == spp)] <- as.data.frame(ARU_counts_SPM)[,spp]
  if (spp %in% colnames(PC_counts)) WT_matrix[which(WT_surveyinfo$Survey_Type == "Point_Count"),which(colnames(WT_matrix) == spp)] <- as.data.frame(PC_counts)[,spp]

}

# Remove species that were never detected
WT_matrix <- WT_matrix[,-which(colSums(WT_matrix)==0)]

# Fix particular species manually
#WT_matrix[,"YRWA"] <- WT_matrix[,"YRWA"] + WT_matrix[,"MYWA"]
#WT_matrix <- WT_matrix[,-which(colnames(WT_matrix) == "MYWA")]

# *********************************************************
# *********************************************************
# Save file
# *********************************************************
# *********************************************************

WT_dat <- list(WT_surveyinfo = WT_surveyinfo,
               WT_matrix = WT_matrix)
saveRDS(WT_dat, file = "../Data_Cleaned/WildTrax/WT_dat.rds")
