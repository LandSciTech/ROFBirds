# ******************************************************************************
# ******************************************************************************
# Implementation of David Isle's Landscape Distribution Modeling ECCC workflow
# for Northern Ontario: Data preparation
# ******************************************************************************
# ******************************************************************************
# Includes downloading data, standardizing format, determining offsets
# ******************************************************************************
# Notable changes of method:
# * Using the wildRtrax PC version of data collected by ARUs. This means that the
# conversion from ARU format to PC format is done on wildRtrax whereas David did
# this in his code.
# * Only using the point count data from NatureCounts, not clear what David did
# since non PC steps are commented out
# * Using location to determine timezone since it does not seem to be recorded
# ******************************************************************************


# run settings #=============================
test_mode <- FALSE
if(test_mode){
  warning("Running in test mode with small data set.")
}


# Load packages and set globals #===============================================
my_packs = c('tidyverse',
             'openxlsx',
             'RColorBrewer',
             'viridis',
             'ggrepel',
             'scales',
             'wildrtrax',
             'lubridate',
             'sf',
             'naturecounts',
             'terra',
             'factoextra',
             'exactextractr',
             'napops',
             'ebirdst',
             'INLA',
             'inlabru',
             'ggtext')

# Note INLA not on CRAN
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE, type = "binary")

lapply(my_packs, library, character.only = TRUE)
rm(my_packs)

# currently using main but might need development branch
#remotes::install_github("ABbiodiversity/wildRtrax@development")
# install.packages("naturecounts",
#                  repos = c(birdscanada = 'https://birdscanada.r-universe.dev',
#                            CRAN = 'https://cloud.r-project.org'))

# load custom functions
devtools::load_all(".")

wt_auth() # Need your wildtrax username

# use reproducible package to Cache long running code
options(reproducible.cachePath = "analysis/data/cache_data")

# ggplot theme
CustomTheme <- theme_set(theme_bw())
CustomTheme <- theme_update(legend.key = element_rect(colour = NA),
                            legend.key.height = unit(1.2, "line"),
                            panel.grid.major = element_line(colour = 'gray95'),
                            #panel.grid.minor = element_line(colour = 'transparent'),
                            panel.border = element_rect(linetype = "solid",
                                                        colour = "black",
                                                        linewidth = 1, fill = NA),
                            axis.line = element_line(colour = "black"),
                            strip.text = element_text(size = 12, colour = "black"),
                            strip.background = element_rect(colour = "black",
                                                            fill = "lightblue2",
                                                            linetype = "solid"),
                            axis.title.y = element_text(margin = margin(0,10,0,0)),
                            axis.title.x = element_text(margin = margin(10,0,0,0)),
                            panel.background = element_rect(fill = "white"))


# Species list #================================================================
# Reconcile species lists between Birds Canada and NatureCounts

BSC_species <- naturecounts::search_species_code() %>%
  # drops duplicate species_ids keeping the first.
  distinct(species_id, .keep_all = TRUE) %>%
  rename(BSC_spcd = BSCDATA, species_scientific_name = scientific_name) %>%
  dplyr::select(species_id, BSC_spcd,species_scientific_name,english_name,french_name) %>%
  unique()
  # mutate(index = 1:nrow(.))

WT_species <- wildrtrax::wt_get_species() %>%
  # assuming we only want birds
  filter(species_class == "AVES") %>%
  rename(WT_spcd = species_code) %>%
  dplyr::select(WT_spcd,species_scientific_name) %>%
  unique()

all_species <- full_join(BSC_species %>% select(-species_id),
                         WT_species,
                         by = join_by(species_scientific_name)) %>%
  relocate(BSC_spcd, WT_spcd) %>%
  # select(-index) %>%
  unique()

all_species <- subset(all_species, species_scientific_name %!in% c("NULL NULL"," "))

rm(WT_species)


# Data Prep #===================================================================
{
## WildTrax Download and Cleaning #==============================================
# Download WildTrax datasets within a study area boundary

# Polygon delineating study area boundary

if(!file.exists("analysis/data/raw_data/BCR_Terrestrial/BCR_Terrestrial_master.shp")){
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

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

Study_Area <- Study_Area %>%
  st_transform(st_crs(AEA_proj))

Study_Area_bound <- Study_Area %>% st_union()

write_sf(Study_Area_bound, "analysis/data/derived_data/INLA_data/study_area.gpkg")

# Identify projects that are part of the Study Area


# Use wildRtrax conversion from ARU to PC
# ARU_projects <- wt_get_download_summary(sensor_id = 'ARU') %>% subset(sensor == "ARU") %>% arrange(project)
PC_projects <- wt_get_download_summary(sensor_id = 'PC') %>% arrange(project)

# Subset projects for testing
if(test_mode){
  PC_projects <- PC_projects %>% filter(status == "Published - Public", tasks <1000) %>%
    filter(str_detect(project, "Ontario")) %>%
    slice(1:10)
}


# Download projects with 'PC' data (includes projects converted from ARU)
wt_dl_loc_year <- safely(wt_dl_loc_year, quiet = FALSE)
# Cache not really working as I would hope. Still does the download again...
PC_fulldat <- map(PC_projects$project_id, ~wt_dl_loc_year(.x, sens_id = "PC"))

# There was one error for a bats data set
PC_fulldat <- PC_fulldat %>% map("result") %>% bind_rows()

write.csv(PC_fulldat, file = "analysis/data/interim_data/WildTrax_PC_locations.csv", row.names = FALSE)


# Subset to data within Study Area Boundary


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

PC_summary  <- read.csv(file = "analysis/data/interim_data/keep/WildTrax_PC_summary.csv")


# Download species records and recording info from each Point Count survey


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
# preserve time format when writing to csv
write.csv(mutate(PC_counts, survey_date=format(survey_date, "%FT%H:%M:%S%z")),
                 file = "analysis/data/interim_data/WildTrax_PC_counts.csv",row.names = FALSE)

PC_counts  <- read.csv(file = "analysis/data/interim_data/keep/WildTrax_PC_counts.csv") %>%
  mutate(survey_date = lubridate::ymd_hms(survey_date))

# Convert WildTrax species codes to Birds Canada species codes when they disagree
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

# Dataframe to track information associated with each survey
# 1 row per survey
PC_surveys <- PC_counts %>%
  select(-detection_distance,-detection_time,-species_code,-species_common_name,-species_scientific_name,-individual_count,-detection_heard,-detection_seen,-detection_comments) %>%
  unique()

# Date of each survey

# Three that have no time are made NA by this so I am not converting it is already POSIX
# PC_surveys$survey_date <- lubridate::ymd_hms(PC_surveys$survey_date)

# Append total duration of each survey
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
                            c(survey_duration_method = "0-1min", Survey_Duration_Minutes = 1),
                            c(survey_duration_method = "0-4.5min", Survey_Duration_Minutes = 5),
                            c(survey_duration_method = "0-4.97min", Survey_Duration_Minutes = 5),
                            c(survey_duration_method = "0-2.38min", Survey_Duration_Minutes = 2),
                            c(survey_duration_method = "0-0.69min", Survey_Duration_Minutes = 0.5)) %>%
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

# Append maximum distance of each survey

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

# Remove records outside the study area boundary
PC_surveys <- PC_surveys %>%
  subset(!is.na(latitude) & !is.na(longitude)) %>%
  st_as_sf(coords=c("longitude","latitude"),crs=4326, remove = FALSE) %>%
  st_transform(crs = st_crs(Study_Area)) %>%
  st_set_agr("constant") %>%
  st_intersection(st_set_agr(Study_Area, "constant"))

# Remove resampled surveys
PC_surveys <- PC_surveys %>%
  filter(str_detect(project, "Resample", negate = TRUE))

# Create a matrix of counts for each species
PC_counts <- subset(PC_counts, survey_id %in% PC_surveys$survey_id)
PC_counts$survey_id <- factor(PC_counts$survey_id, levels = PC_surveys$survey_id)

# Total number of individuals detected per survey (sum of individual counts)
PC_counts <- PC_counts %>%
  group_by(species_code,survey_id) %>%
  summarize(total_count = sum(as.numeric(individual_count)), .groups = "drop") %>%
  as.data.frame() %>%
  pivot_wider(names_from = species_code,
              values_from = total_count,
              values_fill = 0,
              id_expand = TRUE,
              names_expand = TRUE) %>%
  arrange(survey_id)

# Same ordering as PC_surveys
stopifnot(mean(PC_counts$survey_id == PC_surveys$survey_id) == 1)


# Format names
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

# don't need to combine so use existing wide table
WT_matrix <- PC_counts %>% column_to_rownames("survey_id") %>% as.matrix()

# Save file wildTrax data
WT_dat <- list(WT_surveyinfo = PC_surveys,
               WT_matrix = WT_matrix)
saveRDS(WT_dat, file = "analysis/data/interim_data/WT_dat.rds")


# NatureCounts Download and Cleaning  #=========================================

# ______________________________________________________________________________
# NatureCounts data cannot be downloaded using the R package for ON breeding
# bird atlas 3. I was able to request access and recieved a download link
# ______________________________________________________________________________

###--------> Last downloaded/run on Feb 28, 2024
if(!file.exists("analysis/data/interim_data/NatureCounts_data_download.csv")){
  nc_dat_tbl <- nc_requests(username = "sarah.endicott@ec.gc.ca") %>%
    filter(str_detect(collection, "ONATLAS3"))

  nc_data <- nc_data_dl(request_id = nc_dat_tbl$request_id[1],
                        username = "sarah.endicott@ec.gc.ca")

  write.csv(nc_data, "analysis/data/interim_data/NatureCounts_data_download.csv",
            row.names = FALSE)
}
nc_data <- read.csv("analysis/data/interim_data/NatureCounts_data_download.csv")

nc_data2 <- read.csv("analysis/data/raw_data/OBBA2_naturecounts_download.csv")

nc_data <- bind_rows(nc_data, nc_data2)

nc_data_cleaned <- nc_data %>%
  # Remove empty columns
  select(where(~sum(!is.na(.x)) > 0)) %>%
  # rename(species_id = SpeciesCode) %>%
  left_join(BSC_species %>% select(species_id, BSC_spcd),
            by = join_by(species_id)) %>%
  rename(spcd = BSC_spcd)

missing_spcd <- nc_data_cleaned %>% filter(is.na(spcd)) %>%
  count(species_id)

nc_species <- meta_species_codes() %>%
  tidyr::spread("authority", "species_code") %>%
  dplyr::select(-"rank") %>%
  dplyr::distinct()

nc_tax <- meta_species_taxonomy()

# all the species match a concept, some are subspecies, some are / as in could be
# one of 2 species and some are sp.

# missing_spcd %>% left_join(nc_tax) %>% View()

# TODO add a look up by scientific name after dropping subspecies to fix those ones

# for now just using hard coded changes but they don't apply well to ON it seems
nc_data_cleaned %<>% mutate(spcd=case_when(species_id==32671 ~ "GHOW",
                                  species_id==40183 ~ "RNEP",
                                  species_id==46051 ~ "EUCD",
                                  species_id==40309 ~ "RTHA",
                                  species_id==10485 ~ "NOFL",
                                  species_id==40320 ~ "MERL",
                                  species_id==45168 ~ "DEJU",
                                  species_id==16563 ~ "YEWA",
                                  species_id==2551 ~ "GBHE",
                                  species_id==40667 ~ "BARS",
                                  species_id==16801 ~ "PAWA",
                                  species_id==40845 ~ "FOSP",
                                  species_id==40182 ~ "COME",
                                  species_id==41568 ~ "NOGO",
                                  species_id==0 ~ "NONE",
                                  .default = spcd)) %>%
  filter(!is.na(species_id)) %>% filter(species_id!=12285) %>%
  filter(!is.na(spcd))

# only using point count for now because I don't know how/if the others should
# be used
nc_data_pc <- nc_data_cleaned %>% filter(collection %>% str_detect("PC$")) %>%
  filter(AllSpeciesReported == "Yes") %>%
  # Remove empty columns or all the same
  select(where(~sum(!is.na(.x)) > 0 & n_distinct(.x) != 1)) %>%
  rename(surveyID = SamplingEventIdentifier) %>%
  mutate(DurationInMinutes = round(DurationInHours * 60),
         ObservationCount = replace_na(ObservationCount, 0),
         ObservationDate = lubridate::ymd(paste0(survey_year, "-", survey_month,
                                                 "-", survey_day))) %>%
  # remove if date and time is not present. Note this removes 13967 rows, all
  # from the Point Count 6 Interval protocol because they have no time data
  filter(!is.na(TimeCollected), !is.na(ObservationDate)) %>%
  # remove if long/lat are NA, removes 3631
  filter(!is.na(longitude)) %>%
  # remove rows where only catalog num is different
  select(-GlobalUniqueIdentifier, -CatalogNumber, -record_id, -any_of("X")) %>%
  distinct() %>%
  arrange(surveyID) %>%
  st_as_sf(.,coords = c("longitude","latitude"), crs = st_crs(4326), remove = FALSE) %>%
  st_transform(st_crs(Study_Area_bound)) %>%
  st_filter(Study_Area_bound)

# Construct a dataset that contains JUST survey information (not counts)
nc_pc_surveyinfo <- nc_data_pc %>% st_drop_geometry() %>%
  dplyr::select(surveyID, Locality, latitude, longitude, ObservationDate,
                TimeCollected, DurationInMinutes, ProtocolType) %>%
  distinct() %>%
  rename(sq_id = Locality) %>%
  st_as_sf(.,coords = c("longitude","latitude"), crs = st_crs(4326), remove = FALSE)

# Counts for each species on each survey are stored in a matrix
# Note: I checked and the result is the same as the for loop
nc_pc_matrix <- nc_data_pc %>% st_drop_geometry() %>%
  select(surveyID, spcd, ObservationCount) %>%
  pivot_wider(names_from = spcd, values_from = ObservationCount, values_fill = 0) %>%
  column_to_rownames("surveyID") %>%
  as.matrix()

# Remove columns (species) that were never observed to save storage space
nc_pc_matrix <- nc_pc_matrix[,-which(colSums(nc_pc_matrix,na.rm = T)==0)]

# Same ordering as survey info
stopifnot(mean(rownames(nc_pc_matrix) == nc_pc_surveyinfo$surveyID) == 1)

# Ensure column names are consistent WildTrax dataset
nc_pc_surveyinfo <- nc_pc_surveyinfo %>%

  rename(Survey_Type = ProtocolType,
         survey_ID = surveyID,
         Latitude = latitude,
         Longitude = longitude,
         Survey_Duration_Minutes = DurationInMinutes) %>%

  mutate(Max_Distance_Metres = Inf,
         Date_Time = ymd(ObservationDate) + hours(floor(TimeCollected)) + minutes(round(60*(TimeCollected-floor(TimeCollected)))),
         Data_Source = "NatureCounts",
         Project_Name = "OBBA3") %>%

  dplyr::select(Data_Source,
                Project_Name,
                Survey_Type,
                survey_ID,
                Latitude,
                Longitude,
                Date_Time,
                Survey_Duration_Minutes,
                Max_Distance_Metres)

# Combine NC and WT Data #======================================================

# ************************************************************
# ************************************************************
# CHECK FOR OVERLAP BETWEEN NATURECOUNTS AND WILDTRAX POINT COUNT DATASETS
#
# - assume ARU data in WildTrax data is 'authoritative' because it can be reviewed/updated more often
#
# - assume Human Point counts in NatureCounts are authoritative because they have been thoroughly QA/QC'd
# ************************************************************
# ************************************************************

# load WT data

WT_dat <- readRDS(file = "analysis/data/interim_data/keep/WT_dat.rds")
WT_surveyinfo <- WT_dat$WT_surveyinfo %>% st_transform(crs = AEA_proj)
WT_matrix <- WT_dat$WT_matrix

nc_pc_surveyinfo <- nc_pc_surveyinfo %>% st_transform(crs = AEA_proj)

# Identify NatureCounts point counts that are within 10 m of WildTrax locations
WT_buff <- st_buffer(WT_surveyinfo,10) %>% st_union()

nc_pc_surveyinfo$to_evaluate <- FALSE
nc_pc_surveyinfo$to_evaluate[which(as.matrix(st_intersects(nc_pc_surveyinfo,WT_buff)))] <- TRUE

distance_matrix <- st_distance(subset(nc_pc_surveyinfo, to_evaluate),WT_surveyinfo)

# Flags for which observations to remove
nc_pc_surveyinfo$Obs_Index <- 1:nrow(nc_pc_surveyinfo)
WT_surveyinfo$Obs_Index <- 1:nrow(WT_surveyinfo)
nc_pc_surveyinfo$to_remove <- FALSE
WT_surveyinfo$to_remove <- FALSE

# Check for and flag overlap
for (i in 1:nrow(subset(nc_pc_surveyinfo, to_evaluate))){

  obs_to_evaluate <- subset(nc_pc_surveyinfo, to_evaluate)[i,]
  dists <- distance_matrix[i,]
  # TODO: why is this 5m and above intersection is 10m?
  WT_within_5m <- WT_surveyinfo[which(dists <= units::set_units(5,"m")),]
  if(nrow(WT_within_5m) == 0){
    next
  }
  time_diff <- abs(difftime(obs_to_evaluate$Date_Time, WT_within_5m$Date_Time, units='mins')) %>% as.numeric()

  # If the survey was ARU-based, use WildTrax as authoritative
  if (min(time_diff)<=10 & obs_to_evaluate$Survey_Type != "Point_Count") nc_pc_surveyinfo$to_remove[which(nc_pc_surveyinfo$Obs_Index == obs_to_evaluate$Obs_Index)] <- TRUE

  # If the survey was Human-based, use NatureCounts as authoritative
  if (min(time_diff)<=10 & obs_to_evaluate$Survey_Type == "Point_Count") WT_surveyinfo$to_remove[which(WT_surveyinfo$survey_ID %in% WT_within_5m[which(time_diff<=10),]$survey_ID)] <- TRUE
}

# Remove duplicated records from NatureCounts dataset
PC_removed <- subset(nc_pc_surveyinfo,to_remove)
nc_pc_surveyinfo <- subset(nc_pc_surveyinfo,!to_remove)
nc_pc_matrix <- nc_pc_matrix[nc_pc_surveyinfo$Obs_Index,]
nc_pc_surveyinfo <- nc_pc_surveyinfo %>% dplyr::select(-Obs_Index,-to_evaluate, -to_remove)

# Remove duplicated records from WildTrax dataset
WT_removed <- subset(WT_surveyinfo,to_remove)
WT_surveyinfo <- subset(WT_surveyinfo,!to_remove)
WT_matrix <- WT_matrix[WT_surveyinfo$Obs_Index,]
WT_surveyinfo <- WT_surveyinfo %>% dplyr::select(-Obs_Index, -to_remove)


# COMBINE ALL DATASETS INTO A SINGLE FILE

WT_surveyinfo$survey_ID <- as.character(WT_surveyinfo$survey_ID)
nc_pc_surveyinfo$survey_ID <- as.character(nc_pc_surveyinfo$survey_ID)

all_surveys <- bind_rows(WT_surveyinfo %>% st_transform(crs = AEA_proj) %>%
                           mutate(Survey_Class = "WT"),
                         nc_pc_surveyinfo %>% st_transform(crs = AEA_proj) %>%
                           mutate(Survey_Class = "NC_PC"))

# Set time zone
# this forces this timezone which is not local for all of Ontario. Need to check
# tz(all_surveys$Date_Time) <- "Canada/Eastern"

# lubridate forces the same timezone for all the rows in the column so not sure this is helpful
all_surveys <- all_surveys %>%
  mutate(timez = lutz::tz_lookup_coords(Latitude, Longitude, method = "accurate")) %>%
  group_by(timez) %>%
  mutate(Date_Time = force_tz(Date_Time, unique(timez))) %>%
  ungroup()

# TODO double check time zones and sunrise to ensure they are correct

# Hours since sunrise
all_surveys$Sunrise <- suntools::sunriset(crds = st_transform(all_surveys,crs = 4326),
                                          dateTime = all_surveys$Date_Time,
                                          direction = c("sunrise"),
                                          POSIXct.out = TRUE)$time
all_surveys$Hours_Since_Sunrise <- with(all_surveys, difftime(Date_Time,Sunrise,units="hours") )

# might be worth checking some of these are correct
all_surveys$Hours_Since_Sunrise %>% range()

# Combine species count matrices

NC_species <- colnames(nc_pc_matrix)
WT_species <- colnames(WT_matrix)
#DO_species <- colnames(DO_matrix)

# NC_species only detected in NatureCounts point counts
NC_only <- NC_species[NC_species %!in% WT_species]

# NC_species only detected in WildTrax
WT_only <- WT_species[WT_species %!in% NC_species ]


subset(all_species,BSC_spcd %in% NC_only)
subset(all_species,BSC_spcd %in% WT_only)

# need to check why there are still mammals in WT species...

full_count_matrix <- bind_rows(as.data.frame(WT_matrix), as.data.frame(nc_pc_matrix)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>% as.matrix()

analysis_data <- list(all_surveys = all_surveys,
                      full_count_matrix = full_count_matrix)

saveRDS(analysis_data, file = "analysis/data/interim_data/analysis_data.rds")

# Plot Data Availability #======================================================

# ONBoundary <- ON_BCR %>% st_union()

all_surveys$Survey_Class <- factor(all_surveys$Survey_Class, levels = c("WT","NC_PC")) # "DO"
ggplot()+
  geom_sf(data = Study_Area, fill = "gray95", col = "transparent") +
  geom_sf(data = all_surveys, aes(col = year(Date_Time) <2010), size = 0.2)+
  # geom_sf(data = ONBoundary, fill = "transparent", col = "black") +
  # scale_color_manual(values=c("orangered","black","dodgerblue"),name = "Data Source")+
  facet_grid(.~Survey_Class)+
  ggtitle("Data availability")

# Tidy up data prep objects
rm('counts', 'distance_matrix', 'dists', 'i', 'method_definitions',
   'missing_spcd', 'nc_data', 'nc_data2','nc_data_cleaned', 'nc_data_pc', 'NC_only',
   'nc_pc_matrix', 'nc_pc_surveyinfo', 'nc_species', 'NC_species', 'nc_tax',
   'obs_to_evaluate', 'PC_counts', 'PC_fulldat', 'PC_projects', 'PC_removed',
   'PC_sf', 'PC_summary', 'PC_surveys', 'PID', 'project_name', 'species',
   'species_to_fix', 'time_diff', 'WT_buff', 'WT_dat', 'WT_matrix', 'WT_only',
   'WT_removed', 'WT_species', 'WT_surveyinfo', 'WT_within_5m')
}

# Select Analysis Data #========================================================

# Select data meeting criteria for inclusion (dates, time since sunrise, etc)
yr_day_start <- yday(ymd("2022-05-15"))
yr_day_end <- yday(ymd("2022-07-15"))
yr_start <- 2000
yr_end <- 2025
hr_ssr_start <- -2
hr_ssr_end <- 4
# survey duration min and max are different for point and stationary counts

analysis_data <- readRDS(file = "analysis/data/interim_data/analysis_data.rds")

all_surveys <- analysis_data$all_surveys %>% ungroup() %>%
  mutate(Obs_Index = 1:n())
full_count_matrix <- analysis_data$full_count_matrix
all_surveys$Survey_Type[all_surveys$Survey_Type == "Point Count"] <- "Point_Count"

# Select Point Counts / ARUs to use
PC_to_use <- subset(all_surveys,
                    Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM") &

                      Survey_Duration_Minutes > 1 &
                      Survey_Duration_Minutes <= 10 &

                      Hours_Since_Sunrise >= hr_ssr_start &
                      Hours_Since_Sunrise <= hr_ssr_end &

                      yday(Date_Time) >= yr_day_start &
                      yday(Date_Time) <= yr_day_end &

                      year(Date_Time) >= yr_start &
                      year(Date_Time) <= yr_end
)


# Select STATIONARY COUNT data to use
SC_to_use <- subset(all_surveys,

                    Survey_Type %in% c("Stationary Count") &

                      Hours_Since_Sunrise >= -2 &
                      Hours_Since_Sunrise <= 4 &

                      Survey_Duration_Minutes >= 1 &
                      Survey_Duration_Minutes <= 120 &

                      yday(Date_Time) >= yr_day_start &
                      yday(Date_Time) <= yr_day_end &

                      year(Date_Time) >= yr_start &
                      year(Date_Time) <= yr_end)

# Subset matrices
surveys_to_use <- c(PC_to_use$Obs_Index) # , SC_to_use$Obs_Index , LT_to_use$Obs_Index
all_surveys <- subset(all_surveys, Obs_Index %in% surveys_to_use)
full_count_matrix <- full_count_matrix[surveys_to_use,]
all_surveys$Obs_Index <- 1:nrow(all_surveys)

rm(PC_to_use, SC_to_use, hr_ssr_end, hr_ssr_start, yr_day_end, yr_day_start,
   yr_end, yr_start, surveys_to_use)

# Atlas Squares # 10 km x 10 km grid
ONSquares <- st_make_grid(
  Study_Area_bound,
  cellsize = units::set_units(10*10,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE)%>%
  st_as_sf() %>%
  st_intersection(Study_Area_bound) %>%
  na.omit() %>%
  mutate(sq_id = 1:nrow(.)) %>%
  rename(geometry = x)

all_surveys <- all_surveys %>%
  mutate(Obs_Index = 1:nrow(.)) %>% st_intersection(ONSquares)
full_count_matrix <- full_count_matrix[all_surveys$Obs_Index,]

all_surveys$Obs_Index <- 1:nrow(all_surveys)

# Prepare Covariates #===========================================================
{
## Download and crop #==========================================================

# Dataframe to store covariates
all_surveys_covariates <- all_surveys %>% dplyr::select(Obs_Index)

# 1 km x 1 km grid
# TODO check this transformation makes sense
ONGrid <- reproducible::Cache(st_make_grid(
  Study_Area_bound,
  cellsize = units::set_units(10*10,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE)%>%
  st_as_sf() %>%
  st_intersection(Study_Area_bound) %>%
  na.omit())

ONGrid$point_id <- 1:nrow(ONGrid)
ONGrid_centroid <- st_centroid(ONGrid)

# Add 20 km buffer so that covariates can extend slightly outside province if possible
SABoundary_buffer <- Study_Area_bound %>% st_buffer(20000)

# Download covariate layers
url_nalcms_2020_Can <- "http://www.cec.org/files/atlas_layers/1_terrestrial_ecosystems/1_01_0_land_cover_2020_30m/can_land_cover_2020_30m_tif.zip"
nalcms_2020_Can <- reproducible::prepInputs(
  url = url_nalcms_2020_Can,
  destinationPath = "analysis/data/raw_data/",
  targetFile = "analysis/data/raw_data/CAN_NALCMS_landcover_2020_30m.tif",
  alsoExtract = NA,
  fun = terra::rast,
  studyArea = SABoundary_buffer,
  writeTo = "analysis/data/derived_data/covariates/NALCMS_NorthON.tif",
  useCache = TRUE
  )

nalcms_2020_SA <- terra::rast("analysis/data/derived_data/covariates/NALCMS_NorthON.tif") %>%


url_adaptwest_norm <- "https://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6v73/normals/Normal_1961_1990_bioclim.zip"

clim_norm <- reproducible::prepInputs(
  url = url_adaptwest_norm,
  destinationPath = "analysis/data/raw_data/",
  targetFile = "Normal_1961_1990/Normal_1961_1990_bioclim/Normal_1961_1990_MAT.tif"
)

mat_SA <- terra::crop(clim_norm,
                      terra::vect(SABoundary_buffer %>% st_transform(st_crs(clim_norm))),
                      mask = TRUE, overwrite = TRUE,
                      filename = "analysis/data/derived_data/covariates/mat_NorthON.tif")

# Path to spatial covariates
covar_folder <- "analysis/data/derived_data/covariates/"

# TODO: get these variables?
# # Annual mean temperature
# AMT <- rast(paste0(covar_folder,"National/AnnualMeanTemperature/wc2.1_30s_bio_1.tif")) %>%
#   crop(st_transform(ONBoundary_buffer,crs(.))) %>%
#   project(st_as_sf(ONBoundary_buffer), res = 250)
#
# # Land cover of Canada 2020 (not reprojected... faster to reproject grid)
# lcc2020 <- rast(paste0(covar_folder,"National/LandCoverCanada2020/landcover-2020-classification.tif")) %>%
#   crop(st_transform(ONBoundary_buffer,crs(.)))
#
# # Stand canopy closure
# SCC <- rast(paste0(covar_folder,"National/NationalForestInventory/NFI_MODIS250m_2011_kNN_Structure_Stand_CrownClosure_v1.tif")) %>%
#   crop(st_transform(ONBoundary_buffer,crs(.))) %>%
#   project(st_as_sf(ONBoundary_buffer), res = 250)


# Extract survey covariates within 1 km of location #===============

# Create a stack of covariate rasters
covars <- covar_folder %>% list.files(full.names = TRUE) %>%
  set_names(list.files(covar_folder) %>% str_remove("_NorthON\\.tif")) %>%
  map(rast)

covars$NALCMS <- covars$NALCMS %>% resample(covars$mat, method = "near")

covars <- rast(covars)

terra::writeRaster(covars, "analysis/data/derived_data/INLA_data/covariate_stack.grd")

# Note Using terra very slow for fractional coverage

# Continuous covariates
all_surveys_1km <- all_surveys_covariates %>%
  st_buffer(1000) %>%
  mutate(#elevation_1km = exact_extract(elevation, ., 'mean'),
    AMT_1km = exact_extract(covars$mat, ., 'mean'))#,
    #SCC_1km = exact_extract(SCC, ., 'mean'))


# Proportion of each land cover class
prop_LCC_1km <- exact_extract(covars$NALCMS,all_surveys_1km,"frac") %>%
  suppressWarnings()
names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km <- prop_LCC_1km %>% dplyr::select(sort(names(.)))

# TODO: this is not very safe could easily miss one in a new area, would be
# better as look up table

# Fix land cover class names
prop_LCC_1km <- prop_LCC_1km %>%
  mutate(
    Needleleaf_forest_1km = LCC_1_1km + LCC_2_1km,
    Mixed_forest_1km = LCC_5_1km + LCC_6_1km,
    Grass_shrub_1km = LCC_8_1km + LCC_10_1km,
    Barren_lichen_1km = LCC_13_1km + LCC_16_1km,
    Crop_1km = LCC_15_1km,
    Urban_1km = LCC_17_1km,
    Wetland_1km = LCC_14_1km,
    Water_1km = LCC_18_1km) %>%
  dplyr::select(Needleleaf_forest_1km:Water_1km)

all_surveys_1km <- bind_cols(all_surveys_1km, prop_LCC_1km) %>%
  st_drop_geometry()

# Join with survey dataset top get back to points
all_surveys_covariates <- left_join(all_surveys_covariates, all_surveys_1km,
                                    by = join_by(Obs_Index))


# Extract grid covariates within 1 km of location #=============================

# TODO: consider whether it is really helpful to have this as vector or if it
# could be a raster?
# Continuous covariates
ONGrid_1km <- ONGrid_centroid %>%
  st_buffer(1000) %>%
  mutate(#elevation_1km = exact_extract(elevation, ., 'mean'),
    AMT_1km = exact_extract(covars$mat, ., 'mean'))#,
    # SCC_1km = exact_extract(SCC, ., 'mean'))

# TODO: avoid this repetition
# Proportion of each land cover class
prop_LCC_1km <- exact_extract(covars$NALCMS, ONGrid_1km,"frac") %>% suppressWarnings()
names(prop_LCC_1km) <- paste0(str_replace(names(prop_LCC_1km),"frac","LCC"),"_1km")
prop_LCC_1km[setdiff(paste0("LCC_",seq(1,18),"_1km"),names(prop_LCC_1km))] <- 0
prop_LCC_1km <- prop_LCC_1km %>% dplyr::select(sort(names(.)))

# Fix land cover class names
prop_LCC_1km <- prop_LCC_1km %>%
  mutate(
    Needleleaf_forest_1km = LCC_1_1km + LCC_2_1km,
    Mixed_forest_1km = LCC_5_1km + LCC_6_1km,
    Grass_shrub_1km = LCC_8_1km + LCC_10_1km,
    Barren_lichen_1km = LCC_13_1km + LCC_16_1km,
    Crop_1km = LCC_15_1km,
    Urban_1km = LCC_17_1km,
    Wetland_1km = LCC_14_1km,
    Water_1km = LCC_18_1km) %>%
  dplyr::select(Needleleaf_forest_1km:Water_1km)

ONGrid_1km <- bind_cols(ONGrid_1km, prop_LCC_1km)

# Join with ONGrid
ONGrid <- ONGrid %>% left_join(st_drop_geometry(ONGrid_1km), by = join_by(point_id))

# Covariate PCA #===============================================================

# Conduct Principal Components Analysis on covariates to identify axes of major
# variation in habitat

# Remove surveys with no covariate information
surveys_to_remove <- st_drop_geometry(all_surveys_covariates)
surveys_to_remove <- which(is.na(rowSums(surveys_to_remove)))
all_surveys_covariates <- all_surveys_covariates[-surveys_to_remove,]
all_surveys <- all_surveys[-surveys_to_remove,]
full_count_matrix <- full_count_matrix[-surveys_to_remove,]

covars_for_PCA <- all_surveys_covariates %>%
  st_drop_geometry() %>%
  dplyr::select(AMT_1km:Water_1km)

pca <- prcomp(covars_for_PCA, scale = TRUE)

# Interpretation of specific axes (e.g., axes 1 and 2)

summary(pca)   # Proportion variance explaind by axes
fviz_eig(pca)  # Scree plot (first 5 axes explain 85% of variation in habitat between sites)
pca            # Variable loadings

fviz_pca_var(pca,
             axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = viridis(10),
             repel = TRUE     # Avoid text overlapping
)

# Predict PCA values for each survey location and standardize (mean = 0, sd = 1)

all_surveys_PCA <- predict(pca, newdata = st_drop_geometry(all_surveys_covariates)[,names(pca$center)])
ONGrid_PCA <- predict(pca, newdata = st_drop_geometry(ONGrid)[,names(pca$center)])

# TODO: use across and scale
for (covar in colnames(all_surveys_PCA)){

  covar_mean <- mean(all_surveys_PCA[,covar],na.rm = TRUE)
  covar_sd <- sd(all_surveys_PCA[,covar],na.rm = TRUE)

  all_surveys_PCA[,covar] <- (as.data.frame(all_surveys_PCA)[,covar] - covar_mean)/covar_sd
  ONGrid_PCA[,covar] <- (as.data.frame(ONGrid_PCA)[,covar] - covar_mean)/covar_sd

}

all_surveys_covariates <- all_surveys_covariates %>%
  st_drop_geometry() %>%
  bind_cols(all_surveys_PCA)


all_surveys <- full_join(all_surveys,all_surveys_covariates)
ONGrid <- bind_cols(ONGrid,ONGrid_PCA)
}

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

# TODO: this is not really working right because sq_id was missing using survey
# id instead but that is not quite the same
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

  ebirdst_download(species_name, pattern = "range")

  path <- get_species_path(species_name)

  range <- load_ranges(path, resolution = "lr")

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

  pca = pca,

  # Contains covariates on ON-wide grid
  ONGrid = ONGrid,

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

