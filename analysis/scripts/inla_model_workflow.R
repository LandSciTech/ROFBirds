# ------------------------------------------------------------------------------
# ******************************************************************************
# Implementation of David Isle's Landscape Distribution Modeling ECCC workflow
# for Northern Ontario
# ******************************************************************************
# ------------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Notable changes of method:
# * Using the wildRtrax PC version of data collected by ARUs. This means that the
# conversion from ARU format to PC format is done on wildRtrax whereas David did
# this in his code.
# * Only using the point count data from NatureCounts, not clear what David did
# since non PC steps are commented out
# * Using location to determine timezone since it does not seem to be recorded
#----------------------------------------------------------------------------


# run settings #=============================
test_mode <- TRUE
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
             'wildRtrax',
             'lubridate',
             'sf',
             'naturecounts')

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

WT_species <- wildRtrax::wt_get_species() %>%
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

# wildTrax Download and Cleaning #==============================================
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


PC_fulldat <- reproducible::Cache(map(PC_projects$project_id, ~wt_dl_loc_year(.x, sens_id = "PC")) %>%
  bind_rows())

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

PC_summary  <- read.csv(file = "analysis/data/interim_data/WildTrax_PC_summary.csv")


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
PC_counts  <- read.csv(file = "analysis/data/interim_data/WildTrax_PC_counts.csv") %>%
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

# Create a matrix of counts for each species

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
stopifnot(mean(PC_counts$survey_id == PC_surveys$survey_id) == 1)


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
  filter(!is.na(species_id)) %>% filter(species_id!=12285)

# only using point count for now because I don't know how/if the others should
# be used
nc_data_pc <- nc_data_cleaned %>% filter(collection == "ONATLAS3PC") %>%
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
  # remove rows where only catalog num is different
  select(-GlobalUniqueIdentifier, -CatalogNumber, -record_id, -any_of("X")) %>%
  distinct() %>%
  arrange(surveyID)

# Construct a dataset that contains JUST survey information (not counts)
nc_pc_surveyinfo <- nc_data_pc %>%
  dplyr::select(surveyID, Locality, latitude, longitude, ObservationDate,
                TimeCollected, DurationInMinutes, ProtocolType) %>%
  distinct() %>%
  rename(sq_id = Locality) %>%
  st_as_sf(.,coords = c("longitude","latitude"), crs = st_crs(4326), remove = FALSE)

# Counts for each species on each survey are stored in a matrix
# Note: I checked and the result is the same as the for loop
nc_pc_matrix <- nc_data_pc %>% select(surveyID, spcd, ObservationCount) %>%
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

# combine NC and WT data #======================================================

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
AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

WT_dat <- readRDS(file = "analysis/data/interim_data/WT_dat.rds")
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
tz(all_surveys$Date_Time) <- "Canada/Eastern"

# lubridate forces the same timezone for all the rows in the column so not sure this is helpful
all_surveys <- all_surveys %>%
  mutate(timez = lutz::tz_lookup_coords(Latitude, Longitude, method = "accurate")) %>%
  group_by(timez) %>%
  mutate(Date_Time = force_tz(Date_Time, unique(timez)))

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

Study_Area <- Study_Area %>%
  st_transform(st_crs(AEA_proj))

# ONBoundary <- ON_BCR %>% st_union()

all_surveys$Survey_Class <- factor(all_surveys$Survey_Class, levels = c("WT","NC_PC")) # "DO"
ggplot()+
  geom_sf(data = Study_Area, fill = "gray95", col = "transparent") +
  geom_sf(data = subset(all_surveys, year(Date_Time) > 2010),aes(col = Survey_Class), size = 0.2)+
  # geom_sf(data = ONBoundary, fill = "transparent", col = "black") +
  scale_color_manual(values=c("orangered","black","dodgerblue"),name = "Data Source")+
  facet_grid(.~Survey_Class)+
  ggtitle("Data availability")
