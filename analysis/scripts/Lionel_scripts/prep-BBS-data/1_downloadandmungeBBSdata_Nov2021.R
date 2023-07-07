#The purpose of this script is to take raw North American
#Breeding Bird survey data and process it to a format that 
#is ready to use in models like generalized linear models
#or mixed effects models. In most analyses of BBS data, the
#replicate of the analysis is either individual BBS routes
#or strata containing multiple routes. In contrast, this 
#R project is set up to prepare individual BBS stops in a 
#given year as the unit of analysis. The reason for analyzing 
#individual stops rather than whole routes is if there are
#analyses where you want to combine BBS survey data with non-
#BBS survey data, as in studies conducted using methods developed
#by the Boreal Avian Modeling Project. In the final format, there will
#be one row per route stop per year and one column per species
#showing raw total abundance counted per stop per year. There 
#will also be columns containing detection offsets (Solymos et
#al. 2013) with one offset column per species. Detection offsets
#are estimated based on availability of birds for detection (singing
#rate, affected by time of day and season), detection probability
#given availability (probability of detection affected by surrounding
#vegetation and depending on survey methods, distance of individual
#birds), and survey methods.
#of processing.

#First the raw data must be downloaded. There are two data sources that are needed:
#the BBS data itself and the coordinates for individual stops as opposed to the coordinates
#of the first stop on each route (included in the BBS release).

#1. Breeding Bird Survey Data

#The R project folder "0_data/raw/" contains the BBS 2020 release obtained
#from the website https://www.pwrc.usgs.gov/bbs/rawdata/
#using the Science Base option at 
#https://www.pwrc.usgs.gov/bbs/rawdata/Choose-Method.cfm

#Doing so allows you to potentially alter this script to download data then filter
#it to a different jurisdiction than Ontario, although latitude and longitude for 
#individual stops is currently only available from the Canadian Wildlife Service
#(U.S. data have the coordinates for the first stop within BBS routes).

#Choose the most recent release at 
#https://www.sciencebase.gov/catalog/item/52b1dfa8e4b0d9b325230cd9

#The release contains folders or files for:
#North American Breeding Bird Survey Dataset
#A completeness PDF report
#A migrant nonbreeder zip file
#A zip file of the routes
#A PDF description of run type (survey method differences among routes)
#A zip file of state data. Get the values for "StateNum" for the jurisdictions you want to
#filter to from this folder (e.g., Ontario: StateNum=68).
#A zip file of vehicle data
#A zip file of weather data for the routes
#A text file for the species list and another for the singing-species.
#A zip file of raw count data for every stop
#along every route in each route's survey years
#(50-StopData)

#All of these files and folders in the BBS release should be copied to "0_data/raw/"
#and uncompressed.

#2. BBS Routes&CurrentStops&Discontinued.gdb

#This is a geodatabase of BBS point count coordinates for individual stops (ALL stops,
#not just stop 1 on each route). These coordinates are currently available for 
#BBS routes in Canada but not yet in the United States. Assuming that you want to 
#update the stops to get new routes, you would get this geodatabase from:

#Veronica Isabelle Aponte
#National Coordinator of the North American Breeding Bird Survey
#Canadian Wildlife Service / Environment and Climate Change Canada / Government of Canada
#ec.RON-BBS.ec@canada.ca
#Coordonnatrice nationale du Relevé des oiseaux nicheurs de l’Amérique du Nord
#Service canadien de la faune / Environnement et Changement climatique Canada / Gouvernement du Canada
#ec.RON-BBS.ec@canada.ca

#In ArcGIS, you would export two attribute tables from the geodatabase as CSV files:
#(CurrentCanadaBBSStops, DiscontinuedCanadaBBSStops) and store them in "0_data/raw/".


memory.limit(size=56000)

#This script for downloading and munging BBS data from 1991 to 2019 (the 2020 release)
#doesn't use bbsAssistant.

#This script uses geographic coordinates for individual stations provided (currently
#just by the Canadian Wildlife Service but hopefully eventually by the USGS as 
#well). 

#For now the alternative is to use the station coordinates for stations that are
#already in the BAM database.

#By looking at files for individual provinces and states I can get province number
#and state number.
#I would then use the country and state numbers and years to filter
#the 50-stop data
library(dplyr)
library(forcats)
library(lubridate)
library(tidyr)
library(tidyverse)
library(suncalc)

Fifty1<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty1/fifty1.csv",header=TRUE)
Fifty2<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty2/fifty2.csv",header=TRUE)
Fifty3<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty3/fifty3.csv",header=TRUE)
Fifty4<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty4/fifty4.csv",header=TRUE)
Fifty5<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty5/fifty5.csv",header=TRUE)
Fifty6<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty6/fifty6.csv",header=TRUE)
Fifty7<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty7/fifty7.csv",header=TRUE)
Fifty8<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty8/fifty8.csv",header=TRUE)
Fifty9<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty9/fifty9.csv",header=TRUE)
Fifty10<-read.csv("0_data/raw/50-StopData/50-StopData/1997ToPresent_SurveyWide/Fifty10/fifty10.csv",header=TRUE)

#Combine the files
All<-rbind(Fifty1,Fifty2,Fifty3,Fifty4,Fifty5,Fifty6,Fifty7,Fifty8,Fifty9,Fifty10)
nrow(All)#3786794
range(All$Year)#1967 2019
#Filter the data to years 1991 or more recent
All.recent<-All[All$Year>1990,]
nrow(All.recent)#3785777
range(All.recent$Year)#1991 2019

All.recent$StateNum.f<-as.factor(All.recent$StateNum)
#New Brunswick=Countrynum 124 StateNum 56
#Nova Scotia=Countrynum 124 StateNum 65
#Ontario=Countrynum 124 StateNum 68
#Prince Edward Island=Countrynum 124 StateNum 75
#Quebec=Countrynum 124 StateNum 76
#Connecticut=Countrynum 840 StateNum 18 (high priority)
#Delaware=Countrynum 840 StateNum 21 (high priority)
#Illinois=Countrynum 840 StateNum 34
#Indiana=Countrynum 840 StateNum 35
#Maine=Countrynum 840 StateNum 44 
#Maryland=Countrynum 840 StateNum 46 (high priority?)
#Massachussetts=Countrynum 840 StateNum 47 (high priority)
#Michigan=Countrynum 840 StateNum 49 
#Minnesota=Countrynum 840 StateNum 50  
#New Hampshire=Countrynum 840 StateNum 58 (high priority)
#New Jersey=Countrynum 840 StateNum 59 (high priority)
#New York=Countrynum 840 StateNum 61 (high priority)
#Ohio=Countrynum 840 StateNum 66 (high priority)
#Pennsylvania=Countrynum 840 StateNum 72 (high priority)
#Rhode Island=Countrynum 840 StateNum 77 (high priority)
#Vermont=Countrynum 840 StateNum 87 (high priority)
#West Virginia=Countrynum 840 StateNum 90 (high priority)
#Wisconsin=Countrynum 840 StateNum 91 
All.recent.studyarea<-All.recent[All.recent$StateNum.f %in% c("68"),]
#currently, only the "state number" for Ontario is included, but you can put whatever numbers
#you want if you want to filter to additional or other provinces and states.

nrow(All.recent.studyarea)#1398266
write.csv(All.recent.studyarea, file="0_data/processed/1_stopsfilteredbyprovstate.csv")
#some times are improperly input. Currently start time is for the first stop and not 
#for individual stops. Since a BBS route can take all morning or all day, we need times
#of individual route stop surveys, since survey time affects availability of singing birds
#for detection.

#The 50-stop data identifies species using AOU numbers, which take up less memory
#than a whole name.
#AOU contains the species numbers. There is a text file for species list that could be
#joined to the data. 
SpeciesList<-read.csv("0_data/raw/SpeciesList.csv", header=TRUE)
str(SpeciesList)#AOU is the variable to use for merging 
range(SpeciesList$AOU)#10 52021

#There is a "singingspecies" CSV file we want to join to SpeciesList first, before
#joining the species list to the 50-stop data, because "singingspecies" contains alpha
#code names for each species. These alpha codes are shorter than species names so
#they're nicer for analyses, but they are easier to remember for each species than AOU
#numbers.
singingspecies<-read.csv("0_data/raw/singing-species.csv", header=TRUE)
SpeciesListB<-merge(SpeciesList, singingspecies, by=c("English_Common_Name"))

#Detection covariates
#The weather CSV file contains information about the suitability of conditions along individual
#BBS routes for detecting birds. These aren't the actual covariates used for calculating detection
#offsets, but they may be used for selectively removing survey-route years with suboptimal
#conditions for counting birds, so those survey route-years don't influence analyses.
weather<-read.csv("0_data/raw/Weather/weather.csv",header=TRUE)
str(weather)
weather$StateNum.f<-as.factor(weather$StateNum)
weather.studyarea<-weather[weather$StateNum.f %in% c("68"),]
#Remember, if you are filtering data to additional/other provinces and states,
#"%in% c("68")" will need to be edited.

nrow(weather.studyarea)#46036
write.csv(weather.studyarea, file="0_data/processed/weather.StudyArea.csv")
#weather data filtered to just the data for Ontario.

weather.studyarea$JURIS<-as.factor(weather.studyarea$StateNum)
weather.studyareaJ<-weather.studyarea %>%  dplyr::mutate(JURIS = fct_recode(JURIS,
                                                                            ON = "68"
))
str(weather.studyareaJ)

#BBS Station Coordinates for Canada (from Canadian Wildlife Service)
BBSstop.Canada.current<-read.csv("0_data/raw/CurrentCanadaBBSStops.csv", header=TRUE)
BBSstop.Canada.discontin<-read.csv("0_data/raw/DiscontinuedCanadaBBSStops.csv", header=TRUE)
BBSstop.Canada.discontin$Year<-as.integer(BBSstop.Canada.discontin$Year)
BBSstop.Canada.all<-bind_rows(BBSstop.Canada.current, BBSstop.Canada.discontin)
names(BBSstop.Canada.all)
# [1] "OBJECTID"   "route"      "POINT_X"    "POINT_Y"    "Active"    
#[6] "Province"   "Province_R" "Stop"       "Nbr_FullNa" "ProvRoute_"
#[11] "Year"
BBSstop.Canada.all$Province
BBSstop.Canada.all$JURIS<-as.factor(as.character(BBSstop.Canada.all$Province))
BBSstop.Canada.all.studyareaJ<-BBSstop.Canada.all %>%  dplyr::mutate(JURIS = fct_recode(JURIS,
                                                                                  ON = "68"))
#Remember, if you are filtering data to additional/other provinces (coordinates unavailable for
#U.S. states), "%in% c("68")" will need to be edited.

levels(BBSstop.Canada.all.studyareaJ$JURIS)
names(BBSstop.Canada.all.studyareaJ)
# [1] "OBJECTID"   "route"      "POINT_X"    "POINT_Y"    "Active"    
#[6] "Province"   "Province_R" "Stop"       "Nbr_FullNa" "ProvRoute_"
#[11] "Year"

#Create the stop ID "SS" ("SS" has been used as a point count
#station ID in data collated by the Boreal Avian Modeling Project,
#to which the BBS station-level data may be ultimately combined)
BBSstop.Canada.all.studyareaJ$SS<-paste0("BBS",
                                         BBSstop.Canada.all.studyareaJ$JURIS,
                                         ":",
                                         BBSstop.Canada.all.studyareaJ$route,
                                         ":","Stop",
                                         BBSstop.Canada.all.studyareaJ$Stop)

BBSstop.Canada.all.studyareaJ$Route<-BBSstop.Canada.all.studyareaJ$route
BBSstop.Canada.all.studyareaJ$Stop<-paste0("Stop",BBSstop.Canada.all.studyareaJ$Stop)

#This is the raw 50-stop data that was previously filtered to just Ontario
All.recent.studyarea<-read.csv("0_data/processed/1_stopsfilteredbyprovstate.csv", header=TRUE)
#coordinates for BBS stops in Ontario alone
All.recent.studyarea$JURIS<-as.factor(All.recent.studyarea$StateNum)

names(All.recent.studyarea)
# [1] "X"           "RouteDataID" "CountryNum"  "StateNum"    "Route"      
# [6] "RPID"        "Year"        "AOU"         "Stop1"       "Stop2"      
# [11] "Stop3"       "Stop4"       "Stop5"       "Stop6"       "Stop7"      
# [16] "Stop8"       "Stop9"       "Stop10"      "Stop11"      "Stop12"     
# [21] "Stop13"      "Stop14"      "Stop15"      "Stop16"      "Stop17"     
# [26] "Stop18"      "Stop19"      "Stop20"      "Stop21"      "Stop22"     
# [31] "Stop23"      "Stop24"      "Stop25"      "Stop26"      "Stop27"     
# [36] "Stop28"      "Stop29"      "Stop30"      "Stop31"      "Stop32"     
# [41] "Stop33"      "Stop34"      "Stop35"      "Stop36"      "Stop37"     
# [46] "Stop38"      "Stop39"      "Stop40"      "Stop41"      "Stop42"     
# [51] "Stop43"      "Stop44"      "Stop45"      "Stop46"      "Stop47"     
# [56] "Stop48"      "Stop49"      "Stop50"      "StateNum.f"  "JURIS"

All.recent.studyareaJ<-All.recent.studyarea %>%  dplyr::mutate(JURIS = fct_recode(JURIS, 
                                     ON = "68"
))

All.recent.studyareaJ$SpecificRoute<-paste0("BBS",All.recent.studyareaJ$JURIS,":",All.recent.studyareaJ$Route)
All.recent.studyareaJ<-All.recent.studyareaJ[All.recent.studyareaJ$RPID=="101",]
#Limits BBS survey routes to those that were conducted in adequate conditions for counting birds.

levels(as.factor(All.recent.studyareaJ$RPID))
nrow(All.recent.studyareaJ)#149388
nlevels(as.factor(All.recent.studyareaJ$SpecificRoute))#216
levels(as.factor(All.recent.studyareaJ$SpecificRoute))


All.recent.studyareaJ.Canada<-All.recent.studyareaJ[All.recent.studyareaJ$JURIS %in% c("ON"),]
nrow(All.recent.studyareaJ.Canada)#149388
names(All.recent.studyareaJ.Canada)


prov<-c("ON")
#Remember, if you are filtering data to additional/other provinces (coordinates unavailable for
#U.S. states), "prov<-c("68")" will need to be edited.

for (i in prov){
  All.recent.studyareaJ.summ<-All.recent.studyareaJ[All.recent.studyareaJ$JURIS==i,]
  #Use gather function to convert wide format (1 column
  #for each stop) to long format (1 column for count of a species at each stop plus
  #one column indicating stop ID, with separate rows for each species).
  nrow(All.recent.studyareaJ.summ)#e.g. 149388 for Ontario
  Studyarea.long<-All.recent.studyareaJ.summ%>%
    gather(Stop1:Stop50, key="Stop", value="Abund")
  nrow(Studyarea.long)#e.g.7469400 for Ontario; 69913300 if we use a file with all prov/state in study area
  range(Studyarea.long$AOU)#e.g.20 7660 for Ontario; 10 22860 for all jurisdictions
  save(Studyarea.long, file="0_data/processed/2_gatherfilterstopsintolong.RData")
  
  
  Studyarea.long.B<-merge(Studyarea.long, SpeciesListB, by=c("AOU"), all.x=TRUE)
  str(Studyarea.long.B)
  
  #get rid of columns we don't need
  Studyarea.long.B$Species<-NULL
  Studyarea.long.B$French_Common_Name<-NULL
  Studyarea.long.B$Spanish_Common_Name<-NULL
  Studyarea.long.B$Genus<-NULL
  Studyarea.long.B$ORDER<-NULL
  Studyarea.long.B$Family<-NULL
  Studyarea.long.B$Seq<-NULL
  save(Studyarea.long.B, file="0_data/processed/3_gatheredfilteredstopsplusspeciesname.RData")
  
  #Create same stop ID as in other files prior to merging those files.
  Studyarea.long.B$SS<-paste0("BBS",Studyarea.long.B$JURIS,":",Studyarea.long.B$Route,":",Studyarea.long.B$Stop)
  print(paste0("Species English names added to AOU codes for species counts on points in ",i))
  
  m1<-merge(Studyarea.long.B,BBSstop.Canada.all.studyareaJ, by=c("SS","JURIS","Route","Stop"))
  str(m1)
  print(paste0("Counts and individual stop coordinates merged for ",i))
  save(m1, file="0_data/processed/4_addedindividualstopcoordinates.RData")
  
  m2<-m1[m1$JURIS %in% c("ON"),]
  m2$Year<-m2$Year.x
  m2$Year.x<-NULL
  m2$Year.y<-NULL
  
  m3<-merge(m2, weather.studyareaJ, by=c("JURIS","Route","Year","RouteDataID","CountryNum","StateNum","StateNum.f","RPID"))
  print(paste0("Counts, individual stop coordinates, and route-level detection covariates stored for ",i))
  #str(m3)
  n_last <- 2                                # Specify number of characters to extract
  m3$StartTime.Hour<-substr(m3$StartTime, 1, nchar(m3$StartTime) - n_last) # Extract last two characters
  m3$StartTime.Hour<-ifelse(as.numeric(m3$StartTime.Hour)<10,paste0("0",m3$StartTime.Hour),m3$StartTime.Hour)
  m3$StartTime.Min<-substr(m3$StartTime, nchar(m3$StartTime) - n_last + 1, nchar(m3$StartTime)) # Extract last two characters
  #m3$StartTime.Min<-ifelse(as.numeric(m3$StartTime.Min)<10,paste0("0",m3$StartTime.Min),m3$StartTime.Min)
  m3$EndTime.Hour<-substr(m3$EndTime, 1, nchar(m3$EndTime) - n_last) # Extract last two characters
  m3$EndTime.Hour<-ifelse(as.numeric(m3$EndTime.Hour)<10,paste0("0",m3$EndTime.Hour),m3$EndTime.Hour)
  m3$EndTime.Min<-substr(m3$EndTime, nchar(m3$EndTime) - n_last + 1, nchar(m3$EndTime)) # Extract last two characters
  #m3$EndTime.Min<-ifelse(as.numeric(m3$EndTime.Min)<10,paste0("0",m3$EndTime.Min),m3$EndTime.Min)
  
  m3$Month<-ifelse(m3$Month<10,paste0("0",m3$Month),m3$Month)
  m3$Day<-ifelse(m3$Day<10,paste0("0",m3$Day),m3$Day)
  m3$date<-as.Date(paste0(m3$Year,"-",m3$Month,"-",m3$Day))
  
  m3$start.timeL<-paste0(m3$date," ",paste0(m3$StartTime.Hour,":",m3$StartTime.Min,":00"))
  m3$starttime.posix<-as.POSIXct(m3$start.timeL, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
  m3$end.timeL<-paste0(m3$date," ",paste0(m3$EndTime.Hour,":",m3$EndTime.Min,":00"))
  m3$endtime.posix<-as.POSIXct(m3$end.timeL, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
  
  m3$fullrec.interval <- m3$starttime.posix %--% m3$endtime.posix
  m3$fullrec.interval
  
  m3$fullrec.duration <- as.duration(m3$fullrec.interval)
  m3$fullrec.duration 
  m3$fullrec.seconds<-m3$fullrec.duration@.Data  
  
  #2. For each BBS route in each year (routeXyear combination) create a counter that
  #gives us the number of stops per visit. 
  m3$routeyear<-paste0(m3$JURIS,"-",m3$Route,"-",m3$Year)
  
  #3. Store the number of stops per routeXyear combination. All observations from the 
  #same stop should have the same counter value
  m3$Stopnumber <- as.numeric(gsub("Stop", "", m3$Stop))
  m4 <- m3 %>% 
    group_by(routeyear) %>% 
    arrange(routeyear,Stopnumber) %>%
    mutate(counter=match(Stopnumber, unique(Stopnumber))) %>%
    mutate(maxnumstops=max(counter))
  
  save(m4, file="0_data/processed/5_addedweatherstartandstoptimes.RData")
  #write.csv(m4, file="m4.csv")
  
  #4. Divide Time Difference by number of Stops to get Time Increment
  m4$timeincrement<-m4$fullrec.seconds/m4$maxnumstops
  
  #5. Visit Time = Start Time + Time Increment*(Number of Stops Counter)
  m4$visit.time<-m4$starttime.posix+(m4$counter-1)*m4$timeincrement
  
  #6. Hour/Minute/Second = Based on Visit Time
  m4$visit.hour<-hour(m4$visit.time)
  m4$visit.minute<-minute(m4$visit.time)
  
  #7. TSSR based on Hour/Minute/Second
  data<-data.frame(date=m4$date, lat=m4$POINT_Y, lon=m4$POINT_X)
  time.sunrise<-getSunlightTimes(data = data,
                                 keep = c("sunrise", "sunriseEnd", "dawn"))#, tz = "UTC"
  time.sunrise$date<-NULL
  m5<-cbind(m4,time.sunrise)
  rm(m1,m2,m3,m4)
  gc()
  
  save(m5, file=paste0("0_data/processed/6_BBSdata_.",i,".indiv.stoptimesadded.Nov2021.RData"))
  print(paste0("R Data file saved for ",i," BBS region."))
  rm(m5)
  gc()
}
