#To use the "recurring" package on Github to generate offsets, I need a single data frame as 
#input, in wide format with a column for spp, dt (YYYY-MM-DD format,
#zero-padded), tm (hh:mm format, zero-padded), lon, lat, dur, dis
#
#If I extract existing estimates based on most recent version of QPAD (2016) and use functions
#in "recurring" package:
#I'll start out by taking separate data frames from each project in my Wood Thrush study in long
#format, each containing one row per species/distance interval used/time interval used/site visit
#and summarizing abundance of each species per visit. "dur" and "dis" in each long data frame are
#categories that will be converted to numeric variables (maximum distance and length of point count)
#the column for Wood Thrush will become the "spp" column while all other species will be discarded.
#Output for each data source will be the following: PKEY, PCODE, SS, spp, dt, tm, lon, lat, DURMETH,
#DISMETH (DURMETH and DISMETH will be used to reclassify/recalculate dis and dur).

#PKEY (and PKEY_V6) identify unique surveys (stationXtimeXdate) in
#the Boreal Avian Modelling Project, along with SS and location_name_V6 
#(unique station ID). PKEY (and PKEY_V6) and SS (and location_name_V6) 
#are created within the BBS data because it is anticipated that 
#at some point, station-level BBS data will be combined with non-BBS data 
#from the Boreal Avian Modeling Project in a single analysis.
memory.limit(size=56000)

library(dplyr)
library(tidyr)
library(lubridate)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

prov<-c("ON")
#Remember, if you are filtering data to additional/other provinces (coordinates unavailable for
#U.S. states), "prov<-c("ON")" will need to be edited.

for (i in prov){
  load(paste0("0_data/processed/6_BBSdata_.",i,".indiv.stoptimesadded.Nov2021.RData"))
  BBS<-m5
  names(BBS)
  BBS<-BBS[BBS$RunType==1,]#Conditions good enough for surveys
  BBS<-BBS[BBS$Stopnumber<46,]#Drop the last 5 stops in case they were omitted from surveys
  BBS$dataset_code<-"BBSON"#formerly PCODE
  BBS$dataset_name<-"Breeding Bird Survey - Ontario"#formerly PCODE
  BBS$PKEY_V6<-paste0(BBS$SS,":",BBS$Year)#formerly PKEY
  BBS$location_name_V6<-BBS$SS#formerly SS
  BBS$survey_year<-BBS$Year
  BBS$dt<-BBS$date
  BBS$Hour<-as.character(ifelse(BBS$visit.hour<10,paste0("0",BBS$visit.hour),BBS$visit.hour))
  BBS$Minute<-as.character(ifelse(BBS$visit.minute<10,paste0("0",BBS$visit.minute),BBS$visit.minute))
  BBS$tm<-paste0(BBS$Hour,":",BBS$Minute)
  BBS$DISMETH<-"D"
  BBS$DURMETH<-"D"
  BBS$dis<-"Inf"
  BBS$dur<-3
  
  BBS.spp.sum<-BBS %>% 
    setNames(make.names(names(.), unique = TRUE)) %>%
    group_by(PKEY_V6,Species_ID,DISMETH,DURMETH,dis,dur,dataset_code,location_name_V6,lat,lon,dt) %>% #ProtocolCode,YearCollected,MonthCollected,DayCollected,JulianDay,
    arrange(tm) %>%
    summarise(Abund= sum(Abund))
  nrow(BBS.spp.sum)#6020839
  #
  
  BBS.spp.wide<-BBS.spp.sum %>%
    pivot_wider(names_from = Species_ID, values_from = Abund)
  nrow(BBS.spp.wide)#103389 obs. more plausible given 45 visits per route X #routes in Ontario X #years/route
  str(BBS.spp.wide)
  write.csv(BBS.spp.wide, file=paste0("0_data/processed/7_BBS.",i,".specieswide.summ.csv"))
  
  BBS.spp.sumB<-BBS %>% 
    setNames(make.names(names(.), unique = TRUE)) %>%
    group_by(PKEY_V6,DISMETH,DURMETH,dis,dur,dataset_code,location_name_V6,lat,lon,dt) %>% #ProtocolCode,YearCollected,MonthCollected,DayCollected,JulianDay,Species_ID
    arrange(tm) %>%
    summarise(tm= min(tm)) 
  nrow(BBS.spp.sumB)#103389 obs. same number as wide-formatted tibble with species counts
  
  BBS.spp.sumC<-merge(BBS.spp.sumB, BBS.spp.wide, by=c("PKEY_V6","dataset_code","location_name_V6","lat","lon","dt","dis","dur","DISMETH","DURMETH"))
  nrow(BBS.spp.sumC)#103389 obs.
  names(BBS.spp.sumC)
  
  spp_count<-BBS.spp.sumC%>%
    mutate_at(vars(AMCR:WEKI), ~replace(., is.na(.), 0))
  
  write.csv(spp_count, file=paste0("0_data/processed/8_BBS.",i,".readygetoffsets.csv"))
  save(spp_count, file=paste0("0_data/processed/8_BBS.",i,".readygetoffsets.RData"))
  
}
  
