
#' Extract data from GEE
#'
#' Extract data from GEE for buffered point locations for a set of covariates
#' described in a table
#'
#' @param meth_path path to table of extraction methods and sources. Will be
#'   updated with download progress. Must match BAM National Models format.
#'   https://docs.google.com/spreadsheets/d/1XATuq8BOYC2KkbJObaturaf4NFD2ufvn/edit?usp=sharing&ouid=104837701987164094932&rtpof=true&sd=true
#' @param loc.yr sf object with a unique row id and a year column identifying
#'   each unique site and year
#' @param out_save directory where outputs should be saved
#' @param projection.trans EPSG code for projection of output
#' @param local.buffer Buffer for smaller radius TODO: improve so this is only
#'   set in spreadsheet
#' @param neighbour.buffer Buffer for larger radius
#' @param do.test.run Should a test run of downloading 20 points from GEE be
#'   attempted first?
#' @param year_get Default is NA which means match the year column. If "max"
#'   then only the most recent year is extracted
#'
#' @return Files are saved to the out_save directory. And the spreadsheet at
#'   meth_path is updated with Complete = 1 and Running = 0 for variables that
#'   were downloaded successfully
#' @export
#'
#' @examples
gee_points_extract <- function(
    meth_path, loc.yr, out_save,
    projection.trans = 5072, local.buffer = 200,
    neighbour.buffer = 2000, do.test.run = TRUE, year_get = NA){

  meth <- read.csv(meth_path,
                   na.strings=c("NA",""))

  meth <- filter(meth, Complete == 0)

  if(!dir.exists(out_save)) dir.create(out_save)

  loc.n <- loc.yr %>%
    dplyr::mutate(year=as.integer(year)) %>%
    st_transform(crs=projection.trans)

  loc.buff <- st_buffer(loc.n, local.buffer)
  loc.buff2 <- st_buffer(loc.n, neighbour.buffer)

  if(do.test.run){
    loc.buff.small<-loc.buff %>% slice(1:20)
    gee_data <-ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")
    # Not working because SpatialResolution is a character
    # gee_data1<-ee_extract(gee_data, loc.buff.small, fun = ee$Reducer$mean(),
    #   scale = as.integer(meth.gee$SpatialResolution[3]))

    # Since National file says to use native I am extracting that from the object and using it.
    nominal_scale <- gee_data$projection()$nominalScale()$getInfo()

    # set scale to default for now
    gee_data1 <- ee_extract(gee_data, loc.buff.small, fun = ee$Reducer$mean(),
                            scale = nominal_scale)
    message("Test run successful")

  }



  #2.Extract Temporally Static #===============================================

  ##2.1 Get list of static layers to run
  meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine", Use==1,
                            TemporalResolution=="static")

  if(!nrow(meth.gee) == 0){
    ##2.3. Create a plain dataframe
    loc.gee.static <- loc.n %>% sf::st_drop_geometry()

    ##2.4. Make method loop
    for(i in 1:nrow(meth.gee)){
      print(i)
      if(meth.gee$GEEtype[i]=="image"){
        if(!is.na(meth.gee$GEEBand[i]))
          img.i <- ee$Image(meth.gee$Link[i])$
            select(meth.gee$GEEBand[i])
        else img.i <- ee$Image(meth.gee$Link[i])
      }
      if(meth.gee$GEEtype[i]=="imagecollection"){
        img.i <- ee$ImageCollection(meth.gee$Link[i])$
          select(meth.gee$GEEBand[i])$
          toBands()
      }

      loc.gee.i <- switch(as.character(meth.gee$RadiusExtent[i]),
                          "200"=loc.buff, #buffer width!
                          "2000"=loc.buff2, #buffer width!
                          "NA"=loc.n)

      n = 1000 #number of calls per query/
      #GEE only accepts 5000 at time, but payload are restricted then 1000 was chosen
      loc.gee.i<-loc.gee.i%>%
        select(id)%>%
        group_by(row_number() %/% n) %>% group_map(~ .x)


      gee.data.static<-data.frame()

      for (j in 1:length(loc.gee.i)){
        # Issue with ee_extract so I created a modified version. See Issue:
        # https://github.com/r-spatial/rgee/issues/367
        print("check mean or cv")
        if(meth.gee$RadiusFunction[i]=="mean"){
          gee.data.static.extract <- ee_extract(img.i, loc.gee.i[[j]],
                                                fun = ee$Reducer$mean(),
                                                scale = as.integer(meth.gee$GEEScale[i]))
        }

        if(meth.gee$RadiusFunction[i]=="cv"){
          gee.data.static.extract <- ee_extract(img.i, loc.gee.i[[j]],
                                                fun= ee$Reducer$stdDev(),
                                                scale =as.integer(meth.gee$GEEScale[i]))
        }

        if(ncol(gee.data.static.extract)<2){ #Corrects for extraction when all data NA
          gee.data.static.extract<-  gee.data.static.extract %>%
            mutate(namecol=NA)
          names(gee.data.static.extract)[2]<-meth.gee$GEEBand[i]
          message("all data is NA for: ", meth.gee$Link[i])
        }

        gee.data.static<-rbind(gee.data.static, gee.data.static.extract)
        print(paste0(j, " of ", length(loc.gee.i)))
      }
      print("left join")
      names(gee.data.static)[2]<-meth.gee$Label[i]
      loc.gee.static <- left_join(loc.gee.static, gee.data.static)

      print(paste0("Finished ", i, " of ", nrow(meth.gee)))

      # Update meth to track progress
      meth$Complete[which(meth$Label == meth.gee$Label[i])] <- 1
      meth$Running[which(meth$Label == meth.gee$Label[i])] <- 0
    }

    #Zerofill
    zerocols <- meth.gee %>%
      dplyr::filter(Zerofill==1)
    loc.gee.static <- loc.gee.static%>%
      mutate_at(vars(zerocols$Label), ~replace(., is.na(.), 0))


    write.csv(loc.gee.static, file = file.path(out_save, "data_covariates_GEE_static.csv"),
              row.names = FALSE)

    write.csv(meth, file = meth_path, row.names = FALSE)
  }else{
    message("All static variables have already been completed")
  }


  #3. Extract Temporally matched #==============================================

  ##3.1. Get list of temporally matched layers to run
  meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine",
                            TemporalResolution=="match")


  if(!nrow(meth.gee) == 0){
    ##3.2. Plain dataframe for joining to output
    loc.gee.match <- loc.n %>% sf::st_drop_geometry()


    ##3.3. Set up to loop

    for(i in 1:nrow(meth.gee)){

      if(year_get == "max"){
        loc.gee.match$year <- meth.gee$GEEYearMax[i]
        meth.gee$GEEYearMin[i] <- meth.gee$GEEYearMax[i]
      }

      #Identify years of imagery
      years.gee <- seq(meth.gee$GEEYearMin[i], meth.gee$GEEYearMax[i])
      #considers having all years

      #Match year of data to year of data
      dt = data.table::data.table(year=years.gee, val=years.gee)
      data.table::setattr(dt, "sorted", "year")
      data.table::setkey(dt, year)

      loc.n.i <- switch(as.character(meth.gee$RadiusExtent[i]),
                        "200"= loc.buff, #buffer width!
                        "2000"=loc.buff2, #buffer width!
                        "NA"=loc.n)

      loc.n.i$yearrd <- dt[J(loc.n.i$year), roll = "nearest"]$val

      #Set up to loop through years
      loc.j <- data.frame()
      for(j in 1:length(years.gee)){

        loc.n.yr <- dplyr::filter(loc.n.i, yearrd==years.gee[j]) %>%
          select (id)

        if(nrow(loc.n.yr) > 0){

          #Set start & end date for image filtering---
          start.k <- paste0(years.gee[j]+meth.gee$YearMatch[i], "-",
                            meth.gee$GEEMonthMin[i], "-01")

          if(meth.gee$GEEMonthMax[i] > meth.gee$GEEMonthMin[i]){
            end.k <- paste0(years.gee[j]+meth.gee$YearMatch[i], "-",
                            meth.gee$GEEMonthMax[i], "-28")
          }
          if(meth.gee$GEEMonthMax[i] < meth.gee$GEEMonthMin[i]){
            end.k <- paste0(years.gee[j], "-", meth.gee$GEEMonthMax[i], "-28")
          }

          #Get the image
          img.i <- ee$ImageCollection(meth.gee$Link[i])$
            filter(ee$Filter$date(start.k, end.k))$select(meth.gee$GEEBand[i])$mean()

          #Extract
          if (nrow(loc.n.yr)>1000) { #5000 entries is the limit for google earth engine
            #Some maps exceed the payload size, therefore we selected 1000
            n = 1000
            split <- loc.n.yr %>% group_by(row_number() %/% n) %>% group_map(~ .x)
          }

          if (nrow(loc.n.yr)<=1000){
            split <- list (loc.n.yr)
          }

          gee.data<-data.frame()
          for (h in 1:length(split)){
            if(meth.gee$Extraction[i]=="point"){
              gee.ext.data <- ee_extract(img.i, split[[h]],
                                         scale =as.integer(meth.gee$GEEScale[i]))
            }

            if(meth.gee$Extraction[i]=="radius"){

              if(meth.gee$RadiusFunction[i]=="mean"){
                gee.ext.data <- ee_extract(img.i, split[[h]],
                                           fun= ee$Reducer$mean(),
                                           scale =as.integer(meth.gee$GEEScale[i]))
              }

              if(meth.gee$RadiusFunction[i]=="mode"){
                gee.ext.data <- ee_extract(img.i, split[[h]],
                                           fun= ee$Reducer$mode(),
                                           scale =as.integer(meth.gee$GEEScale[i]))
              }
            }

            if(ncol(gee.ext.data)<2){ #Corrects for extraction when all data NA
              gee.ext.data<-  gee.ext.data %>%
                mutate(namecol=NA)
              names(gee.ext.data)[2]<-meth.gee$GEEBand[i]

              message("all data is NA for: ", meth.gee$Link[i], " and year ", years.gee[j])
            }
            gee.data <- rbind(gee.data, gee.ext.data)
          }

          #Collapse data across years
          names(gee.data)[2]<-meth.gee$Label[i]
          loc.j<-rbind(loc.j, gee.data)
        }

        print(paste0("Finished year ", j, " of ", length(years.gee)))

      }
      loc.gee.match <- left_join(loc.gee.match, loc.j)
      print(paste0("Finished ", i, " of ", nrow(meth.gee)))

      meth$Complete[which(meth$Label == meth.gee$Label[i])] <- 1
      meth$Running[which(meth$Label == meth.gee$Label[i])] <- 0
    }

    ##3.4. Save
    write.csv(loc.gee.match, file = file.path(out_save, "data_covariates_GEE_match.csv"),
              row.names = FALSE)

    write.csv(meth, file = meth_path, row.names = FALSE)
  }else{
    message("All temporally matched variables have already been completed")
  }
}

