# functions for inla model workflow

#' Not in
`%!in%` <- Negate(`%in%`)

#' Download location and year information from WildTrax
#'
#' @param proj_id proect_id in wildTrax
#'
#' @return
#' @export
#'
#' @examples
wt_dl_loc_year <- function(proj_id, sens_id){
  rec_id <- ifelse(sens_id == "ARU", "recording", "point_count")
  date_col <- ifelse(sens_id == "ARU", "recording_date_time", "survey_date")

  loc <- wt_download_report(project_id = proj_id, sensor_id = sens_id,
                            reports = c("location"), weather_cols = FALSE)

  if (is.null(nrow(loc))) return(NULL)

  loc <- loc %>%
    select(location_id,latitude,longitude) %>%
    unique()

  rec <- wt_download_report(project_id = proj_id, sensor_id = sens_id,
                            reports = c(rec_id), weather_cols = FALSE) %>%
    mutate(year = as.numeric(year(.data[[date_col]]))) %>%
    select(project_id,location_id,year) %>%
    unique()

  dat <- full_join(loc, rec, by = join_by(location_id))
  dat
}


#' Rasterize a series of spatial predictions (needed for plotting)

cut.fn <- function(df = NA,
                   target_raster = NA,
                   column_name = NA,
                   lower_bound = NA,
                   upper_bound = NA){

  max_val <- upper_bound
  max_val <- ifelse(is.na(max_val), 0, max_val)
  max_lev <- ifelse(max_val > 1.6, 4,ifelse(max_val > 0.8, 4, 3))

  cut_levs <- signif(max_val/(2^((max_lev-1):0)), 2)
  cut_levs <- unique(cut_levs)
  cut_levs <- ifelse(is.na(cut_levs), 0, cut_levs)

  if (lower_bound %in% cut_levs) cut_levs <- cut_levs[-which(cut_levs == lower_bound)]
  if (lower_bound > min(cut_levs)) cut_levs = cut_levs[-which(cut_levs < lower_bound)]

  max_lev <- length(cut_levs)

  cut_levs_labs <- c(paste0("0-",lower_bound),
                     paste(lower_bound, cut_levs[1], sep="-"),
                     paste(cut_levs[-max_lev], cut_levs[-1], sep="-"),
                     paste(cut_levs[max_lev], "+"))

  cut_levs <- c(-1, lower_bound, cut_levs, 1000) %>% unique()

  df <- mutate(df, levs = cut(.data[[column_name]], cut_levs, labels = cut_levs_labs))
# TODO: change this to use terra. It is faster now and easier I think
  tgt <- stars::st_as_stars(target_raster)
  tmp = stars::st_rasterize(df %>% dplyr::select(levs, geometry),
                            nx = dim(tgt)[1],ny = dim(tgt)[2])

  return(list(raster = tmp,cut_levs = cut_levs))
}


#' Fit boosted regression tree models
#'
fit_brt <- function(model_data,
                    response_column = NA,
                    covariate_columns = NA){
  mod_brt <- NULL

  ntrees <- 50
  tcomplexity <- 5
  lrate <- 0.01
  m <- 0

  while(is.null(mod_brt)){

    m <- m + 1
    if(m < 11){
      ntrees <- 50
      lrate <- 0.01
    } else if(m < 21){
      lrate <- 0.001
    } else if(m < 31){
      ntrees <- 25
      lrate <- 0.001
    } else if(m < 41){
      ntrees <- 25
      lrate <- 0.0001
    } else{
      break
    }

    ptm <- proc.time()
    if(inherits(try(
      mod_brt <- dismo::gbm.step(data = model_data,
                                 gbm.x = covariate_columns,
                                 gbm.y = response_column,
                                 offset = model_data$log_QPAD_offset,
                                 family = "poisson",
                                 tree.complexity = tcomplexity,
                                 learning.rate = lrate,
                                 n.trees = ntrees,
                                 n.folds = 5,
                                 max.trees = 10000)
    ), "try-error")){
      cat("Couldn't fit model", n, "in the iteration", m, "\n")
    }
    t <- proc.time() - ptm
  }
  if(is.null(mod_brt)){
    next
  }
  return(mod_brt)
}

