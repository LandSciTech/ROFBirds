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

# Plot INLA results maps

do_res_plot <- function(pred_rast, title, subtitle, subsubtitle = "", samp_grid,
                        col_pal_fn, species_label, levs_nm, file_nm){
  res_plot <- ggplot() +

    stars::geom_stars(data = pred_rast) +
    scale_fill_manual(name = paste0("<span style='font-size:13pt'>", title,
                                    "</span><br><span style='font-size:7pt'>", subtitle,
                                    "</span><br><span style='font-size:7pt'>", subsubtitle,
                                    "</span>"),
                      values = col_pal_fn(length(levels(pred_rast[[levs_nm]]))),
                      drop = FALSE, na.translate = FALSE)+

    # BCR boundaries
    geom_sf(data = BCR_PROV, fill = "transparent", col = "gray20", linewidth = 0.5)+

    # Point count detections and surveyed squares
    geom_sf(data = subset(samp_grid, !is.na(PC_detected)),
            aes(col = as.factor(PC_detected)), size = 0.5, stroke = 0, shape = 16)+
    scale_colour_discrete(type = c("gray70", "black"), guide = NULL)+

    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(5,10,5,-20),
          legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = 1),
          legend.justification = c(1,1),
          legend.position = "inside",
          legend.position.inside = c(0.99,0.95))+
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))+

    annotate(geom="text",x=1700000,y=1930000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=2155000,y=530000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray60")+

    guides(fill = guide_legend(order = 1),
           size = guide_legend(order = 2))

  png(file_nm, width=10, height=6.5, units="in", res=1000, type="cairo")
  print(res_plot)
  dev.off()

  return(res_plot)
}

prep_sp_dat <- function(analysis_data, proj_use){
  sp_dat <- analysis_data$all_surveys %>%
    mutate(count = analysis_data$full_count_matrix[,sp_code]) %>%

    # Only select point counts
    subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>%
    st_transform(proj_use)


  # Extract ebird range for this species (if it exists)

  if (sp_code %in% names(analysis_data$species_ranges)){

    range <- analysis_data$species_ranges[[sp_code]] %>% st_transform(st_crs(sp_dat))

    # Identify distance of each survey to the edge of species range (in km)
    sp_dat$distance_from_range <- ((st_distance(sp_dat, range)) %>% as.numeric())/1000
  } else{
    sp_dat$distance_from_range <- 0
  }

  # ---
  # Generate QPAD offsets for each survey (assumes unlimited distance point counts)
  # ---

  species_offsets <- subset(analysis_data$species_to_model, Species_Code_BSC == sp_code)

  if (species_offsets$offset_exists == FALSE) sp_dat$log_QPAD_offset <- 0

  if (species_offsets$offset_exists == TRUE){
    # Calculate offset for duration of survey from species overall offset value
    # how come this doesn't include time since sunrise? Maybe because it is in the model directly?
    A_metres <- pi*species_offsets$EDR^2
    p <- 1-exp(-sp_dat$Survey_Duration_Minutes*species_offsets$cue_rate)
    sp_dat$log_QPAD_offset <- log(A_metres * p)
  }

  sp_dat
}

fit_inla <- function(sp_code, analysis_data, proj_use, study_poly){
  message("starting model for: ", sp_code)

  model_file <- paste0("analysis/data/derived_data/INLA_results/models/", sp_code, "_mod.rds")

  if(!dir.exists("analysis/data/derived_data/INLA_results/models")){
    dir.create("analysis/data/derived_data/INLA_results/models")
  }

  # Skip this species if already run
  #if (file.exists(map_file)) next

  # Prepare data for this species
  sp_dat <- prep_sp_dat(analysis_data, proj_use)

  # n_det_sq <- subset(sp_dat,count>0) %>%
  #   as.data.frame() %>%
  #
  #   summarize(n_sq = length(unique(sq_id)),
  #             n_det = sum(count>0))
  #
  # if (n_det_sq$n_sq < 30 | n_det_sq$n_det < 60) next

  # FIT MODEL WITH INLA

  covariates_to_include <- paste0("PC",1:8)

  # ---
  # Create a spatial mesh, which is used to fit the residual spatial field
  # ---

  # make a two extension hulls and mesh for spatial model
  hull <- fm_extensions(
    study_poly,
    convex = c(50000, 200000),
    concave = c(350000, 500000)
  )
  mesh_spatial <- fm_mesh_2d_inla(
    boundary = hull,
    max.edge = c(70000, 100000), # km inside and outside
    cutoff = 30000,
    crs = fm_crs(sp_dat)
  ) # cutoff is min edge
  mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
  dim(mesh_locs)
  #plot(mesh_spatial)
  #lines(BCR_PROV)


  # # Note: mesh developed using tutorial at: https://rpubs.com/jafet089/886687
  # max.edge = 50000 #diff(range(st_coordinates(sp_dat)[,1]))/30
  # bound.outer = diff(range(st_coordinates(sp_dat)[,1]))/3
  # cutoff = max.edge/5
  #
  # mesh_spatial <- inla.mesh.2d(loc = st_coordinates(sp_dat),
  #                              cutoff = max.edge/2,
  #                              max.edge = c(1,2)*max.edge,
  #                              offset=c(max.edge, bound.outer))
  # mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
  # dim(mesh_locs)
  # plot(mesh_spatial)
  # lines(BCR_PROV)

  prior_range <- c(300000,0.1) # 10% chance range is smaller than 300000
  prior_sigma <- c(0.5,0.1) # 10% chance sd is larger than 0.5
  matern_coarse <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = prior_range,
                                       prior.sigma = prior_sigma
  )

  # ---
  # Create mesh to model effect of time since sunrise (TSS)
  # ---
  sp_dat$Hours_Since_Sunrise <- as.numeric(sp_dat$Hours_Since_Sunrise)
  TSS_range <- range(sp_dat$Hours_Since_Sunrise)
  TSS_meshpoints <- seq(TSS_range[1]-0.1,TSS_range[2]+0.1,length.out = 11)
  TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
  TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
                                 prior.range = c(6,0.1),
                                 prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1

  # ---
  # Model formulas
  # ---

  sd_linear <- 1
  prec_linear <-  c(1/sd_linear^2,1/(sd_linear/2)^2)
  model_components = as.formula(paste0('~
            Intercept_PC(1)+
            range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
            spde_coarse(main = coordinates, model = matern_coarse) +',
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[1],')', collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = ', prec_linear[2],')', collapse = " + ")))


  model_formula_PC = as.formula(paste0('count ~
                  Intercept_PC +
                  log_QPAD_offset +
                  TSS +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       '+',
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))


  # ---
  # Fit with nbinomial error
  # ---

  PC_sp <- sp_dat %>% as('Spatial')
  start <- Sys.time()
  fit_INLA <- NULL
  while(is.null(fit_INLA)){
    fit_INLA <- bru(components = model_components,

                    like(family = "nbinomial",
                         formula = model_formula_PC,
                         data = PC_sp),

                    options = list(
                      control.compute = list(waic = FALSE, cpo = FALSE),
                      bru_verbose = 4))
    if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
  }

  end <- Sys.time()
  runtime_INLA <- difftime( end,start, units="mins") %>% round(2)
  message(paste0(sp_code," - ",runtime_INLA," min to fit model"))

  saveRDS(fit_INLA, model_file)

  return(fit_INLA)
}

predict_inla <- function(analysis_data, mod){

  if (sp_code %in% names(analysis_data$species_ranges)){

    range <- analysis_data$species_ranges[[sp_code]] %>% st_transform(st_crs(ONGrid))
  }

  # For every pixel on landscape, extract distance from eBird range limit
  distance_from_range = (st_centroid(analysis_data$ONGrid) %>%
                           st_distance( . , range) %>%
                           as.numeric())/1000

  ONGrid_species <- analysis_data$ONGrid %>% mutate(distance_from_range = distance_from_range)
  # TODO: Not ideal, dropping several areas that are complex due to BCR intersection
  # Needed for conversion to sp
  ONGrid_species <- ONGrid_species %>%
    st_set_geometry("geometry") %>%
    filter(st_geometry_type(geometry) == "POLYGON") %>%
    st_cast()

  # ---
  # QPAD offsets associated with a 5-minute unlimited distance survey
  # ---

  species_offsets <- subset(analysis_data$species_to_model, Species_Code_BSC == sp_code)

  log_offset_5min <- 0
  if (species_offsets$offset_exists == TRUE) log_offset_5min <- species_offsets$log_offset_5min

  # ---
  # Generate predictions on ONGrid_species raster (1 km x 1 km pixels)
  # ---
  covariates_to_include <- paste0("PC",1:8)

  pred_formula_PC = as.formula(paste0(' ~
                  Intercept_PC +
                  log_offset_5min +
                  range_effect * distance_from_range +
                  spde_coarse +',
                                      paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                      '+',
                                      paste0("Beta2_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))

  # Predictions are on log scale, and do not include variance components
  start2 <- Sys.time()
  pred <- NULL
  pred <- generate(mod,
                   as(ONGrid_species,'Spatial'),
                   formula =  pred_formula_PC,
                   n.samples = 1000)

  pred <- exp(pred)

  # Median and upper/lower credible intervals (90% CRI)
  prediction_quantiles = apply(pred,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  ONGrid_species$pred_q05 <- prediction_quantiles[1,]
  ONGrid_species$pred_q50 <- prediction_quantiles[2,]
  ONGrid_species$pred_q95 <- prediction_quantiles[3,]
  ONGrid_species$pred_CI_width_90 <- prediction_quantiles[3,] - prediction_quantiles[1,]
  ONGrid_species$CV <- apply(pred,1,function(x) sd(x,na.rm = TRUE)/mean(x,na.rm = TRUE))

  # Probability of observing species in 5-minute point count
  size <- mod$summary.hyperpar$'0.5quant'[1] # parameter of negative binomial

  # Probability of detecting species in a 5-minute point count
  prob_zero_PC <- dnbinom(0,mu=prediction_quantiles[2,],size=size)
  ONGrid_species$pObs_5min <- 1-prob_zero_PC

  end2 <- Sys.time()
  runtime_pred <- difftime( end2,start2, units="mins") %>% round(2)
  message(paste0(sp_code," - ",runtime_pred," min to generate predictions"))
  return(ONGrid_species)
}

map_inla_preds <- function(sp_code, analysis_data, preds, proj_use, atlas_squares){
  map_file <- paste0("analysis/data/derived_data/INLA_results/maps/", sp_code, "_q50.png")

  if(!dir.exists("analysis/data/derived_data/INLA_results/maps")){
    dir.create("analysis/data/derived_data/INLA_results/maps")
  }


  # Prepare data for this species
  sp_dat <- prep_sp_dat(analysis_data, proj_use)

  # Summarize atlas_squares where species was detected
  PC_detected <- sp_dat %>%
    subset(Survey_Type %in% c("Point_Count","ARU_SPT","ARU_SPM")) %>%
    st_intersection(atlas_squares)%>%
    as.data.frame() %>%
    group_by(sq_id) %>%
    summarize(PC_detected = as.numeric(sum(count)>0),
              PC_mean_count = mean(count) %>% round(2))

  # CL_detected <-sp_dat %>%
  #   subset(Survey_Type %in% c("Breeding Bird Atlas","Linear transect")) %>%
  #   as.data.frame() %>%
  #   group_by(sq_id) %>%
  #   summarize(CL_detected = as.numeric(sum(count)>0),
  #             CL_mean_count = mean(count))

  atlas_squares_species <- atlas_squares %>%
    relocate(geometry,.after = last_col()) %>%
    left_join(PC_detected, by = join_by(sq_id)) #%>% left_join(CL_detected)

  atlas_squares_centroids <- st_centroid(atlas_squares_species)

  # Label for figure and ebird range limit

  species_name = analysis_data$ON_spcd$CommonName[which(analysis_data$ON_spcd$spcd == sp_code)]
  species_label = analysis_data$ON_spcd$Label[which(analysis_data$ON_spcd$spcd == sp_code)]

  # sf object for ebird range limit (optional - not needed for plotting)

  range <- NA
  if (sp_code %in% names(analysis_data$species_ranges)){
    range <- analysis_data$species_ranges[[sp_code]]  %>%
      st_transform(st_crs(Study_Area_bound)) %>%
      st_intersection(Study_Area_bound)
  }

  # Plot median prediction

  colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
                        "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
  colpal_relabund <- colorRampPalette(colscale_relabund)

  lower_bound <- 0.01
  upper_bound <- quantile(preds$pred_q50,0.99,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)

  sp_cut <- cut.fn(df = preds,
                   target_raster = target_raster,
                   column_name = "pred_q50",
                   lower_bound = lower_bound,
                   upper_bound = upper_bound)

  raster_q50 <- sp_cut$raster

  # Median of posterior
  plot_q50 <- do_res_plot(raster_q50, "Relative Abundance",
                          "Per 5-minute point count", "(Posterior Median)",
                          atlas_squares_centroids,
                          colpal_relabund, species_label,
                          "levs",
                          map_file)

  # Plot uncertainty in prediction (width of 90% CRI)

  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)

  lower_bound <- 0.01
  upper_bound <- quantile(preds$pred_CI_width_90,0.99,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)

  raster_CI_width_90 <- cut.fn(df = preds,
                               target_raster = target_raster,
                               column_name = "pred_CI_width_90",
                               lower_bound = lower_bound,
                               upper_bound = upper_bound)$raster

  plot_CI_width_90 <- do_res_plot(raster_CI_width_90, "Relative Uncertainty",
                                  "Per 5-minute point count", "Width of 90% CI",
                                  atlas_squares_centroids,
                                  colpal_uncertainty, species_label,
                                  "levs",
                                  map_file %>% str_replace("_q50", "_CI_width_90"))

  # Plot uncertainty in prediction (coefficient of variation)

  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)

  cut_levs <- c(-0.1,0.25,0.5,1,2,5,2000)
  cut_levs_labs <- c("0 to 0.25",
                     "0.25 to 0.5",
                     "0.5 to 1",
                     "1 to 2",
                     "2 to 5",
                     "> 5")

  preds$CV_levs <- cut(as.data.frame(preds)[,"CV"],
                                cut_levs,labels=cut_levs_labs)
  raster_CV <- stars::st_rasterize(preds %>% dplyr::select(CV_levs, geometry),
                                   nx = dim(raster_q50)[1],ny = dim(raster_q50)[2])

  plot_CV <- do_res_plot(raster_CV, "Coef. of Variation",
                         "Per 5-minute point count", "", atlas_squares_centroids,
                         colpal_uncertainty, species_label,
                         "CV_levs",
                         map_file %>% str_replace("_q50", "_CV"))

  # Plot probability of observing species in a 5-minute point count

  colscale_pObs <- c("#FEFEFE",RColorBrewer::brewer.pal(5,"BuGn")[2:5])
  colpal_pObs <- colorRampPalette(colscale_pObs)

  cut_levs <- c(-0.1,0.01,0.05,0.125,0.5,1)
  cut_levs_labs <- c("0 to 0.01",
                     "0.01 to 0.05",
                     "0.05 to 0.125",
                     "0.125 to 0.50",
                     "0.50 to 1")

  preds$pObs_levs <- cut(as.data.frame(preds)[,"pObs_5min"],
                                  cut_levs,labels=cut_levs_labs)

  raster_pObs = stars::st_rasterize(preds %>% dplyr::select(pObs_levs, geometry),
                                    nx = dim(raster_q50)[1], ny = dim(raster_q50)[2])

  plot_pObs <- do_res_plot(raster_pObs, "Prob. of Observation",
                           "Per 5-minute point count",
                           "(Posterior Median)", atlas_squares_centroids,
                           colpal_pObs, species_label,
                           "pObs_levs",
                           map_file %>% str_replace("_q50", "_PObs"))

  # Density estimate (per m2) - subtract detectability offset
  species_offsets <- subset(analysis_data$species_to_model, Species_Code_BSC == sp_code)

  log_offset_5min <- 0
  if (species_offsets$offset_exists == TRUE) log_offset_5min <- species_offsets$log_offset_5min

  if (log_offset_5min != 0){

    preds$density_per_ha_q50 <- preds$pred_q50 / exp(log_offset_5min) * 10000

    colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
                          "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
    colpal_relabund <- colorRampPalette(colscale_relabund)

    lower_bound <- 0.01
    upper_bound <- quantile(preds$density_per_ha_q50,0.99,na.rm = TRUE) %>% signif(2)
    if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)

    sp_cut <- cut.fn(df = preds,
                     target_raster = target_raster,
                     column_name = "density_per_ha_q50",
                     lower_bound = lower_bound,
                     upper_bound = upper_bound)

    raster_dens <- sp_cut$raster

    # Median of posterior
    plot_dens <- do_res_plot(raster_dens, "Density", "Males per hectare",
                             "(Posterior Median)", atlas_squares_centroids,
                             colpal_relabund, species_label, "levs",
                             map_file %>% str_replace("_q50", "_density"))
  }
}


