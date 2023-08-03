# Modified versions of napops functions
# TODO: replace by adding to napops package through a PR if approved by maintainers

#
# cue_rate2 <- function (species = NULL, model = NULL, od = NULL, tssr = NULL,
#           pairwise = FALSE, quantiles = NULL, samples = 1000)
# {
#   rem_vcv_list <- NULL
#   rm(rem_vcv_list)
#   napops:::check_data_exists()
#   if (!is.null(species)) {
#     napops:::check_valid_species(species = species, mod = "rem")
#   }
#   if (!is.null(model)) {
#     napops:::check_valid_model(model = model, mod = "rem")
#   }
#   if (pairwise) {
#     if (length(od) != length(tssr)) {
#       stop("Pairwise set to TRUE but OD and TSSR are not the same length.")
#     }
#   }
#   sp_covars <- covariates_removal(species = species)
#   if (!is.null(od)) {
#     od_range <- range(sp_covars$OD)
#     if (any(od < od_range[1]) || any(od > od_range[2])) {
#       warning(paste0("You are providing some OD values that are outside the training values of [",
#                      od_range[1], ",", od_range[2], "] for species ",
#                      species))
#     }
#   }
#   if (!is.null(tssr)) {
#     tssr_range <- range(sp_covars$TSSR)
#     if (any(tssr < tssr_range[1]) || any(tssr > tssr_range[2])) {
#       warning(paste0("You are providing some TSSR values that are outside the training values of [",
#                      tssr_range[1], ",", tssr_range[2], "] for species ",
#                      species))
#     }
#   }
#   if (isFALSE(pairwise)) {
#     tssr_values <- rep(tssr, each = length(od))
#     sim_data <- data.frame(Intercept = rep(1, times = length(tssr_values)),
#                            TSSR = tssr_values, OD = rep(od, length(tssr)))
#   }
#   else {
#     sim_data <- data.frame(Intercept = rep(1, times = length(tssr)),
#                            TSSR = tssr, OD = od)
#   }
#   design <- sim_data
#   tssr_median <- stats::median(sp_covars$TSSR)
#   design$TSSR <- (design$TSSR - tssr_median)/24
#   design$TSSR2 <- design$TSSR^2
#   od_sp_median <- stats::median(sp_covars$OD)
#   design$OD <- (design$OD - od_sp_median)/365
#   design$OD2 <- design$OD^2
#   design <- design[, c("Intercept", "TSSR", "TSSR2", "OD",
#                        "OD2")]
#   coefficients <- coef_removal(species = species, model = model)
#   coefficients <- as.numeric(coefficients[, c("Intercept",
#                                               "TSSR", "TSSR2", "OD", "OD2")])
#   if (is.null(quantiles)) {
#     coefficients[which(is.na(coefficients))] <- 0
#     phi <- exp(as.matrix(design) %*% coefficients)
#     sim_data$CR_est <- phi
#     sim_data <- sim_data[, c("TSSR", "OD", "CR_est")]
#     return(sim_data)
#   }
#   else {
#     load(paste0(rappdirs::app_dir(appname = "napops")$data(),
#                 "/rem_vcv.rda"))
#     vcv <- rem_vcv_list[[model]][[species]]
#     bootstrap_df <- bootstrap(vcv = vcv, coefficients = coefficients,
#                               design = design, quantiles = quantiles, samples = samples,
#                               model = "rem")
#     return(cbind(sim_data[, c("TSSR", "OD")], bootstrap_df))
#   }
# }

# Version that uses loaded table to get species covariates
# also returns dataframe with all NA if species not valid
# and changes result so the CR_est is a numeric not an array
# Note: not tested with quantiles
cue_rate3 <- function (species = NULL, model = NULL, od = NULL, tssr = NULL,
                       pairwise = FALSE, quantiles = NULL, samples = 1000)
{
  rem_vcv_list <- NULL
  rm(rem_vcv_list)
  napops:::check_data_exists()
  if (!is.null(species)) {
    val <- check_valid_species(species = species, mod = "rem")
    if(!val){
      return(data.frame(TSSR = NA_real_, OD = NA_real_, CR_est = NA_real_))
    }
  }
  if (!is.null(model)) {
    napops:::check_valid_model(model = model, mod = "rem")
  }
  if (pairwise) {
    if (length(od) != length(tssr)) {
      stop("Pairwise set to TRUE but OD and TSSR are not the same length.")
    }
  }
  sp_covars <- all_cov_tbl[all_cov_tbl$Species == species, ]
  if (!is.null(od)) {
    if (any(od < sp_covars$OD_min) || any(od >  sp_covars$OD_max)) {
      warning(paste0("You are providing some OD values that are outside the training values of [",
                     sp_covars$OD_min, ",",  sp_covars$OD_max, "] for species ",
                     species))
    }
  }
  if (!is.null(tssr)) {
    if (any(tssr <  sp_covars$TSSR_min) || any(tssr > sp_covars$TSSR_max)) {
      warning(paste0("You are providing some TSSR values that are outside the training values of [",
                     sp_covars$TSSR_min, ",", sp_covars$TSSR_max, "] for species ",
                     species))
    }
  }
  if (isFALSE(pairwise)) {
    tssr_values <- rep(tssr, each = length(od))
    sim_data <- data.frame(Intercept = rep(1, times = length(tssr_values)),
                           TSSR = tssr_values, OD = rep(od, length(tssr)))
  }
  else {
    sim_data <- data.frame(Intercept = rep(1, times = length(tssr)),
                           TSSR = tssr, OD = od)
  }
  design <- sim_data
  tssr_median <- sp_covars$TSSR_median
  design$TSSR <- (design$TSSR - tssr_median)/24
  design$TSSR2 <- design$TSSR^2
  od_sp_median <- sp_covars$OD_median
  design$OD <- (design$OD - od_sp_median)/365
  design$OD2 <- design$OD^2
  design <- design[, c("Intercept", "TSSR", "TSSR2", "OD",
                       "OD2")]
  coefficients <- coef_removal(species = species, model = model)
  coefficients <- as.numeric(coefficients[, c("Intercept",
                                              "TSSR", "TSSR2", "OD", "OD2")])
  if (is.null(quantiles)) {
    coefficients[which(is.na(coefficients))] <- 0
    phi <- exp(as.matrix(design) %*% coefficients)
    sim_data <- cbind(sim_data, CR_est = phi)
    sim_data <- sim_data[, c("TSSR", "OD", "CR_est")]
    return(sim_data)
  }
  else {
    load(paste0(rappdirs::app_dir(appname = "napops")$data(),
                "/rem_vcv.rda"))
    vcv <- rem_vcv_list[[model]][[species]]
    bootstrap_df <- bootstrap(vcv = vcv, coefficients = coefficients,
                              design = design, quantiles = quantiles, samples = samples,
                              model = "rem")
    return(cbind(sim_data[, c("TSSR", "OD")], bootstrap_df))
  }
}

# Version using loaded table
# also returns dataframe with all NA if species not valid
# and changes result so the CR_est is a numeric not an array
# Note: not tested with quantiles
edr2 <- function (species = NULL, model = NULL, road = NULL, forest = NULL,
                  pairwise = FALSE, quantiles = NULL, samples = 1000) {
  dis_vcv_list <- NULL
  rm(dis_vcv_list)
  napops_dir <- NULL
  rm(napops_dir)
  roadside <- NULL
  rm(roadside)
  napops:::check_data_exists()
  if (!is.null(species)) {
    val <- check_valid_species(species = species, mod = "dis")
    if(!val){
      return(data.frame(Road = NA_real_, Forest = NA_real_, EDR_est = NA_real_))
    }
  }
  if (!is.null(model)) {
    napops:::check_valid_model(model = model, mod = "dis")
  }
  if (pairwise) {
    if (length(road) != length(forest)) {
      stop("Pairwise set to TRUE but Road and Forest are not the same length.")
    }
  }
  sp_covars <- all_cov_tbl[all_cov_tbl$Species == species, ]
  if (!is.null(forest)) {
    if (any(forest > 1) || any(forest < 0)) {
      stop("Forest coverage values must be between 0 and 1.")
    }
    if (any(forest < sp_covars$Forest_min) || any(forest > sp_covars$Forest_max)) {
      warning(paste0("You are providing some Forest Coverage values that are outside the training values of [",
                     sp_covars$Forest_min, ",", sp_covars$Forest_max, "] for species ",
                     species))
    }
  }
  if (isFALSE(pairwise)) {
    forest_values <- rep(forest, each = length(road))
    sim_data <- data.frame(Intercept = rep(1, times = length(forest_values)),
                           Forest = forest_values, Road = as.integer(rep(road,
                                                                         length(forest))))
  }
  else {
    sim_data <- data.frame(Intercept = rep(1, times = length(forest)),
                           Forest = forest, Road = as.integer(road))
  }
  sim_data$RoadForest <- sim_data$Road * sim_data$Forest
  design <- sim_data
  coefficients <- coef_distance(species = species, model = model)
  coefficients <- as.numeric(coefficients[, c("Intercept",
                                              "Forest", "Road", "RoadForest")])
  if (is.null(quantiles)) {
    coefficients[which(is.na(coefficients))] <- 0
    tau <- exp(as.matrix(design) %*% coefficients)
    sim_data <- cbind(sim_data, EDR_est = tau)
    sim_data <- sim_data[, c("Road", "Forest", "EDR_est")]
    return(sim_data)
  }
  else {
    load(paste0(rappdirs::app_dir(appname = "napops")$data(),
                "/dis_vcv.rda"))
    vcv <- dis_vcv_list[[model]][[species]]
    bootstrap_df <- bootstrap(vcv = vcv, coefficients = coefficients,
                              design = design, quantiles = quantiles, samples = samples,
                              model = "dis")
    return(cbind(sim_data[, c("Road", "Forest")], bootstrap_df))
  }
  if (is.null(species)) {
    stop("No argument passed for species\n")
  }
  distance <- read.csv(file = paste0(napops_dir$data(), "/distance.csv"))
  species <- toupper(species)
  if (isFALSE(species %in% distance$Species)) {
    stop(paste0("Species ", species, " does not exist in NA-POPS database.\n"))
  }
  dis_sp <- distance[which(distance$Species == species), ]
  model_number <- NULL
  if (model == "best") {
    model_number <- which(dis_sp$aic == min(dis_sp$aic))
    dis_sp <- dis_sp[which(dis_sp$model == model_number),
    ]
  }
  else if (model == "full") {
    model_number <- 5
    dis_sp <- dis_sp[which(dis_sp$model == 5), ]
  }
  else if (is.numeric(model)) {
    if (model < 1 || model > 5) {
      stop(paste0(model, " is an invalid model."))
    }
    else {
      model_number <- model
      dis_sp <- dis_sp[which(dis_sp$model == model), ]
    }
  }
  else {
    stop(paste0(model, " is an invalid model."))
  }
  if (model_number == 1) {
    return(exp(dis_sp$intercept))
  }
  else if (model_number == 2) {
    if (is.null(roadside)) {
      stop("No argument supplied for roadside status.")
    }
    else if (isFALSE(is.logical(roadside))) {
      stop("Invalid argument for roadside status. Must be TRUE or FALSE.")
    }
    return(exp(dis_sp$intercept + as.numeric(roadside) *
                 dis_sp$road))
  }
  else if (model_number == 3) {
    if (is.null(forest)) {
      stop("No argument supplied for forest coverage.")
    }
    else if (forest < 0 || forest > 1) {
      stop("Invalid argument for forest. Must be a proportion between 0 and 1.")
    }
    return(exp(dis_sp$intercept + forest * dis_sp$forest))
  }
  else if (model_number == 4) {
    if (is.null(roadside)) {
      stop("No argument supplied for roadside status.")
    }
    else if (isFALSE(is.logical(roadside))) {
      stop("Invalid argument for roadside status. Must be TRUE or FALSE.")
    }
    if (is.null(forest)) {
      stop("No argument supplied for forest coverage.")
    }
    else if (forest < 0 || forest > 1) {
      stop("Invalid argument for forest. Must be a proportion between 0 and 1.")
    }
    return(exp(dis_sp$intercept + as.numeric(roadside) *
                 dis_sp$road + forest * dis_sp$forest))
  }
  else if (model_number == 5) {
    if (is.null(roadside)) {
      stop("No argument supplied for roadside status.")
    }
    else if (isFALSE(is.logical(roadside))) {
      stop("Invalid argument for roadside status. Must be TRUE or FALSE.")
    }
    if (is.null(forest)) {
      stop("No argument supplied for forest coverage.")
    }
    else if (forest < 0 || forest > 1) {
      stop("Invalid argument for forest. Must be a proportion between 0 and 1.")
    }
    return(exp(dis_sp$intercept + as.numeric(roadside) *
                 dis_sp$road + forest * dis_sp$forest + (as.numeric(roadside) *
                                                           forest * dis_sp$roadforest)))
  }
}

# would be better if invalid species returned a consistent result that indicated that.
check_valid_species <- function (species = NULL, mod = NULL) {
  if(mod == "rem"){
    mod_cov_tbl <- all_cov_tbl[all_cov_tbl$Removal == 1,]
  } else if(mod == "dis"){
    mod_cov_tbl <- all_cov_tbl[all_cov_tbl$Distance == 1,]
  } else {
    stop("mod must be 'rem' or 'dis', not ", mod)
  }

  not_in <- setdiff(species, mod_cov_tbl$Species)
  if (length(not_in) > 0) {
    warning(paste0("The following species do not exist for the selected model type in NA-POPS:\n",
                   not_in), call. = FALSE)
    return(FALSE)
  } else {
    return(TRUE)
  }
}

