
#' Run a BRT with cross validation to determine the best ntrees and important
#' variables
#'
#' @param spp species name
#' @param pred_data model data
#' @param learning_rate learning rate
#' @param ... other arguments passed to gbm.step
#'
#'
brt_cross_val <- function(spp, pred_data, learning_rate = 0.001, ...) {
  out <- try(
    dismo::gbm.step(pred_data,
             gbm.y = 1,
             gbm.x = 3:ncol(pred_data),
             offset = pred_data$offset,
             family = "poisson",
             tree.complexity = 3,
             learning.rate = learning_rate,
             bag.fraction = 0.5,
             ...)
  )
  # if (!inherits(out, "try-error")) {
  #   out$rof_settings <- list(learning_rate = learning_rate, spp = spp)
  # }
  out
}

#' Make dataframe of response and predictors for one species
#'
#' @param y sparse matrix with response variable
#' @param off offsets
#' @param spp species names
#' @param samp_row_ind index of rows to use. Based on resampling
#' @param pred_vars predictor variables
#' @param pred_var_nms predictor variable names
#'
make_pred_data <- function(spp, samp_row_ind, y, off, pred_vars, pred_var_nms) {
  if (sum(y[samp_row_ind, spp]) < 1) {
    return(structure(sprintf("0 detections for %s", spp), class = "try-error"))
  }
  data.frame(
    count = as.numeric(y[samp_row_ind, spp]),
    offset = off[samp_row_ind, spp],
    pred_vars[samp_row_ind, pred_var_nms]
  )
}



#' Plot marginal effect plots for the top 12 variables
#'
#' @param res model object
#' @param rel_inf table of relative influence
#' @param spp species name

plot_fun <- function(res, rel_inf, spp) {
  xx <- do.call(rbind, lapply(1:12, .plot_fun, res, rel_inf))
  p <- ggplot(xx, aes(x=x, y=y)) +
    geom_line() +
    facet_wrap(vars(var), scales="free_x") +
    ylab(paste(spp, "density (males/ha)")) +
    xlab("Predictor values") +
    theme_minimal()
}
.plot_fun <- function(i, res, rel_inf) {
  j <- as.character(rel_inf$var[i])
  x <- gbm::plot.gbm(res, j,
                     n.trees = res$n.trees,
                     return.grid=TRUE,
                     type="response")
  colnames(x) <- c("x", "y")
  x$var <- paste0(j, " (", round(rel_inf$rel.inf[i], 2), "%)")
  attr(x, "out.attrs") <- NULL
  x
}

#' Do cross validation step
#'
#' make data, fit model, save model, plot response, return ntrees, relative
#' influence in table
#'
#' @param spp
#' @param samp_row_ind
#' @param samp_col_ind
#' @param y
#' @param off
#' @param pred_vars
#' @param pred_var_nms
#' @param learning_rate
#' @param save_dir
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
run_cross_val <- function(spp, samp_row_ind, samp_col_ind, y, off, pred_vars,
                          pred_var_nms, learning_rate = 0.001, save_dir, ...) {

  dat <- make_pred_data(spp = spp, samp_row_ind = samp_row_ind, y = y,
                        off = off, pred_vars = pred_vars,
                        pred_var_nms = pred_var_nms)

  mod <- brt_cross_val(spp, dat, learning_rate, ...)

  saveRDS(mod, file.path(save_dir, paste0(spp,"_brt_xv_samp_",
                                          samp_col_ind, ".rds")))

  rel_inf <- summary(mod, plotit = FALSE) %>% as_tibble()

  # make plots
  p <- plot_fun(mod, rel_inf, spp)

  ggsave(file.path(save_dir,
                   paste0(spp, "_samp_", samp_col_ind, "_effects12.png")),
         p, height = 7, width = 7)

  tibble(spp = spp, sample = samp_col_ind, ntrees = mod$n.trees,
         rel_inf = list(rel_inf))
}

brt_boot <- function(spp, ntrees, rel_inf, pred_data, learning_rate = 0.001,
                     n_cores = 1, ...) {

  out <- try(gbm::gbm(pred_data$count ~ . + offset(pred_data$offset),
                      data = pred_data[,-(1:2)],
                      n.trees = ntrees,
                      interaction.depth = 3,
                      shrinkage = learning_rate,
                      bag.fraction = 0.5,
                      distribution = "poisson",
                      keep.data = FALSE,
                      verbose = FALSE,
                      n.cores = n_cores))
  out
}

run_boot <- function(spp, ntrees, rel_inf, resamps,
                     y, off, pred_vars, learning_rate = 0.001, save_dir,
                     n_cores = 1, ...){

  pred_var_nms <- filter(rel_inf, rel.inf > 0) %>% pull(var)

  dat_lst <- map(resamps,
                 ~make_pred_data(spp = spp,
                                 samp_row_ind = .x, y = y,
                                 off = off, pred_vars = pred_vars,
                                 pred_var_nms = pred_var_nms))

  mod_lst <- map(dat_lst, ~brt_boot(spp, ntrees, rel_inf, .x,
                                    learning_rate, n_cores, ...))

  fl_nm <- map(1:length(mod_lst),
                 ~file.path(save_dir, paste0(spp, "_brt_", ntrees,
                                             "_trees_samp_", .x, ".rds")))

  walk2(mod_lst, fl_nm, ~saveRDS(.x, file = .y))

  fl_nm
}
