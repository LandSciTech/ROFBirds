#' Crop and then transform
#'
#' Transforms bbox to raster CRS, then crops raster, then transforms raster to
#' given CRS
#'
#' @param ras SpatRaster
#' @param bbx sf polygon of bbox
#' @param crs_out crs to use. As returned by `terra::crs`
#'
#' @return cropped and transformed SpatRaster
crop_transform <- function(ras, bbx, crs_out){
  bbx1 <- st_transform(bbx, st_crs(ras)) %>% terra::vect()

  bbx2 <- st_transform(bbx, crs_out) %>% terra::vect()

  ras <- terra::crop(ras, bbx1, snap = "out")

  ras <- terra::project(ras, crs_out)

  ras <- terra::crop(ras, bbx2, snap = "out")

  return(ras)
}


#' Transform points then extract
#'
#' @param ras SpatRaster
#' @param pts SpatVector of points that will be projected to match the raster
#'   and then used to extract raster values
#'
#' @return data.frame

extract_transform <- function(ras, pts){
  pts <- terra::project(pts, terra::crs(ras))

  terra::extract(ras, pts)
}


#' Compare species across models
#'
#' @param sp6 6 letter species common name code
#' @param sp4 4 letter species common name code
#' @param samp_pts SpatVector of points to compare across
#' @param in_dat_pth root folder containing all models, will be searched
#'   recursively to match the species code
#'
#' @return Prints species code, maps and comparison graphs

do_sp_compare <- function(sp6, sp4 = "XXXX", samp_pts, val_dat, in_dat_pth){

  cat(paste0("\n\n## ", sp6, "\n\n"))

  pats <- list(bateman = paste0(sp6, ".*_breeding.*suitability.tif$"),
               bam_rof = paste0(sp4, ".tif$"),
               bam_nat = paste0(sp4, ".*CAN-Mean.tif$"),
               ebird_st = paste0(sp6, "_abundance.*.tif"))

  fls <- map(pats, \(x) list.files(in_dat_pth, x, recursive = TRUE)) %>%
    # drop empty
    compact()

  rasts <- map(fls, \(x) rast(file.path(in_dat_pth, x)))

  # for some species the ebird file has multiple layers, keep breeding or resident
  lyr_use <- which(names(rasts$ebird_st) %in% c("resident", "breeding"))
  if(length(lyr_use) == 0){
    rasts$ebird_st <- NULL
  } else {
    rasts$ebird_st <- rasts$ebird_st[[lyr_use]]
  }

  bbx_use <- samp_pts %>% st_as_sf() %>% st_bbox() %>% st_as_sfc()
  crs_use <- terra::crs(rasts$bateman)
  rasts_crop <- map(rasts, \(x) crop_transform(x, bbx = bbx_use, crs_out = crs_use))

  oldpar <- par(mfrow = c(2,2))
  iwalk(rasts_crop, \(x, nm) terra::plot(x, main = nm, background = "grey70",
                                         col = viridis::viridis(50)))
  par(oldpar)


  mod_samps <- map(rasts, \(x) extract_transform(x, samp_pts)[,2]) %>% bind_cols() %>%
    # Bateman models were multiplied by 10000 to store as integer
    mutate(ID = 1:n(), bateman = bateman/10000) %>%
    pivot_longer(-ID, names_to = "model", values_to = "prediction") %>%
    group_by(model) %>%
    mutate(scaled_pred = scale(prediction, center = FALSE)[,1]) %>%
    # need to remove NAs before ranking
    group_by(ID) %>%
    mutate(anyNA = ifelse(any(is.na(prediction)), NA, 1)) %>%
    group_by(model) %>%
    mutate(rank_pred = rank(prediction*anyNA, na.last = "keep",
                            ties.method = "average")) %>%
    group_by(ID) %>%
    mutate(rank_sd = sd(rank_pred, na.rm = T)) %>%
    ungroup() %>%
    mutate(sp4 = sp4, sp6 = sp6) %>%
    select(-anyNA)

  mod_samps_l <- mod_samps %>%
    pivot_longer(matches("pred|rank"), names_to = "metric") %>%
    mutate(metric = metric %>% factor(levels = c("prediction", "scaled_pred",
                                                 "rank_pred", "rank_sd")))

  print(ggplot(mod_samps_l, aes(ID, value, col = model))+
          geom_point()+
          facet_wrap(~metric, nrow = 2, scales = "free_y"))

  # calculating mean abundance at each site to compare to predictions
  val_dat <- val_dat %>% select(location_id, sp4) %>% rename(observed = sp4)
  val_dat2 <- val_dat %>% group_by(location_id) %>%
    summarise(observed = mean(observed))

  val_samps <- map(rasts, \(x) extract_transform(x, terra::vect(val_dat2))[,2]) %>% bind_cols() %>%
    bind_cols(val_dat2 %>% st_drop_geometry())

  val_samps_l <- val_samps %>%
    # Bateman models were multiplied by 10000 to store as integer
    mutate(bateman = bateman/10000) %>%
    pivot_longer(-c(location_id, observed),
                 names_to = "model", values_to = "prediction")

  print(ggplot(val_samps_l, aes(observed, prediction))+
          geom_point()+
          geom_smooth(method = "lm", formula = y ~ x)+
          facet_wrap(~model)+
          coord_equal(xlim = c(0,1), ylim = c(0,1)))

  perf_cont <- val_samps_l %>% nest_by(model) %>%
    mutate(corr = cor(data$observed, data$prediction, use = "na.or.complete",
                      method = "pearson"),
           lm = list(lm(observed ~ prediction, data = data)),
           coefs = broom::tidy(lm) %>% select(term, estimate) %>%
                          pivot_wider(names_from = term, values_from = estimate) %>%
             set_names(c("intercept", "slope")),
           perf = broom::glance(lm)) %>%
    select(-c(data, lm)) %>% unnest(cols = c(coefs, perf))

  # connect binary obs with predictions
  perf_bi <- val_samps %>% select(-observed) %>%
    left_join(val_dat %>% st_drop_geometry(), by = "location_id") %>%
    # Bateman models were multiplied by 10000 to store as integer
    mutate(bateman = bateman/10000) %>%
    pivot_longer(-c(location_id, observed),
                 names_to = "model", values_to = "prediction") %>%
    group_by(model) %>%
    summarise(auc = Metrics::auc(observed, prediction),
              rmse = Metrics::rmse(observed[which(!is.na(prediction))],
                                   prediction[which(!is.na(prediction))]))

  mod_perf <- left_join(perf_cont, perf_bi, by = "model") %>%
    mutate(sp4 = sp4, sp6 = sp6)

  return(lst(mod_samps, mod_perf))

}

