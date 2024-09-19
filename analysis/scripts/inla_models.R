# ******************************************************************************
# ******************************************************************************
# Implementation of David Isle's Landscape Distribution Modeling ECCC workflow
# for Northern Ontario: Running models
# ******************************************************************************
# ******************************************************************************
# Includes Fitting models, generating predictions, mapping predictions and model
# validation through cross-validation
# ******************************************************************************
# Notable changes of method:
# ******************************************************************************

library(INLA)
library(inlabru)
library(tidyverse)
library(sf)
library(ggtext)


devtools::load_all(".")

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)

Study_Area_bound <- read_sf("analysis/data/derived_data/INLA_data/study_area.gpkg")

analysis_data <- readRDS("analysis/data/derived_data/INLA_data/analysis_data_package.rds")
# attach(analysis_data)

covars <- terra::rast("analysis/data/derived_data/INLA_data/covariate_stack.grd")

# Raster with target properties
target_raster <- covars$mat

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

BCR_PROV <- st_read("analysis/data/raw_data/BCR_Terrestrial/BCR_Terrestrial_master.shp")  %>%
  subset(PROVINCE_S == "ONTARIO") %>%
  st_make_valid() %>%
  group_by(BCR,PROVINCE_S)%>%
  summarize(geometry = st_union(geometry)) %>%
  st_transform(st_crs(Study_Area_bound))

# Loop through species, fit models, generate maps
species_to_model <- analysis_data$species_to_model %>%
  arrange(desc(n_squares),desc(n_detections))

# adding a more southern dispersed species to compare
species_to_fit <- species_to_model %>%
  subset(Species_Code_BSC %in% c("SOSA","LEYE","GRYE", "DEJU"))

# This is a functionified version of what David had in 4_Analysis_INLA
for (sp_code in species_to_fit$Species_Code_BSC){

  sp_mod <- fit_inla(sp_code, analysis_data, st_crs(Study_Area_bound), Study_Area_bound)

  sp_pred <- predict_inla(analysis_data$ONGrid, analysis_data, sp_mod, sp_code)

  map_inla_preds(sp_code, analysis_data, sp_pred, st_crs(Study_Area_bound), ONSquares)

}

# Evaluate model performance by splitting data into modeling and test sets
# Splits of interest include:
# # train on old data (eg atlas 2) test on far north
# # Spatial cross-validation to test overall performance David's
#   5_Model_Comparison_Crossval makes grid and assigns to n folds randomly

# Evaluate across spatial cross-validation folds #===============================
# Create spatial folds
set.seed(999)
n_folds <- 5
Crossval_Grid <- st_make_grid(
  Study_Area_bound,
  cellsize = units::set_units(10*10,km^2),
  what = "polygons",
  square = TRUE,
  flat_topped = FALSE
) %>%
  st_as_sf() %>%
  na.omit() %>%
  mutate(id = sample(1:nrow(.),nrow(.))) %>%
  mutate(Crossval_Fold = cut(id,breaks = seq(0,max(id)+1,length.out = n_folds+1)) %>% as.factor() %>% as.numeric())%>%
  dplyr::select(Crossval_Fold,x)

analysis_data$all_surveys <- analysis_data$all_surveys %>% st_intersection(Crossval_Grid)

model_performance <- data.frame()
for (sp_code in species_to_fit$Species_Code_BSC){
  for (n in unique(analysis_data$all_surveys$Crossval_Fold)) {
    sp_mod <- fit_inla(sp_code, analysis_data, st_crs(Study_Area_bound), Study_Area_bound,
                       train_dat_filter = paste0("Crossval_Fold != ", n),
                       save_mod = TRUE, file_name_bit = n)

    sp_pred <- predict_inla(analysis_data$all_surveys %>% filter(Crossval_Fold == n),
                            analysis_data, sp_mod, sp_code)

    sp_perf <- evaluate_preds(sp_pred, sp_mod, sp_code, analysis_data)

    model_performance <- bind_rows(model_performance, sp_perf)
  }
}

model_performance %>% group_by(species) %>%
  summarise(across(-c(fold), lst(mean, min, max))) %>%
  pivot_longer(-species, names_to = c("variable", "sum_var"), names_sep = "_m") %>%
  mutate(sum_var = paste0("m", sum_var)) %>%
  pivot_wider(names_from = "sum_var", values_from = "value") %>%
  mutate(range = max - min,
         across(-c(species, variable), \(x)round(x, 4))) %>%
  View()

# TODO: I added the eval metric that NEON uses CRPS but I am wondering if it is
# a good one when we have so many zeros, I noticed when looking at the data it
# seems like the CRPS very low (aka good) score whenever the
# observed value is zero, regardless of the CI width around 0.
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011393

# Train on old data test on new data #===============================
analysis_data$all_surveys <- analysis_data$all_surveys %>%
  mutate(Crossval_Fold = ifelse(year(Date_Time) > 2015, "new", "old"))
model_performance2 <- data.frame()
for (sp_code in species_to_fit$Species_Code_BSC){
  for (n in unique(analysis_data$all_surveys$Crossval_Fold)) {
    sp_mod <- fit_inla(sp_code, analysis_data, st_crs(Study_Area_bound), Study_Area_bound,
                       train_dat_filter = paste0("Crossval_Fold == '", n, "'"),
                       save_mod = TRUE, file_name_bit = n)

    sp_pred <- predict_inla(analysis_data$all_surveys %>% filter(Crossval_Fold != n),
                            analysis_data, sp_mod, sp_code)

    sp_perf <- evaluate_preds(sp_pred, sp_mod, sp_code, analysis_data)

    model_performance2 <- bind_rows(model_performance2, sp_perf)
  }
}

model_performance2 %>% group_by(species) %>%
  summarise(across(-c(fold), lst(mean, min, max)))

# use each model to make maps and compare
for (sp_code in species_to_fit$Species_Code_BSC){
  for (n in unique(analysis_data$all_surveys$Crossval_Fold)) {
    sp_mod <- readRDS(paste0("analysis/data/derived_data/INLA_results/models/",
                             sp_code,"_", n, "_mod.rds"))
    sp_pred <- predict_inla(analysis_data$ONGrid, analysis_data, sp_mod, sp_code)

    map_inla_preds(sp_code, analysis_data, sp_pred, st_crs(Study_Area_bound),
                   ONSquares, file_name_bit = n,
                   train_dat_filter = paste0("Crossval_Fold == '", n, "'"))
  }
  map_list <- list.files("analysis/data/derived_data/INLA_results/maps",
                         pattern = paste0(sp_code, ".*png$"),
                         full.names = TRUE)

  name_list <- list.files("analysis/data/derived_data/INLA_results/maps",
                          pattern = paste0(sp_code, "_new.*png$|", sp_code, "_old.*png$"))

  # WIP need to order by these so maps can be compared
  map_vars <- name_list %>%
    str_remove(paste0(sp_code, "_new_|", sp_code, "_old_|", sp_code, "_all_")) %>%
    unique()

  map_list <- map(map_vars, \(x){
    str_subset(map_list, x)
  }) %>% unlist()
  dat_list <- str_extract(map_list, "new|old|all") %>% str_to_sentence()

  plots <- map2(map_list, dat_list,
                  function(x, y){
    img <- as.raster(png::readPNG(x))
    gridExtra::arrangeGrob(grid::textGrob(y),
                           grid::rasterGrob(img, interpolate = FALSE),
                           nrow = 11, ncol = 1,
                           layout_matrix = matrix(c(1, rep(2, 10)), ncol = 1))

  })
  png(paste0("analysis/data/derived_data/INLA_results/maps/", sp_code,"_compare%03d.png"), width=10, height=6.5, units="in", res=1000, type="cairo")
  print(gridExtra::marrangeGrob(grobs = plots, nrow=3, ncol=1,top=NULL))
  dev.off()

}


