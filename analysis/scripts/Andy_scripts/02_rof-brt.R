# R code for fitting BRT models and calculating variable importance for birds in Far North Ontario 

library(mefa4)
library(gbm)
library(dismo)
library(tidyverse)

load("0_data/processed/BAMv6_RoFpackage_2022-04.RData")


## run BRT with xv

run_brt_xv <- function(spp, RATE=0.001) {
  i <- 1
  si <- BB[,i]
  if (sum(y[si, spp]) < 1)
    return(structure(sprintf("0 detections for %s", spp), class="try-error"))
  xi <- data.frame(
    count=as.numeric(y[si, spp]),
    offset=off[si, spp],
    ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
    xx2[si, cn2])
  out <- try(gbm.step(xi,
                      gbm.y = 1,
                      gbm.x = 3:ncol(xi),
                      offset = xi$offset,
                      family = "poisson",
                      tree.complexity = 3,
                      learning.rate = RATE,
                      bag.fraction = 0.5))
  if (!inherits(out, "try-error"))
    out$rof_settings <- list(RATE=RATE, spp=spp, i=i)
  out
}


if(!dir.exists("2_pipeline/store/brt2-xv"))
  (dir.create("2_pipeline/store/brt2-xv"))
system.time({
  for (spp in SPP) {
    cat("\n\n------------------------------", spp, "------------------------------\n\n")
    res <- run_brt_xv(spp)
    save(res, file=paste0("2_pipeline/store/brt2-xv/", spp, ".RData"))
  }
})


## check variable importance

inherits(res, "try-error")

library(mefa4)
library(gbm)
library(dismo)
library(ggplot2)

rel_inf <- function(res) {
  rel.inf <- relative.influence(res, res$n.trees)
  rel.inf[rel.inf < 0] <- 0
  i <- order(-rel.inf)
  rel.inf <- 100 * rel.inf/sum(rel.inf)
  out <- data.frame(var = res$var.names[i], rel.inf = rel.inf[i])
  attr(out, "n.trees") <- res$n.trees
  out
}

.plot_fun <- function(i, res, u) {
  j <- as.character(u$var[i])
  x <- plot.gbm(res, j,
                n.trees = res$n.trees,
                return.grid=TRUE,
                type="response",
                ylab=paste(res$rof_settings$spp, "density (males/ha)"),
                xlab=paste0(j, " (", round(u$rel.inf[i], 2), "%)"))
  colnames(x) <- c("x", "y")
  x$var <- paste0(j, " (", round(u$rel.inf[i], 2), "%)")
  attr(x, "out.attrs") <- NULL
  x
}

plot_fun <- function(res) {
  u <- rel_inf(res)
  xx <- do.call(rbind, lapply(1:12, .plot_fun, res, u))
  p <- ggplot(xx, aes(x=x, y=y)) +
    geom_line() +
    facet_wrap(vars(var), scales="free_x") +
    ylab(paste(res$rof_settings$spp, "density (males/ha)")) +
    xlab("Predictor values") +
    theme_minimal()
}


if(!dir.exists("2_pipeline/store/brt2-xv-pred-mosaic")){dir.create("2_pipeline/store/brt2-xv-pred-mosaic")}


mod.fit <- NULL
for (spp in SPP) {
  cat(spp, "\n")
  load(paste0("2_pipeline/store/brt2-xv/", spp, ".RData"))
  if (inherits(res, "gbm")) {
    cvstats <- as.data.frame(res$cv.statistics[c(1,3)])
    cvstats$deviance.null <- res$self.statistics$mean.null
    cvstats$deviance.exp <- (cvstats$deviance.null-cvstats$deviance.mean)/cvstats$deviance.null
    
    mod.fit <- rbind(mod.fit, data.frame(spp = spp, cvstats))
  }
}
mod.fit$deviance.exp[which(mod.fit$deviance.exp < 0)] <- 0
saveRDS(mod.fit, "2_pipeline/store/model_fit.rds")


# See best variables
RIall <-NULL
for (spp in SPP) {
  cat(spp, "\n")
  load(paste0("2_pipeline/store/brt2-xv/", spp, ".RData"))
  if (inherits(res, "gbm")) {
    u <- rel_inf(res)
    u$spp <- spp
    RIall <- rbind(RIall, u)
    p <- plot_fun(res)
    ggsave(sprintf("2_pipeline/store/brt2-xv-pred-mosaic/%s-effects12.png", spp), p)
  }
}
write.csv(RIall, row.names=FALSE, file="3_outputs/tables/SppBRTVarImp_v2.csv")

RIall$n0 <- ifelse(RIall$rel.inf > 0, 1, 0)
z <- xtabs(~var+n0,RIall)
z <- z[,"1"]/rowSums(z)
sort(z)
z["eskerpoint"]
which(names(sort(z)) == "eskerpoint")
""

RIall <- read.csv(file="3_outputs/tables/SppBRTVarImp_v2.csv")
varImp <- data.frame(RIall %>% group_by(var) %>% summarize(med = median(rel.inf)) %>% arrange(desc(med)))
colnames(varImp) <- c("var", "median_value")

cn2x <- c("ecozone", cn2)

cn2x_var <- c("Ecozone", "Agriculture", "Bedrock", "Biomass (2015)", "Bog", "Communities", "Conifer treed", 
             "Deciduous treed", "Disturbance", "Elevation", "Esker", "Fen", "Vegetated", "Vegetated (non-treed)", 
             "Vegetated (treed)", "Abies balsamia", "Acer negundo", "Acer pensylvanicum", "Acer rubrum", 
             "Acer saccharum", "Acer saccharinum", "Acer spicatum", "Acer spp", "Alnus spp", "Betula alleghaniensis", 
             "Betula papyrifera", "Betula populifolia", "Fagus grandifolia", "Fraxinus americana", "Fraxinus nigra", 
             "Fraxinus pennsylvanica", "Conifer spp", "Hardwood spp", "Larix laricina", "Picea abies", "Picea glauca", 
             "Picea mariana", "Picea rubens", "Pinus banksiana", "Pinus contorta", "Pinus resinosa", "Pinus strobus", 
             "Populus balsamifera", "Populus grandifolia", "Populus spp", "Populus tremuloides", "Prunus pensylvanica", 
             "Quercus macrocarpa", "Quercus rubra", "Salix spp", "Sorbus americana", "Thuja occidentalis", 
             "Tilia americana", "Tsuga canadensis", "Ulmus americana", "Broadleaf spp", "Needleleaf spp", "Unknown spp", 
             "Total dead biomass", "Stand age", "Total biomass", "Heath", "Height", "Lidar height", "Marsh", "Mixed treed", 
             "Mudflat", "Open water", "Road", "Slope", "Sparsely treed", "Swamp", "Topopgraphic position index", 
             "Tree cover", "Turbid water", "Volume (2015)")


cn2x_scale <- c(250, 750, 750, 250, 750, 750, 750, 750, 750, 250, 250, rep(750, 51), 250, 250, rep(750, 4), "none", 250, 750, 750, 250, 250, 750, 250)

var_table <- data.frame(variable = cn2x, variable_name = cn2x_var, scale = cn2x_scale)

varImp_table <- data.frame(var_table[match(varImp$var, var_table$variable), ], med_influence = varImp$med)

saveRDS(varImp_table, "2_pipeline/store/varImp_table.rds")

saveRDS(var_table, file = "2_pipeline/store/var_table.rds")



# Retrieve the deviance and correlation statistics
cv.stats <- NULL
for (spp in SPP) {
  cat(spp, "\n")
  load(paste0("2_pipeline/store/brt2-xv/", spp, ".RData"))
  if (inherits(res, "gbm")) {
    p <- data.frame(spp = spp, deviance = res$cv.statistics$deviance.mean, cor = res$cv.statistics$correlation.mean)
    cv.stats <- rbind(cv.stats, p)
  }
}
write.csv(cv.stats, row.names=FALSE, file="3_outputs/tables/cv_stats.csv")
saveRDS(cv.stats, "2_pipeline/store/cv_stats.rds")

