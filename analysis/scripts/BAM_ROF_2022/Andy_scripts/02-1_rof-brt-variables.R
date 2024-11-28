


## check variable importance

inherits(res, "try-error")

library(mefa4)
library(gbm)
library(dismo)
library(ggplot2)

SPP <- list.files("2_pipeline/store/brt2-xv", full.names=FALSE)
SPP <- sapply(1:length(SPP), function(x) strsplit(SPP[x], "[.]")[[1]][1])

var_table <- readRDS("2_pipeline/store/var_table.rds")

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
  j1 <- var_table$variable_name[which(var_table$variable == j)]
  x$var <- paste0(j1, " (", round(u$rel.inf[i], 2), "%)")
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

d2 <- "2_pipeline/store/brt2-xv-pred-mosaic"
if(!dir.exists(d2)){dir.create(d2)}


for (spp in SPP) {
  cat(spp, "\n")
  load(paste0("2_pipeline/store/brt2-xv/", spp, ".RData"))
  if (inherits(res, "gbm")) {
    p <- plot_fun(res)
    ggsave(sprintf("2_pipeline/store/brt2-xv-pred-mosaic/%s-effects12.png", spp), p)
  }
}


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

# Retrieve the deviance and correlation statistics
cv.stats <- NULL
for (spp in SPP) {
  cat(spp, "\n")
  load(paste0("../../2_pipeline/store/brt2-xv/", spp, ".RData"))
  if (inherits(res, "gbm")) {
    p <- data.frame(spp = spp, deviance = res$cv.statistics$deviance.mean, cor = res$cv.statistics$correlation.mean)
    cv.stats <- rbind(cv.stats, p)
  }
}
write.csv(cv.stats, row.names=FALSE, file="../../3_outputs/tables/cv_stats.csv")
saveRDS(cv.stats, "2_pipeline/store/cv_stats.rds")
