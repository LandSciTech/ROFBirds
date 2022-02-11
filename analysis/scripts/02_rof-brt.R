# R code for fitting BRT models and calculating variable importance for birds in Far North Ontario 

library(mefa4)
library(gbm)
library(dismo)


#load("d:/bam/2021/rof/BAMv6_RoFpackage.RData")
load("0_data/processed/BAMv6_RoFpackage_2022-01.RData")

cn2 <- c("eskerpoint",
         "agriculture_G750.O", "bedrock_G750.O", "biomass2015.ntems",
         "bog_G750.O", "communities_G750.O", "coniftreed_G750.O", "decidtreed_G750.O",
         "disturbance_G750.O", "elev", "fen_G750.O", "G750LandCover_Veg_v1.grd",
         "G750LandCover_VegNonTreed_v1.grd", "G750LandCover_VegTreed_v1.grd",
         "G750Species_Abie_Bal_v1.grd", "G750Species_Acer_Neg_v1.grd",
         "G750Species_Acer_Pen_v1.grd", "G750Species_Acer_Rub_v1.grd",
         "G750Species_Acer_Sac_v1.grd", "G750Species_Acer_Sah_v1.grd",
         "G750Species_Acer_Spi_v1.grd", "G750Species_Acer_Spp_v1.grd",
         "G750Species_Alnu_Spp_v1.grd", "G750Species_Betu_All_v1.grd",
         "G750Species_Betu_Pap_v1.grd", "G750Species_Betu_Pop_v1.grd",
         "G750Species_Fagu_Gra_v1.grd", "G750Species_Frax_Ame_v1.grd",
         "G750Species_Frax_Nig_v1.grd", "G750Species_Frax_Pen_v1.grd",
         "G750Species_Genc_Spp_v1.grd", "G750Species_Genh_Spp_v1.grd",
         "G750Species_Lari_Lar_v1.grd", "G750Species_Pice_Abi_v1.grd",
         "G750Species_Pice_Gla_v1.grd", "G750Species_Pice_Mar_v1.grd",
         "G750Species_Pice_Rub_v1.grd", "G750Species_Pinu_Ban_v1.grd",
         "G750Species_Pinu_Con_v1.grd", "G750Species_Pinu_Res_v1.grd",
         "G750Species_Pinu_Str_v1.grd", "G750Species_Popu_Bal_v1.grd",
         "G750Species_Popu_Gra_v1.grd", "G750Species_Popu_Spp_v1.grd",
         "G750Species_Popu_Tre_v1.grd", "G750Species_Prun_Pen_v1.grd",
         "G750Species_Quer_Mac_v1.grd", "G750Species_Quer_Rub_v1.grd",
         "G750Species_Sali_Spp_v1.grd", "G750Species_Sorb_Ame_v1.grd",
         "G750Species_Thuj_Occ_v1.grd", "G750Species_Tili_Ame_v1.grd",
         "G750Species_Tsug_Can_v1.grd", "G750Species_Ulmu_Ame_v1.grd",
         "G750SpeciesGroups_Broadleaf_Spp_v1.grd", "G750SpeciesGroups_Needleleaf_Spp_v1.grd",
         "G750SpeciesGroups_Unknown_Spp_v1.grd", "G750Structure_Biomass_TotalDead_v1.grd",
         "G750Structure_Stand_Age_v1.grd", "G750Structure_Volume_Total_v1.grd",
         "heath_G750.O", "height2015.ntems", "LIDARheight", "marsh_G750.O",
         "mixedtreed_G750.O", "mudflat_G750.O", "openwater_G750.O", "road_yesno",
         "slope", "sparsetreed_G750.O", "swamp_G750.O", "TPI", "treecover",
         "turbidwater_G750.O", "volume2015.ntems")

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

RIall$n0 <- ifelse(RIall$rel.inf > 0, 1, 0)
z <- xtabs(~var+n0,RIall)
z <- z[,"1"]/rowSums(z)
sort(z)
z["eskerpoint"]
which(names(sort(z)) == "eskerpoint")



