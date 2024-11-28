#There are data products that can be used to estimate the age and 
#structure of Canadian forests as of the year 2015. These estimates
#allow for more accurate estimates of predictor values in bird models
#in parts of Canada where there have been many disturbances more recent
#than what is estimated for forests in 2011 (from the raster layers
#in Beaudoin et al. [2014] in the previous R script).

#These data products include estimates of the year of disturbance by
#fire and harvest, estimates of forest height, volume, and biomass.

#These rasters are derived from a series of 30-m resolution
#national-scale (Canada-wide) Landsat remote-sensing products.
#Because of the extent and fine-scale spatial resolution of these
#products, processing them within R, even clipping them takes up
#a lot of memory and processing power. Initially, I tried to do
#all spatial processing in R but it took at least several hours to
#just do a clip to Ontario and some results looked like there were
#errors.

#I ended up deciding to clip and resample to 250-m resolution all 
#of these derived raster layers in ArcGIS since it was much faster
#than doing the clips in R and results looked more accurate. 

#Due to their size, original data products are not included in this
#R project.
#The original raster products can all be obtained from:
#https://opendata.nfis.org/mapserver/nfis-change_eng.html

#fire mask (Hermosilla et al. 2015)
#fire year (Hermosilla et al. 2015)
#fire severity (Hermosilla et al. 2015)
#harvest mask (Hermosilla et al. 2015)
#harvest year (Hermosilla et al. 2015)
#height2015 (Matasci et al. 2018)
#biomass2015 (Matasci et al. 2018)
#volume2015 (Matasci et al. 2018)

