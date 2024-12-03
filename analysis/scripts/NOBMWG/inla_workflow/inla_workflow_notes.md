# Sarah's notes on David Ilse's INLA bird modelling workflow


## Bird data download

TODO add notes on gaining data access and download process as well as questions
about species codes and data that is on both naturecounts and WildTrax.

## Covariate extraction

Waiting on final list of variables to use. In the mean time I have incorporated
BAM's method for getting variables from Google Earth Engine (GEE) with David's
existing covariate extraction from rasters.

### Issues/suggestions
  - For categorical variables BAM uses the mode to get value but in David's we
  use proportion of the polygon with each land cover class. Differences will be
  especially pronounced for larger radius. It looks like the
  ee.Reducer.frequencyHistogram() function might give the frequency of the
  different classes which could be converted to proportion. See this post
  https://gis.stackexchange.com/questions/415445/calculate-coverage-of-class-types-in-gee

  - Don't know how BAM gets their prediction rasters. David currently makes
  prediction grids in vector form. Could modify BAM GEE extraction code to download whole
  rasters. I have used David's method of getting values from prediction grid centre point and
  combined that with BAM's method for getting raster values in a radius around a
  point and modified it to use the most recent year for predictions. Currently
  radius is set in the spreadsheet but not sure if that makes sense if we use
  those prediction for mapping

  - Various conversions for covariates are hard coded into the covariate
  assembly step. It would be better if columns were added to the spreadsheet to
  define these post extraction steps so the code could run without all the
  variables.

