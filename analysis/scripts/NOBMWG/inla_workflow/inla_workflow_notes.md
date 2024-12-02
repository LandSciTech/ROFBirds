# Sarah's notes on David Ilse's INLA bird modelling workflow


## Bird data download

TODO add notes on gaining data access and download process as well as questions
about species codes and data that is on both naturecounts and WildTrax.

## Covariate extraction

Waiting on final list of variables to use. In the mean time I have incorporated
BAM's method for getting variables from Google Earth Engine (GEE) with David's
existing covariate extraction from rasters.

### Issues
  - For categorical variables BAM uses the mode to get value but in David's we
  use proportion of the polygon with each land cover class. Differences will be
  especially pronounced for larger radius. It looks like the
  ee.Reducer.frequencyHistogram() function might give the frequency of the
  different classes which could be converted to proportion. See this post
  https://gis.stackexchange.com/questions/415445/calculate-coverage-of-class-types-in-gee

  - Don't know how BAM gets their prediction rasters. David currently makes
  prediction grids in vector form. Could modify BAM code to download whole
  rasters or use David's method of getting values from grid centre point and
  combine that with BAM's method for getting raster values in a radius around a
  point. The later would require new code to just get the most recent year.
