Simplified code to run the model fitting and data preparation procedure
developed by BAM for birds in the Ring of Fire. 


prep_BBS_data: partially re-created Lionel's scripts but in a more reproducible
way including downloading data and offset scripts. Didn't redo all of it but
added a description of what I understand Lionel did.

prepare_covariate_data: re-created process for extracting raster data to points

fit_models: re-created process for fitting bootstrapped models of bird density, 
calculates mean and sd from bootstrapped models

make_figure: create a figure for Ring of Fire project report 

Lionel_scripts: contains original scripts from Lionel Leston at BAM. This
includes code for the downloading and calculating offsets for the BBS data  

Andy_scripts: contains original scripts from Andy Crosby at BAM. This uses
inputs from Lionel's code and runs brt models and calculates density by land
cover class and conservation scores. Only the brt models are implemented in my
scripts
