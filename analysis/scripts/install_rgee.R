# Installation of Google Earth Engine R package and dependencies
#
# Needed for covariate extraction

# Documentation from the rgee package has instructions for this but they did not
# work so I followed this Youtube video instead:
# https://www.youtube.com/watch?v=olNYYynSJfI

# Installation #======================================================
# Only needs to be run once
install.packages("rgee")
library(rgee)


install.packages("reticulate")
library(reticulate)
reticulate::py_available()
reticulate::py_discover_config()
# copy the python path returned

rgee::ee_install_set_pyenv(
  py_path = "C:/Users/EndicottS/AppData/Local/r-miniconda/envs/r-reticulate/python.exe",
  py_env = "rgee"
)

# gives an error so using parts below
# rgee::ee_check()

earthengine_python <- Sys.getenv("EARTHENGINE_PYTHON", unset = NA)
if (!is.na(earthengine_python))
  Sys.setenv(RETICULATE_PYTHON = earthengine_python)

# ee_check_python caused the error because not character so doing separately
py_version <- py_discover_config()[["version"]] %>% as.character()
if (utils::compareVersion(py_version, "3.5") == -1) {
  stop("rgee needs Python 3.5 >=")
}

ee_check_python_packages()

reticulate::py_install('earthengine-api==0.1.370', envname='rgee')

# Initialize rgee session #=====================================================

# rgee function doesn't work
# rgee::ee_Initialize(user = "sarah.endicott@canada.ca", project = 'ee-sarahendicott-eccc')

# Can use this instead at the start of scripts using rgee
# recommended in this issue https://github.com/r-spatial/rgee/issues/355
ee$Authenticate(auth_mode='notebook')
ee$Initialize(project='ee-sarahendicott-eccc')

ee$String('Hello from the Earth Engine servers!')$getInfo()


