# STEP 1 - Load project data
# -----------------------------------------------------------------------------
# 
# This workflow necessitates data from the experiment and the DRomics pipeline !
# This script will load the following data :
#   - DRomics workflow results : 
#           *object of class "drcfit" from the dose-response modelling of responsive transcripts (holds the background transcript list from the experiment)
#           *object of class "bmdboot" from the computation of CI on benchmark doses by bootstrap
#
# The "bmdboot" object is used throughout the workflow.
# The "drcfit" object is used to obtain the background transcript list and the tested doses for the "curves2pdf" function.
#
#
# All two files are stored in `data/derived-data/`.

# Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = "data/raw-data/zebra_fitted_05")

# Load DRomics bmdboot object
b <- readRDS(file = "data/raw-data/zebra_boots_05")

# We filter the bmdboot result by selecting only items with a defined confidence interval around the BMD
BMDres_definedCI <- bmdfilter(b$res, BMDfilter = "definedCI")
