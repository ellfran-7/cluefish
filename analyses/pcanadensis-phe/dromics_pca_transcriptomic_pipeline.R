#> =============================================================================
#> DRomics pipeline: Dose-response analysis on published data (Gréau et al. 2025)
#> =============================================================================
#> 
#> This script applies the DRomics pipeline to an external dataset:
#> 
#> Dataset background information:
#> - Authors: Gréaut L, Blaudez D, Le Jean M, Gallois N, Paysant Le-Roux C, Huguet S, Beguiristain T, Billoir E, Cébron A 
#> - Title: Transcriptomics highlights dose‑dependent response of poplar to a phenanthrene contamination
#> - GEO accession: GSE263776
#> - Publish date: April 08, 2025
#> - Associated paper: Transcriptomics highlights dose‑dependent response of poplar to a phenanthrene contamination
#> - Journal: Environmental Science and Pollution Research
#> - DOI: 10.1007/s11356-025-36002-5 (https://doi.org/10.1007/s11356-025-36002-5)
#>
#> Dataset specifics:
#> - Organism: Populus canadensis
#> - Contaminant: Phenanthrene (PHE) 
#> - Experiment type: Expression profiling by high throughput sequencing
#> - Number of doses: control + 7
#> - Number of replicates: 4


# If you're trying to reproduce the Cluefish workflow for the paper by Franklin et al. (submitted),
# follow these steps carefully. It's also recommended to visit the DRomics vignette for a detailed 
# understanding of each step (and if you want to simply use the tool): https://lbbe-software.github.io/DRomics/articles/DRomics_vignette.html.


#> 0. Load and list the packages needed
#> ---------------------------------
require(DRomics)
require(ggplot2)
require(stringr)
require(readr)
require(dplyr)
require(GEOquery)


#> 1. Download the Processed Dataset from GEO
#> ------------------------------------------

# Download the count dataset from the GEO repository.
# This can be done manually or as follows:
download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/name-of-file",
  destfile = "data/raw-data/name-of-file-and-extension.gz"
)

# Uncompress the .gz file
R.utils::gunzip("data/raw-data/name-of-file-and-extension.gz", 
                remove = FALSE)



#> 2. Load the Datasets into R
#> ---------------------------

# Load the count dataset. This will be referred to as `pop_counts_df`.
pop_counts_df <- read.table( # or read.csv2()
  file = "data/raw-data/name-of-file-and-extension",
  header = FALSE)

# Check the structure of the dataframe to ensure it is loaded correctly.
str(pop_counts_df)

# =================== TEMPORARY =======================================

require(readxl)

pop_roots_counts_df <- read_excel(path = "data/raw-data/Results_NGS2021_12_EndOMIX_Raw_COUNTS_F_R.xlsx", 
                            sheet = "Roots")

# =====================================================================

pop_roots_counts_df <- as.data.frame(pop_roots_counts_df)

# Check the structure of the dataframe to ensure it is appropriate
str(pop_roots_counts_df)


# Convert the column names from "R_2000PHE_rep2" to "2000"
colnames <- names(pop_roots_counts_df)

for (i in 1:length(colnames)) {
  if (grepl("^R_\\d+PHE", colnames[i])) {
    colnames[i] <- gsub("R_(\\d+)PHE.*", "\\1", colnames[i])
  } else {
    colnames[i] <- "item"
  }
}

# Add the modified column names as the first row of data
pop_roots_counts_df <- rbind(colnames, pop_roots_counts_df)

# Rename all columns to a sequential pattern V1, V2, ..., VX
names(pop_roots_counts_df) <- paste0("V", 1:ncol(pop_roots_counts_df))

# Reset row names
rownames(pop_roots_counts_df) <- NULL

# Check the structure of the dataframe to ensure it is appropriate
str(pop_roots_counts_df)

# The structure of the dataframe should look like this for the first five rows:

# 'data.frame':	34700 obs. of  33 variables:
# $ V1 : chr  "item" "Potri.001G000400" "Potri.001G000450" "Potri.001G000500" ...
# $ V2 : chr  "0" "730" "31" "0" ...
# $ V3 : chr  "0" "796" "18" "0" ...
# $ V4 : chr  "0" "941" "47" "0" ...
# $ V5 : chr  "0" "652" "27" "0" ...
# $ V6 : chr  "100" "742" "30" "0" ...



#> 4. Save the prepared dataframe
#> ---------------------

write.table(pop_roots_counts_df, 
            file = "data/derived-data/pop_phe_counts_df.txt")

# Reload the saved prepared text file to ensure it was saved correctly.
pop_roots_counts_df <- read.table("data/derived-data/pop_phe_counts_df.txt")





# ============ SELECTION OF RESPONSIVE TRANSCRIPTS AND MODELLING ===============

#> -------------------------------------------------------------
#> Importation, checking and normalization of RNAseq count data 
#> -------------------------------------------------------------

set.seed(1234) # Fixing the seed to reproduce the results

(o <- DRomics::RNAseqdata(pop_roots_counts_df, # or use the file directory
                          transfo.method = 'rlog', 
                          round.counts = FALSE)) 

#> output : object of class "RNAseqdata", a list with 7 components
head(o) 
str(o) 

# Plot the distribution of signal across all transcripts for each sample before and after normalization.
plot(o, cex.main = 0.8, col = "green")



#> Perform a Principal Component Analysis (PCA) 
#> --------------------------------------------

# Prepare metadata for PCA plot.
data4PCA <- list(dose = c(rep(0, 4), rep(100, 4), rep(400, 4), rep(700, 4), 
                          rep(1000, 4), rep(1500, 4), rep(2000, 4)),
                 replicate = as.factor(rep(c("rep1", "rep2", "rep3", "rep4"), times = 8)))

# Generate the PCA plot
DRomics::PCAdataplot(o, batch = data4PCA$replicate) + 
  ggplot2::theme_bw()

# After examining the PCA plot, everything seems to be in order



#> --------------------------------------------
#> Selection of significantly responsive items
#> --------------------------------------------

#> input : object of class "RNAseqdata" (or "microarraydata", "metabolomicdata",
#                                       "continuousanchoringdata")

#> parameter details :
# - select.method : "quadratic" for a quadratic trend test on dose ranks 
#                 "linear" for a linear trend test on dose ranks
#                 "ANOVA" for an ANOVA-type test
# - FDR : the threshold in term of False Discovery Rate for selection responsive items

#> parameter choice : 
# - select.method = "quadratic" : detecting monotonic and biphasic trends 
# - FDR = 0.05 (default)

(s <- DRomics::itemselect(o, select.method = "quadratic", FDR = 0.05))

#> output : object of class "itemselect", a list of 5 components
head(s)
str(s)



#> ---------------------------------------------
#> Dose response modelling for responsive items
#> ---------------------------------------------

#> input : object of class "itemselect" returned by function 'itemselect'

#> parameter details :
# - information.criterion : the info criterion used to select the best fit model,
#                       "AICc" ; the corrected version of the AIC 
# - parallel : type of parallel operation to be used, "snow" (Windows/Linux) or 
#           "multicore" (Linux)
# - ncpus : number of processes to be used in the parallel operation: typically
#         fixed to the number of available CPUs on the computer (apply detectCores()
#         to know the number of CPUs on computer)
# 

#> parameter choice : 
# - information.criterion = "AICc" : recommended and default choice for small samples
# - parallel = "snow" : allow parallel computing for high number of responsive items (snow because my computer is Windows)
# - ncpus = 4 : my computer has 8 cpus, so anything less is appropriate

parallel::detectCores() # number of computer CPUs

system.time(f <- DRomics::drcfit(s, parallel = "snow", ncpus = 4, information.criterion = "AICc")) 

#> output : object of class "drcfit", a list with 4 components
head(f)
str(f)

#> We can now save the resulting drcfit object to a file using the 'saveRDS'
#> (for futur representations of the raw data of target items with fitted curves
#> if available)
saveRDS(f, file = "data/derived-data/fitres_pop_root_phe.rds")




#> ----------------------------------------------------
#> Computation of benchmark doses for responsive items
#> ----------------------------------------------------

#> input : object of class "drcfit" returned by the function 'drcfit'

#> parameter details :
# - z : value of z defining the BMD-zSD as the dose at which the response is reaching
#     y0 +/- z*SD, with y0 the level at the control given by the dose-response 
#     fitted model and SD the residual standard deviation of the dose-response fitted model
# - x : value of x given as a percentage and defining the BMD-xfold as the dose at
#     which the response is reaching y0 +/- (x/100) * y0, with y0 the level at the
#     control given by the dose-response fitted model

#> parameter choice : 
# - z = 1 : default (detailed by EFSA in 2017)
# - x = 10 : default (detailed by EFSA in 2017))

# Load the saved fitres
f <- readRDS(file = "data/derived-data/fitres_pop_root_phe.rds")

(r <- DRomics::bmdcalc(f, z = 1, x = 10)) 

#> output : object of class "bmdcalc", a list with 4 components
head(r)
str(r)
nrow(r$res)
bmdplot(r$res, colorby = "trend")




#> ----------------------------------------------------------------
#> Computation of confidence interval benchmark doses by bootstrap
#> ----------------------------------------------------------------

#> input : object of class "bmdcalc" returned by 'bmdcalc'

#> parameter details : 
# - niter : the number of samples drawn by bootstrap
# - parallel : same as for 'drcfit' function
# - ncpus : same as for 'drcfit' function

#> parameter choice
# - niter = 5000 : high number of samples drawn by bootstrap
# - parallel = "snow" : same as for 'drcfit' function
# - ncpus = 6 : this step is more computationally intensive, so using more cpus is required

system.time(b <- DRomics::bmdboot(r, niter = 5000, parallel = "snow", ncpus = 4))
# This can take quite some time depending on your computer performance...
# so grab a coffee or tea :)

#> output : object of class "bmdboot", a list with 3 components
head(b)
str(b)

#> We can now save the resulting bmdboot object to a file using the 'saveRDS'
#> (this is the data used for the workflow and for many representations at the end)
saveRDS(b, file = "data/derived-data/bootres_pop_root_phe_seed1234_5000iter.rds")

#> Let's see which items have a calculated BMD, confidence interval etc.
b <- readRDS(file = "data/derived-data/bootres_pop_root_phe_seed1234_5000iter.rds")

BMDres <- b$res
nrow(BMDres)

# We can filter the `bmdboot` results to retain only those transcripts with a defined confidence interval around the BMD.
# The `DRomics::bmdfilter` function provides two other filtering options based on the desired stringency:
# 
# - **"finiteCI"**: Retains transcripts where both point and interval estimates of the BMD were successfully calculated and fall within the range of tested or observed doses.
# - **"definedBMD"**: Retains transcripts where the point estimate of the BMD falls within the range of tested or observed doses.
# 
# Choose the appropriate filter based on the level of stringency required for your analysis.
BMDres_definedCI <- DRomics::bmdfilter(BMDres, BMDfilter = "definedCI")
nrow(BMDres_definedCI)

saveRDS(BMDres_definedCI, file = "data/derived-data/bootres_pop_root_phe_seed1234_5000iter_definedCI.rds")
