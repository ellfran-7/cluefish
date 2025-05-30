#> =============================================================================
#> DRomics pipeline: Dose-response analysis on published data (Gwinn et al. 2020)
#> =============================================================================
#> 
#> This script applies the DRomics pipeline to an external dataset:
#> 
#> Dataset background information:
#> - Authors: Gwinn WM, Auerbach SS, Parham F, Stout MD et al. 
#> - Title: Evaluation of 5-Day In Vivo Rat Liver and Kidney with High-Throughput Transcriptomics for Estimating Benchmark Doses of Apical Outcomes
#> - GEO accession: GSE132815
#> - Publish date: May 12, 2020
#> - Associated paper: Evaluation of 5-day In Vivo Rat Liver and Kidney With High-throughput Transcriptomics for Estimating Benchmark Doses of Apical Outcomes
#> - Journal: Toxicol Sci
#> - DOI: 10.1093/toxsci/kfaa081 (https://doi.org/10.1093/toxsci/kfaa081)
#> - PMID: 32492150 (https://www.ncbi.nlm.nih.gov/pubmed/32492150)
#>
#> Dataset specifics:
#> - Organism: Rattus norvegicus
#> - Contaminant: Perfluoro-octanoic acid (PFOA) 
#> - Experiment type: Expression profiling by high throughput sequencing
#> - Number of doses: control + 8
#> - Number of replicates: 4


# To reproduce the Cluefish workflow results as presented in Franklin et al. (submitted),
# please follow the steps below carefully. For a more detailed understanding of each phase—
# or if you're simply looking to apply the tool in your own work—refer to the DRomics vignette: 
# https://lbbe-software.github.io/DRomics/articles/DRomics_vignette.html.


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
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE147072&format=file&file=GSE147072%5FSRA%5F20%5FChemical%5Fmatrix%2Etxt%2Egz",
  destfile = "data/raw-data/GSE147072_SRA_20_Chemical_matrix.txt.gz"
)

# Uncompress the .gz file
R.utils::gunzip("data/raw-data/GSE147072_SRA_20_Chemical_matrix.txt.gz", 
                remove = FALSE)

# Download the series matrix file from the GEO repository.
gse_sm <- GEOquery::getGEO(GEO = "GSE147072") # Sometimes this codes download fails due to connection instability, so don't hesitate to try again

# Format the dataset information for appropriate sample distinction.
# Select the data specifics dataframe from the large list
gse147072_data <- gse_sm$GSE147072_series_matrix.txt.gz@phenoData@data

# Select columns and remove unecessary characters 
gse147072_info <- gse147072_data |> 
  dplyr::mutate(
    code = str_remove(title, pattern = "HSD male rat liver tissue |HSD male rat kidney tissue "),
    id = str_remove(characteristics_ch1.3, pattern = "animal id: "),
    chem = str_remove(characteristics_ch1.5, pattern = "chemical: "),
    dose = substr(characteristics_ch1.6, 14, stop = 25),
    organ = str_remove(characteristics_ch1.4, pattern = "tissue: ")) |> 
  select(code, id, chem, dose, organ)

# Save the GSE data information to CSV
write.csv(gse147072_info, file = "data/raw-data/GSE147072_series_info.csv", row.names = FALSE)




#> 2. Load the Datasets into R
#> ---------------------------

# Load the count dataset. This will be referred to as `rat_counts_df`.
rat_counts_df <- read.table(
  file = "data/raw-data/GSE147072_SRA_20_Chemical_matrix.txt",
  header = FALSE)

# Check the structure of the dataframe to ensure it is loaded correctly.
str(rat_counts_df)

# Load the information data for correspondance between sample information and count data.
rat_info_df <- read.csv(
  file = "data/raw-data/GSE147072_series_info.csv"
)

# Check the structure of the dataframe to ensure it is loaded correctly.
str(rat_info_df)



#> 3. Select the tissue, pollutant and prepare the data for DRomics 
#> ------------------------------------------------------------

# The GSExxxx dataset consists of over 20 tested chemicals on two different organs. 
# In our case, we aim to test one combination of chemical and organ, this being PFOA on the liver.
# We therefor need to subset the data and create the appropriate input for the DRomics pipeline.

# Filter out samples with "Plate2" in the code
rat_info_df <- rat_info_df |>
  dplyr::filter(!stringr::str_detect(code, "Plate2"))

# Display the contingency table of chemicals and organs
table(rat_info_df$chem, rat_info_df$organ)

# Define the contaminant of interest: PFOA
contaminant_select <- "PFOA" 

# Define the organ of interest: liver
organ_select <- "liver" 

# Print the number of unique doses for the selected contaminant
length(unique(rat_info_df$dose[rat_info_df$chem == contaminant_select]))

# Filter the information data to include only the selected contaminant and organ
rat_liver_pfoa_info_df <- rat_info_df |>
  dplyr::filter(chem == contaminant_select, organ == organ_select) |>
  dplyr::arrange(code)

# Get sample names from the first row
sample_names <- as.character(rat_counts_df[1, ])

# Create a logical vector of which columns to keep
columns_to_keep <- c(TRUE, sample_names[-1] %in% rat_liver_pfoa_info_df$code)

# Select the columns
rat_liver_pfoa_counts_df <- rat_counts_df[, columns_to_keep]

# Get the codes from the first row of filtered dataframe
selected_codes <- as.character(rat_liver_pfoa_counts_df[1, -1])

# Look up the corresponding doses
doses <- rat_liver_pfoa_info_df$dose[match(selected_codes, rat_liver_pfoa_info_df$code)]

# Replace the codes with doses in the first row
rat_liver_pfoa_counts_df[1, -1] <- doses

# Rename columns to be sequential
colnames(rat_liver_pfoa_counts_df) <- paste0("V", 1:ncol(rat_liver_pfoa_counts_df))

# Rename the first cell to "Item"
rat_liver_pfoa_counts_df[1,1] <- "item"

# Check the structure of the dataframe to ensure it is appropriate
str(rat_liver_pfoa_counts_df)

# The structure of the dataframe should look like this for the first five rows:

# 'data.frame':	2655 obs. of  37 variables:
# $ V1 : chr  "item" "A1BG_7930" "A2M_7932" "AADAC_7934" ...
# $ V2 : chr  "0" "1" "1524" "1249" ...
# $ V3 : chr  "0.3125" "0" "4" "396" ...
# $ V4 : chr  "1.25" "0" "6" "336" ...
# $ V5 : chr  "5" "1" "96" "1351" ...



#> 4. Save the prepared dataframe
#> ---------------------

write.table(rat_liver_pfoa_counts_df, 
            file = "data/derived-data/rat_pfoa_counts_df.txt")

# Reload the saved prepared text file to ensure it was saved correctly.
rat_liver_pfoa_counts_df <- read.table("data/derived-data/rat_pfoa_counts_df.txt")





# ============ SELECTION OF RESPONSIVE TRANSCRIPTS AND MODELLING ===============

#> -------------------------------------------------------------
#> Importation, checking and normalization of RNAseq count data 
#> -------------------------------------------------------------

set.seed(1234) # Fixing the seed to reproduce the results

(o <- DRomics::RNAseqdata(rat_liver_pfoa_counts_df, # or use the file directory
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
data4PCA <- list(dose = c(rep(0, 4), rep(0.156, 4), rep(0.3125, 4), rep(0.6250, 4), 
                          rep(1.25, 4), rep(2.5, 4), rep(5, 4), rep(10, 4), rep(20, 4)),
                 replicate = as.factor(rep(c("rep1", "rep2", "rep3", "rep4"), times = 9)))

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
saveRDS(f, file = "data/derived-data/fitres_rat_liver_pfoa.rds")




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
f <- readRDS(file = "data/derived-data/fitres_rat_liver_pfoa.rds")

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

system.time(b <- DRomics::bmdboot(r, niter = 5000, parallel = "snow", ncpus = 6))
# This can take quite some time depending on your computer performance...
# so grab a coffee or tea :)

#> output : object of class "bmdboot", a list with 3 components
head(b)
str(b)

#> We can now save the resulting bmdboot object to a file using the 'saveRDS'
#> (this is the data used for the workflow and for many representations at the end)
saveRDS(b, file = "data/derived-data/bootres_rat_liver_pfoa_seed1234_5000iter.rds")

#> Let's see which items have a calculated BMD, confidence interval etc.
b <- readRDS(file = "data/derived-data/bootres_rat_liver_pfoa_seed1234_5000iter.rds")

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

saveRDS(BMDres_definedCI, file = "data/derived-data/bootres_rat_liver_pfoa_seed1234_5000iter_definedCI.rds")
