#> ======================================================================
#> DRomics pipeline: Reproducing results from Franklin et al. (submitted)
#> ======================================================================
#> 
#> Case : dose-response transcriptomic data !!!!!
#> 
#> This script applies the DRomics pipeline to the dataset associated with the publication.
#> 
#> Dataset background information:
#> - Authors: Prud'homme SM, Sohm B, Cossu-Leguille C
#> - Title: 	Dose-dependent effects of dibutyl phthalate on embryonic development in zebrafish (Danio rerio): The value of integrating concentration gradients in transcriptomic analyses
#> - GEO accession: GSE283957
#> - Publish date: Dec 10, 2024
#> - Associated paper: (Submitted) Cluefish: mining the dark matter of transcriptional data series with over-representation analysis enhanced by aggregated biological prior knowledge 
#> - Journal: (Submitted) NAR Genomics and Bioinformatics
#> - DOI: bioRxiv > https://doi.org/10.1101/2024.12.18.627334
#>
#> Dataset specifics:
#> - Organism: Danio rerio
#> - Contaminant: di-butyl phthalate (DBP)
#> - Experiment type: Expression profiling by high throughput sequencing
#> - Number of doses: control + 4
#> - Number of replicates: 3


# If you're trying to reproduce the Cluefish workflow for the paper by Franklin et al. (submitted),
# follow these steps carefully. It's also recommended to visit the DRomics vignette for a detailed 
# understanding of each step (and if you want to simply use the tool): https://lbbe-software.github.io/DRomics/articles/DRomics_vignette.html.


#> 0. Load and list the packages needed
#> ---------------------------------
require(DRomics)
require(ggplot2)



#> 1. Download the Processed Dataset from GEO
#> ------------------------------------------

# Download the "GSE283957_raw_counts_All_Samples.csv" dataset from the GEO repository.
# This can be done manually or as follows:
download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE283957&format=file&file=GSE283957%5Fraw%5Fcounts%5FAll%5FSamples%2Ecsv%2Egz",
  destfile = "data/raw-data/raw_counts_All_Samples.csv.gz"
)

# Uncompress the .gz file
R.utils::gunzip("data/raw-data/raw_counts_All_Samples.csv.gz", 
                remove = FALSE)




#> 2. Load the Dataset into R
#> ---------------------------

# Load the downloaded CSV file into an R dataframe. This will be referred to as `zebra_dbp_counts_df`.

zebra_dbp_counts_df <- read.csv2(file = "data/raw-data/raw_counts_All_Samples.csv",
                          header = FALSE)

# Check the structure of the dataframe to ensure it is loaded correctly.
str(zebra_dbp_counts_df)


# The structure of the dataframe should look like this for the first five rows:

# 'data.frame':	46489 obs. of  16 variables:
# $ V1 : chr  "Item" "ENSDART00000000004.5" "ENSDART00000000005.7" "ENSDART00000000042.11" ...
# $ V2 : chr  "0B" "119" "385" "0" ...
# $ V3 : chr  "0D" "117" "553" "0" ...
# $ V4 : chr  "0E" "134" "572" "2.644" ...
# $ V5 : chr  "5A" "107" "548" "0" ...



#> 3. Uniformize the Dose Names (Remove Non-Numeric Characters) 
#> ------------------------------------------------------------

# Each dose has a replicate distinction, which GEO requires to be identifiable.
# However, for the DRomics pipeline, we are only interested in the dose levels themselves (e.g., 0, 5, 10, etc.), and we don't need the replicate-specific information (e.g., "0B", "5A", "5B", etc.).
# We need to remove these replicate-specific characters (like 'A', 'B', 'C') from the dose labels.

# Remove non-numeric characters from the first row (except the first column).
zebra_dbp_counts_df[1, -1] <- gsub("[^0-9]", "", zebra_dbp_counts_df[1, -1])

# Check the structure of the dataframe after the transformation.
str(zebra_dbp_counts_df)



#> 4. Save the dataframe
#> ---------------------

write.table(zebra_dbp_counts_df, 
            file = "data/derived-data/zebra_dbp_counts_df.txt")

# Reload the saved file to ensure it was saved correctly.
zebra_dbp_counts_df <- read.table("data/derived-data/zebra_dbp_counts_df.txt")



# ============ SELECTION OF RESPONSIVE TRANSCRIPTS AND MODELLING ===============


#> -------------------------------------------------------------
#> Importation, checking and normalization of RNAseq count data 
#> -------------------------------------------------------------

#> input : raw count data matrix with :
#   -> first row = dose 
#   -> first column = item identification (in our case Ensembl Transcript IDs)

#> parameters details :
# - transfo.method : the method chosen to transform raw counts in log2 scale using
#                  DESeq2: "rlog" or "vst" (preferably "rlog" as default)
# - round.counts : put if TRUE if the counts come from Kallisto or Salmon to round
#                them before treatment 

#> parameter choice : 
# - transfo.method = "rlog" : performs shrinkage estimation for dispersions 
#                           to improve stability and interpretability of estimates
# - round.counts = TRUE : counts come from Salmon, so round counts

set.seed(1234) # Fixing the seed to reproduce the results

(o <- DRomics::RNAseqdata(zebra_dbp_counts_df, # or use the file directory
                          transfo.method = 'rlog', 
                          round.counts = TRUE)) 

#> output : object of class "RNAseqdata", a list with 7 components
head(o) 
str(o) 

# Plot the distribution of signal across all transcripts for each sample before and after normalization.
plot(o, cex.main = 0.8, col = "green")



#> Perform a Principal Component Analysis (PCA) 
#> --------------------------------------------

# Prepare metadata for PCA plot.
data4PCA <- list(dose = c(0, 0, 0, 5, 5, 5, 10, 10, 10, 50, 50, 50, 100, 100, 100),
                 replicate = as.factor(rep(c("rep1", "rep2", "rep3"), times = 5)))

# Generate the PCA plot
DRomics::PCAdataplot(o, batch = data4PCA$replicate) + 
  ggplot2::theme_bw()

# After examining the PCA plot, we observe that the control dose (0D) appears to be an outlier, which may affect the analysis. 
# To address this, we can remove the 0D sample from the dataset and re-run the DRomics steps.




#> Removing outliers and re-performing DRomics step 1
#> --------------------------------------------------

# Remove the 0B replicate (the second column) from the dataframe.
zebra_dbp_counts_df_rm0 <- zebra_dbp_counts_df[, -2]

# Reset column names to V1, V2, V3, etc. after removing the column.
colnames(zebra_dbp_counts_df_rm0) <- paste0("V", seq_along(colnames(zebra_dbp_counts_df_rm0)))

(o <- DRomics::RNAseqdata(zebra_dbp_counts_df_rm0,  
                          transfo.method = 'rlog', 
                          round.counts = TRUE)
)

head(o)
str(o)

plot(o, cex.main = 0.8, col = "green")

# Perform PCA again with the updated dataset to check the impact of removing the outlier.
data4PCA <- list(
  dose = c(0, 0, 5, 5, 5, 10, 10, 10, 50, 50, 50, 100, 100, 100),  # Updated dose information
  replicate = as.factor(c("rep1", "rep2", rep(c("rep1", "rep2", "rep3"), times = 4)))  # Updated replicate information
)

# Generate the PCA plot after outlier removal.
DRomics::PCAdataplot(o, batch = data4PCA$replicate) + ggplot2::theme_bw()

# Now that the outlier is removed, there doesn't seem to be any irregularities when viewing the positions of the control points for the two replicates.




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
saveRDS(f, file = "data/derived-data/fitres_zebrafish_phtalate.rds")




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
f <- readRDS(file = "data/derived-data/fitres_zebrafish_phtalate.rds")

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
saveRDS(b, file = "data/derived-data/bootres_zebrafish_phtalate_UF_seed3_5000iter.rds")

#> Let's see which items have a calculated BMD, confidence interval etc.
b <- readRDS(file = "data/derived-data/bootres_zebrafish_phtalate_UF_seed3_5000iter.rds")

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

saveRDS(BMDres_definedCI, file = "data/derived-data/bootres_zebrafish_phtalate_UF_seed3_5000iter_definedCI.rds")