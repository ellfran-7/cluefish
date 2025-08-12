#> ======================================================================
#> DRomics pipeline: Reproducing results from Franklin et al. (2025)
#> ======================================================================
#> 
#> Case : dose-response anchoring data !!!!!
#> 
#> This script applies the DRomics pipeline to the individual endpoint data associated with the publication.


# To reproduce the Cluefish workflow results as presented in Franklin et al. (2025),
# please follow the steps below carefully. For a more detailed understanding of each phase—
# or if you're simply looking to apply the tool in your own work—refer to the DRomics vignette: 
# https://lbbe-software.github.io/DRomics/articles/DRomics_vignette.html.


#> 0. Load and list the packages needed
#> ---------------------------------
require(DRomics)
require(ggplot2)



#> 1. Save and load the excel spreadsheet supplementary file associated with the submitted paper
#> -------------------------------------------------------------------------------------

# This can be done manually, but let's suppose you have the csv saved, under the name "zebra_anchoring_dbp_data.csv" to the "data/raw-data/" directory then, you can load the file as such : 

zebra_anchoring_dbp_data <- read.csv2(file = "data/raw-data/apical/zebra_anchoring_dbp_data.csv")

# Check the structure of the dataframe to ensure it is loaded correctly.
str(zebra_anchoring_dbp_data)

# The dataframe should look like this for the first four rows and three columns :

# Condition	Eye surface (mm²)	Body length (mm)
# 0	0,061	3,38
# 0	0,058	3,18
# 0	0,063	3,17





#> 2. Format the data for DRomics
#> ------------------------------

# First we can rename the columns to avoid errors:
colnames(zebra_anchoring_dbp_data) <- c("condition", "eye_surface_mm2", "body_length_mm")

# Format the data for DRomics:

# Create a first column for the dataframe with the endpoints measured
col1 <- data.frame(V1 = as.factor(c("endpoint", colnames(zebra_anchoring_dbp_data)[-1])))

# Create the other columns for the dataframe consisting of the conditions (and replicates), the measurements for each endpoint
othercolumns <- as.data.frame(as.matrix(rbind(zebra_anchoring_dbp_data[,1], zebra_anchoring_dbp_data[,2], zebra_anchoring_dbp_data[,3])))

# Combine the two 
zebra_anchoring_4DRomics <- cbind(column1, othercolumns)

# Readjust the rownames, in case it is now sequential
rownames(zebra_anchoring_4DRomics) <- seq_len(nrow(zebra_anchoring_4DRomics))

# The formatted df should look like this now: 

#        V1    V2    V3    V4   V5
# 1    endpoint 0.000 0.000 0.000 0.00
# 2    eye_area 0.061 0.058 0.063 0.06
# 3 body_length 3.380 3.180 3.170 3.35




#> 3. Save the dataframe
#> ---------------------

write.table(zebra_anchoring_4DRomics, 
            file = "data/raw-data/zebra_anchoring_dbp_df.txt")

# Reload the saved file to ensure it was saved correctly.
zebra_anchoring_4DRomics <- read.table("data/raw-data/zebra_anchoring_dbp_df.txt")





# ============ SELECTION OF RESPONSIVE ENDPOINTS AND MODELLING ===============


#> --------------------------------------------------------------------
#> Importation, checking and normalization of endpoint measurement data 
#> --------------------------------------------------------------------

#> input : raw count data matrix with :
#   -> first row = dose 
#   -> first column = endpoint identification 

#> parameters details :
# - transfo.method : the method chosen to transform raw counts in log2 scale using
#                  DESeq2: "rlog" or "vst" (preferably "rlog" as default)
# - round.counts : put if TRUE if the counts come from Kallisto or Salmon to round
#                them before treatment 

#> parameter choice : 
# - transfo.method = "rlog" : performs shrinkage estimation for dispersions 
#                           to improve stability and interpretability of estimates
# - round.counts = TRUE : counts come from Salmon, so round counts


set.seed(1234) # Fixing the seed to reproduce te results

o <- DRomics::continuousanchoringdata(zebra_anchoring_4DRomics, backgrounddose = 0, check = TRUE)

print(o)
plot(o) + theme_bw()

#> output : object of class "continuousanchoringdata", a list with 7 components
head(o) # check the first 5 rows of the output
str(o) # check the structure of the output





#> -----------------------------------------------
#> Selection of significantly responsive endpoints
#> -----------------------------------------------

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

(s <- itemselect(o, select.method = "quadratic", FDR = 0.05))

#> output : object of class "itemselect", a list of 5 components
head(s)
str(s)




#> ------------------------------------------------
#> Dose response modelling for responsive endpoints
#> ------------------------------------------------

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
# - parallel = "snow" : my computer is Window
# - ncpus = 4 : my computer has 4 cpus

parallel::detectCores() # number of computer CPUs

system.time(f <- DRomics::drcfit(s, parallel = "snow", ncpus = 8, information.criterion = "AICc")) 

#> output : object of class "drcfit", a list with 4 components
head(f)
str(f)

#> We can now save the resulting drcfit object to a file using the 'saveRDS'
#> (for futur representations of the raw data of target items with fitted curves
#> if available)
saveRDS(f, file = "data/derived-data/apical/fitres_apical_zebra_dbp.rds")




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
# - z = 1 : défault (detailed by EFSA in 2017)
# - x = 10 : défault (detailed by EFSA in 2017))

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
# - niter = 1000 : default

set.seed(3) # Fixing seed, but not compatible with parallel computation
system.time(b <- DRomics::bmdboot(r, niter = 5000, progressbar = TRUE))

#> output : object of class "bmdboot", a list with 3 components
head(b)
str(b)

#> We can now save the resulting bmdboot object to a file using the 'saveRDS'
#> (this is the data used for the workflow and for many representations at the end)
saveRDS(b, file = "data/derived-data/apical/bootres_apical_zebra_dbp_UF_seed3_5000iter.rds")