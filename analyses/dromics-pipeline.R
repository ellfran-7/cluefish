# Prior to commencing the workflows, it is imperative to  characterize 
# the dose-response transcriptomic data. In order to accomplish this, we use the 
# DRomics turn-key tool, which offers responsive item selection and modeling. 

#>> Load and list the packages needed
require(DRomics)
require(ggplot2)


# ---------------- SELECTION OF RESPONSIVE ITEMS AND MODELING --------------


#>> Importation, checking and normalization of RNAseq count data 
# ---------------------------------------------------

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

set.seed(1234) # Fixing the seed to reproduce te results

(o <- RNAseqdata(here::here("data", "raw-data", "Sal_RawCounts.txt"), 
                 transfo.method = 'rlog', 
                 round.counts = TRUE)) 

#> output : object of class "RNAseqdata", a list with 7 components
head(o) # check the first 5 rows of the output
str(o) # check the structure of the output

# One can plot of the output shows the distribution of the signal on all the 
# transcripts, for each sample, before and after the normalization and transformation 
# of data.
plot(o, cex.main = 0.8, col = "green")

# Just getting a global overview of the data : determine intra- and intergroup Sample variability   
# + allows us to see if there is a batch effect, shown by the position of control points (if grouped = good)

data4PCA <- list(dose = c(0, 0, 5, 5, 5, 10, 10, 10, 50, 50, 50, 100, 100, 100),
                 batch = as.factor(c("A", "B", rep(c("A", "B", "C"), times = 4))))

DRomics::PCAdataplot(o, batch = l$batch) + ggplot2::theme_bw()
# There doesn't seem to be any irregularities when viewing the positions of the control points for the two batches



#>> Selection of significantly responsive items
# ---------------------------------------------------

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


#>> Dose response modelling for responsive items
# ---------------------------------------------------

#> input : object of class "intemselect" returned by function 'itemselect'

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

system.time(f <- drcfit(s, parallel = "snow", ncpus = 4, information.criterion = "AICc")) 

#> output : object of class "drcfit", a list with 4 components
head(f)
str(f)

#> We can now save the resulting drcfit object to a file using the 'saveRDS'
#> (for futur representations of the raw data of target items with fitted curves
#> if available)
saveRDS(f, file = here::here("data", "raw-data", "fitres_zebrafish_phtalate.rds"))


#>> Computation of benchmark doses for responsive items
# ---------------------------------------------------

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

# Load the saved fitres
f <- readRDS(file = here::here("data", "raw-data", "fitres_zebrafish_phtalate.rds"))

(r <- bmdcalc(f, z = 1, x = 10)) 

#> output : object of class "bmdcalc", a list with 4 components
head(r)
str(r)
nrow(r$res)
bmdplot(r$res, colorby = "trend")



#>> Computation of confidence interval benchmark doses by bootstrap
# ------------------------------------------------------

#> input : object of class "bmdcalc" returned by 'bmdcalc'

#> parameter details : 
# - niter : the number of samples drawn by bootstrap
# - parallel : same as for 'drcfit' function
# - ncpus : same as for 'drcfit' function

#> parameter choice
# - niter = 1000 : default

set.seed(3) # Fixing graine, but not compatible with parallel computation
system.time(b <- bmdboot(r, niter = 5000, progressbar = TRUE))
# 100 iter -> 247 s = 4 min sur portable qui ne faisait que cela
# 1000 iter -> 5400 s pendant CSO LBBE sur portable
# 1000 iter -> 3100 s sur UF
# 5000 iter -> 15000 =  4h sur UF

#> output : object of class "bmdboot", a list with 3 components
head(b)
str(b)

#> We can now save the resulting bmdboot object to a file using the 'saveRDS'
#> (this is the data used for the workflow and for many representations at the end)
saveRDS(b, file = here::here("data", "raw-data", "bootres_zebrafish_phtalate_UF_seed3_5000iter.rds"))

#> Let's see which items have a calculated BMD, confidence interval etc.
b <- readRDS(file = here::here("data", "raw-data", "bootres_zebrafish_phtalate_UF_seed3_5000iter.rds"))

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

saveRDS(BMDres_definedCI, file = here::here("data", "raw-data", "bootres_zebrafish_phtalate_UF_seed3_5000iter_definedCI.rds"))