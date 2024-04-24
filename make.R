#' bio_int_workflow: A workflow to alleviate biological interpretation of the dose-response (DR) modeling results of DR transcriptomic data
#' 
#' @description 
#' This project harbors the construction of a biological interpretation workflow for dose-response transcriptomic data. It incorporates various steps to enhance the understanding of dose-response modeling results. The main concept revolves around combining biological function annotations, gene regulation status, Protein-Protein Interaction Network (PPIN) analysis, cluster enrichment, cluster fusion, and lonely gene fishing to create a holistic view of the functional implications of omic data.
#' 
#' @author Ellis Franklin \email{ellis.franklin@univ-lorraine.fr}
#' 
#' @date 2024/02/29


## Install Dependencies (listed in DESCRIPTION) ----

devtools::install_deps(upgrade = "never")



## Load Project Addins (R Functions and Packages) ----

devtools::load_all(here::here())



## Run Project using separate files ----

# source(here::here("analyses", "proposed-approach", "0-download-data.R"))
# source(here::here("analyses", "proposed-approach", "1-load-data.R"))
# source(here::here("analyses", "proposed-approach", "2-get-ids.R"))
# source(here::here("analyses", "proposed-approach", "3-get-reg-annot.R"))
# source(here::here("analyses", "proposed-approach", "4-create-data-4-string.R"))
# source(here::here("analyses", "proposed-approach", "5-get-clustr-data.R"))
# source(here::here("analyses", "proposed-approach", "6-clustr-size-filtr.R"))
# source(here::here("analyses", "proposed-approach", "7-clustr-enrich.R"))
# source(here::here("analyses", "proposed-approach", "8-clustr-fusion.R"))
# source(here::here("analyses", "proposed-approach", "9-lonely-gene-fishing.R"))
# source(here::here("analyses", "proposed-approach", "10-generate-summary-table.R"))
# source(here::here("analyses", "proposed-approach", "11-curves2pdf.R"))




## Run project in the one file ----


#>> STEP 0 - Download TF and CoTF data 
#>------------------------------------

dl_regulation_data(
  url_tf = "https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Danio_rerio_TF",
  url_cof = "https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Danio_rerio_Cof",
  path = "data/derived-data/",
  filename_tf = "Danio_rerio_TF.txt",
  filename_cof = "Danio_rerio_Cof.txt",
  overwrite = TRUE
)







#>> STEP 1 - Load and filter the DRomics results
#>----------------------------------------------
# This workflow necessitates data from the experiment and the DRomics pipeline !
# This script will load the following data :
#   - DRomics workflow results : 
#           *object of class "drcfit" from the dose-response modelling of responsive transcripts (holds the background transcript list from the experiment)
#           *object of class "bmdboot" from the computation of CI on benchmark doses by bootstrap
#
# The "bmdboot" object is used throughout the workflow.
# The "drcfit" object is used to obtain the background transcript list and the tested doses for the "curves_to_pdf" function.
#
#
# All two files are to be created beforehand and are stored in `data/raw-data/`.

# Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = "data/raw-data/fitres_zebrafish_phtalate.rds")

# Load DRomics bmdboot object
b <- readRDS(file = "data/raw-data/bootres_zebrafish_phtalate_UF_seed3_5000iter.rds")

# We filter the bmdboot result by selecting only transcripts with a defined confidence interval around the BMD
BMDres_definedCI <- DRomics::bmdfilter(b$res, BMDfilter = "definedCI")







#>> STEP 2 - Retrieve deregulated gene identifiers from Ensembl
#>-------------------------------------------------------------

bg_t_ids <- getids(
  id_query = f$omicdata$item, 
  species_dataset = "drerio_gene_ensembl",
  id_filter = "ensembl_transcript_id_version", 
  id_attribut = c("ensembl_transcript_id_version",
                  "ensembl_gene_id", 
                  "external_gene_name")
  )

# Save the "time-consuming" data if already created
write.table(bg_t_ids, "outputs/bg_t_ids_2024_04_22.txt")
# Load the "time-consuming" data if already created
bg_t_ids <- read.table("outputs/bg_t_ids_2024_04_22.txt")







#>> STEP 3 - Retrieve regulatory status of deregulated genes
#>----------------------------------------------------------

# We only need the annotations for the deregulated transcripts/genes derived from the DRomics pipeline :
dr_t_ids <- bg_t_ids[bg_t_ids$ensembl_transcript_id_version %in% BMDres_definedCI$id,]

dr_t_regs <- getregs(
  getids_data = dr_t_ids,
  regulator_file = "data/derived-data/Danio_rerio_TF.txt",
  coregulator_file = "data/derived-data/Danio_rerio_Cof.txt")







#>> STEP 4 - Create and retrieve the clustered PPIN data from the StringApp in Cytoscape
#>--------------------------------------------------------------------------------------

# Create the data to be exported
DR_output4string <- merge(BMDres_definedCI, dr_t_regs, 
                          by.x = "id", by.y = "ensembl_transcript_id_version")

# Save the data 
write.table(DR_output4string, file = "outputs/DR_output4string_2024_04_22.txt", row.names = FALSE, sep = "\t")

# Once the clustered network is created, the resulting *.csv* files need to be stored in `outputs/`.

# The Ensembl gene IDs of all the responsive transcripts are manually exported into the StringApp plug-in found within the Cytoscape platform. They are queried with the parameters configured to include all types of interactions, a confidence score of 0.9 and no additional interactions. This creates an interaction network based on the known and predicted protein-protein associations data in the STRING database. Subsequently, to group the interaction network, we use clusterMaker2 to run Markov Chain Clustering (MCL). The inflation parameter is set as the default value (4.0) to reduce cluster size. To extract the relevant cluster information, the 'node table' containing the clustered elements is manually exported from the StringApp as a comma-separated values (CSV) file into the outputs/ folder. This exported table encompasses all the transcripts found within a cluster, their identifications, along with an appended cluster ID column. 

# This table, found in the "outputs/" folder, can then be imported back into the Rstudio environment in order to pursue the workflow : 

dr_t_clustrs <- getclustrs(
  getregs_data = dr_t_regs,
  path = "outputs/",
  nodetable_filename = "Resp_PPIN_clustered_220424.csv"
)







#>> STEP 5 - Filter clusters based on their gene set size
#>-------------------------------------------------------

dr_t_clustrs_filtr <- clustrfiltr(
  getclustrs_data = dr_t_clustrs,
  size_filtr = 3
)







#>> STEP 6 - Functional enrichment by clusters and annotation retrieval
#>---------------------------------------------------------------------

clustr_enrichres <- clustrenrich(
  clustrfiltr_data = dr_t_clustrs_filtr,
  dr_genes = dr_t_regs$ensembl_gene_id,
  bg_genes = bg_t_ids$ensembl_gene_id,,
  organism = "drerio",
  user_threshold = 0.05,
  correction_method = "fdr",
  enrich_size_filtr = TRUE,
  min_term_size = 5,
  max_term_size = 500,
  only_highlighted_GO = TRUE,
  ngenes_enrich_filtr = 3,
  path = "outputs/",
  output_filename = "clustr_enrichres_2024_04_22.rds",
  overwrite = FALSE
)







#>> STEP 7 - Fusion clusters based on shared cluster enrichment
#>-------------------------------------------------------------

clustr_fusionres <- clustrfusion(
  clustrenrich_data = clustr_enrichres
)







#>> STEP 8 - Fishing lonely genes sharing annotations with clusters enrichment
#>----------------------------------------------------------------------------

lonely_fishres <- lonelyfishing(
  dr_data = dr_t_regs,
  clustrenrich_data = clustr_enrichres,
  clustrfusion_data = clustr_fusionres,
  friendly_limit = 0,
  path = "outputs/",
  output_filename = "lonely_fishres_2024_04_22.rds",
  overwrite = TRUE
)







#>> STEP 9 - Generate a summary dataframe of the workflow
#>---------------------------------------------------

# With the workflow fulfilled, we can create a concise summary dataframe capturing the key details from the final results post the lonely gene fishing step. It encompasses all essential information for result exploration, striking a balance by avoiding an overwhelming amount of data that might hinder ease of exploration. 

results_to_csv(
  lonelyfishing_data = lonely_fishres,
  bmdboot_data = BMDres_definedCI,
  path = "outputs/",
  output_filename = "summary_workflow_2024_04_22.csv",
  overwrite = TRUE
)







#>> STEP 10 - Generate cluster-level curvesplots to a PDF file
#>------------------------------------------------------------

# Finally, this last step consists of generating the output PDF file containing a plot of dose-response curves for each cluster of genes, with each plot labeled with the cluster ID and the number of transcripts in that cluster. The curves are color-coded according to whether the trend is increasing, decreasing, U-shaped, or bell-shaped. The plot axes are labeled with "Dose (µg/L)" and "Signal", and the y-axis is scaled to be the same across all plots.
require(ggplot2)

curves_to_pdf(
  lonelyfishing_data = lonely_fishres,
  bmdboot_data = BMDres_definedCI, 
  clustrfusion_data = clustr_fusionres,
  tested_doses = unique(f$omicdata$dose), 
  colorby = "trend",
  addBMD = TRUE,
  scaling = TRUE,
  npoints= 100,
  free.y.scales = FALSE,
  xmin = 0.01, 
  xmax = 100, 
  dose_log_transfo = TRUE, 
  line.size = 0.7, 
  line.alpha = 0.4, 
  point.size = 2, 
  point.alpha = 0.4,
  xunit = "µg/L",
  xtitle = "Dose (µg/L)",
  ytitle = "Signal",
  colors = c("inc" = "#1B9E77", "dec" = "#D95F02", "U" = "#7570B3", "bell" = "#E7298A"),
  path = "outputs/",
  output_filename = "workflow_curvesplots_2024_04_22.pdf",
  overwrite = TRUE
)









