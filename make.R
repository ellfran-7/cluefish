#' bio_int_workflow: A workflow to alleviate biological interpretation of the dose-response (DR) modeling results of DR transcriptomic data
#' 
#' @description 
#' This project harbors the construction of a biological interpretation workflow for dose-response transcriptomic data. It incorporates various steps to enhance the understanding of dose-response modeling results. The main concept revolves around combining biological function annotations, gene regulation status, Protein-Protein Interaction Network (PPIN) analysis, cluster enrichment, cluster fusion, and lonely gene fishing to create a holistic view of the functional implications of omic data. Each step in the workflow builds on the results of the previous steps, although some steps can be performed independently. Additionally, in some cases, certain steps are not specifically required to proceed to the next phase.
#' 
#' @author Ellis Franklin \email{ellis.franklin@univ-lorraine.fr}
#' 
#' @date 2024/07/07



## Install Dependencies (listed in DESCRIPTION) ----

devtools::install_deps(upgrade = "never")



## Load Project Addins (R Functions and Packages) ----

devtools::load_all(here::here())


## State the Time Variable for file saving and reading

file_date <- "2024-07-07"



## Run Project Workflow  ----


#>> STEP 0 - Download TF and CoTF Data 
#>------------------------------------

dl_regulation_data(
  url_tf = "https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Danio_rerio_TF",
  url_cof = "https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Danio_rerio_Cof",
  path = "data/derived-data/",
  filename_tf = paste0("Danio_rerio_TF_", file_date, ".txt"),
  filename_cof = paste0("Danio_rerio_Cof_", file_date, ".txt"),
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

# First, find the names of the BioMart services that Ensembl is currently providing.
# Use the listEnsembl() (or one of the following) function to display all available Ensembl BioMart web services.
biomaRt::listMarts()
biomaRt::listEnsembl()

# Annotation from the vertebrate genomes is provided by the main Ensembl project across taxonomic space, with separate BioMart interfaces for Protists, Plants, Metazoa and Fungi.
biomaRt::listEnsemblGenomes()

# Next, connect to the desired BioMart database using the useMart(), useEnsembl(), or useEnsemblGenomes() function.
# The biomart argument should be a valid name from the output of one of the previous functions.
ensembl <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL")

# BioMart databases can contain several datasets. For instance, within the Ensembl genes mart, each species is a different dataset.
biomaRt::listDatasets(mart = ensembl)

# Now, run the getids() function with the correct input based on the organism of the study.
# This retrieves IDs using the specified parameters for the query.bg_t_ids <- getids(
  id_query = f$omicdata$item, 
  biomart_db = "ENSEMBL_MART_ENSEMBL",
  species_dataset = "drerio_gene_ensembl",
  transcript_id = "ensembl_transcript_id_version",
  gene_id = "ensembl_gene_id",
  gene_name = "external_gene_name"
)

# Note: The output dataframe must include identifiers supporting the organism:
#               *in the STRING db to create Protein-Protein Interaction Networks
#               *in the g:profiler db to perform functional enrichment
#       For example, with Danio rerio, the "ensembl_gene_id" identifier is supported in STRING and g:profiler.

# Save the "time-consuming" data if already created
write.table(bg_t_ids, paste0("outputs/bg_t_ids_", file_date, ".txt"))
# Load the "time-consuming" data if already created
bg_t_ids <- read.table(paste0("outputs/bg_t_ids_", file_date, ".txt"))

# The "gene_id" from the background gene list (bg_t_ids) is only needed for function enrichment. However, the "gene_id" from the deregulated transcripts (DRomics pipeline) is needed for the whole workflow, including creating a STRING PPI network and function enrichment. Therefore, we need to subset the bg_t_ids dataframe.
dr_t_ids <- bg_t_ids[bg_t_ids$transcript_id %in% BMDres_definedCI$id,]





#>> STEP 3 - Retrieve regulatory status of deregulated genes
#>----------------------------------------------------------

dr_t_regs <- getregs(
  getids_data = dr_t_ids,
  regulator_file = paste0("data/derived-data/Danio_rerio_TF_", file_date, ".txt"),
  coregulator_file = paste0("data/derived-data/Danio_rerio_Cof_", file_date, ".txt"))







#>> STEP 4 - Create and retrieve the clustered PPIN data from the StringApp in Cytoscape
#>--------------------------------------------------------------------------------------

# Create the data to be exported into Cytoscape
DR_output4string <- merge(BMDres_definedCI, dr_t_regs, 
                          by.x = "id", by.y = "transcript_id")

# Save the data 
write.table(DR_output4string, file = paste0("outputs/DR_output4string_", file_date, ".txt"), row.names = FALSE, sep = "\t")

# Once the clustered network is created, the resulting *.csv* files need to be stored in `outputs/`.

# The Ensembl gene IDs of all the responsive transcripts are manually exported into the StringApp plug-in found within the Cytoscape platform. They are queried with the parameters configured to include all types of interactions, a confidence score of 0.9 and no additional interactions. This creates an interaction network based on the known and predicted protein-protein associations data in the STRING database. Subsequently, to group the interaction network, we use clusterMaker2 to run Markov Chain Clustering (MCL). The inflation parameter is set as the default value (4.0) to reduce cluster size. To extract the relevant cluster information, the 'node table' containing the clustered elements is manually exported from the StringApp as a comma-separated values (CSV) file into the outputs/ folder. This exported table encompasses all the transcripts found within a cluster, their identifications, along with an appended cluster ID column. 


# This table, found in the "outputs/" folder, can then be imported back into the Rstudio environment in order to pursue the workflow : 

dr_t_clustrs <- getclustrs(
  gene_data = dr_t_regs,
  colname_for_merge = "gene_id",
  path = "outputs/cytoscape-files/",
  nodetable_filename = paste0("Resp_PPIN_clustered_cs09_mcl4_", file_date, ".csv")
)





#>> STEP 5 - Filter clusters based on their gene set size
#>-------------------------------------------------------

dr_t_clustrs_filtr <- clustrfiltr(
  getclustrs_data = dr_t_clustrs,
  size_filtr = 4
)






#>> STEP 6 - Functional enrichment by clusters and annotation retrieval
#>---------------------------------------------------------------------

clustr_enrichres <- clustrenrich(
  clustrfiltr_data = dr_t_clustrs_filtr,
  dr_genes = dr_t_regs$gene_id,
  bg_genes = bg_t_ids$gene_id,
  bg_type = "custom_annotated",
  sources = c("GO:BP", "KEGG", "WP"), 
  organism = "drerio",
  user_threshold = 0.05,
  correction_method = "fdr",
  exclude_iea = FALSE, 
  enrich_size_filtr = TRUE,
  min_term_size = 5,
  max_term_size = 500,
  only_highlighted_GO = TRUE,
  ngenes_enrich_filtr = 3,
  path = "outputs/cs09-cf4/",
  output_filename = paste0("clustr_enrichres_cs09_cf4_", file_date, ".rds"),
  overwrite = FALSE
)






#>> STEP 7 - Fusion clusters based on shared cluster enrichment
#>-------------------------------------------------------------

clustr_fusionres <- clustrfusion(
  clustrenrich_data = clustr_enrichres,
  monoterm_fusion = FALSE
)







#>> STEP 8 - Fishing lonely genes sharing annotations with clusters enrichment
#>----------------------------------------------------------------------------

lonely_fishres <- lonelyfishing(
  dr_data = dr_t_regs,
  clustrenrich_data = clustr_enrichres,
  clustrfusion_data = clustr_fusionres,
  friendly_limit = 0,
  path = "outputs/cs09-cf4/",
  output_filename = paste0("lonely_fishres_cs09_cf4_", file_date, ".rds"), 
  overwrite = FALSE
)




#>> STEP 9 - Generate a summary dataframe of the workflow
#>---------------------------------------------------

# With the workflow fulfilled, we can create a concise summary dataframe capturing the key details from the final results post the lonely gene fishing step. It encompasses all essential information for result exploration, striking a balance by avoiding an overwhelming amount of data that might hinder ease of exploration. 

results_to_csv(
  lonelyfishing_data = lonely_fishres,
  bmdboot_data = BMDres_definedCI,
  path = "outputs/cs09-cf4/",
  output_filename = paste0("summary_workflow_cs09_cf4_", file_date, ".csv"),
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
  annot_order = c("GO:BP", "KEGG", "WP"),
  colorby = "trend",
  addBMD = TRUE,
  scaling = TRUE,
  npoints= 100,
  free.y.scales = FALSE,
  xmin = 0.1, 
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
  path = "outputs/cs09_cf4/",
  output_filename = paste0("workflow_curvesplots_cs09_cf4_", file_date, ".pdf"),
  overwrite = TRUE
)





#>> STEP 11 - Generate the quarto report 
#>--------------------------------------

# # Define the date, corresponding to the date that the files of the workflow are saved. This only functions if all the outputs have the same date identifiers in the filename (e.g. lonely_fishres_2024-07-07)(you can get todays date dynamically, e.g., from file_date)
# new_date <- "2024-07-07"  
# 
# # Render and preview the html report in the Viewer panel
# quarto::quarto_render(
#   input = here::here("analyses", "quarto", "workflow_results_report.qmd"),
#   execute_params = list(`file-date` = new_date),
#   output_file = "workflow_results_report_sc09_cf4_2024-07-07.html",
#   output_format = "html"
# )



#>> Additional steps 
#> -----------------

# # Characterisation of the lonely cluster: basic functional enrichment 
# source(here::here("analyses", "lonely_cluster_analysis", "lonely_cluster_analysis.R"))
# 
# # Render and preview the lonely_results_report html report contextualising the lonely cluster
# quarto::quarto_render(
#   input = here::here("analyses", "quarto", "lonely_results_report.qmd"),
#   execute_params = list(`file-date` = new_date),
#   output_file = "lonely_results_report_sc09_cf4_2024-07-07.html",
#   output_format = "html"
# )
# 
# # Basic enrichment of the deregulated transcripts genes from the DRomics workflow, for comparison with the proposed workflow
# source(here::here("analyses", "standard_approach", "standard_pipeline.R"))
# 
# # Render and preview the comparison_results_report html report comparing the results of both approaches on the same data
# quarto::quarto_render(
#   input = here::here("analyses", "quarto", "comparison_results_report.qmd"),
#   execute_params = list(`file-date` = new_date),
#   output_file = "comparison_results_report_sc09_cf4_2024-07-07.html",
#   output_format = "html"
# )



