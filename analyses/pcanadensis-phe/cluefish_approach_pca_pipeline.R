#> =============================================================================
#> Cluefish: A workflow to optimise biological interpretation of transcriptomic data series
#> =============================================================================

#> This script applies the cluefish approach pipeline for biological interpretation 
#> of results of the DRomics analysis of the poplar root transcriptomic dataset
#> 
#> NOTE: This pipeline is specifically configured for the poplar root dataset.
#> 
#> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


## Install packages listed in DESCRIPTION (and/or R and Rmd files) ----
renv::install()

# OR --

## Restore the local environment to ensure reproducibility of the analysis with the same package versions and states ---
renv::restore()


## Load Project Addins (R Functions and Packages) ----
devtools::load_all(here::here())


## State the Time Variable for file saving and reading ---
file_date <- "2025-04-06"


## Create directory for saving Cluefish output files (if it doesn't already exist) ---
dir_path <- paste0("outputs/", file_date)

if (!dir.exists(dir_path)) { # Check if the directory path exists
  
  dir.create(dir_path) # Create it if not
  
}


## Create directory for saving Cytoscape files (if it doesn't already exist) ---
dir_path <- paste0("outputs/", file_date, "/cytoscape-files")

if (!dir.exists(dir_path)) { # Check if the directory path exists
  
  dir.create(dir_path) # Create it if not
  
}


## Run Project Workflow  ----


#>> STEP 1 - Download TF and CoTF Data 
#>------------------------------------

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Step 1 is not performed as the species is not referenced in the knowledgebase
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#>> STEP 2 - Load the DRomics results
#>----------------------------------------------

# This workflow requires data from both the experimental results and the DRomics pipeline.
#
# The following code loads the necessary data:
#   - **DRomics workflow results:**
#     - An object of class `"drcfit"` from the dose-response modeling of responsive transcripts. 
#       °This object contains:
#       - The background transcript list (accessible via `f$omicdata$item`), which is used by the `getids()` and `clustrenrich()` functions.
#       - The tested doses (accessible via `f$omicdata$dose`), which are required for the `curves_to_pdf()` function.
#     - The dataframe of results provided by "bmdboot" (res). TThis dataframe may have been filtered based on the quality of BMD estimates using DRomics::bmdfilter(b$res, BMDfilter = "definedCI") for example.
#       °This dataframe provides the deregulated transcript data, which is used throughout the workflow.
#
# Both files must be created in advance and can be stored in the `data/derived-data/` directory.
# If these files are not already available, you can generate them using the `dromics_dre_transcriptomic_pipeline.R` script found in the `analyses/` folder.


# Load DRomics "drcfit" object
f <- readRDS(file = "data/derived-data/fitres_pop_root_phe.rds")

# Load DRomics "bmdboot" results, filtered to include only transcripts with a defined confidence interval around the BMD.
b_definedCI <- readRDS(file = "data/derived-data/bootres_pop_root_phe_seed1234_5000iter_definedCI.rds")






#>> STEP 3 - Retrieve deregulated gene identifiers from Ensembl
#>-------------------------------------------------------------

# During pre-processing in this study, reads were mapped to the "Populus trichocarpa" genome (V4, Phytozome). Therefor we are searching for "Populus trichocarpa" gene identifiers.

# Nevertheless, there element missing from the identifiers. The Ensembl stable gene IDs of Populus trichocarpa end in ".v4.1" suffix to the gene IDs. We need to add this back.

# Create a vector of all the genes from the experiment but adding the suffixe ".v4.1"
f_id_mod <- paste0(f$omicdata$item, ".v4.1")

# Create an additional column in the DRomics output that holds the deregulated genes with the added suffixe ".v4.1"
b_definedCI_mod <- b_definedCI |> 
  dplyr::mutate(id_mod = paste0(id, ".v4.1"))

# Create the only_IDs version for merging
b_definedCI_mod_onlyids <- b_definedCI_mod |> 
  dplyr::select(id, id_mod)

# First, find the names of the BioMart services that Ensembl is currently providing.
# Use the listEnsembl() (or one of the following) function to display all available Ensembl BioMart web services.
biomaRt::listMarts()
biomaRt::listEnsembl()

# Annotation from the vertebrate genomes is provided by the main Ensembl project across taxonomic space, with separate BioMart interfaces for Protists, Plants, Metazoa and Fungi.
biomaRt::listEnsemblGenomes()

# Next, connect to the desired BioMart database using the useMart(), useEnsembl(), or useEnsemblGenomes() function.
# The biomart argument should be a valid name from the output of one of the previous functions.
ensembl <- biomaRt::useEnsemblGenomes(biomart = "plants_mart")

# BioMart databases can contain several datasets. For instance, within the Ensembl genes mart, each species is a different dataset.
biomaRt::listDatasets(mart = ensembl)

# Now, run the getids() function with the correct input based on the organism of the study.
# This retrieves IDs using the specified parameters for the query.
bg_t_ids <- getids(
  id_query = f_id_mod, 
  biomart_db = "plants_mart",
  species_dataset = "ptrichocarpa_eg_gene",
  # version = "111",
  transcript_id = "ensembl_gene_id",
  gene_id = "ensembl_gene_id",
  other_ids = "uniprotsptrembl"
)

# Note: The output dataframe must include identifiers supporting the organism:
#               *in the STRING db to create Protein-Protein Interaction Networks
#               *in the g:profiler db to perform functional enrichment
#       For example, with Danio rerio, the "ensembl_gene_id" identifier is supported in STRING and g:profiler.

# Readjust the ensembl_gene_id.1 column to gene_id, as we query it twice, as the transcript_id is also ensembl_gene_id
bg_t_ids <- bg_t_ids |> 
  dplyr::rename(gene_id = ensembl_gene_id.1)

# Save the "time-consuming" data if already created
write.table(bg_t_ids, paste0("outputs/", file_date, "/bg_t_ids_", file_date, ".txt"))

# Load the "time-consuming" data if already created
bg_t_ids <- read.table(paste0("outputs/", file_date, "/bg_t_ids_", file_date, ".txt"))

# The "gene_id" from the background gene list (bg_t_ids) is only needed for functional enrichment. However, the "gene_id" from the deregulated transcripts (DRomics transcriptomics pipeline) is needed for the whole workflow, including creating a STRING PPI network and function enrichment. Therefore, we need to subset the deregulated transcript data from the bg_t_ids dataframe.
dr_t_ids <- bg_t_ids[bg_t_ids$transcript_id %in% b_definedCI_mod_onlyids$id_mod,]

# We need the original ids from the study (e.g. A2M_7932)n which are those used in DRomics. This especially needed for plotting later on 
dr_t_ids <- merge(dr_t_ids, b_definedCI_mod_onlyids, by.x = "transcript_id", by.y = "id_mod")




#>> STEP 4 - Retrieve regulatory status of deregulated genes
#>----------------------------------------------------------

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Step 4 is not performed as Step 1 is not performed
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




#>> STEP 5 - Create and retrieve the clustered PPIN data from the StringApp in Cytoscape
#>--------------------------------------------------------------------------------------

# Create the data to be exported into Cytoscape
# DR_output4string <- merge(dr_t_regs, b_definedCI_mod,
#                           by.x = "transcript_id", by.y = "id_mod")

# As the dr_t_regs is not created, no need to re-merge
DR_output4string <- dr_t_ids

# Save the data 
write.table(DR_output4string, file = paste0("outputs/", file_date, "/DR_output4string_", file_date, ".txt"), row.names = FALSE, sep = "\t")

# Once the clustered network is created, the resulting *.csv* files need to be stored in `outputs/`.

# The Ensembl gene IDs of all the responsive transcripts are manually exported into the StringApp plug-in found within the Cytoscape platform. They are queried with the parameters configured to include all types of interactions, a confidence score of 0.9 and no additional interactions. This creates an interaction network based on the known and predicted protein-protein associations data in the STRING database. Subsequently, to group the interaction network, we use clusterMaker2 to run Markov Chain Clustering (MCL). The inflation parameter is set as the default value (4.0) to reduce cluster size. To extract the relevant cluster information, the 'node table' containing the clustered elements is manually exported from the StringApp as a comma-separated values (CSV) file into the outputs/ folder. This exported table encompasses all the transcripts found within a cluster, their identifications, along with an appended cluster ID column. 


# This table, found in the "outputs/" folder, can then be imported back into the Rstudio environment in order to pursue the workflow : 

dr_t_clustrs <- getclustrs(
  gene_data = dr_t_ids,
  colname_for_merge = "uniprotsptrembl",
  path = paste0("outputs/", file_date, "/cytoscape-files/"),
  nodetable_filename = paste0("string_clusteredppin_cs06_i2_", file_date, ".csv")
)





#>> STEP 6 - Filter clusters based on their gene set size
#>-------------------------------------------------------

dr_t_clustrs_filtr <- clustrfiltr(
  getclustrs_data = dr_t_clustrs,
  size_filtr = 4
)





#>> STEP 7 - Functional enrichment by clusters and annotation retrieval
#>---------------------------------------------------------------------

# You can access different versions of g:Profiler !

# Check the current g:Profiler version
gprofiler2::get_base_url()

# You can set the beta version of g:Profiler  
gprofiler2::set_base_url("http://biit.cs.ut.ee/gprofiler_beta")

# You can use an archived version for reproducibility (e.g., release e94_eg41_p11).
gprofiler2::set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e111_eg58_p18_f463989d")

# You can simply use the stable version of g:Profiler 
gprofiler2::set_base_url("http://biit.cs.ut.ee/gprofiler")
# This corresponds to the latest stable release version 'e111_eg58_p18_f463989d' at the time of this run ["2024-08-07T13:01:09.277618+00:00"]

# Make sure that you are using the appropriate version of the g:Profiler database !!!

# Now, we can run the clustrenrich() to perform cluster-wise functional enrichment. 
clustr_enrichres <- clustrenrich(
  clustrfiltr_data = dr_t_clustrs_filtr,
  dr_genes = dr_t_ids$gene_id,
  bg_genes = bg_t_ids$gene_id,
  bg_type = "custom_annotated",
  sources = c("GO:BP", "KEGG"), # No WP
  organism = "ptrichocarpa",
  user_threshold = 0.05,
  correction_method = "fdr",
  exclude_iea = FALSE, 
  min_term_size = 5,
  max_term_size = 800,
  only_highlighted_GO = TRUE,
  ngenes_enrich_filtr = 2,
  path = paste0("outputs/", file_date, "/"),
  output_filename = paste0("clustr_enrichres_", file_date, "test.rds"),
  overwrite = FALSE
)


# Selecting appropriate minimum and maximum term sizes is challenging and depends on various factors, such as the organism, the type of transcriptomic data, and the definition of "generalistic terms."
# To understand the relationship between gene set size and the generality of terms, the following code prints all terms associated with the deregulated genes in decreasing order of size. This allows us to check if the chosen size limits are reasonable based on the data.

clustr_enrichres$dr_g_a_whole |> 
  dplyr::group_by(term_name) |> 
  dplyr::summarise(count = dplyr::n()) |> 
  dplyr::arrange(desc(count)) |> 
  print(n = 100) # Number of rows to print, adjustable based on the study

# OR #

# View the summarized data in the viewer pane as an interactive datatable. This allows for easier exploration and filtering of terms associated with the deregulated genes to assess the appropriateness of the chosen size limits.

clustr_enrichres$dr_g_a_whole |> 
  dplyr::group_by(term_name) |> 
  dplyr::summarise(count = dplyr::n()) |>
  DT::datatable(
    options = list(pageLength = 10
    ),
    filter = 'top',
    class = c("compact")
  )





#>> STEP 8 - Fusion clusters based on shared cluster enrichment
#>-------------------------------------------------------------

clustr_fusionres <- clustrfusion(
  clustrenrich_data = clustr_enrichres,
  monoterm_fusion = FALSE
)





#>> STEP 9 - Fishing lonely genes sharing annotations with clusters enrichment
#>----------------------------------------------------------------------------

lonely_fishres <- lonelyfishing(
  dr_data = dr_t_ids,
  clustrenrich_data = clustr_enrichres,
  clustrfusion_data = clustr_fusionres,
  friendly_limit = 2,
  path = paste0("outputs/", file_date, "/"),
  output_filename = paste0("lonely_fishres_", file_date, ".rds"), 
  overwrite = FALSE
)



#>> STEP 10 - Characterize the lonely cluster by simple functional enrichment
#>---------------------------------------------------------------------------

# Only select the transcripts part of the lonely cluster

lonelycluster_data <- lonely_fishres$dr_t_c_a_fishing |> 
  dplyr::filter(new_clustr == "Lonely")

# Perform a simple functional enrichment using the simplenrich() homemade function 

lonely_clustr_analysis_res <- simplenrich(
  input_genes = lonelycluster_data$gene_id,
  bg_genes = bg_t_ids$gene_id,
  bg_type = "custom_annotated",
  sources = c("GO:BP", "KEGG"), # No WP
  organism = "ptrichocarpa",
  user_threshold = 0.05,
  correction_method = "fdr",
  only_highlighted_GO = TRUE,
  min_term_size = 5,
  max_term_size = 800,
  ngenes_enrich_filtr = 2,
  path = paste0("outputs/", file_date, "/"),
  output_filename = paste0("lonely_clustr_analysis_res_", file_date, ".rds"),
  overwrite = FALSE
)




#>> Step 11 - Generate a summary dataframe of the workflow -------

# With the workflow fulfilled, we can create a concise summary dataframe capturing the key details from the final results post the lonely gene fishing step. It encompasses all essential information for result exploration, striking a balance by avoiding an overwhelming amount of data that might hinder ease of exploration. 

results_to_csv(
  lonelyfishing_data = lonely_fishres,
  bmdboot_data = b_definedCI_mod,
  path = paste0("outputs/", file_date, "/"),
  output_filename = paste0("summary_workflow_", file_date, ".csv"),
  overwrite = TRUE
)


#>> Generate cluster-level curvesplots to a PDF file -------

# Finally, this last step consists of generating the output PDF file containing a plot of dose-response curves for each cluster of genes, with each plot labelled with the cluster ID and the number of transcripts in that cluster. The curves are color-coded according to whether the trend is increasing, decreasing, U-shaped, or bell-shaped. The plot axes are labelled with "Dose (µg/L)" and "Signal", and the y-axis is scaled to be the same across all plots.

require(ggplot2)

curves_to_pdf(
  lonelyfishing_data = lonely_fishres,
  bmdboot_data = b_definedCI, 
  clustrfusion_data = clustr_fusionres,
  id_col_for_curves = "id",
  tested_doses = unique(f$omicdata$dose), 
  annot_order = c("GO:BP", "KEGG"),
  colorby = "trend",
  addBMD = TRUE,
  scaling = TRUE,
  npoints= 100,
  free.y.scales = FALSE,
  xmin = 0.1, 
  xmax = max(f$omicdata$dose), 
  dose_log_transfo = TRUE, 
  line.size = 0.7, 
  line.alpha = 0.4, 
  point.size = 2, 
  point.alpha = 0.4,
  xunit = "mg/kg",
  xtitle = "Dose (mg/kg)",
  ytitle = "Signal",
  colors = c("inc" = "#1B9E77", "dec" = "#D95F02", "U" = "#7570B3", "bell" = "#E7298A"),
  path = paste0("outputs/", file_date, "/"),
  output_filename = paste0("workflow_curvesplots_", file_date, ".pdf"),
  overwrite = TRUE
)


#>> ADDITIONAL STEPS: Generation Quarto reports
#> --------------------------------------------

# IMPORTANT: Before executing this section, please note that the .qmd files referenced below 
# need to be thoroughly checked and evaluated prior to report generation. If you are using 
# a different dataset or your own custom dataset, these .qmd files will require  
# adaptation and adjustment based on your specific results and data structure.
#
# DATASET COMPATIBILITY WARNING: These .qmd reports are currently adapted specifically 
# for the Zebrafish dataset and were not originally generated for the Rat liver or 
# Poplar root datasets. Users working with alternative datasets should expect to modify 
# the report templates accordingly to match their data characteristics and analysis outputs.

#>> Generate the cluefish quarto report -------

# Note: To ensure that Quarto reports are generated automatically and correctly, follow these guidelines:
# 
# 1. Consistent Filenames:
#    Ensure that output files use the exact filenames specified in your Quarto (.qmd) files.
#    This ensures that Quarto can correctly locate and reference your generated outputs.

# 2. Matching Dates:
#    Output files should:
#    - Share the same date in their filenames, reflecting that they were produced simultaneously.
#    - If necessary, rename files to align with the expected filename format and date consistency.

# 3. Verify and Update the `load-fun-data` Chunk:
#    The `load-fun-data` chunk is the initial code block that loads functions and data required for the report.
#    - Review and update this chunk to ensure it is correctly configured for your specific datasets.
#
# For example: `clustr_enrichres_2024-07-07.rds` and `lonely_fishres_2024-07-07.rds`
# Today's date can be dynamically generated using the "Sys.date()" function.
file_date <- "2025-04-06"

# Now, render and preview the html report. The output is moved to the directory chosen in 'output_path'.
#
# Note: The automation of visualizations for each cluster is currently not implemented. 
#   The challenge lies in dynamically adjusting the visualization code chunks based on the number of clusters 
#   retained after running Cluefish. Each cluster should have associated interactive visualizations, including:
#   - an interactive curveplot
#   - an interactive BMDplot
#   - an interactive table
#   Currently, there is no automated solution for this; manual adjustments to the report are needed to 
#   include or exclude visualization code chunks according to the final number of clusters.

render_qmd(
  input_file = here::here("report_workflow_results.qmd"), 
  output_file = paste0("report_workflow_results_", file_date),
  file_ext = "html",
  output_path = here::here("outputs", file_date, "quarto-outputs"), 
  execute_params = list(`file-date` = file_date)
)



#>> Generate the lonely cluster analysis quarto report -------

# Render and preview the lonely_results_report html report contextualising the lonely cluster
render_qmd(
  input_file = here::here("report_lonely_results.qmd"), 
  output_file = paste0("report_lonely_results_", file_date),
  file_ext = "html",
  output_path = here::here("outputs", file_date, "quarto-outputs"), 
  execute_params = list(`file-date` = file_date)
)



#>> Generate the comparison of cluefish and standard workflow quarto report -------

# Basic enrichment of the deregulated transcripts genes from the DRomics workflow, for comparison with the cluefish workflow
source(here::here("analyses", "standard_approach_dre_pipeline.R"))

# Render and preview the comparison_results_report html report comparing the results of both approaches on the same data
render_qmd(
  input_file = here::here("report_comparison_results.qmd"), 
  output_file = paste0("report_comparison_results_", file_date),
  file_ext = "html",
  output_path = here::here("outputs", file_date, "quarto-outputs"), 
  execute_params = list(`file-date` = file_date)
)


