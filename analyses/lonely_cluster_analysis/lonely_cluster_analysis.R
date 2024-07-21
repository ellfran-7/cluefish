## Characterizing the lonely cluster

# This script serves to characterize the lonely cluster derived from the Lonely Fishing step within the proposed workflow. This step represents the final phase of data manipulation. While the road towards interpretations for the clustered data is already established, this analysis focuses on understanding/characterizing the content of the lonely cluster by conducting functional enrichment one last time. This allows us to pool even more transcripts into the interpretation phase.
## ============================================================================

## State the Time Variable for file saving and reading
file_date <- "2024-07-07"

## Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = "data/raw-data/fitres_zebrafish_phtalate.rds")

## Load DRomics bmdboot object
b <- readRDS(file = "data/raw-data/bootres_zebrafish_phtalate_UF_seed3_5000iter.rds")

# We filter the bmdboot result by selecting only transcripts with a defined confidence interval around the BMD
BMDres_definedCI <- DRomics::bmdfilter(b$res, BMDfilter = "definedCI")

## Retrieve deregulated gene identifiers from Ensembl
bg_t_ids <- getids(
  id_query = f$omicdata$item, 
  biomart_db = "ENSEMBL_MART_ENSEMBL",
  species_dataset = "drerio_gene_ensembl",
  transcript_id = "ensembl_transcript_id_version",
  gene_id = "ensembl_gene_id",
  gene_name = "external_gene_name"
)

# OR

bg_t_ids <- read.table(paste0("outputs/bg_t_ids_", file_date, ".txt"))

## Load the workflow results 
lonelyfishing_data <- readRDS(paste0("outputs/cs09-cf4/lonely_fishres_cs09_cf4_", file_date, ".rds"))

# Only select the transcripts part of the lonely cluster
lonelycluster_data <- lonelyfishing_data$dr_t_c_a_fishing |> 
  dplyr::filter(new_clustr == "Lonely")


## In a second part, perform a simple functional enrichment using the simplenrich() homemade function 
## ----------------------------------------------------------------------------

lonely_cluster_analysis_res <- simplenrich(
  input_genes = lonelycluster_data$gene_id,
  bg_genes = bg_t_ids$gene_id,
  bg_type = "custom_annotated",
  sources = c("GO:BP", "KEGG", "WP"), 
  organism = "drerio",
  user_threshold = 0.05,
  correction_method = "fdr",
  min_term_size = 5,
  max_term_size = 500,
  only_highlighted_GO = TRUE,
  ngenes_enrich_filtr = 3,
  path = "analyses/lonely_cluster_analysis/results",
  output_filename = paste0("lonely_cluster_analysis_res_", file_date, ".rds"),
  overwrite = TRUE
)
