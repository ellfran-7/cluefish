#> =============================================================================
#> Performing the standard pipeline in order to derive biological interpretations of dose-response modelling results
#> =============================================================================

#> This script applies the cluefish approach pipeline for biological interpretation 
#> of results of the DRomics analysis of the rat liver pfoa transcriptomic dataset

#> Before continuing you will need to follow the Cluefish workflow in the 
#> "make.R" file at the root up to Step 3, in order to retrieve the bg_t_ids
#> dataset holding the backgrounds identifiers retrieved from Ensembl
#> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


## State the Time Variable for file saving and reading
file_date <- "2024-12-04"

# Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = "data/derived-data/fitres_zebrafish_phtalate.rds")

# Load DRomics "bmdboot" results, filtered to include only transcripts with a defined confidence interval around the BMD.
b_definedCI <- readRDS(file = "data/derived-data/bootres_zebrafish_phtalate_UF_seed3_5000iter_definedCI.rds")

## Retrieve the getids() output holding the identifiers for the deregulated transcripts
bg_t_ids <- read.table(paste0("outputs/", file_date, "/bg_t_ids_", file_date, ".txt"))

## Create a subset of the items of interests (e.g. deregulated transcripts derived from the DRomics workflow)
dr_t_ids <- bg_t_ids[bg_t_ids$transcript_id %in% b_definedCI$id,]


#> In a second part, perform a simple functional enrichment using the simplenrich() homemade function
#> ----------------------------------------------------------------------------

standard_pipeline_res <- simplenrich(
  input_genes = dr_t_ids$gene_id,
  bg_genes = bg_t_ids$gene_id,
  bg_type = "custom_annotated",
  sources = c("GO:BP", "KEGG", "WP"), 
  organism = "drerio",
  user_threshold = 0.05,
  correction_method = "fdr",
  only_highlighted_GO = TRUE,
  min_term_size = 5,
  max_term_size = 500,
  ngenes_enrich_filtr = 3,
  path = paste0("outputs/", file_date, "/"),
  output_filename = paste0("standard_pipeline_res_", file_date, ".rds"),
  overwrite = FALSE
)
