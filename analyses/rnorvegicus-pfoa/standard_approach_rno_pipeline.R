#> =============================================================================
#> Performing the standard pipeline in order to derive biological interpretations of dose-response modelling results
#> =============================================================================

#> This script applies the standard approach pipeline for biological interpretation 
#> of results of the DRomics analysis of the main zebrafish dbp transcriptomic dataset

#> Before continuing you will need to follow the Cluefish workflow in the 
#> "cluefish_approach_dre_pipeline.R" file at the root up to Step 3, in order 
#> to retrieve the bg_t_ids dataset holding the backgrounds identifiers 
#> retrieved from Ensembl
#> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# In a first, the standard approach begins similarly to the Cluefish approach

## State the Time Variable for file saving and reading
file_date <- "2025-04-01"

# Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = "data/raw-data/fitres_rat_liver_pfoa.rds")

# Load DRomics bmdboot object
b_definedCI <- readRDS(file = "data/raw-data/bootres_rat_liver_pfoa_seed1234_5000iter_definedCI.rds")

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
  organism = "rnorvegicus",
  user_threshold = 0.05,
  correction_method = "fdr",
  only_highlighted_GO = TRUE,
  min_term_size = 5,
  max_term_size = 500,
  ngenes_enrich_filtr = 3,
  path = paste0("outputs/", file_date, "/"),
  output_filename = paste0("standard_pipeline_rat_liver_pfoa_res_", file_date, ".rds"),
  overwrite = FALSE
)
