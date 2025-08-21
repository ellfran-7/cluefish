#> =============================================================================
#> Performing the standard pipeline in order to derive biological interpretations of dose-response modelling results
#> =============================================================================

#> This script applies the standard approach pipeline for biological interpretation 
#> of results of the DRomics analysis of the poplar root transcriptomic dataset

#> PREREQUISITES: Before continuing you will need to follow the Cluefish workflow in the 
#> "cluefish_approach_pca_pipeline.R" file at the root up to Step 3, in order 
#> to retrieve the bg_t_ids dataset holding the backgrounds identifiers 
#> retrieved from Ensembl
#> 
#> NOTE: This pipeline is specifically configured for the poplar root dataset.
#> 
#> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# In a first, the standard approach begins similarly to the Cluefish approach

## State the Time Variable for file saving and reading
file_date <- "2025-04-06"

# Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = "data/derived-data/fitres_pop_root_phe.rds")

# Load DRomics bmdboot object
b_definedCI <- readRDS(file = "data/derived-data/bootres_pop_root_phe_seed1234_5000iter_definedCI.rds")

# Create a vector of all the genes from the experiment but adding the suffixe ".v4.1"
f_id_mod <- paste0(f$omicdata$item, ".v4.1")

# Create an additional column in the DRomics output that holds the deregulated genes with the added suffixe ".v4.1"
b_definedCI_mod <- b_definedCI |> 
  dplyr::mutate(id_mod = paste0(id, ".v4.1"))

## Retrieve the getids() output holding the identifiers for the deregulated transcripts
bg_t_ids <- read.table(paste0("outputs/", file_date, "/bg_t_ids_", file_date, ".txt"))

# The "gene_id" from the background gene list (bg_t_ids) is only needed for function enrichment. However, the "gene_id" from the deregulated transcripts (DRomics transcriptomics pipeline) is needed for the whole workflow, including creating a STRING PPI network and function enrichment. Therefore, we need to subset the deregulated transcript data from the bg_t_ids dataframe.
dr_t_ids <- bg_t_ids[bg_t_ids$transcript_id %in% b_definedCI_mod$id_mod,]

#> In a second part, perform a simple functional enrichment using the simplenrich() homemade function
#> ----------------------------------------------------------------------------

standard_pipeline_res <- simplenrich(
  input_genes = dr_t_ids$gene_id,
  bg_genes = bg_t_ids$gene_id,
  bg_type = "custom_annotated",
  sources = c("GO:BP", "KEGG"),  # No WP
  organism = "ptrichocarpa",
  user_threshold = 0.05,
  correction_method = "fdr",
  only_highlighted_GO = TRUE,
  min_term_size = 5,
  max_term_size = 800,
  ngenes_enrich_filtr = 2,
  path = paste0("outputs/", file_date, "/"),
  output_filename = paste0("standard_pipeline_res_", file_date, ".rds"),
  overwrite = FALSE
)
