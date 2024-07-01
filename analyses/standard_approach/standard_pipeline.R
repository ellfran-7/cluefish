## Performing the standard pipeline in order to derive biological interpretations of dose-response modelling results
## ============================================================================

## First and foremost, as the proposed workflow, the standard one begins as such : 
## ----------------------------------------------------------------------------

# Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = "data/raw-data/fitres_zebrafish_phtalate.rds")

# Load DRomics bmdboot object
b <- readRDS(file = "data/raw-data/bootres_zebrafish_phtalate_UF_seed3_5000iter.rds")

# Filter the bmdboot result by selecting only transcripts with a defined confidence interval around the BMD
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

bg_t_ids <- read.table("outputs/bg_t_ids_2024-07-01.txt")

## Create a subset of the items of interests (e.g. deregulated transcripts derived from the DRomics workflow)
dr_t_ids <- bg_t_ids[bg_t_ids$transcript_id %in% BMDres_definedCI$id,]


## In a second part, perform a simple functional enrichment using the simplenrich() homemade function 
## ----------------------------------------------------------------------------

standard_pipeline_res <- simplenrich(
  input_genes = dr_t_ids$gene_id,
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
  path = "analyses/standard_approach/results",
  output_filename = paste0("standard_pipeline_res_", Sys.Date(), ".rds"),
  overwrite = TRUE
)
