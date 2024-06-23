## Performing the standard pipeline in order to derive biological interpretations of dose-response modelling results

## First and foremost, as the proposed workflow, the standard one necessitates : 
# ------------------------------------------------------------------------------

# ** Load DRomics drcfit object (which holds the background transcript list) **
f <- readRDS(file = "data/raw-data/fitres_zebrafish_phtalate.rds")

# ** Load DRomics bmdboot object **
b <- readRDS(file = "data/raw-data/bootres_zebrafish_phtalate_UF_seed3_5000iter.rds")

# We filter the bmdboot result by selecting only transcripts with a defined confidence interval around the BMD
BMDres_definedCI <- DRomics::bmdfilter(b$res, BMDfilter = "definedCI")



# ** mapping the experiment transcript identifiers to corresponding gene identifiers allowing us to perform functional enrichment **

ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl") 
# enables connection to ensembl database and species dataset within

bg_t_ids <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id_version"), 
                                filters = "ensembl_transcript_id_version",
                                values = f$omicdata$item,
                                mart = ensembl)

# Only select the selected gene transcripts from the DRomics workflow 
dr_t_ids <- bg_t_ids[bg_t_ids$ensembl_transcript_id_version %in% BMDres_definedCI$id,]


## Secondly, we can perform functional enrichment analysis using the gprofiler2::gost function : 
## --------------------------------------------------------------------------------------------
standard_gostres <- gprofiler2::gost(
  query = dr_t_ids$ensembl_gene_id, 
  organism = "drerio", 
  ordered_query = FALSE, 
  multi_query = FALSE, 
  significant = TRUE, 
  exclude_iea = FALSE, 
  measure_underrepresentation = FALSE, 
  evcodes = TRUE, 
  user_threshold = 0.05, 
  correction_method = "fdr", 
  domain_scope = "custom_annotated",
  custom_bg = bg_t_ids$ensembl_gene_id, 
  numeric_ns = "", 
  sources = c("GO:BP", "KEGG", "WP"), 
  as_short_link = FALSE, 
  highlight = TRUE 
) 

# Filter the results by gene set size and only keep GO terms when they are highlighted
dr_c_a_standard <- standard_gostres$result |> 
  dplyr::filter(
    ((grepl("GO", source) & highlighted == TRUE) | (!(grepl("GO", source))))
  ) 

dr_c_a_standard_filter <- dr_c_a_standard |> 
  dplyr::filter(
    5 <= term_size & term_size <= 500
  ) 

# Reset the row numbers
rownames(dr_c_a_standard) <- NULL

# Create a list containing the enrichment results, annotations, and the trace of biological function filtering
dr_c_a_standard_gostres <- list(
  unfiltered = dr_c_a_standard,
  filtered = dr_c_a_standard_filter
)

# Save the output
saveRDS(dr_c_a_standard_gostres, file.path("analyses/standard-approach/standard_gostres.rds"))

## The output is a named list where 'unfiltered' is the gost 'result' dataframe before filtering and 'filtered' is the gost 'results dataframe after filtering. 