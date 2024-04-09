# STEP 8 - PER-CLUSTER FUNCTIONAL ENRICHMENT 
# -----------------------------------------------------------------------------
#
# This script performs over-representation analysis independently for each cluster, using the 'clustrenrich()' function. It harbors additional features for keeping only driver GO terms and set geneset size lower/upper limits for enriched terms. Intermediate steps derived annotations for all the deregulated genes and a number of biological function summary for each cluster for each biological after each filter.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '5-clustrenrich.R'.

clustr_enrichres <- clustrenrich(
  clustrfiltr_data = dr_t_clustrs_filtr,
  dr_genes = dr_t_regs$ensembl_gene_id,
  bg_genes = bg_t_ids$ensembl_gene_id,
  organism = "drerio",
  user_threshold = 0.05,
  correction_method = "fdr",
  sources = c("GO:BP", "KEGG", "WP"),
  enrich_size_filtr = TRUE,
  n_genes_filter = 3,
  min_term_size = 5,
  max_term_size = 500,
  only_highlighted_GO = TRUE,
  path = "outputs/",
  output_filename = "clustr_enrichres.rds",
  overwrite = TRUE
)



