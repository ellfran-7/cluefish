# STEP 8 - PER-CLUSTER FUNCTIONAL ENRICHMENT 
# -----------------------------------------------------------------------------
#
# This script performs over-representation analysis independently for each cluster, using the 'clustrenrich()' function. It harbors additional features for keeping only driver GO terms and set geneset size lower/upper limits for enriched terms. Intermediate steps derived annotations for all the deregulated genes and a number of biological function summary for each cluster for each biological after each filter.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '5-clustrenrich.R'.

clustr_enrichres <- clustrenrich(
  clustr.data = dr_t_clustrs_filtr$kept,
  responsiv.genes = dr_t_regs$ensembl_gene_id,
  background = bg_t_ids$ensembl_gene_id,
  org = "drerio",
  pcutoff = 0.05,
  padjmethod = "fdr",
  databases = c("GO:BP", "KEGG", "WP"),
  background.type = "custom_annotated",
  exclude.iea = FALSE,
  enrich.size.filtr = TRUE,
  min.term.size = 5,
  max.term.size = 500,
  only.highlighted.GO = TRUE,
  path = "outputs/",
  output.filename = "clustr_enrichres",
  overwrite = TRUE
  )



