#' Workflow Validation Script: Exploring Functionality Step-by-Step 
#' 
#' @description 
#' To become acquainted with the workflow's functionality, this script serves as a test to explore and validate each step and function within it. This script is specifically dedicated to the later steps 5 to 9, which are longer, more complex, and therefore prone to potential uncertainties.
#' 
#' @author Ellis Franklin \email{ellis.franklin@univ-lorraine.fr}
#' 
#' @date 2024/04/19


# Load the sample data needed to test the workflow from steps 5 to 8. These can be found in the *data/derived-data/tests/* folder.

# A sample data from Step 4, holding 100 unique Ensembl genes and 104 unique Ensembl transcripts (ENSDARG00000041619 and ENSDARG00000077291 are each associated to two transcripts). 

sample_5_clustrs <- read.table(file = "data/derived-data/tests/sample_input_5_clusters.txt")
head(sample_5_clustrs)

# A sample background dataframe containing Ensembl gene and transcript identifiers, alongside an additional external gene name representing the all the gene transcripts from a theoretical experience.

sample_bg <- read.table(file = "data/derived-data/tests/sample_background.txt")
head(sample_bg)

# A sample deregulated dataframe containing Ensembl gene and transcript identifiers, alongside an additional external gene name representing the selected gene transcripts from the DRomics::drcfit() step.

sample_dr <- read.table(file = "data/derived-data/tests/sample_deregulated.txt")
head(sample_dr)




#>> STEP 5 TEST - Filter clusters based on their gene set size
#>-------------------------------------------------------
sample_dr_t_clustrs_filtr <- clustrfiltr(
  getclustrs_data = sample_5_clustrs,
  size_filtr = 3
)
# Remark:
# All the clusters have passed the filter, as they harbor each at-least 3 Ensembl genes.




#>> STEP 6 TEST - Functional enrichment by clusters and annotation retrieval
#>---------------------------------------------------------------------
sample_clustr_enrichres <- clustrenrich(
  clustrfiltr_data = sample_dr_t_clustrs_filtr,
  dr_genes = sample_dr_t_regs$ensembl_gene_id,
  bg_genes = bg_t_ids$ensembl_gene_id,
  organism = "drerio",
  user_threshold = 0.05,
  correction_method = "fdr",
  enrich_size_filtr = TRUE,
  min_term_size = 5,
  max_term_size = 500,
  only_highlighted_GO = TRUE,
  ngenes_enrich_filtr = 3,
  path = "outputs/tests/",
  output_filename = "sample_clustr_enrichres_2024_04_19.rds",
  overwrite = TRUE
)
# Remark: 
# All four clusters contribute to the enrichment process, resulting in 8 out of 13 enriched pathways retained after term filtering. Two filters are applied: one based on the number of genes involved in the enrichment, and the other based on the size of the term gene set. In this scenario, the analysis retains 8 pathways following the first filter, indicating that they are enriched by at least 3 genes.




#>> STEP 7 TEST - Fusion clusters based on shared cluster enrichment
#>-------------------------------------------------------------
sample_clustr_fusionres <- clustrfusion(
  clustrenrich_data = sample_clustr_enrichres
)
# Remark:
# After the fusion process, only 2 out of the initial 4 clusters remain. Reviewing the fusion log reveals that cluster 6 merges with cluster 2 due to a shared unique GO term. The resulting cluster 2, now merged with cluster 6, is then further merged with cluster 1 based on a shared KEGG term. Cluster 5 remains unaffected.




#>> STEP 8 TEST - Fishing lonely genes sharing annotations with clusters enrichment
#>----------------------------------------------------------------------------
sample_lonely_fishres <- lonelyfishing(
  dr_data = sample_dr_t_regs,
  clustrenrich_data = sample_clustr_enrichres,
  clustrfusion_data = sample_clustr_fusionres,
  friendly_limit = 0,
  path = "outputs/tests/",
  output_filename = "sample_lonely_fishres_2024_04_19.rds",
  overwrite = TRUE
)
# Remark:
# The lonely fishing process identified 75 annotated lonely genes, with 13 of them fished into clusters. Following fusion, cluster 1 contained 83 unique Ensembl genes, which increased to 95 after the lonely fishing. Cluster 5 gained 1 gene, expanding from 17 to 18 Ensembl genes.

