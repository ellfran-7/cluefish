# STEP 2 - ID MAPPING 
# -----------------------------------------------------------------------------
#
# This script performs ID mapping using the 'getids()' function.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '1-getids.R'.

bg_t_ids <- getids(
  id.query = f$omicdata$item, 
  species.dataset = "drerio_gene_ensembl",
  id.filter = "ensembl_transcript_id_version", 
  id.attribut = c("ensembl_gene_id", "ensembl_transcript_id_version", 
                  "external_gene_name", "entrezgene_id"))
