# STEP 2 - ID MAPPING 
# -----------------------------------------------------------------------------
#
# This script performs ID mapping using the 'getids()' function.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '1-getids.R'.

bg_t_ids <- getids(
  id_query = f$omicdata$item, 
  species_dataset = "drerio_gene_ensembl",
  id_filter = "ensembl_transcript_id_version", 
  id_attribut = c("ensembl_transcript_id_version",
                  "ensembl_gene_id", 
                  "external_gene_name"
  )
  
)