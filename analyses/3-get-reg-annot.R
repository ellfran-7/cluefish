# STEP 4 - REGULATORY ANNOTATION RETRIEVAL 
# -----------------------------------------------------------------------------
#
# This script retrieves regulatory status of the deregulated genes using the 'getregs()' function.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '2-getregs.R'.

# We only need the annotations for the deregulated transcripts/genes derived from the DRomics pipeline :
dr_t_ids <- bg_t_ids[bg_t_ids$ensembl_transcript_id_version %in% BMDres_definedCI$id,]

dr_t_regs <- getregs(
  responsiv_ids = dr_t_ids,
  regulator_file = "data/derived-data/Danio_rerio_TF.txt",
  coregulator_file = "data/derived-data/Danio_rerio_Cof.txt")
