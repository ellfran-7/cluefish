# STEP 4 - REGULATORY ANNOTATION RETRIEVAL 
# -----------------------------------------------------------------------------
#
# This script retrieves regulatory status of the deregulated genes using the 'getregs()' function.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '2-getregs.R'.

# We only need the annotations for the deregulated transcripts/genes from the DRomics workflow :
dr_t_ids <- bg_t_ids[bg_t_ids$ensembl_transcript_id_version %in% BMDres_definedCI$id,]

dr_t_regs <- getregs(
  responsiv.ids = dr_t_ids,
  regulator.file = "data/derived-data/Danio_rerio_TF.txt",
  coregulator.file = "data/derived-data/Danio_rerio_Cof.txt")
