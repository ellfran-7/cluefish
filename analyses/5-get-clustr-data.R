# STEP 6 -  RETRIEVE THE CLUSTERED NETWORK NODE TABLE 
# -----------------------------------------------------------------------------
#
# This script retrieves, and formats the node table created in the Cytoscape app, using the 'getclustrs()' function.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '3-getclustrs.R'.

dr_t_clustrs <- getclustrs(
  responsiv_data = dr_t_regs,
  string_clustr_file = "outputs/Resp_PPIN_clustered_051023.csv"
)

