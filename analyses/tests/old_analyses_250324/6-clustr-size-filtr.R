# STEP 7 - CLUSTER SIZE FILTERING
# -----------------------------------------------------------------------------
#
# This script filters out clusters with a chosen size limit using the 'clustrfiltr' function.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '4-clustrfiltr.R'.

dr_t_clustrs_filtr <- clustrfiltr(
  clustr.data = dr_t_clustrs,
  size.filtr = 3
  )

