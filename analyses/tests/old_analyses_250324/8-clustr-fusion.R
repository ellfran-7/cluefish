# STEP 9 -  CLUSTER FUSION
# -----------------------------------------------------------------------------
#
# This script merges clusters sharing singular biological function enrichment in at-least one annotation database (GO, KEGG...), using the 'clustrfusion()' function. 
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '6-clustrfusion.R'.

clustr_fusionres <- clustrfusion(
  clustr.enrich.results = clustr_enrichres
  )

