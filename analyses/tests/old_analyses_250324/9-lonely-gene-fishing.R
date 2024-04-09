# STEP 10 -  LONELY GENE FISHING
# -----------------------------------------------------------------------------
#
# This script fishes genes not part of a cluster up to now, know as lonely genes, into clusters sharing a same biological function enrichment(s) to the lonely gene's annotation(s). This is done using the 'lonelyfishing()' function. Any gene not fished is kept in a "Lonely" cluster. The function harbors an "friendly" metric that informs on the belongingness of genes between clusters, and the user can set a overfriendly limit, associating any gene over this limit to a "Friendly" cluster.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '7-lonelyfishing.R'.

lonely_fishres <- lonelyfishing(
  responsiv.data = dr_t_regs,
  clustr.enrich.results = clustr_enrichres,
  clustr.fusion.data = clustr_fusionres$dr_g_a_fusion,
  friendly.limit = 3,
  path = "outputs/",
  output.filename = "lonely_fishres.rds",
  overwrite = TRUE)
