# STEP 11 -  GENERATING A GENE-LEVEL SUMMARY TABLE
# -----------------------------------------------------------------------------
#
# This script creates a concise gene-level summary table capturing the key details from the final results post the lonely gene fishing step. It encompasses all essential information for result exploration, striking a balance by avoiding an overwhelming amount of data that might hinder ease of exploration.
#
# The output should be saved to the "outputs" file
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name 'results_to_csv.R'.
#

summary_to_csv(
  lonelyfishing_data = lonely_fishres,
  bmdboot_data = BMDres_definedCI,
  path = "outputs/",
  output_filename = "summary_workflow.csv",
  overwrite = TRUE
)