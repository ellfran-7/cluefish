# Characterizing the lonely cluster

# This script serves to characterize the lonely cluster derived from the Lonely Fishing step within the proposed workflow. This step represents the final phase of data manipulation. While the road towards interpretations for the clustered data is already established, this analysis focuses on understanding the content of the lonely cluster by conducting functional enrichment one last time. This allows us to pool even more transcripts into the interpretation phase.

## Load data -----------

# Load DRomics bmdboot object
b <- readRDS(file = here::here("data", "raw-data", "bootres_zebrafish_phtalate_UF_seed3_5000iter.rds"))

# We filter the bmdboot result by selecting only transcripts with a defined confidence interval around the BMD
BMDres_definedCI <- DRomics::bmdfilter(b$res, BMDfilter = "definedCI")

##Load the workflow results 
lonelyfishing_data <- readRDS(here::here("outputs", "cs09-cf4", "lonely_fishres_cs09_cf4_2024-05-13.rds"))


## Format the data for performing function analysis  -----------

# Combine DRomics info and the workflow info
b_lonely_fishres <- merge(lonelyfishing_data$dr_t_c_a_fishing, BMDres_definedCI,  by = "id")

# Only select the transcripts part of the lonely cluster
b_only_lonely_fishres <- b_lonely_fishres |> 
  dplyr::filter(new_clustr == "Lonely")

# Perfom functional enrichment 
lonely_gostres <- gprofiler2::gost(
  query = b_only_lonely_fishres$ensembl_gene_id, 
  organism = "drerio", 
  ordered_query = FALSE, 
  multi_query = FALSE, 
  significant = TRUE, 
  exclude_iea = FALSE, 
  measure_underrepresentation = FALSE, 
  evcodes = TRUE, 
  user_threshold = 0.05, 
  correction_method = "fdr", 
  domain_scope = "custom_annotated",
  custom_bg = bg_t_ids$ensembl_gene_id, 
  numeric_ns = "", 
  sources = c("GO:BP", "KEGG", "WP"), 
  as_short_link = FALSE, 
  highlight = TRUE 
) 

# Transform the dataframe from "cluster per row" to "gene per row":
dr_c_a_lonely <- lonely_gostres$result |> 
  dplyr::filter(
    ((grepl("GO", source) & highlighted == TRUE) | (!(grepl("GO", source))))
  ) |> 
  dplyr::filter(
    5 <= term_size & term_size <= 500
  ) 
dplyr::select(term_name, term_size, intersection_size, p_value, source) 



# Remove duplicate rows
dr_c_a_lonely <- unique(dr_c_a_lonely)

# Reset the row numbers
rownames(dr_c_a_lonely) <- NULL