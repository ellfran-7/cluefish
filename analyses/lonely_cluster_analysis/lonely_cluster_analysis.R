# Characterizing the lonely cluster

# This script serves to characterize the lonely cluster derived from the Lonely Fishing step within the proposed workflow. This step represents the final phase of data manipulation. While the road towards interpretations for the clustered data is already established, this analysis focuses on understanding/characterizing the content of the lonely cluster by conducting functional enrichment one last time. This allows us to pool even more transcripts into the interpretation phase.

## Load data -----------

# Load DRomics bmdboot object
b <- readRDS(file = "data/raw-data/bootres_zebrafish_phtalate_UF_seed3_5000iter.rds")

# We filter the bmdboot result by selecting only transcripts with a defined confidence interval around the BMD
BMDres_definedCI <- DRomics::bmdfilter(b$res, BMDfilter = "definedCI")

# Retrieve the getids() results: the identifiers of the background transcript list
bg_t_ids <- read.table("outputs/bg_t_ids_2024-07-01.txt")

# Load the workflow results 
lonelyfishing_data <- readRDS("outputs/cs09-cf4/lonely_fishres_cs09_cf4_2024-07-01.rds")


## Format the data for performing function analysis  -----------

# Modify the "transcript_id" column to "id" for easier merge and integration to DRomics visualisations
names(lonelyfishing_data$dr_t_c_a_fishing)[names(lonelyfishing_data$dr_t_c_a_fishing) == "transcript_id"] <- "id"

# Combine DRomics and the workflow results
b_lonely_fishres <- merge(lonelyfishing_data$dr_t_c_a_fishing, BMDres_definedCI,  by = "id")

# Only select the transcripts part of the lonely cluster
b_only_lonely_fishres <- b_lonely_fishres |> 
  dplyr::filter(new_clustr == "Lonely")

# Perfom functional enrichment 
lonely_gostres <- gprofiler2::gost(
  query = b_only_lonely_fishres$gene_id, 
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
  custom_bg = bg_t_ids$gene_id, 
  numeric_ns = "", 
  sources = c("GO:BP", "KEGG", "WP"), 
  as_short_link = FALSE, 
  highlight = TRUE 
) 

## Create dataframe in "gene x annotation per row" (g_a) format for the unfiltered category. 
## The "annotation per row" (a) is already the output of the gprofiler2::gost() function.
dr_g_a_lonely_unfiltered <- lonely_gostres$result |> 
  dplyr::select(intersection, term_name, term_size, query_size, 
                intersection_size, effective_domain_size, p_value, source) |> 
  tidyr::separate_rows(intersection, sep = ",") |> 
  dplyr::rename(gene_id = intersection) |> 
  dplyr::distinct()

## Create dataframes in "gene x annotation per row" (g_a) and "annotation per row" (a) formats for the filtered category

# Filter the results by gene set size and only keep GO terms when they are highlighted
dr_a_lonely_filtered <- lonely_gostres$result |> 
  dplyr::filter(
    ((grepl("GO", source) & highlighted == TRUE) | (!(grepl("GO", source))))) |> 
  dplyr::filter(5 <= term_size & term_size <= 500)

# Keep only rows where the number of IDs in the intersection column is 3 or more
dr_a_lonely_filtered <- dr_a_lonely_filtered |> 
  dplyr::mutate(num_genes = sapply(strsplit(intersection, ","), length)) |> 
  dplyr::filter(num_genes >= 3) |> 
  dplyr::select(-num_genes)  # Remove the temporary column

# Format the dataframe from "annotation per row" (a) to "gene x annotation per row" (g_a)
dr_g_a_lonely_filtered <- dr_a_lonely_filtered |> 
  dplyr::select(intersection, term_name, term_size, query_size, 
                intersection_size, effective_domain_size, p_value, source) |> 
  tidyr::separate_rows(intersection, sep = ",") |> 
  dplyr::rename(gene_id = intersection) |> 
  dplyr::distinct()

# Reset the row numbers
rownames(dr_g_a_lonely_unfiltered) <- NULL
rownames(dr_g_a_lonely_filtered) <- NULL
rownames(dr_a_lonely_filtered) <- NULL

# Create a list containing the enrichment results, annotations, and the trace of biological function filtering
lonely_cluster_analysis_res <- list(
  unfiltered = list(dr_g_a_lonely_unfiltered = dr_g_a_lonely_unfiltered,
                    lonely_gostres = lonely_gostres),
  filtered = list(dr_g_a_lonely_filtered = dr_g_a_lonely_filtered,
                  dr_a_lonely_filtered = dr_a_lonely_filtered)
  )

# Save the output
saveRDS(lonely_cluster_analysis_res, file.path(paste0("analyses/lonely_cluster_analysis/lonely_cluster_analysis_res_", Sys.Date(), ".rds")))

## The output is a named list where 'unfiltered' is the gost 'result' dataframe before filtering and 'filtered' is a sub-list of the gost 'results dataframe after filtering with 'dr_g_a_lonely_filtered' being the reuslts with each row representing a gene and 'dr_a_lonely_filtered' being the results with each row representing a cluster.