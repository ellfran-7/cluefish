#> ==========================================================================
#> Script to generate the supplementary tables in Franklin et al. (submitted)
#> ==========================================================================
#> Case: drerio - RNA-seq ****
#> Run : 12/04/24 ****



#> Load and set everything necessary for the script (e.g., functions, variables)
#> -----------------------------------------------------------------------------

#> Load the project functions --------------------

devtools::load_all(here::here())

#> Load packages --------------------

require(ggplot2)
require(patchwork)
require(here)

#> Load fonts --------------------

require(extrafont)
# First, import the fonts installed on the system, if not already done
# font_import()

# Second, register the fonts with the PDF output device (only necessary in session where you ran font_import())
loadfonts()



#> Set the file_date variable --------------------
file_date = "2024-12-04"



#> Create directory for saving supp. data (if it doesn't already exist) --------
dir_path <- "figures/for-supp"

if (!dir.exists(dir_path)) { # Check if the directory path exists
  
  dir.create(dir_path) # Create it if not
  
}



#> Load the necessary data --------------------

## DRomics pipeline - Transcriptomics data ----

# Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = here::here("data", "derived-data", "fitres_zebrafish_phtalate.rds"))

# Load DRomics "bmdboot" results filtered with only transcripts with a defined confidence interval around the BMD
b_definedCI <- readRDS(file = here::here("data", "derived-data", "bootres_zebrafish_phtalate_UF_seed3_5000iter_definedCI.rds"))



## Cluefish workflow ----

# Load the clustrenrich() result
clustr_enrichres <- readRDS(here::here("outputs", file_date, paste0("clustr_enrichres_", file_date, ".rds")))

# Retrieve the clustrfusion() result
clustr_fusionres <- clustrfusion(
  clustrenrich_data = clustr_enrichres
)

# Load the lonelyfishing() result
lonelyfishing_data <- readRDS(here::here("outputs", file_date, paste0("lonely_fishres_", file_date, ".rds")))

# Modify the "ensembl_transcript_id_version" column to "id" for easier merge and integration to DRomics visualisations
names(lonelyfishing_data$dr_t_c_a_fishing)[names(lonelyfishing_data$dr_t_c_a_fishing) == "transcript_id"] <- "id"

# Combine DRomics info and the workflow info
b_lonely_fishres <- merge(lonelyfishing_data$dr_t_c_a_fishing, b_definedCI,  by = "id")

# Load the lonely cluster analysis simplenrich() result
lonely_cluster_analysis_res <- readRDS(here::here("outputs", file_date, paste0("lonely_clustr_analysis_res_", file_date, ".rds")))

## Standard workflow ----

# Load the standard workflow results
stand_res <- readRDS(here::here("outputs", file_date, paste0("standard_pipeline_res_", file_date, ".rds")))

# ------------------------------------------------------------------------------




#> Resulting table of all clusters BEFORE MERGING AND FISHING with:
#>  - cluster number
#>  - function name
#>  - function ID
#>  - source
#>  - p-value
#>  - highlighted status
#>  - term_size (number of genes associated with the term)
#>  - gene count in query
#>  - precision and recall of term in the gene list
#>  - enrichment ratio (overlap/expect)
#>  - domain scope (effective domain size)
#> ----------------------------------------

# Display the structure of the initial cluster enrichment results
dplyr::glimpse(clustr_enrichres$gostres$result)

# Select relevant columns for enrichment data, ensuring each row represents a unique term and its enrichment data
enrich_data <- clustr_enrichres$gostres$result |> 
  dplyr::select(query, term_name, term_id, source, p_value, term_size, query_size, intersection_size, precision, recall, enrichment_ratio, effective_domain_size) |> 
  dplyr::mutate(p_value = signif(p_value, 2),
                precision = signif(precision, 2),
                recall = signif(recall, 2),
                enrichment_ratio = signif(enrichment_ratio, 2)) |> 
  dplyr::distinct()

# Uncomment below lines if you want to inspect the dataset interactively
# dplyr::glimpse(enrich_data)  # View structure of selected enrichment data
# DT::datatable(enrich_data)   # Show data in an interactive table format
# head(enrich_data)            # Show first few rows of enrichment data
# dim(enrich_data)             # Check dimensions of enrichment data

# # Write the selected enrichment data to a CSV file for supplementary material
# write.csv2(enrich_data, 
#            paste0("figures/for-supp/supp-table-original-gprofiler-res", Sys.Date(), ".csv"),
#            row.names = FALSE)

#> ----------------------------------------



#> Resulting table of all clusters AFTER MERGING AND FISHING (end of workflow) with:
#>  - original cluster number
#>  - new cluster number (post-fusion with related clusters)
#>  - function names 
#>  - function size
#>  - function id
#>  - data sources for functions
#> ----------------------------------------

# Check structure of post-merging data for fishing clusters
dplyr::glimpse(lonelyfishing_data$dr_c_a_fishing)

# Select and filter merged clusters that are highlighted as significant, keeping relevant columns for supplementary results
lonelyfish_data <- lonelyfishing_data$dr_c_a_fishing |> 
  dplyr::select(old_clustr, new_clustr, term_name, 
                term_size, term_id, source) |> 
  dplyr::filter(term_name %in% clustr_fusionres$dr_g_a_fusion$term_name) |> 
  dplyr::distinct()

# Export merged and filtered cluster data to CSV for supplementary materials
write.csv2(lonelyfish_data,
           paste0("figures/for-supp/supp-table-lonelyfish-res", Sys.Date(), ".csv"),
           row.names = FALSE)

#> ----------------------------------------




#> Resulting table of lonely cluster enrichment analysis (after filtering) with:
#>  - function name and ID
#>  - data source
#>  - p-value for enrichment
#>  - highlighted status (indicates driver terms)
#>  - term_size and gene count within each term
#>  - precision and recall scores for enrichment
#>  - enrichment ratio (overlap/expect)
#>  - domain scope
#> ----------------------------------------

# Display the structure of the filtered enrichment analysis results for lonely clusters
dplyr::glimpse(lonely_cluster_analysis_res$filtered)

# Select columns from the lonely cluster analysis, post-filtering
lonely_enrich_data <- lonely_cluster_analysis_res$filtered$dr_a |> 
  dplyr::select(term_name, term_id, source, p_value, term_size, query_size, intersection_size, precision, recall, enrichment_ratio, effective_domain_size) |> 
  dplyr::mutate(p_value = signif(p_value, 2),
                precision = signif(precision, 2),
                recall = signif(recall, 2),
                enrichment_ratio = signif(enrichment_ratio, 2)) |> 
  dplyr::distinct()

# # Save filtered lonely cluster enrichment data to CSV
# write.csv2(lonely_enrich_data, 
#            paste0("figures/for-supp/supp-table-lonelycluster-gprofiler-res", Sys.Date(), ".csv"),
#            row.names = FALSE)

#> ----------------------------------------




#> Resulting table of the standard enrichment analysis (after filtering) with:
#>  - function name and ID
#>  - data source
#>  - adjusted p-value for significance
#>  - highlighted status (indicates terms with strong signals)
#>  - term_size and query gene counts
#>  - precision and recall metrics
#>  - enrichment ratio (overlap/expect)
#>  - domain scope for the term
#> ----------------------------------------

# Glimpse of the filtered standard enrichment analysis results
dplyr::glimpse(stand_res$filtered$dr_a)

# Selecting columns of interest for the standard enrichment analysis results
stand_enrich_data <- stand_res$filtered$dr_a |> 
  dplyr::select(term_name, term_id, source, p_value, term_size, query_size, intersection_size, precision, recall, enrichment_ratio, effective_domain_size) |> 
  dplyr::mutate(p_value = signif(p_value, 2),
                precision = signif(precision, 2),
                recall = signif(recall, 2),
                enrichment_ratio = signif(enrichment_ratio, 2)) |> 
  dplyr::distinct()

# # Save filtered standard enrichment data to CSV for supplemental output
# write.csv2(stand_enrich_data, 
#            paste0("figures/for-supp/supp-table-standard-gprofiler-res", Sys.Date(), ".csv"),
#            row.names = FALSE)

#> ----------------------------------------





#> Resulting table of the summarized BMD values for each cluster
#>  - cluster number
#>  - BMD first quantile
#>  - BMD median
#>  - BMD third quantile
#> ----------------------------------------

# View structure of the BMD (Benchmark Dose) values before summarization
dplyr::glimpse(b_lonely_fishres)

# Summarize BMD quantiles by cluster, grouping by new cluster number
summary_bmd_clusters <- b_lonely_fishres |> 
  dplyr::group_by(new_clustr) |> 
  dplyr::summarize(
    first_quantile = signif(quantile(BMD.zSD, 0.25, na.rm = TRUE), 2),  # 25th percentile
    median = signif(quantile(BMD.zSD, 0.50, na.rm = TRUE), 2),          # 50th percentile (median)
    third_quantile = signif(quantile(BMD.zSD, 0.75, na.rm = TRUE), 2)   # 75th percentile
  ) |> 
  # Sort clusters numerically, including mixed text/numeric values, for consistent ordering in output
  dplyr::mutate(new_clustr = gtools::mixedsort(new_clustr)) |> 
  dplyr::ungroup()  # Remove grouping structure for clarity in output

# Check the structure of the summarized BMD clusters
dplyr::glimpse(summary_bmd_clusters)

# # Write the summarized BMD values by cluster to a CSV file
# write.csv2(summary_bmd_clusters, 
#            paste0("figures/for-supp/supp-table-cluster-summarized-bmd", Sys.Date(), ".csv"),
#            row.names = FALSE)

#> ----------------------------------------






#> Hox genes and their transcript modelling metrics
#>  - gene name
#>  - BMD.zSD
#>  - trend
#>  - typology
#> ----------------------------------------

# Catch a glimpse of the data
dplyr::glimpse(b_lonely_fishres)

# Create the dataframe
hox_gene_data <- b_lonely_fishres |> 
  dplyr::filter(grepl("^hox", gene_name)) |> 
  dplyr::select(gene_name, BMD.zSD, trend, TF) |> 
  dplyr::mutate(gene_name = gsub("_t[0-9]+", "", gene_name)) |>
  dplyr::mutate(BMD.zSD = signif(BMD.zSD, 2)) |> 
  dplyr::distinct() |> 
  dplyr::arrange(gene_name)

hox_gene_data

# # Write the summarized BMD values by cluster to a CSV file
# write.csv2(hox_gene_data, 
#            paste0("figures/for-supp/supp-table-hox-gene-data", Sys.Date(), ".csv"),
#            row.names = FALSE)

#> ----------------------------------------






#> Crystiallin genes and their transcript modelling metrics
#>  - gene name
#>  - BMD.zSD
#>  - trend
#>  - typology
#> ----------------------------------------

# Catch a glimpse of the data
dplyr::glimpse(b_lonely_fishres)

# Create the dataframe
cry_gene_data <- b_lonely_fishres |> 
  dplyr::filter(grepl("^cry", gene_name)) |> 
  dplyr::select(gene_name, gene_id, id, BMD.zSD, trend, TF) |> 
  dplyr::mutate(BMD.zSD = signif(BMD.zSD, 2)) |> 
  dplyr::distinct() |> 
  dplyr::arrange(gene_name)

cry_gene_data

# # Write the summarized BMD values by cluster to a CSV file
# write.csv2(cry_gene_data, 
#            paste0("figures/for-supp/supp-table-cry-gene-data", Sys.Date(), ".csv"),
#            row.names = FALSE)

#> ----------------------------------------






#> We can export all dataframes to separate sheets within a same excel file:
#> ----------------------------------------

library(openxlsx)

#define sheet names for each data frame
dataset_names <- list('enrich_data' = enrich_data, 
                      'lonelyfish_data' = lonelyfish_data, 
                      'lonely_enrich_data' = lonely_enrich_data,
                      'stand_enrich_data' = stand_enrich_data, 
                      'summary_bmd_clusters' = summary_bmd_clusters, 
                      'hox_gene_data' = hox_gene_data,
                      'cry_gene_data' = cry_gene_data)

#export each data frame to separate sheets in same Excel file
openxlsx::write.xlsx(dataset_names, 
                     file = paste0("figures/for-supp/supp-table-data-drerio-", Sys.Date(), ".xlsx"),
                     rowNames = FALSE)

#> ----------------------------------------