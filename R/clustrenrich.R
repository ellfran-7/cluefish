#' Cluster-specific functional enrichment and whole annotation retrieval
#' 
#' @description
#' This function performs Over-Representation Analysis (ORA) on clusters to identify enriched biological functions using the clustrenrich() function. It leverages the gprofiler2::gost() function and offers customization options, including the choice of background gene list, background type (e.g., “custom” or “custom_annotated”), database sources (e.g., GO, KEGG, WP), adjusted p-value correction methods, and the option to exclude IEA (Electronically Inferred Annotations) GO terms. The function is adaptable to various organisms and biological annotation sources.
#' 
#' Users can filter terms/pathways based on gene set size (min_term_size and max_term_size) and the number of genes enriched (ngenes_enrich_filtr). For example, terms with fewer than the minimum required genes or more than the maximum allowed genes are excluded, and terms enriched by fewer than the specified number of genes are filtered out. 
#' 
#' Additionally, users can choose to retain only highlighted/driver GO terms to reduce redundancy and focus on key biological functions. A secondary gprofiler2::gost() run with significant = FALSE retrieves annotations for all deregulated genes, which is utilized later in the lonelyfishing() function. Throughout the process, a dataframe tracks the number of biological functions linked to each cluster after each filtering step, categorized by source. All main parameters used are saved in the output for transparency and reproducibility.
#'
#' @param clustrfiltr_data The named `list` output from the `clustrfiltr()` function. 
#' @param dr_genes The character vector of deregulated genes that can correspond to the `gene_id` column in the output of the `getids()` or  `getregs()`  function. The `gprofiler2::gost()` function handles mixed types of gene IDs and even duplicates by treating them as a single unique occurrence of the identifier, disregarding any duplication.
#' @param bg_genes The character vector of background  genes (preferably from the experiment) that typically corresponds to the `gene_id` column in the output of the `getids()` function.
#' @param bg_type The background type, i.e. the statistical domain, that can be one of "annotated", "known", "custom" or "custom_annotated"
#' @param sources A vector of data sources to use. Currently, these are set at GO:BP, KEGG and WP.
#' @param organism Organism ID defined for the chosen sources (e.g. zebrafish = "drerio")
#' @param user_threshold Adjusted p-value cutoff for Over-Representation analysis (default at 0.05 in `gost()` function)
#' @param correction_method P-value adjustment method: one of “gSCS” ,“fdr” and “bonferroni (default set at "fdr")
#' @param exclude_iea Option to exclude GO electronic annotations (IEA)
#' @param only_highlighted_GO Whether to retain only highlighted driver GO terms in the results. Default is set to TRUE.
#' @param min_term_size Minimum size of gene sets to be included in the analysis. If NULL (default), no filtering by size is applied.
#' @param max_term_size Maximum size of gene sets to be included in the analysis. If NULL (default), no filtering by size is applied.
#' @param ngenes_enrich_filtr Minimum number of genes in a cluster needed for a gene set to be considered enriched. If NULL (default), no filtering by gene count is applied.
#' @param path Destination folder for the output data results.
#' @param output_filename Output enrichment result filename.
#' @param overwrite If `TRUE`, the function overwrites existing output files; otherwise, it reads the existing file. (default is set to `FALSE`).
#'
#' @return A named `list` holding 4 components, where :
#'      -`dr_g_a_enrich` is a dataframe of type *g_a* holding the enrichment results with each row being a combination of gene and biological function annotation
#'      -`gostres` is a named list where 'result' contains the data frame with enrichment analysis results, and 'meta' contains metadata necessary for creating a Manhattan plot. This is the original output of a gprofiler2::gost(), with added enrichment ratios in the 'result' dataframe.
#'      -`dr_g_a_whole` is a dataframe of type *g_a* holding all the biological function annotations found in the g:profiler database for all the deregulated genes.
#'      -`c_simplifylog` is a dataframe tracing the number of biological functions enriched per cluster before and after each filtering step for each source
#'      -`params` is a list of the main parameters used
#'
#' @export
#'
#' @examples

clustrenrich <- function(
    clustrfiltr_data,
    dr_genes,
    bg_genes,  
    bg_type = "custom_annotated",
    sources = c("GO:BP", "KEGG", "WP"), 
    organism, 
    user_threshold = 0.05,  
    correction_method = "fdr",
    exclude_iea = FALSE, 
    enrich_size_filtr = TRUE, 
    only_highlighted_GO = TRUE,
    min_term_size = NULL,
    max_term_size = NULL, 
    ngenes_enrich_filtr = NULL, 
    path, 
    output_filename, 
    overwrite = FALSE 
)

{
  
  # Check if the directory exists, and create it if it does not
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    # `recursive = TRUE` creates intermediate directories as needed
    
    # Message stating that the directory is created
    message("Directory '", path, "' created.")
  }
  
  # Check if the output file already exists locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, output_filename)) && !overwrite) {
    
    # Message stating that the results file already exists and will be read
    message("The enrichment results file already exists. Reading this file. Use 'overwrite = TRUE' to replace it.")
    
    # Read the existing file and return it 
    readRDS(file.path(path, output_filename))
    
  } else {
    
    # Select the necessary columns from kept clusters data: Ensembl gene ID and cluster ID, removing duplicates
    dr_g_clustrfiltr_data <- clustrfiltr_data$kept |> 
      dplyr::select(gene_id, clustr) |> 
      dplyr::group_by(clustr) 
    
    # Set the “clustr” column values to factors for proper list names
    dr_g_clustrfiltr_data$clustr <- factor(dr_g_clustrfiltr_data$clustr,
                                           levels = unique(dr_g_clustrfiltr_data$clustr))
    
    # Convert the “gene_id” column values to characters
    dr_g_clustrfiltr_data$gene_id <- as.character(dr_g_clustrfiltr_data$gene_id)
    
    # Create a list by dividing gene set data (“gene_id”) into groups (“clustr”). This is reassembled in the form of a list of vectors containing the values for the groups. This list format is input for gprofiler2::gost() function, allowing for gost() analysis on each cluster
    cluster_list <- split(dr_g_clustrfiltr_data$gene_id, 
                          dr_g_clustrfiltr_data$clustr)
    
    
    # Perform Over-Representation Analysis (ORA) using gprofiler2::gost()
    multi_gostres <- gprofiler2::gost(
      query = cluster_list, 
      organism = organism, 
      ordered_query = FALSE, 
      multi_query = FALSE, 
      significant = TRUE, 
      exclude_iea = exclude_iea, 
      measure_underrepresentation = FALSE, 
      evcodes = TRUE,
      user_threshold = user_threshold, 
      correction_method = correction_method, 
      domain_scope = bg_type,
      custom_bg = bg_genes, 
      numeric_ns = "", 
      sources = sources, 
      as_short_link = FALSE, 
      highlight = TRUE 
    ) 
    
    # Format the results for convenient exploration
    multi_gostres$result <- multi_gostres$result |> 
      dplyr::mutate(query = as.numeric(query)) |> 
      dplyr::arrange(query)
    
    # Compute enrichment ratios, corresponding to the overlap/expected. This allows the us to quantify the magnitude of a pathway overrepresentation.
    multi_gostres$result <- multi_gostres$result |> 
      dplyr::mutate(enrichment_ratio = (intersection_size/query_size) / (term_size/effective_domain_size))
    
    # ---------------
    # using the "multi_gostres$result" the function gradually creates one of the components the output of clustrenrich(), which is the $c_simplifylog. It holds the number of biological functions enriched per cluster, before and after each filtering step for each source. In this first block, we create the base dataframe with no filters performed yet.
    
    # Define the mapping of original column names to new names based on sources
    column_mapping <- setNames(paste0("all_", gsub(":", "_", sources)), sources)
    
    # Create a dataframe holding the number of terms (term_name) enriched by clusters before and after each condition. This will be updated after each sub-step of this function.
    dr_c_a_gostres <- multi_gostres$result |> 
      dplyr::select(query, term_name, source)
    
    # Group the data by clustr ('query') and biological function database ('source'), count the number of distinct term names per combination of query and source, turn the source column categories into columns associated to their proper count by cluster
    dr_c_a_termcount <- dr_c_a_gostres |>
      dplyr::group_by(query, source) |>
      dplyr::summarize(count = dplyr::n_distinct(term_name), .groups = "drop") |>
      tidyr::pivot_wider(names_from = source, values_from = count, values_fill = 0) |> 
      as.data.frame()
    
    # Rename the columns based on the mapping if they exist
    for (original_name in names(column_mapping)) {
      if (original_name %in% colnames(dr_c_a_termcount)) {
        new_name <- column_mapping[[original_name]]
        dr_c_a_termcount <- dplyr::rename(dr_c_a_termcount, !!new_name := !!rlang::sym(original_name))
      }
    }
    
    # Rename the 'query' column to 'clustr'
    dr_c_a_termcount <- dplyr::rename(dr_c_a_termcount, clustr = query)
    # ---------------
    
    # Conditionally filter out non-highlighted GO terms
    if (only_highlighted_GO == TRUE) {
      
      multi_gostres$result <- multi_gostres$result |> 
        dplyr::filter(
          ((grepl("GO", source) & highlighted == TRUE) | (!(grepl("GO", source))))
        )
      
      # ----------------
      # When we filter out non-highlighted GO terms, we want to know how many highlighted GO terms are left, and add it to our "dr_c_a_termcount" dataframe
      
      # Select only the cluster ('query'), term_name and source columns.
      dr_c_a_high_go <- multi_gostres$result |> 
        dplyr::select(query, term_name, source) |> 
        dplyr::filter(grepl("GO", source))
      
      # Generate column mappings for GO categories
      column_mapping <- setNames(paste0("driver_", gsub(":", "_", dr_c_a_high_go$source)), dr_c_a_high_go$source)
      
      # Group the data by clustr ('query') and biological function database ('source'), count the number of distinct term names per combination of query and source, turn the source column categories into columns associated to their proper count by cluster
      driver_dr_c_a_termcount <- dr_c_a_high_go |>
        dplyr::group_by(query, source) |>
        dplyr::summarize(count = dplyr::n_distinct(term_name), .groups = "drop") |>
        tidyr::pivot_wider(names_from = source, values_from = count, values_fill = 0) |> 
        dplyr::rename(clustr = query) |> 
        as.data.frame()
      
      # Rename the columns based on the mapping if they exist
      for (original_name in names(column_mapping)) {
        if (original_name %in% colnames(driver_dr_c_a_termcount)) {
          new_name <- column_mapping[[original_name]]
          driver_dr_c_a_termcount <- dplyr::rename(driver_dr_c_a_termcount, !!new_name := !!rlang::sym(original_name))
        }
      }
      
      # Merge the previous dr_c_a_termcount with the dataframe also holding driver_GO term occurrence counting, arranging by cluster and keeping all rows from both dataframes
      dr_c_a_termcount <- merge(dr_c_a_termcount, driver_dr_c_a_termcount, by = 'clustr', all = TRUE)
      
      # Replace all NA values with 0 to handle cases where enrichment filtering removed clusters
      dr_c_a_termcount <- dr_c_a_termcount |> 
        dplyr::mutate_all(~replace(., is.na(.), 0))
      # ----------------
      
      cat(paste0("Only highlighted GO terms are kept \n"))
      
    } else {
      
      cat(paste0("All GO terms are kept \n"))
      
    }
    
    # Create a filtered version of the `multi_gostres` data for further analysis
    multi_gostres_filtr <- multi_gostres
    
    # Check if both 'min_term_size' and 'max_term_size' are NULL
    if (is.null(min_term_size) && is.null(max_term_size)) {
      
      # Both parameters are NULL, so no filtering is applied
      cat("Both `min_term_size` and `max_term_size` are NULL. No gene set size filtering \n")
      
    } else {
      
      # If 'min_term_size' is provided (not NULL) and 'max_term_size' is NULL
      if (!is.null(min_term_size) && is.null(max_term_size)) {
        
        # Filter the data to include only gene sets with size greater than or equal to `min_term_size`
        multi_gostres_filtr$result <- multi_gostres_filtr$result |>
          dplyr::filter(term_size >= min_term_size)
        
        cat(paste0("Filtered gene sets sizes for at least: ", min_term_size, " genes \n"))
        
      }
      
      # If 'max_term_size' is provided (not NULL) and 'min_term_size' is NULL
      if (is.null(min_term_size) && !is.null(max_term_size)) {
        
        # Filter the data to include only gene sets with size less than or equal to `max_term_size`
        multi_gostres_filtr$result <- multi_gostres_filtr$result |>
          dplyr::filter(term_size <= max_term_size)
        
        cat(paste0("Filtered gene sets sizes for at most: ", max_term_size, "genes \n"))
        
      }
      
      # If both 'min_term_size' and 'max_term_size' are provided (not NULL)
      if (!is.null(min_term_size) && !is.null(max_term_size)) {
        
        # Filter the data to include only gene sets with size between `min_term_size` and `max_term_size`
        multi_gostres_filtr$result <- multi_gostres_filtr$result |>
          dplyr::filter(term_size >= min_term_size & term_size <= max_term_size)
        
        cat(paste0("Filtered gene sets sizes for: ", min_term_size, " to ", max_term_size, " genes \n"))
        
      }
    }
    
    # Transform the dataframe from c_a to g_a format. The intersection column holds the gene ids that intersect between the query and term. The term_name a,d term_id are the respective names and ids for each term. The source represents the databases from which the term is from.
    dr_g_a_gostres <- multi_gostres_filtr$result |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::select(intersection, query, term_name, term_id, source) |> 
      dplyr::rename(gene_id = intersection, 
                    clustr = query)
    
    # Conditionally remove biological functions that are not sufficiently enriched by a cluster
    if (!is.null(ngenes_enrich_filtr)) {
      
      # Group the data by cluster and term name, count the number of unique Ensembl gene IDs per combination of cluster and term name, then ungroup the data to remove groupings      
      dr_g_a_termcount <- dr_g_a_gostres  |> 
        dplyr::group_by(clustr, term_name) |> 
        dplyr::summarise(n_genes = length(unique(gene_id)), .groups = "drop") |> 
        dplyr::ungroup()  
      
      # Filter out biological functions enriched by less than the specified number of genes
      dr_g_a_termfiltr <- subset(dr_g_a_termcount, n_genes >= ngenes_enrich_filtr)
      
      # Merge the filtered biological functions with the gprofiler2::gost() results
      dr_g_a_termkept <- merge(dr_g_a_gostres, dr_g_a_termfiltr, by = c("clustr", "term_name"))
      
      # Remove the "n_genes" column
      dr_g_a_termkept$n_genes <- NULL
      
      # Remove repeated rows
      dr_g_a_termkept <- unique(dr_g_a_termkept)
      
      # Order the data by cluster and keep the selected columns
      dr_g_a_termkept <- dr_g_a_termkept |> 
        dplyr::select(gene_id, clustr, term_name, term_id, source)
      
      # Transform the cluster column to numeric to order the data by cluster
      dr_g_a_termkept$clustr <- as.numeric(dr_g_a_termkept$clustr)
      dr_g_a_termkept <- dr_g_a_termkept[order(dr_g_a_termkept$clustr), ]
      
      # Create a copy of the dr_g_gostres df in order to avoid squashing
      dr_g_a_gostres_og <- dr_g_a_gostres
      
      # Update the dr_g_a_gostres df with the filtered and sorted data
      dr_g_a_gostres <- dr_g_a_termkept
      
      cat(paste0("Filtered gene sets enriched by at least: ", ngenes_enrich_filtr, " genes \n"))
      cat("--- \n")
      
    } else { 
      
      # Transform the clustr column to numeric to order the data by clustr
      dr_g_a_gostres$clustr <- as.numeric(dr_g_a_gostres$clustr)
      dr_g_a_gostres <- dr_g_a_gostres[order(dr_g_a_gostres$clustr), ]
      
      # The parameter is NULL, so no filtering is applied
      cat("`ngenes_enrich_filtr` is NULL. No gene set enrichment size filtering \n")
      cat("--- \n")
    }
    
    
    ## Print out a summary of the enrichment results with or without filtering :
    # Print the ratio of clusters participating in enrichment
    cat(length(unique(dr_g_a_gostres_og$clustr)), "/",  length(unique(clustrfiltr_data$kept$clustr)), "clusters participating in enrichment", "\n")
    
    # Print the ratio of terms removed because of gene set filtering
    cat(length(unique(dr_g_a_gostres_og$term_name)), "/",  length(unique(multi_gostres$result$term_name)), "enriched terms kept after gene set size filters", "\n")
    
    # Print the ratio of terms removed because of filtering
    cat(length(unique(dr_g_a_gostres$term_name)), "/",  length(unique(dr_g_a_gostres_og$term_name)), "enriched terms kept after enrichment size filter", "\n")
    
    
    # --------------
    # After filtering the gene set sizes and the enrichment sizes of biological functions, we also want to get the number of occurences for each cluster and source combination.
    
    # Define the mapping of original column names to new names based on sources
    if (only_highlighted_GO == TRUE) {
      
      column_mapping <- setNames(paste0("kept_", ifelse(grepl("^GO", sources), "driver_", ""), gsub(":", "_", sources)), sources)
      
    } else {
      
      column_mapping <- setNames(paste0("kept_", gsub(":", "_", sources)), sources)
      
    }
    
    # Group and count the number of distinct term names per combination of clustr and source
    dr_g_a_terms_afterfiltr <- dr_g_a_gostres |> 
      dplyr::select(clustr, term_name, source)
    
    # Group the data by clustr ('query') and biological function database ('source'), count the number of distinct term names per combination of query and source, turn the source column categories into columns associated to their proper count by cluster
    kept_dr_c_a_termcount <- dr_g_a_terms_afterfiltr |>
      dplyr::group_by(clustr, source) |>
      dplyr::summarize(count = dplyr::n_distinct(term_name), .groups = "drop") |>
      tidyr::pivot_wider(names_from = source, values_from = count, values_fill = 0) |> 
      as.data.frame()
    
    # Rename the columns based on the mapping if they exist
    for (original_name in names(column_mapping)) {
      if (original_name %in% colnames(kept_dr_c_a_termcount)) {
        new_name <- column_mapping[[original_name]]
        kept_dr_c_a_termcount <- dplyr::rename(kept_dr_c_a_termcount, !!new_name := !!rlang::sym(original_name))
      }
    }
    
    
    # Merge the previous dr_c_a_termcount with the dataframe also holding driver_GO term occurrence counting, arranging by cluster and keeping all rows from both dataframes
    dr_c_a_termcount <- merge(dr_c_a_termcount, kept_dr_c_a_termcount, by = 'clustr', all = TRUE)
    
    # Replace all NA values with 0 to handle cases where enrichment filtering removed clusters
    dr_c_a_termcount <- dr_c_a_termcount |> 
      dplyr::mutate_all(~replace(., is.na(.), 0))
    # --------------
    
    
    
    # After conducting enrichment analysis, the dataset now exclusively contains genes and data involved in the enrichment of terms. Consequently, genes within a string cluster but not engaged in any enrichment are absent from the direct results. These genes must be reintegrated into the results.
    
    # Retrieve data for Ensembl gene IDs not involved in any enrichment
    dr_g_no_enrich <- clustrfiltr_data$kept[!(clustrfiltr_data$kept$gene_id %in% dr_g_a_gostres$gene_id),]
    
    # Construct a dataframe for the rbind operation with the GO enrichment output
    dr_g_no_enrich_4_rbind <- data.frame(gene_id = dr_g_no_enrich$gene_id,
                                         clustr = dr_g_no_enrich$clustr,
                                         term_name = NA,
                                         term_id = NA,
                                         source = NA)
    
    # Combine both results: enriching and non-enriching genes
    dr_g_all_data <- rbind(dr_g_a_gostres, dr_g_no_enrich_4_rbind)
    
    # Remove duplicate rows
    dr_g_all_data <- unique(dr_g_all_data)
    
    # Order the data by cluster
    dr_g_all_data <- dr_g_all_data[order(dr_g_all_data$clustr),]
    
    
    # STEP 7 bis - Retrieve annotations for the deregulated genes using the gprofiler2::gost() function
    # =================================================================================================
    
    # Retrieve GO:BP, KEGG, and WP annotations using the gprofiler::gost() function:
    annot_gostres <- gprofiler2::gost(
      query = dr_genes, 
      organism = organism, 
      ordered_query = FALSE, 
      multi_query = FALSE, 
      significant = FALSE, 
      exclude_iea = exclude_iea, 
      measure_underrepresentation = FALSE, 
      evcodes = TRUE,
      user_threshold = user_threshold, 
      correction_method = correction_method, 
      domain_scope = bg_type,
      custom_bg = bg_genes, 
      numeric_ns = "", 
      sources = sources, 
      as_short_link = FALSE, 
      highlight = FALSE 
    ) 
    
    # Transform the dataframe from "cluster per row" to "gene per row":
    dr_g_a_annots <- annot_gostres$result |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::select(intersection, term_name, term_id, source) |> 
      dplyr::rename(gene_id = intersection)
    
    # Remove duplicate rows
    dr_g_a_annots <- unique(dr_g_a_annots)
    
    # Reset the row numbers
    rownames(dr_g_all_data) <- NULL
    rownames(dr_g_a_annots) <- NULL
    rownames(dr_c_a_termcount) <- NULL
    
    # Create a list of the results
    clustr_enrichres <- list(
      
      dr_g_a_enrich = dr_g_all_data,
      
      gostres = multi_gostres_filtr,
      
      dr_g_a_whole = as.data.frame(dr_g_a_annots),
      
      c_simplifylog = dr_c_a_termcount,
      
      params = list(
        bg_type = bg_type,
        sources = sources, 
        user_threshold = user_threshold,
        correction_method = correction_method,
        min_term_size = min_term_size,
        max_term_size = max_term_size,
        only_highlighted_GO = only_highlighted_GO,
        ngenes_enrich_filtr = ngenes_enrich_filtr
      )
      
    )
    
    # Define the class of the output
    clustr_enrichres <- structure(clustr_enrichres, class = "clustrenrichres")
    
    # Save the output
    saveRDS(clustr_enrichres, file.path(path, output_filename))
    
    # Return the clustrenrich results
    return(clustr_enrichres)
    
  }
  
}