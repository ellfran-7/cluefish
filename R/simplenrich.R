#' Simple functional enrichment with added filtering
#' 
#' @description
#' This function utilizes gprofiler2::gost() to perform standard functional enrichment analysis on a list of genes of interest. It incorporates filters similar to those in the cluefish workflow, enabling users to set limits on gene set sizes (lower and upper), specify the minimum number of genes involved in enrichment, and restrict results to driver GO terms if requested. The function provides flexibility by applying the same filters and features as the cluefish workflow. The output includes both unfiltered and filtered enrichment results, available in two formats: a combination of gene and annotation per row, or annotation per row only.
#'
#' @param input_genes A character vector of genes of interest. The `gprofiler2::gost()` function handles mixed types of gene IDs and even duplicates by treating them as a single unique occurrence of the identifier, disregarding any duplication.
#' @param bg_genes The vector of background Ensembl genes (preferably from the experiment).
#' @param bg_type The background type, i.e. the statistical domain, that can be one of "annotated", "known", "custom" or "custom_annotated"
#' @param sources A vector of data sources to use. Currently, these are set at GO:BP, KEGG and WP.
#' @param organism Organism ID defined for the chosen sources (e.g. if zebrafish = "drerio")
#' @param user_threshold Adjusted p-value cutoff for Over-Representation analysis (default at 0.05 in `gost()` function)
#' @param correction_method P-value adjustment method: one of “gSCS” ,“fdr” and “bonferroni (default set at "fdr")
#' @param exclude_iea Option to exclude GO electronic annotations (IEA)
#' @param only_highlighted_GO Whether to retain only highlighted driver GO terms in the results. Default is set to TRUE.
#' @param min_term_size Minimum size of gene sets to be included in the analysis. If NULL (default), no filtering by size is applied.
#' @param max_term_size Maximum size of gene sets to be included in the analysis. If NULL (default), no filtering by size is applied.
#' @param ngenes_enrich_filtr Minimum number of genes in the  input gene list needed for a gene set to be considered enriched. If NULL (default), no filtering by gene count is applied.
#' @param path Destination folder for the output data results.
#' @param output_filename Output enrichment result filename.
#' @param overwrite If `TRUE`, the function overwrites existing output files; otherwise, it reads the existing file. (default is set to `FALSE`).
#'
#' @return A named `list` holding two components:
#'      -`unfiltered` is a named `list` holding two sub-components::
#'            - `dr_g_a` is a dataframe of type *g_a* holding the unfiltered enrichment results (e.g. all GO terms, no limits set on gene set size )
#'            - `gostres` is a named list where 'result' contains the data frame with enrichment analysis results, and 'meta' contains metadata necessary for creating a Manhattan plot. This is the original output of a gprofiler2::gost(), with added enrichment ratios in the 'result' dataframe.
#'      -`filtered` is named `list` holding two sub-components:
#'            - `dr_g_a` is a dataframe of type *g_a* holding the filtered enrichment results.
#'            - `dr_a` is a dataframe of type *a* holding the filtered enrichment results.
#'      -`params` is a list of the main parameters used
#' @export
#'
#' @examples

simplenrich <- function(
    input_genes,
    bg_genes,
    bg_type = "custom_annotated",
    sources = c("GO:BP", "KEGG", "WP"), 
    organism, 
    user_threshold = 0.05,  
    correction_method = "fdr",
    exclude_iea = FALSE, 
    only_highlighted_GO = TRUE,
    min_term_size = NULL,
    max_term_size = NULL, 
    ngenes_enrich_filtr = NULL, 
    path, 
    output_filename, 
    overwrite = FALSE)

{
  
  # Check if the directory exists, and create it if it does not
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    # `recursive = TRUE` creates intermediate directories as needed
    
    message("Directory '", path, "' created.")
  }
  
  # Check if the output file already exists locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, output_filename)) && !overwrite) {
    
    message("The simplenrich() result file already exists. Reading this file. Use 'overwrite = TRUE' to replace it.")
    
    # Read the existing file and return it 
    readRDS(file.path(path, output_filename))
    
  } else {
    
    # Perform Over-Representation Analysis (ORA) using gprofiler2::gost()
    gostres <- gprofiler2::gost(
      query = input_genes, 
      organism = organism, 
      ordered_query = FALSE, 
      multi_query = FALSE, 
      significant = TRUE, 
      exclude_iea = FALSE, 
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
    
    # Compute enrichment ratios, corresponding to the overlap/expected. This allows the us to quantify the magnitude of a pathway overrepresentation.
    gostres$result <- gostres$result |> 
      dplyr::mutate(enrichment_ratio = (intersection_size/query_size) / (term_size/effective_domain_size))
    
    ## Create dataframe in "gene x annotation per row" (g_a) format for the unfiltered category
    ## The "annotation per row" (a) is already the output of the gprofiler2::gost() function.
    dr_g_a_unfiltered <- gostres$result |> 
      dplyr::select(intersection, term_name, term_size, query_size, 
                    intersection_size, effective_domain_size, p_value, enrichment_ratio, source) |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::rename(gene_id = intersection) |> 
      dplyr::distinct() # Remove duplicate rows
    
    # Turn the tibble to a dataframe
    dr_g_a_unfiltered <- as.data.frame(dr_g_a_unfiltered)
    
    ## Create dataframes in "gene x annotation per row" (g_a) and "annotation per row" (a) formats for the filtered category
    
    # Filter the results by only keeping highlighted GO terms
    if (only_highlighted_GO == TRUE) {
      
      dr_a_filtered <- gostres$result |> 
        dplyr::filter(
          ((grepl("GO", source) & highlighted == TRUE) | (!(grepl("GO", source))))
        )
      
      cat(paste0("Only highlighted GO terms are kept \n"))
      
    } else {
      
      dr_a_filtered <- gostres$result
      
      cat(paste0("All GO terms are kept \n"))
      
    }
    
    
    # Check if both 'min_term_size' and 'max_term_size' are NULL
    if (is.null(min_term_size) && is.null(max_term_size)) {
      
      dr_a_size_filtered <- dr_a_filtered
      
      # Both parameters are NULL, so no filtering is applied
      cat("Both `min_term_size` and `max_term_size` are NULL. No gene set size filtering \n")
      
    } else {
      
      # If 'min_term_size' is provided (not NULL) and 'max_term_size' is NULL
      if (!is.null(min_term_size) && is.null(max_term_size)) {
        
        # Filter the data to include only gene sets with size greater than or equal to `min_term_size`
        dr_a_size_filtered <- dr_a_filtered |> 
          dplyr::filter(term_size >= min_term_size)
        
        cat(paste0("Filtered gene sets sizes for at least: ", min_term_size, " genes \n"))
        
      }
      
      # If 'max_term_size' is provided (not NULL) and 'min_term_size' is NULL
      if (is.null(min_term_size) && !is.null(max_term_size)) {
        
        # Filter the data to include only gene sets with size less than or equal to `max_term_size`
        dr_a_size_filtered <- dr_a_filtered |> 
          dplyr::filter(term_size <= max_term_size)
        
        cat(paste0("Filtered gene sets sizes for at most: ", max_term_size, "genes \n"))
        
      }
      
      # If both 'min_term_size' and 'max_term_size' are provided (not NULL)
      if (!is.null(min_term_size) && !is.null(max_term_size)) {
        
        # Filter the data to include only gene sets with size between `min_term_size` and `max_term_size`
        dr_a_size_filtered <- dr_a_filtered |> 
          dplyr::filter(term_size >= min_term_size & term_size <= max_term_size)
        
        cat(paste0("Filtered gene sets sizes for: ", min_term_size, " to ", max_term_size, " genes \n"))
        
      }
    }
    
    
    # Conditionally remove biological functions that are not sufficiently enriched by a cluster
    if (!is.null(ngenes_enrich_filtr)) {
      
      # Keep only rows where the number of IDs in the intersection column is 3 or more
      dr_a_enrich_size_filtered <- dr_a_size_filtered |> 
        dplyr::mutate(num_genes = sapply(strsplit(intersection, ","), length)) |> 
        dplyr::filter(num_genes >= ngenes_enrich_filtr) |> 
        dplyr::select(-num_genes)  # Remove the temporary column
      
      cat(paste0("Filtered gene sets enriched by at least: ", ngenes_enrich_filtr, " genes \n"))
      cat("--- \n")
      
    } else {
      
      dr_a_enrich_size_filtered <- dr_a_size_filtered
      
      # The parameter is NULL, so no filtering is applied
      cat("`ngenes_enrich_filtr` is NULL. No gene set enrichment size filtering \n")
      cat("--- \n")
      
    }
    
    # Format the dataframe from "annotation per row" (a) to "gene per row": (g_a)
    dr_g_a_enrich_size_filtered <- dr_a_enrich_size_filtered |> 
      dplyr::select(intersection, term_name, term_size, query_size, 
                    intersection_size, effective_domain_size, p_value, enrichment_ratio, source) |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::rename(gene_id = intersection) |> 
      dplyr::distinct() # Remove duplicate rows
    
    
    # Print the ratio of terms removed because of gene set filtering
    cat(length(unique(dr_a_size_filtered$term_name)), "/",  length(unique(dr_a_filtered$term_name)), "enriched terms kept after gene set size filters", "\n")
    
    # Print the ratio of terms removed because of filtering
    cat(length(unique(dr_g_a_enrich_size_filtered$term_name)), "/",  length(unique(dr_a_size_filtered$term_name)), "enriched terms kept after enrichment size filter", "\n")
    
    # Turn the tibble to a dataframe
    dr_g_a_enrich_size_filtered <- as.data.frame(dr_g_a_enrich_size_filtered)
    
    # Reset the row numbers
    rownames(dr_g_a_unfiltered) <- NULL
    rownames(dr_g_a_enrich_size_filtered) <- NULL
    rownames(dr_a_filtered) <- NULL
    
    # Create a list containing the unfiltered, filtered enrichment results and the parameters
    results <- list(
      
      unfiltered = list(
        dr_g_a = dr_g_a_unfiltered,
        gostres = gostres),
      
      filtered = list(
        dr_g_a = dr_g_a_enrich_size_filtered,
        dr_a = dr_a_enrich_size_filtered),
      
      params = list(
        bg_type = bg_type,
        sources = sources, 
        user_threshold = user_threshold,
        min_term_size = min_term_size,
        max_term_size = max_term_size,
        only_highlighted_GO = only_highlighted_GO,
        ngenes_enrich_filtr = ngenes_enrich_filtr
      )
      
    )
    
    # Define the class of the output
    results <- structure(results, class = "simplenrichres")
    
    # Save the output
    saveRDS(results, file.path(path, output_filename))
    
    # Return the clustrenrich results
    return(results)
    
  }
  
}
