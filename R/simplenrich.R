#' Simple functional enrichment with added filtering
#' 
#' @description
#' This function utilizes gprofiler2::gost() to perform standard functional enrichment analysis on a list of genes of interest. It incorporates filters similar to those in the proposed workflow, enabling users to set limits on gene set sizes (lower and upper), specify the minimum number of genes involved in enrichment, and restrict results to driver GO terms if requested. The function provides flexibility by applying the same filters and features as the cluefish workflow. The output includes both unfiltered and filtered enrichment results, available in two formats: a combination of gene and annotation per row, or annotation per row only.
#'
#' @param input_genes A character vector of genes of interest. The `gprofiler2::gost()` function handles mixed types of gene IDs and even duplicates by treating them as a single unique occurrence of the identifier, disregarding any duplication.
#' @param bg_genes The vector of background Ensembl genes (preferably from the experiment).
#' @param bg_type The background type, i.e. the statistical domain, that can be one of "annotated", "known", "custom" or "custom_annotated"
#' @param sources A vector of data sources to use. Currently, these are set at GO:BP, KEGG and WP.
#' @param organism Organism ID defined for the chosen sources (e.g. if zebrafish = "drerio")
#' @param user_threshold Adjusted p-value cutoff for Over-Representation analysis (default at 0.05 in `gost()` function)
#' @param correction_method P-value adjustment method: one of “gSCS” ,“fdr” and “bonferroni (default set at "fdr")
#' @param exclude_iea Option to exclude GO electronic annotations (IEA)
#' @param min_term_size Minimum gene set size to consider a biological pathway as relevant to the analysis (default set at 5)
#' @param max_term_size Maximum gene set size to consider a biological pathway as relevant to the analysis (default set at 500)
#' @param only_highlighted_GO Option to only keep highlighted driver GO terms in the analysis results (default set at TRUE)
#' @param ngenes_enrich_filtr Minimum number of cluster genes participating in enrichment to consider sufficient for an enriched pathway to be kept in the output (default set)
#' @param path Destination folder for the output data results.
#' @param output_filename Output enrichment result filename.
#' @param overwrite If `TRUE`, the function overwrites existing output files; otherwise, it reads the existing file. (default is set to `FALSE`).
#'
#' @return A named `list` holding two components:
#'      -`unfiltered` is a named `list` holding two sub-components::
#'            - `dr_g_a` is a dataframe of type *g_a* holding the unfiltered enrichment results (e.g. all GO terms, no limits set on gene set size )
#'            - `gostres` is a named list where 'result' contains the data frame with enrichment analysis results, and 'meta' contains metadata necessary for creating a Manhattan plot. This is the original output of a gprofiler2::gost().
#'      -`filtered` is named `list` holding two sub-components:
#'            - `dr_g_a` is a dataframe of type *g_a* holding the filtered enrichment results.
#'            - `dr_a` is a dataframe of type *a* holding the filtered enrichment results.
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
    min_term_size = 5,
    max_term_size = 500, 
    only_highlighted_GO = TRUE,
    ngenes_enrich_filtr = 3, 
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
    
    ## Create dataframe in "gene x annotation per row" (g_a) format for the unfiltered category
    ## The "annotation per row" (a) is already the output of the gprofiler2::gost() function.
    dr_g_a_unfiltered <- gostres$result |> 
      dplyr::select(intersection, term_name, term_size, query_size, 
                    intersection_size, effective_domain_size, p_value, source) |> 
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
    }
    
    # Filter the results by gene set size
    dr_a_filtered <- dr_a_filtered |> 
      dplyr::filter(min_term_size <= term_size & term_size <= max_term_size)
    
    # Keep only rows where the number of IDs in the intersection column is 3 or more
    dr_a_filtered <- dr_a_filtered |> 
      dplyr::mutate(num_genes = sapply(strsplit(intersection, ","), length)) |> 
      dplyr::filter(num_genes >= ngenes_enrich_filtr) |> 
      dplyr::select(-num_genes)  # Remove the temporary column
    
    # Format the dataframe from "annotation per row" (a) to "gene per row": (g_a)
    dr_g_a_filtered <- dr_a_filtered |> 
      dplyr::select(intersection, term_name, term_size, query_size, 
                    intersection_size, effective_domain_size, p_value, source) |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::rename(gene_id = intersection) |> 
      dplyr::distinct() # Remove duplicate rows
    
    # Turn the tibble to a dataframe
    dr_g_a_filtered <- as.data.frame(dr_g_a_filtered)
    
    # Reset the row numbers
    rownames(dr_g_a_unfiltered) <- NULL
    rownames(dr_g_a_filtered) <- NULL
    rownames(dr_a_filtered) <- NULL
    
    # Create a list containing the unfiltered, filtered enrichment results and the parameters
    results <- list(
      
      unfiltered = list(
        dr_g_a = dr_g_a_unfiltered,
        gostres = gostres),
      
      filtered = list(
        dr_g_a = dr_g_a_filtered,
        dr_a = dr_a_filtered),
      
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