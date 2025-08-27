#' Simple functional enrichment with added filtering
#' 
#' @description
#' This function utilizes gprofiler2::gost() to perform over-representation analysis on a list of genes of interest. Users can choose either standard annotation sources (e.g. GO:BP, KEGG, WP) via the `sources` argument, or provide one or more custom GMT files through `gmt_file_paths`. If both are provided, results from all sources are combined. At least one of `sources` or `gmt_file_paths` must be supplied. 
#'
#' Similar to the `clustrenrich()` workflow, the function supports filters on gene 
#' set size (minimum and maximum), on the minimum number of input genes per 
#' enriched set, and on whether to keep only highlighted driver GO terms. 
#' 
#' The output contains both unfiltered and filtered results, in two formats: 
#' "gene × annotation per row" and "annotation per row".
#'
#' @param input_genes A character vector of genes of interest. The `gprofiler2::gost()` function handles mixed types of gene IDs and even duplicates by treating them as a single unique occurrence of the identifier, disregarding any duplication.
#' @param bg_genes The vector of background Ensembl genes (preferably from the experiment).
#' @param bg_type The background type, i.e. the statistical domain, that can be one of "annotated", "known", "custom" or "custom_annotated"
#' @param sources A vector of data sources to use. Currently, these are set at GO:BP, KEGG and WP. Visit the g:GOSt web tool for the comprehensive list and details on incorporated data sources.
#' @param organism Organism ID defined for the chosen sources (e.g. if zebrafish = "drerio")
#' @param user_threshold Adjusted p-value cutoff for Over-Representation analysis (default at 0.05 in `gost()` function)
#' @param correction_method P-value adjustment method: one of “gSCS” ,“fdr” and “bonferroni (default set at "fdr")
#' @param exclude_iea Option to exclude GO electronic annotations (IEA)
#' @param only_highlighted_GO Whether to retain only highlighted driver GO terms in the results. Default is set to TRUE.
#' @param min_term_size Minimum size of gene sets to be included in the analysis. If NULL (default), no filtering by size is applied.
#' @param max_term_size Maximum size of gene sets to be included in the analysis. If NULL (default), no filtering by size is applied.
#' @param ngenes_enrich_filtr Minimum number of genes in the  input gene list needed for a gene set to be considered enriched. If NULL (default), no filtering by gene count is applied.
#' @param gmt_file_paths If provided, these will be uploaded to g:Profiler and included in the enrichment analysis. For guidance on creating and validating GMT files, see the g:Profiler GMT Helper: https://biit.cs.ut.ee/gmt-helper/.
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
    bg_type = c("custom_annotated", "custom", "annotated", "known"),
    sources = c("GO:BP", "KEGG", "WP"), 
    organism, 
    user_threshold = 0.05,  
    correction_method = c("fdr", "g_SCS", "bonferroni", "false_discovery_rate",
                          "analytical"),
    exclude_iea = FALSE, 
    only_highlighted_GO = TRUE,
    min_term_size = NULL,
    max_term_size = NULL, 
    ngenes_enrich_filtr = NULL, 
    gmt_file_paths = NULL,
    path, 
    output_filename, 
    overwrite = FALSE
    )

{
  
  # Match arguments
  bg_type <- match.arg(bg_type)
  correction_method <- match.arg(correction_method)
  
  # Stop function running if no sources selected or custom GMT file chosen 
  if ((is.null(sources) || length(sources) == 0) && 
      (is.null(gmt_file_paths) || length(gmt_file_paths) == 0)) {
    stop(
      "No annotation sources available. Provide either standard g:Profiler sources (e.g. GO:BP, KEGG, WP) or a valid GMT annotation file (gmt_file_paths)."
    )
  }
  
  # Helper function to read GMT file and extract annotations
  read_gmt_annotations <- function(gmt_path) {
    if (!file.exists(gmt_path)) {
      stop("GMT file not found: ", gmt_path)
    }
    
    lines <- readLines(gmt_path, warn = FALSE)
    if (length(lines) == 0) {
      warning("GMT file is empty: ", gmt_path)
      return(data.frame())
    }
    
    gmt_data <- strsplit(lines, "\t")
    
    results <- lapply(gmt_data, function(parts) {
      if (length(parts) < 3) return(NULL)
      genes <- parts[-c(1, 2)]
      if (length(genes) == 0) return(NULL)
      
      data.frame(
        gene_id = genes,
        term_name = parts[2],
        term_id = parts[1],
        source = paste0("GMT:", basename(gmt_path)),
        stringsAsFactors = FALSE
      )
    })
    
    results <- results[!sapply(results, is.null)]
    if (length(results) == 0) {
      warning("No valid gene sets found in GMT file: ", gmt_path)
      return(data.frame())
    }
    
    do.call(rbind, results)
  }
  
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
    return(readRDS(file.path(path, output_filename)))
    
  }
  
  # Ensure gprofiler2::gost() inputs (dr and bg) are duplicate-free. This will increase gprofiler2::gost() speed.
  input_genes <- unique(input_genes)
  bg_genes <- unique(bg_genes)
  
  # Initialize list to store all gost results
  all_gost_results <- list()
  
  message("Performing simple functional enrichment...")
  
  # Perform standard Over-Representation Analysis using gprofiler2::gost()
  if (!is.null(sources) && length(sources) > 0) {
    message(paste("For standard sources:", paste(sources, collapse = ", ")))
    
    gostres <- gprofiler2::gost(
      query = input_genes, 
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
    
    all_gost_results[["standard"]] <- gostres
  } else {
    
    message(paste("No standard sources chosen (sources == NULL or empty)"))
    
    # Initialize empty result structure with proper format
    gostres <- list(
      result = data.frame(
        query = numeric(0),
        significant = logical(0),
        p_value = numeric(0),
        term_size = numeric(0),
        query_size = numeric(0),
        intersection_size = numeric(0),
        source = character(0),
        term_id = character(0),
        term_name = character(0),
        effective_domain_size = numeric(0),
        source_order = numeric(0),
        parents = character(0),
        evidence_codes = character(0),
        intersection = character(0),
        highlighted = logical(0),
        stringsAsFactors = FALSE
      ), 
      meta = list()
    )
  }
  
  # Process GMT files if provided
  gmt_annotations_all <- data.frame()
  
  if (!is.null(gmt_file_paths) && length(gmt_file_paths) > 0) {
    
    message("For custom uploaded GMT files...")
    
    for (i in seq_along(gmt_file_paths)) {
      gmt_path <- gmt_file_paths[i]
      message(paste(" -"))
      message(paste(" Processing GMT file", i, "of", length(gmt_file_paths), ":", basename(gmt_path)))
      
      tryCatch({
        # Upload GMT file
        gmt_id <- gprofiler2::upload_GMT_file(gmt_path)
        message(paste(" GMT file uploaded with ID:", gmt_id))
        
        # Perform enrichment analysis with GMT file
        gmt_gostres <- gprofiler2::gost(
          query = input_genes,
          organism = gmt_id,
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
          sources = NULL,  # Use all sources from GMT
          as_short_link = FALSE,
          highlight = FALSE  # GMT results don't have highlight
        )
        
        if (!is.null(gmt_gostres$result) && nrow(gmt_gostres$result) > 0) {
          # Add source identifier and highlighted column
          gmt_gostres$result <- gmt_gostres$result |>
            dplyr::mutate(
              source = paste0("GMT:", basename(gmt_path)),
              highlighted = FALSE  # GMT files don't have highlighting
            )
          
          all_gost_results[[paste0("gmt_", i)]] <- gmt_gostres
          message(paste("Found", nrow(gmt_gostres$result), "enriched terms from GMT file"))
        } else {
          message("No enriched terms found in GMT file")
        }
        
        # Read GMT annotations for later use
        gmt_annots <- read_gmt_annotations(gmt_path)
        if (nrow(gmt_annots) > 0) {
          gmt_annotations_all <- rbind(gmt_annotations_all, gmt_annots)
        }
        
      }, error = function(e) {
        warning(paste("Error processing GMT file", gmt_path, ":", e$message))
      })
    }
  }
  
  # Merge all results
  if (length(all_gost_results) > 1) {
    
    # Combine all result dataframes
    all_results <- lapply(all_gost_results, function(x) x$result)
    all_results <- all_results[sapply(all_results, function(x) !is.null(x) && nrow(x) > 0)]
    
    if (length(all_results) > 0) {
      # Ensure all results have the same columns
      all_cols <- unique(unlist(lapply(all_results, names)))
      
      # Add missing columns to each dataframe
      all_results <- lapply(all_results, function(df) {
        missing_cols <- setdiff(all_cols, names(df))
        for (col in missing_cols) {
          if (col == "highlighted") {
            df[[col]] <- grepl("^GO", df$source) & rep(FALSE, nrow(df))  # Default to FALSE for GMT
          } else {
            df[[col]] <- NA
          }
        }
        return(df[, all_cols])
      })
      
      # Combine results
      gostres$result <- do.call(rbind, all_results)
      gostres$result <- gostres$result |>
        dplyr::arrange(source)
      
    }
  } else if (length(all_gost_results) == 1) {
    
    gostres <- all_gost_results[[1]]
  }
  
  # Handle case where no results are found at all
  if (is.null(gostres$result) || nrow(gostres$result) == 0) {
    stop(
      "No enrichment results found. Please check your input genes, background, organism, or try running with significant = FALSE."
    )
    
  } else {
    # Add enrichment ratios to all results (both standard and GMT)
    gostres$result <- gostres$result |> 
      dplyr::mutate(enrichment_ratio = (intersection_size/query_size) / (term_size/effective_domain_size))
  }
  
  ## Create dataframe in "gene x annotation per row" (g_a) format for the unfiltered category
  if (nrow(gostres$result) > 0) {
    dr_g_a_unfiltered <- gostres$result |> 
      dplyr::select(intersection, term_name, term_size, query_size, 
                    intersection_size, effective_domain_size, p_value, source) |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::rename(gene_id = intersection) |> 
      dplyr::distinct() # Remove duplicate rows
  } else {
    dr_g_a_unfiltered <- data.frame(
      gene_id = character(0),
      term_name = character(0),
      term_size = numeric(0),
      query_size = numeric(0),
      intersection_size = numeric(0),
      effective_domain_size = numeric(0),
      p_value = numeric(0),
      source = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  # Turn the tibble to a dataframe
  dr_g_a_unfiltered <- as.data.frame(dr_g_a_unfiltered)
  
  message("---")
  
  ## Create dataframes in "gene x annotation per row" (g_a) and "annotation per row" (a) formats for the filtered category
  
  # Filter the results by only keeping highlighted GO terms
  if (only_highlighted_GO == TRUE && nrow(gostres$result) > 0) {
    
    # Filter only standard GO terms (not GMT sources) for highlighting
    dr_a_filtered <- gostres$result |> 
      dplyr::filter(
        ((grepl("^GO", source) & !grepl("^GMT:", source) & highlighted == TRUE) | 
           (!(grepl("^GO", source) & !grepl("^GMT:", source))))
      )
    
    message(paste0("Non-highlighted GO terms are removed "))
    
  } else {
    
    dr_a_filtered <- gostres$result
    
    message(paste0("All GO terms are kept "))
    
  }
  
  # Check if both 'min_term_size' and 'max_term_size' are NULL
  if (is.null(min_term_size) && is.null(max_term_size)) {
    
    dr_a_size_filtered <- dr_a_filtered
    
    # Both parameters are NULL, so no filtering is applied
    message("Both `min_term_size` and `max_term_size` are NULL. No gene set size filtering ")
    
  } else if (nrow(dr_a_filtered) > 0) {
    
    # If 'min_term_size' is provided (not NULL) and 'max_term_size' is NULL
    if (!is.null(min_term_size) && is.null(max_term_size)) {
      
      # Filter the data to include only gene sets with size greater than or equal to `min_term_size`
      dr_a_size_filtered <- dr_a_filtered |> 
        dplyr::filter(term_size >= min_term_size)
      
      message(paste0("Filtered gene sets sizes for at least: ", min_term_size, " genes "))
      
    }
    
    # If 'max_term_size' is provided (not NULL) and 'min_term_size' is NULL
    if (is.null(min_term_size) && !is.null(max_term_size)) {
      
      # Filter the data to include only gene sets with size less than or equal to `max_term_size`
      dr_a_size_filtered <- dr_a_filtered |> 
        dplyr::filter(term_size <= max_term_size)
      
      message(paste0("Filtered gene sets sizes for at most: ", max_term_size, "genes "))
      
    }
    
    # If both 'min_term_size' and 'max_term_size' are provided (not NULL)
    if (!is.null(min_term_size) && !is.null(max_term_size)) {
      
      # Filter the data to include only gene sets with size between `min_term_size` and `max_term_size`
      dr_a_size_filtered <- dr_a_filtered |> 
        dplyr::filter(term_size >= min_term_size & term_size <= max_term_size)
      
      message(paste0("Filtered gene sets sizes for: ", min_term_size, " to ", max_term_size, " genes "))
      
    }
  } else {
    dr_a_size_filtered <- dr_a_filtered
  }
  
  # Conditionally remove biological functions that are not sufficiently enriched
  if (!is.null(ngenes_enrich_filtr) && nrow(dr_a_size_filtered) > 0) {
    
    dr_a_enrich_size_filtered <- dr_a_size_filtered |> 
      dplyr::filter(intersection_size >= ngenes_enrich_filtr)
    
    message(paste0("Filtered gene sets enriched by at least: ", ngenes_enrich_filtr, " genes "))
    message("--- ")
    
  } else {
    
    dr_a_enrich_size_filtered <- dr_a_size_filtered
    
    # The parameter is NULL, so no filtering is applied
    message("`ngenes_enrich_filtr` is NULL. No gene set enrichment size filtering ")
    message("--- ")
    
  }
  
  # Format the dataframe from "annotation per row" (a) to "gene per row": (g_a)
  if (nrow(dr_a_enrich_size_filtered) > 0) {
    dr_g_a_enrich_size_filtered <- dr_a_enrich_size_filtered |> 
      dplyr::select(intersection, term_name, term_size, query_size, 
                    intersection_size, effective_domain_size, p_value, source) |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::rename(gene_id = intersection) |> 
      dplyr::distinct() # Remove duplicate rows
  } else {
    dr_g_a_enrich_size_filtered <- data.frame(
      gene_id = character(0),
      term_name = character(0),
      term_size = numeric(0),
      query_size = numeric(0),
      intersection_size = numeric(0),
      effective_domain_size = numeric(0),
      p_value = numeric(0),
      source = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  ## Print out a summary of the enrichment results with or without filtering :
  
  # Print the ratio of terms removed because of gene set filtering
  original_terms <- if (nrow(gostres$result) > 0) length(unique(gostres$result$term_name)) else 0
  after_highlight_filter <- if (nrow(dr_a_filtered) > 0) length(unique(dr_a_filtered$term_name)) else 0
  after_size_filter <- if (nrow(dr_a_size_filtered) > 0) length(unique(dr_a_size_filtered$term_name)) else 0
  after_enrich_filter <- if (nrow(dr_a_enrich_size_filtered) > 0) length(unique(dr_a_enrich_size_filtered$term_name)) else 0
  
  message(after_size_filter, "/", after_highlight_filter, " enriched terms kept after gene set size filters")
  message(after_enrich_filter, "/", after_size_filter, " enriched terms kept after enrichment size filter")
  
  # Turn the tibble to a dataframe
  dr_g_a_enrich_size_filtered <- as.data.frame(dr_g_a_enrich_size_filtered)
  dr_a_enrich_size_filtered <- as.data.frame(dr_a_enrich_size_filtered)
  
  # Reset the row numbers
  rownames(dr_g_a_unfiltered) <- NULL
  rownames(dr_g_a_enrich_size_filtered) <- NULL
  rownames(dr_a_enrich_size_filtered) <- NULL

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
      correction_method = correction_method,
      min_term_size = min_term_size,
      max_term_size = max_term_size,
      only_highlighted_GO = only_highlighted_GO,
      ngenes_enrich_filtr = ngenes_enrich_filtr,
      gmt_file_paths = gmt_file_paths
    )
    
  )
  
  # Define the class of the output
  results <- structure(results, class = "simplenrichres")
  
  # Save the output
  saveRDS(results, file.path(path, output_filename))
  
  # Return the results
  return(results)
  
}