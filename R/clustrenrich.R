#' Cluster-specific functional enrichment and whole annotation retrieval
#' 
#' @description
#' Perform Over-Representation Analysis (ORA) on gene clusters to identify enriched biological functions.  
#' This function wraps around `gprofiler2::gost()` and provides flexible options, including:  
#' - custom background gene lists and background types (e.g. "custom", "custom_annotated")  
#' - choice of annotation sources (e.g. GO, KEGG, WP)  
#' - multiple p-value correction methods  
#' - exclusion of IEA (Electronically Inferred Annotations) for GO terms
#' 
#' The function supports a wide range of organisms and annotation sources. Users may rely on standard databases available through g:Profiler and/or extend the analysis by supplying custom GMT annotation file(s).
#' 
#' Users can filter terms/pathways based on gene set size (min_term_size and max_term_size) and the number of genes enriched (ngenes_enrich_filtr). For example, terms with fewer than the minimum required genes or more than the maximum allowed genes are excluded, and terms enriched by fewer than the specified number of genes are filtered out. 
#' 
#' Additionally, users can choose to retain only highlighted/driver GO terms to reduce redundancy and focus on key biological functions. A secondary gprofiler2::gost() run with significant = FALSE retrieves annotations for all deregulated genes, which is utilised later in the lonelyfishing() function. Throughout the process, a dataframe tracks the number of biological functions linked to each cluster after each filtering step, categorized by source. All main parameters used are saved in the output for transparency and reproducibility.
#'
#' @param clustrfiltr_data The named `list` output from the `clustrfiltr()` function. 
#' @param dr_genes The character vector of deregulated genes that can correspond to the `gene_id` column in the output of the `getids()` or  `getregs()`  function. The `gprofiler2::gost()` function handles mixed types of gene IDs and even duplicates by treating them as a single unique occurrence of the identifier, disregarding any duplication.
#' @param bg_genes The character vector of background  genes (preferably from the experiment) that typically corresponds to the `gene_id` column in the output of the `getids()` function.
#' @param bg_type The background type, i.e. the statistical domain, that can be one of "annotated", "known", "custom" or "custom_annotated"
#' @param sources A vector of data sources to use. Currently, these are set at GO:BP, KEGG and WP. Visit the g:GOSt web tool for the comprehensive list and details on incorporated data sources.
#' @param organism Organism ID defined for the chosen sources (e.g. zebrafish = "drerio")
#' @param user_threshold Adjusted p-value cutoff for Over-Representation analysis (default at 0.05 in `gost()` function)
#' @param correction_method P-value adjustment method: one of “gSCS” ,“fdr” and “bonferroni (default set at "fdr")
#' @param exclude_iea Option to exclude GO electronic annotations (IEA)
#' @param only_highlighted_GO Whether to retain only highlighted driver GO terms in the results. Default is set to TRUE.
#' @param min_term_size Minimum size of gene sets to be included in the analysis. If NULL (default), no filtering by size is applied.
#' @param max_term_size Maximum size of gene sets to be included in the analysis. If NULL (default), no filtering by size is applied.
#' @param ngenes_enrich_filtr Minimum number of genes in a cluster needed for a gene set to be considered enriched. If NULL (default), no filtering by gene count is applied.
#' @param gmt_file_paths If provided, these will be uploaded to g:Profiler and included in the enrichment analysis. For guidance on creating and validating GMT files, see the g:Profiler GMT Helper: https://biit.cs.ut.ee/gmt-helper/.
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
    bg_type = c("custom_annotated", "custom", "annotated", "known"),
    sources =  c("GO:BP", "KEGG", "WP"),
    organism, 
    user_threshold = 0.05,  
    correction_method = c("fdr", "g_SCS", "bonferroni", "false_discovery_rate",
                          "analytical"),
    exclude_iea = FALSE, 
    enrich_size_filtr = TRUE, 
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
      "No annotation sources available. ",
      "Provide either standard g:Profiler sources (e.g. GO:BP, KEGG, WP) ",
      "or a valid GMT annotation file (gmt_file_paths)."
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
    
    # Message stating that the directory is created
    message("Directory '", path, "' created.")
  }
  
  # Check if the output file already exists locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, output_filename)) && !overwrite) {
    
    # Message stating that the results file already exists and will be read
    message("The enrichment results file already exists. Reading this file. Use 'overwrite = TRUE' to replace it.")
    
    # Read the existing file and return it 
    return(readRDS(file.path(path, output_filename)))
    
  }
  
  # Select the necessary columns from kept clusters data: Ensembl gene ID and cluster ID, removing duplicates
  dr_g_clustrfiltr_data <- clustrfiltr_data$kept |> 
    dplyr::select(gene_id, clustr) |> 
    dplyr::group_by(clustr) 
  
  # Set the "clustr" column values to factors for proper list names
  dr_g_clustrfiltr_data$clustr <- factor(dr_g_clustrfiltr_data$clustr,
                                         levels = unique(dr_g_clustrfiltr_data$clustr))
  
  # Convert the "gene_id" column values to characters
  dr_g_clustrfiltr_data$gene_id <- as.character(dr_g_clustrfiltr_data$gene_id)
  
  # Create a list by dividing gene set data ("gene_id") into groups ("clustr"). This is reassembled in the form of a list of vectors containing the values for the groups. This list format is input for gprofiler2::gost() function, allowing for gost() analysis on each cluster
  cluster_list <- split(dr_g_clustrfiltr_data$gene_id, 
                        dr_g_clustrfiltr_data$clustr)
  
  # Ensure gprofiler2::gost() inputs (dr and bg) are duplicate-free. This will increase gprofiler2::gost() speed.
  cluster_list <- lapply(cluster_list, unique)
  dr_genes <- unique(dr_genes)
  bg_genes <- unique(bg_genes)
  
  # Initialise list to store all gost results
  all_gost_results <- list()
  
  message("Performing cluster-wise functional enrichment...")
  
  # Perform standard Over-Representation Analysis using gprofiler2::gost()
  if (!is.null(sources) && length(sources) > 0) {
    message(paste("For standard sources:", paste(sources, collapse = ", ")))
    
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
    
    # Format the results
    if (!is.null(multi_gostres$result) && nrow(multi_gostres$result) > 0) {
      multi_gostres$result <- multi_gostres$result |> 
        dplyr::mutate(query = as.numeric(query)) |> 
        dplyr::arrange(query)
    }
    
    all_gost_results[["standard"]] <- multi_gostres
  } else {
    
    message(paste("No standard sources chosen (sources == NULL or empty)"))
    
    # Initialise empty result structure with proper format
    multi_gostres <- list(
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
          query = cluster_list,
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
          # Format results and add source identifier
          gmt_gostres$result <- gmt_gostres$result |>
            dplyr::mutate(
              query = as.numeric(query),
              source = paste0("GMT:", basename(gmt_path)),
              highlighted = FALSE  # GMT files don't have highlighting
            ) |>
            dplyr::arrange(query)
          
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
      multi_gostres$result <- do.call(rbind, all_results)
      multi_gostres$result <- multi_gostres$result |>
        dplyr::arrange(query, source)
      
    }
  } else if (length(all_gost_results) == 1) {
    
    multi_gostres <- all_gost_results[[1]]
  }
  
  # Handle case where no results are found at all
  if (is.null(multi_gostres$result) || nrow(multi_gostres$result) == 0) {
    stop(
      "No enrichment results found. Please check your input genes, background, organism, or try running with significant = FALSE."
    )
    
  } else {
    # Add enrichment ratios to all results (both standard and GMT)
    multi_gostres$result <- multi_gostres$result |> 
      dplyr::mutate(enrichment_ratio = (intersection_size/query_size) / (term_size/effective_domain_size))
  }
  
  # ---------------
  # using the "multi_gostres$result" the function gradually creates one of the components the output of clustrenrich(), which is the $c_simplifylog. It holds the number of biological functions enriched per cluster, before and after each filtering step for each source. In this first block, we create the base dataframe with no filters performed yet.
  
  # Update sources list to include GMT sources for tracking
  all_sources <- if (!is.null(sources) && length(sources) > 0) sources else character(0)
  if (nrow(gmt_annotations_all) > 0) {
    gmt_sources <- unique(gmt_annotations_all$source)
    all_sources <- c(all_sources, gmt_sources)
  }
  
  # Only proceed with tracking if we have results
  if (nrow(multi_gostres$result) > 0) {
    # Create initial term count tracking
    column_mapping <- setNames(paste0("all_", gsub(":", "_", all_sources)), all_sources)
    
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
  } else {
    # Create empty tracking dataframe if no results
    dr_c_a_termcount <- data.frame(clustr = unique(clustrfiltr_data$kept$clustr))
  }
  # ---------------
  
  message("---")
  
  # Conditionally filter out non-highlighted GO terms
  if (only_highlighted_GO == TRUE && nrow(multi_gostres$result) > 0) {
    
    multi_gostres$result <- multi_gostres$result |> 
      dplyr::mutate(query = as.numeric(query)) |> 
      dplyr::arrange(query)
    
    # Filter only standard GO terms (not GMT sources) for highlighting
    multi_gostres$result <- multi_gostres$result |> 
      dplyr::filter(
        ((grepl("^GO", source) & !grepl("^GMT:", source) & highlighted == TRUE) | 
           (!(grepl("^GO", source) & !grepl("^GMT:", source))))
      )
    
    # ----------------
    # When we filter out non-highlighted GO terms, we want to know how many highlighted GO terms are left, and add it to our "dr_c_a_termcount" dataframe
    
    # Select only the cluster ('query'), term_name and source columns.
    dr_c_a_high_go <- multi_gostres$result |> 
      dplyr::select(query, term_name, source) |> 
      dplyr::filter(grepl("^GO", source) & !grepl("^GMT:", source))
    
    if (nrow(dr_c_a_high_go) > 0) {
      column_mapping <- setNames(paste0("driver_", gsub(":", "_", dr_c_a_high_go$source)), dr_c_a_high_go$source)
      
      driver_dr_c_a_termcount <- dr_c_a_high_go |>
        dplyr::group_by(query, source) |>
        dplyr::summarize(count = dplyr::n_distinct(term_name), .groups = "drop") |>
        tidyr::pivot_wider(names_from = source, values_from = count, values_fill = 0) |> 
        dplyr::rename(clustr = query) |> 
        as.data.frame()
      
      for (original_name in names(column_mapping)) {
        if (original_name %in% colnames(driver_dr_c_a_termcount)) {
          new_name <- column_mapping[[original_name]]
          driver_dr_c_a_termcount <- dplyr::rename(driver_dr_c_a_termcount, !!new_name := !!rlang::sym(original_name))
        }
      }
      
      dr_c_a_termcount <- merge(dr_c_a_termcount, driver_dr_c_a_termcount, by = 'clustr', all = TRUE)
      dr_c_a_termcount <- dr_c_a_termcount |> 
        dplyr::mutate_all(~replace(., is.na(.), 0))
    }
    # ----------------
    
    message(paste0("Non-highlighted GO terms are removed "))
    
  } else {
    
    message(paste0("All GO terms are kept "))
    
  }
  
  # Create a filtered version of the `multi_gostres` data for further analysis
  multi_gostres_filtr <- multi_gostres
  
  # Check if both 'min_term_size' and 'max_term_size' are NULL
  if (is.null(min_term_size) && is.null(max_term_size)) {
    
    # Both parameters are NULL, so no filtering is applied
    message("Both `min_term_size` and `max_term_size` are NULL. No gene set size filtering ")
    
  } else if (nrow(multi_gostres_filtr$result) > 0) {
    
    # If 'min_term_size' is provided (not NULL) and 'max_term_size' is NULL
    if (!is.null(min_term_size) && is.null(max_term_size)) {
      
      # Filter the data to include only gene sets with size greater than or equal to `min_term_size`
      multi_gostres_filtr$result <- multi_gostres_filtr$result |>
        dplyr::filter(term_size >= min_term_size)
      
      message(paste0("Filtered gene sets sizes for at least: ", min_term_size, " genes "))
      
    }
    
    # If 'max_term_size' is provided (not NULL) and 'min_term_size' is NULL
    if (is.null(min_term_size) && !is.null(max_term_size)) {
      
      # Filter the data to include only gene sets with size less than or equal to `max_term_size`
      multi_gostres_filtr$result <- multi_gostres_filtr$result |>
        dplyr::filter(term_size <= max_term_size)
      
      message(paste0("Filtered gene sets sizes for at most: ", max_term_size, "genes "))
      
    }
    
    # If both 'min_term_size' and 'max_term_size' are provided (not NULL)
    if (!is.null(min_term_size) && !is.null(max_term_size)) {
      
      # Filter the data to include only gene sets with size between `min_term_size` and `max_term_size`
      multi_gostres_filtr$result <- multi_gostres_filtr$result |>
        dplyr::filter(term_size >= min_term_size & term_size <= max_term_size)
      
      message(paste0("Filtered gene sets sizes for: ", min_term_size, " to ", max_term_size, " genes "))
      
    }
  }
  
  # Handle case where no results exist or all results were filtered out
  if (is.null(multi_gostres_filtr$result) || nrow(multi_gostres_filtr$result) == 0) {
    
    # Create empty g_a dataframe with proper structure
    dr_g_a_gostres <- data.frame(
      gene_id = character(0),
      clustr = numeric(0),
      term_name = character(0),
      term_id = character(0),
      source = character(0),
      stringsAsFactors = FALSE
    )
    
  } else {
    
    # Transform the dataframe from c_a to g_a format. The intersection column holds the gene ids that intersect between the query and term. The term_name and term_id are the respective names and ids for each term. The source represents the databases from which the term is from.
    dr_g_a_gostres <- multi_gostres_filtr$result |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::select(intersection, query, term_name, term_id, source) |> 
      dplyr::rename(gene_id = intersection, 
                    clustr = query)
  }
  
  # Save original gost results for summary purposes
  dr_g_a_gostres_og <- dr_g_a_gostres
  
  # Conditionally remove biological functions that are not sufficiently enriched by a cluster
  if (!is.null(ngenes_enrich_filtr) && nrow(dr_g_a_gostres) > 0) {
    
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
    
    message(paste0("Filtered gene sets enriched by at least: ", ngenes_enrich_filtr, " genes "))
    message("--- ")
    
  } else { 
    
    # Transform the clustr column to numeric to order the data by clustr
    if (nrow(dr_g_a_gostres) > 0) {
      dr_g_a_gostres$clustr <- as.numeric(dr_g_a_gostres$clustr)
      dr_g_a_gostres <- dr_g_a_gostres[order(dr_g_a_gostres$clustr), ]
    }
    
    # The parameter is NULL, so no filtering is applied
    message("`ngenes_enrich_filtr` is NULL. No gene set enrichment size filtering ")
    message("--- ")
  }
  
  
  ## Print out a summary of the enrichment results with or without filtering :
  # Print the ratio of clusters participating in enrichment
  message(length(unique(dr_g_a_gostres_og$clustr)), "/",  length(unique(clustrfiltr_data$kept$clustr)), " clusters participating in enrichment")
  
  # Print the ratio of terms removed because of gene set filtering
  original_terms <- if (nrow(multi_gostres$result) > 0) length(unique(multi_gostres$result$term_name)) else 0
  after_size_filter <- if (nrow(dr_g_a_gostres_og) > 0) length(unique(dr_g_a_gostres_og$term_name)) else 0
  message(after_size_filter, "/", original_terms, " enriched terms kept after gene set size filters")
  
  # Print the ratio of terms removed because of filtering
  after_enrich_filter <- if (nrow(dr_g_a_gostres) > 0) length(unique(dr_g_a_gostres$term_name)) else 0
  message(after_enrich_filter, "/", after_size_filter, " enriched terms kept after enrichment size filter")
  
  
  # --------------
  # After filtering the gene set sizes and the enrichment sizes of biological functions, we also want to get the number of occurences for each cluster and source combination.
  
  if (nrow(dr_g_a_gostres) > 0) {
    # Define the mapping of original column names to new names based on sources
    if (only_highlighted_GO == TRUE) {
      
      # Get sources from the actual results instead of just the input sources
      result_sources <- unique(dr_g_a_gostres$source)
      column_mapping <- setNames(paste0("kept_", ifelse(grepl("^GO", result_sources) & !grepl("^GMT:", result_sources), "driver_", ""), gsub(":", "_", result_sources)), result_sources)
      
    } else {
      
      result_sources <- unique(dr_g_a_gostres$source)
      column_mapping <- setNames(paste0("kept_", gsub(":", "_", result_sources)), result_sources)
      
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
  }
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
  
  # Retrieve annotations using the gprofiler::gost() function:
  if (!is.null(sources) && length(sources) > 0) {
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
  } else {
    dr_g_a_annots <- data.frame(
      gene_id = character(0),
      term_name = character(0),
      term_id = character(0),
      source = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  # Add GMT annotations
  if (nrow(gmt_annotations_all) > 0) {
    # Filter GMT annotations to only include genes in dr_genes
    gmt_annots_filtered <- gmt_annotations_all |>
      dplyr::filter(gene_id %in% dr_genes)
    
    if (nrow(gmt_annots_filtered) > 0) {
      dr_g_a_annots <- rbind(dr_g_a_annots, gmt_annots_filtered)
    }
  }
  
  # Remove duplicate rows
  dr_g_a_annots <- unique(dr_g_a_annots)
  
  # Reset the row numbers
  rownames(dr_g_all_data) <- NULL
  rownames(dr_g_a_annots) <- NULL
  rownames(dr_c_a_termcount) <- NULL
  
  # Create a list of the results
  clustr_enrichres <- list(
    
    dr_g_a_enrich = dr_g_all_data,
    
    gostres = multi_gostres,
    
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
      ngenes_enrich_filtr = ngenes_enrich_filtr,
      gmt_file_paths = gmt_file_paths
    )
    
  )
  
  # Define the class of the output
  clustr_enrichres <- structure(clustr_enrichres, class = "clustrenrichres")
  
  # Save the output
  saveRDS(clustr_enrichres, file.path(path, output_filename))
  
  # Return the clustrenrich results
  return(clustr_enrichres)
  
  
}
