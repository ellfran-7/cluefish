#' Fishing lonely genes into clusters based mapping of annotation to enrichment
#'
#' @description
#' This function identifies and integrates lonely genes (genes not found in any cluster after the `clustrfiltr()` function) into existing clusters. Integration relies on shared gene pathway annotations between lonely genes and cluster pathway enrichments. Lonely genes are pinpointed from the input gene data (`dr_data`) and then linked to clusters based on shared biological annotations from pathways.
#' 
#' @param dr_data A `dataframe` of type *t* that typically corresponds to the output of the `getregs()` function. This input holds at least 'gene_id' and 'term_name' columns, respectively containing Ensembl gene identifiers and biological function annotations for the deregulated genes. 
#' @param clustrenrich_data The named `list` output of the `clustrenrich()` function.
#' @param clustrfusion_data The named `list` output of the `clustrfusion()` function. 
#' @param friendly_limit The maximum number of clusters a gene can be part of to be considered "Friendly". Genes exceeding this limit are assigned to a separate "Friendly" cluster. If the limit is set to 0, the "Friendly" cluster isn't created (default is set to 0)
#' @param path Destination folder for the output data results.
#' @param output_filename Output lonelyfishing result filename.
#' @param overwrite If `TRUE`, the function overwrites existing output files; otherwise, it reads the existing file. (default is set to `FALSE`).
#' 
#' @return A named `list` holding 2 components, where :
#'      -`dr_t_c_a_fishing` is a dataframe of type *t_c_a* holding the lonely fishing results.
#'      -`dr_c_a_fishing` is a dataframe of type *c_a* holding the lonely fishing results. It shares a similar structure to the *clustrfusion_data$dr_c_a_fusion* dataframe with each row being a combination of cluster ID and biological function annotation.
#' 
#' @export
#'
#' @examples
#' 

lonelyfishing <- function(
    dr_data,
    clustrenrich_data,
    clustrfusion_data,
    friendly_limit = 0,
    path,
    output_filename,
    overwrite = FALSE
)
{
  # Check if the output file already exists locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, output_filename)) && !overwrite) {
    
    message("The lonely fishing results file already exists. Reading this file. Use 'overwrite = TRUE' to replace it.")
    
    # Read the existing file and return it
    readRDS(file.path(path, output_filename))
    
  } else {
    
    # If the output file does not exist or overwrite is set to TRUE, proceed with the fishing
    # Filtering all genes not part of clusters : the lonely genes
    dr_t_no_clustr <- dr_data[!(dr_data$gene_id %in% clustrfusion_data$dr_g_a_fusion$gene_id),]
    
    # Filtering annotations based on enriched biological functions
    # In the fishing process, the function only considers terms and pathways that are enriched. In this case, Keep rows with term_names present in cluster fusion data
    dr_g_a_annots_filtr <- clustrenrich_data$dr_g_a_whole |>
      dplyr::filter(term_name %in% clustrfusion_data$dr_g_a_fusion$term_name)
    
    # Merge lonely gene data with enriched function annotation data utilizing an inner_join operation. This ensures that only rows existing in both dataframes are retained from the first dataframe. Given that multiple terms are typically linked to one gene, the relationship is designated as many-to-many.
    dr_g_a_no_clustr <- dplyr::inner_join(dr_t_no_clustr, dr_g_a_annots_filtr, by = "gene_id", relationship = "many-to-many")
    
    # Select columns for fishing: "gene_id" to fish and "term_name" for fishing. Remove duplicates.
    dr_g_a_no_clustr <- dr_g_a_no_clustr |>
      dplyr::select(gene_id, term_name) |>
      dplyr::distinct()
    # This is the dataframe of lonely genes with Ensembl gene IDs and associated term names, susceptible to be fished into clusters.
    
    
    ## Expand clusters by fishing lonely genes sharing the same driver-GO and other (KEGG, WP) annotations as a cluster
    # Remove "gene_id" column to avoid conflicts, then remove duplicate rows to pass from "g_a" to "a".
    dr_a_fusion_modif <- clustrfusion_data$dr_g_a_fusion |>
      dplyr::select(-gene_id) |>
      dplyr::distinct()
    
    # Merge modified cluster fusion data and annotated lonely gene data
    dr_g_a_lonely_data <- merge(dr_g_a_no_clustr, dr_a_fusion_modif, by = "term_name")
    
    # Get biologicaly annotated genes without clusters for the following messages
    dr_g_a_no_clustr_annotated <- dr_data |>
      dplyr::filter((!gene_id %in% clustrfusion_data$dr_g_a_fusion$gene_id) &
                      gene_id %in% clustrenrich_data$dr_g_a_whole$gene_id)
    
    # Get genes without clusters for the following messages
    dr_g_a_no_clustr <- dr_data |>
      dplyr::filter((!gene_id %in% clustrfusion_data$dr_g_a_fusion$gene_id))
    
    
    ## Print info about lonely genes fished
    
    # Before fishing summary
    cat("Before Fishing:", "\n")
    
    # Print the number of total genes
    cat("- Total Lonely Genes:", length(unique(dr_g_a_no_clustr$gene_id)), "\n")
    
    # Print the number of annotated genes
    cat("- Annotated Lonely Genes:", length(unique(dr_g_a_no_clustr_annotated$gene_id)), "\n")
    
    # Space
    cat("", "\n")
    
    # Print the number of lonely genes fishes
    cat("Fished", length(unique(dr_g_a_lonely_data$gene_id)), "Lonely Genes !", "\n")
    
    # Space
    cat("", "\n")
    
    # Before fishing summary
    cat("After Fishing:", "\n")
    
    # Print the number of lonely genes remaining
    cat("- Total Lonely Genes:", length(unique(dr_g_a_no_clustr[!dr_g_a_no_clustr$gene_id %in% dr_g_a_lonely_data$gene_id,]$gene_id)), "\n")
    
    # Print the number of annotated lonely genes remaining
    cat("- Annotated Lonely Genes:", length(unique(dr_g_a_no_clustr_annotated[!dr_g_a_no_clustr_annotated$gene_id %in% dr_g_a_lonely_data$gene_id,]$gene_id)), "\n")
    
    # Combine both objects by rows
    dr_g_c_a_lonely_results <- rbind(clustrfusion_data$dr_g_a_fusion, dr_g_a_lonely_data)
    
    ## Dataframe contains all clusters with genes types:
    ## - genes in original cluster but not annotated
    ## - genes in original cluster enriching kept annotations
    ## - genes from different clusters but fused together
    ## - lonely genes not originally in cluster but sharing biological annotation
    
    ## But what are the remaining genes not fished into a cluster?
    
    # Select remaining lonely data
    dr_t_unfished_lonely_data <- dr_data[!(dr_data$gene_id %in% dr_g_c_a_lonely_results$gene_id), ]
    
    # Prepare the unfished gene data for rbind
    dr_g_c_unfished_lonely_4rbind <- data.frame(gene_id = dr_t_unfished_lonely_data$gene_id,
                                                old_clustr = "Lonely",
                                                new_clustr = "Lonely")
    
    # Prepare the lonely result data (what we have fished so far and the cluster fusion results) for rbind
    dr_g_c_a_lonely_results_4rbind <- dr_g_c_a_lonely_results |>
      dplyr::select(gene_id, old_clustr, new_clustr)
    
    # Combine both datasets by row harboring all the deregulated transcript's genes in order to merge with the annotation dataframe
    dr_g_c_a_fishing_4merge <- rbind(dr_g_c_a_lonely_results_4rbind, dr_g_c_unfished_lonely_4rbind)
    
    
    ## Select enriched GO terms, all KEGG pathways, and Wikipathways.
    
    # The function specifically focuses on the driver GO terms to avoid excessive redundancy, as the GO term tree is already highly redundant and present in massive quantities. Including all GO terms would significantly inflate the size of the dataframes and therefore the results to explore.We also want to remove pointless mother root tems such as "KEGG root term" and "WIKIPATHWAYS"
    dr_g_a_annot_data <- clustrenrich_data$dr_g_a_whole |>
      dplyr::filter(term_name %in% clustrfusion_data$dr_g_a_fusion$term_name & source == "GO:BP" |
                      source %in% c("KEGG", "WP") & !(term_name %in% c("KEGG root term", "WIKIPATHWAYS")))
    
    # Combine the lonely fishing results with the driver GO, KEGG and Wikipathways annotations
    dr_g_c_a_fishing <- merge(dr_g_c_a_fishing_4merge, dr_g_a_annot_data, by = "gene_id", all.x = TRUE)
    
    
    ## But what about gene  part of many clusters ? The friendly genes !
    
    # When fishing lonely genes, we retrieve a lot of information, especially considering multiple annotations (GO, KEGG, WP). This results in many genes scattered among multiple clusters, termed "Friendly" genes. We associate these to a "Friendly" cluster based on a limit of cluster membership, unless they are transcription factors.
    
    if (friendly_limit != 0) {
      
      # Group by gene_id, create a "friendliness" column counting unique clusters per gene. If the "friendliness" value exceeds the limit for a specific gene, update the value in new_clustr to "Friendly". Then ungroup.
      dr_g_c_a_fishing <- dr_g_c_a_fishing |>
        dplyr::group_by(gene_id) |>
        dplyr::mutate(friendliness = dplyr::n_distinct(new_clustr)) |>
        dplyr::mutate(new_clustr = dplyr::if_else(friendliness > friendly_limit, "Friendly", as.character(new_clustr))) |>
        dplyr::ungroup()
    } else {
      
      # Group by gene_id, count unique clusters per gene in "friendliness" column. Then ungroup.
      dr_g_c_a_fishing <- dr_g_c_a_fishing |>
        dplyr::group_by(gene_id) |>
        dplyr::mutate(friendliness = dplyr::n_distinct(new_clustr)) |>
        dplyr::ungroup()
    }
    
    
    ## Add remaining information columns: other gene IDs/transcript ID and TF status
    
    # Merge expanded cluster data with getregs() output
    dr_t_c_a_fishing <- merge(dr_g_c_a_fishing, dr_data, by = "gene_id", all = TRUE)
    
    # Order columns and rows by new_clustr, then source, then term_name
    dr_t_c_a_fishing <- dr_t_c_a_fishing |>
      dplyr::select(transcript_id, gene_id, old_clustr, new_clustr, friendliness, everything()) |>
      dplyr::arrange(new_clustr, source, term_name)
    
    # # Position the readable gene name column following the other ids if provided and exists in the data
    # if ("gene_name" %in% names(dr_t_c_a_fishing)) {
    #   dr_t_c_a_fishing <- dr_t_c_a_fishing |>
    #     dplyr::select(transcript_id, gene_id, gene_name, old_clustr, new_clustr, friendliness, everything())
    # }
    # # Position the TF column if provided and exists in the data
    # if ("TF" %in% names(dr_t_c_a_fishing)) {
    #   dr_t_c_a_fishing <- dr_t_c_a_fishing |>
    #     dplyr::select(transcript_id, gene_id, gene_name, TF, old_clustr, new_clustr, friendliness, everything())
    # }
    
    ## Insert columns if they exist in the input data
    optional_columns <- c("gene_name", "description", "TF")
    
    # Use relocate to move optional columns after 'gene_id'
    dr_t_c_a_fishing <- dr_t_c_a_fishing |> 
      dplyr::relocate(dplyr::any_of(optional_columns), .after = gene_id)
    
    ## The results are in a combined gene/biological function annotation per row format, but we also want the results in a combined clustr/biological function annotation per row format.
    
    # Remove "gene_id" column, select columns of interest to create cluster dataset
    dr_c_a_fishing <- dr_g_c_a_fishing |>
      dplyr::select(-gene_id) |>
      dplyr::select(old_clustr, new_clustr, term_name, term_id, source)
    
    # Prepare a dataframe to integrate "term_size" and "highlighted" columns from enrichment results
    gostres_4_merge <- clustrenrich_data$gostres$result |>
      dplyr::select(query, term_name, term_size, highlighted) |>
      dplyr::rename(old_clustr = query)
    
    # Merge dataframes
    dr_c_a_fishing <- merge(dr_c_a_fishing, gostres_4_merge, by = c("old_clustr", "term_name"))
    
    # Remove row repetitions
    dr_c_a_fishing <- unique(dr_c_a_fishing)
    
    # Organize dataframe for better readability: order columns, and rows by "old_clustr"
    dr_c_a_fishing <- dr_c_a_fishing |>
      dplyr::select(old_clustr, new_clustr, term_name,
                    term_size, term_id, source, highlighted) |>
      dplyr::arrange(as.numeric(old_clustr))
    
    # Reset row numbers
    rownames(dr_t_c_a_fishing) <- NULL
    rownames(dr_c_a_fishing) <- NULL
    
    # Create a named list of the lonelyfishing results
    lonelyfishingres <- list(dr_t_c_a_fishing = dr_t_c_a_fishing,
                             dr_c_a_fishing = dr_c_a_fishing)
    
    # Define the class of the output
    lonelyfishingres <- structure(lonelyfishingres, class = "lonelyfishingres")
    
    # Save the output
    saveRDS(lonelyfishingres, file.path(path, output_filename))
    
    # Return the lonelyfishing results
    return(lonelyfishingres)
  }
}
