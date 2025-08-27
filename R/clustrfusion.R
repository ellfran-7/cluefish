#' Fusion of clusters based on shared cluster enrichment
#' 
#' @description
#' This function fusions clusters that share the same biological function enrichment in at-least one of the sources used (e.g. GO, KEGG...).
#'
#' @param clustrenrich_data The named `list` output from the `clustrenrich()` function.
#' @param monoterm_fusion Option to merge clusters only if they either enrich the same total terms (monoterm_fusion = FALSE) or the same individual term (monoterm_fusion = TRUE) in at least one source (default set to FALSE).
#' 
#' @return A named `list` holding 4 components, where :
#'      -`dr_g_a_fusion` is a dataframe of type *g_a* holding the cluster fusion results. It shares a similar structure to the *clustrenrich_data$dr_g_a_enrich* dataframe with each row being a combination of gene and biological function annotation.
#'      -`dr_c_a_fusion` is a dataframe of type *c_a* holding the cluster fusion results with each row being a combination of cluster ID and biological function annotation
#'      -`c_fusionlog` is a dataframe tracing cluster fusion events, indicating the sources from which they originated (e.g. GO, KEGG).
#'      -`params` is a list of the main parameters used; in this case monoterm_fusion
#' 
#' @examples
#' 

clustrfusion <- function(
    clustrenrich_data,
    monoterm_fusion = FALSE
)
{
  
  # Extract the enrichment results for gene data
  dr_g_a_enrich <- clustrenrich_data$dr_g_a_enrich
  
  # Add new columns "old_clustr" and "new_clustr" for fusion process: this is so that we perform the fusion on the 'new_clustr" column, keeping the 'old_clustr' intact
  dr_g_a_enrich <- dr_g_a_enrich |> 
    dplyr::mutate(old_clustr = clustr,
                  new_clustr = clustr) |> 
    dplyr::select(-clustr)
  
  # Check if there are any enrichment results to work with
  if (is.null(clustrenrich_data$gostres$result) || nrow(clustrenrich_data$gostres$result) == 0) {
    stop("No enrichment results found in the input data. Cluster fusion requires enriched terms to identify clusters that share biological functions. Please run clustrenrich() first with data that produces enrichment results.")
  }
  
  # Arrange the gost() results by effective domain size to prioritize sources with a larger panel of biological information about the organism (details in README). Then, select the source column, remove duplicates, and convert to a vector for use in the for loop.
  ordered_term_sources <- clustrenrich_data$gostres$result |> 
    dplyr::arrange(desc(effective_domain_size)) |> 
    dplyr::select(source) |> 
    dplyr::distinct() |> 
    dplyr::pull(source)
  
  
  # Initialize the fusion log dataframe to track cluster fusion events. This dataframe will record the clusters merged for each source evaluated during fusion.
  c_fusionlog <- data.frame(old_clustr = sort(as.numeric(dr_g_a_enrich$old_clustr)))
  
  # Create a list to store fusion information for each source
  fusion_info <- list()
  
  # Loop through each source
  for (i in 1:length(ordered_term_sources)){
    
    # If the cluster fusion is to be done based on full term enrichment: Combine clusters only if they enrich the same total terms in at-least one source
    if (monoterm_fusion == FALSE) {
      
      # Filter the data to keep only the enrichment for each source.
      # Procedure:
      # - Filter rows with the same value of ordered_term_sources in the "source" column.
      # - Select the columns that are use to identify mergings: new_clustr, term_name and source
      # - Remove duplicated rows as we go from type "g_a" to c_a"
      dr_g_a_enriched_terms <- dr_g_a_enrich |> 
        dplyr::filter(source == ordered_term_sources[i]) |> 
        dplyr::select(new_clustr, term_name, source) |> 
        dplyr::distinct()
      
      # Loop through each cluster
      for (j in 1:length(dr_g_a_enriched_terms$new_clustr)) {
        
        # Retrieve the enriched terms associated to the specific cluster 
        clusterterms <- dr_g_a_enriched_terms[dr_g_a_enriched_terms$new_clustr == dr_g_a_enriched_terms$new_clustr[j],]$term_name
        
        # Determine which clusters enrich all the same terms (all or nothing)
        clusterstomerge <- which(sapply(split(dr_g_a_enriched_terms$term_name, dr_g_a_enriched_terms$new_clustr), function(x) setequal(x, clusterterms)))
        
        # Get the numeric cluster IDs of the clusters that are to merge
        clusterstomerge <- as.numeric(names(clusterstomerge))
        
        # If multiple clusters associated with the term, proceed
        if (length(clusterstomerge) > 1){
          
          # Identify the smallest cluster id, therefore the largest cluster in gene set size
          min_cluster <- min(clusterstomerge)
          
          # Merge clusters by assigning the smallest id to all clusters enriching the term
          dr_g_a_enrich[dr_g_a_enrich$new_clustr %in% clusterstomerge,]$new_clustr <- min_cluster
          
        } else {
          
          next
          
        }
      }
      
      # If the cluster fusion is to be done based on mono-term enrichment: Combine clusters only if they enrich the same individual term in at-least one source
      
    } else {
      
      # Filter the data to keep only one unique enrichment for each source.
      # This line of code aims to identify unique combinations of sources. 
      # Procedure:
      # - Filter rows with the same value of ordered_term_sources in the "source" column.
      # - Group the data by "new_clustr" and "source".
      # - Summarize the groups by creating a new column "term_count" that counts unique occurrences of "term_name" values.
      # - Drop the grouping structure.
      # - Remove NA values if any.
      # - Filter rows where the term_count is equal to 1.
      # - Remove the term_count column and convert the tibble into a regular dataframe.
      dr_g_a_unique_enriched_terms <- dr_g_a_enrich |> 
        dplyr::filter(source == ordered_term_sources[i]) |> 
        dplyr::group_by(new_clustr, source) |> 
        dplyr::summarize(term_count = dplyr::n_distinct(term_name), .groups = "drop") |>
        na.omit() |> 
        dplyr::filter(term_count == 1) |> 
        dplyr::select(-term_count) |> 
        as.data.frame()
      
      # Merge the previously created "dr_g_a_unique_enriched_terms" dataset that holds the info of terms unique to a source within a cluster. This means we filter out all rows where the "term_name" isn't alone to be enriched in the same source.
      dr_g_a_enrich_filtr <- merge(dr_g_a_enrich, dr_g_a_unique_enriched_terms, by = c("new_clustr", "source"))
      
      # Reorder by "new_clustr"
      dr_g_a_enrich_filtr <- dr_g_a_enrich_filtr[order(dr_g_a_enrich_filtr$new_clustr),]
      
      # Retrieve the unique non-NA term names
      term_names <- na.omit(unique(dr_g_a_enrich_filtr$term_name))
      
      # Loop through each term name
      for (j in 1:length(term_names)){
        
        # Select data associated with the unique term within the source
        dr_g_a_unique_enriched_term <- dr_g_a_enrich_filtr[dr_g_a_enrich_filtr$term_name %in% term_names[j],]
        
        # Create a vector of clusters enriching the unique term and remove duplicates
        clusterstomerge <- unique(dr_g_a_unique_enriched_term$new_clustr)
        
        # Get the numeric cluster IDs of the clusters that are to merge
        clusterstomerge <- as.numeric(names(clusterstomerge))
        
        # If multiple clusters associated with the term, proceed
        if (length(clusterstomerge) > 1){
          
          # Identify the smallest cluster id, therefore the largest cluster in gene set size
          min_cluster <- min(clusterstomerge)
          
          # Merge clusters by assigning the smallest id to all clusters enriching the term
          dr_g_a_enrich[dr_g_a_enrich$new_clustr %in% clusterstomerge,]$new_clustr <- min_cluster
          
        } else {
          
          next
          
        }
        
        
      }
      
      
    }
    
    # Store fusion information for each source
    fusion_info[[ordered_term_sources[i]]] <- as.numeric(dr_g_a_enrich$new_clustr)
    
  }
  
  # Add fusion information to the fusion log dataframe
  for (source in names(fusion_info)) {
    c_fusionlog[[paste0("after_", source, "_fusion")]] <- fusion_info[[source]]
  }
  
  # Count the number of clusters before and after fusion
  total_clusters_before_fusion <- length(unique(dr_g_a_enrich$old_clustr))
  total_clusters_after_fusion <- length(unique(dr_g_a_enrich$new_clustr))
  
  # Calculate the number of clusters that participated in fusion
  clusters_participated_in_fusion <- total_clusters_before_fusion - total_clusters_after_fusion
  
  # Print the ratio of clusters after/before the fusion
  cat(total_clusters_after_fusion, "/", total_clusters_before_fusion, "clusters left after the fusion process. \n")
  
  # Order columns
  dr_g_a_fusion <- dr_g_a_enrich |> 
    dplyr::select(gene_id, old_clustr, new_clustr, everything())
  
  # Renew rownames
  rownames(dr_g_a_fusion) <- NULL
  
  # Create cluster dataset from gene dataset
  dr_c_a_fusion <- dr_g_a_fusion |> 
    dplyr::select(-gene_id) |> 
    dplyr::select(old_clustr, new_clustr, term_name, term_id, source)
  
  # Prepare dataframe to integrate term_size and highlighted columns from enrichment results
  gostres_4_merge <- clustrenrich_data$gostres$result |> 
    dplyr::select(query, term_name, term_size, highlighted, source) |> 
    dplyr::rename(old_clustr = query)
  
  # Merge dataframes
  dr_c_a_fusion <- merge(dr_c_a_fusion, gostres_4_merge, by = c("old_clustr", "term_name", "source"))
  
  # Remove row repetitions
  dr_c_a_fusion <- unique(dr_c_a_fusion)
  
  # Organize dataframe for readability
  dr_c_a_fusion <- dr_c_a_fusion |> 
    dplyr::select(old_clustr, new_clustr, term_name,
                  term_size, term_id, source, highlighted) |> 
    dplyr::arrange(as.numeric(old_clustr))
  
  # Reset row numbers
  rownames(dr_g_a_fusion) <- NULL
  rownames(dr_c_a_fusion) <- NULL
  rownames(c_fusionlog) <- NULL
  
  # Create a list of the results
  clustr_fusionres <- list(
    
    dr_g_a_fusion = dr_g_a_fusion,
    
    dr_c_a_fusion = dr_c_a_fusion,
    
    c_fusionlog = c_fusionlog,
    
    params = list(
      monoterm_fusion = monoterm_fusion
    )
  )
  
  # Define the class of the output
  clustr_fusionres <- structure(clustr_fusionres, class = "clustrfusion")
  
  
  # Return the clusterfusion results
  return(clustr_fusionres)
  
}