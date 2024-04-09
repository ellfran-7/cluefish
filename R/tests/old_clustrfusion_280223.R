#' Biological function based cluster fusion
#' 
#' @description
#' This function fusions clusters that share a common single biological function enrichment.
#'
#' @param clustr.ora.filtr.data A table with at least the columns `clustr`, and `go_term`, `kegg_pathway` and/or `clustr` 
#' respectively holding KEGG pathway annotations and cluster identifiers for the selected genes.
#' The input should be the output of the `clustrenrich` function.
#' @return A `.txt` file similar to the input but certain genes belong to a different cluster, 
#' based on the fusion of their old cluster to their new cluster.
#'
#' @examples
#' 

clustrfusion <- function(
    clustr.ora.filtr.data
)
{
  
  # Group by cluster and source, count the number of unique biological functions for each cluster and source
  term.count <- clustr.ora.filtr.data |>
    dplyr::group_by(clustr, source) |>
    dplyr::summarize(n_terms = dplyr::n_distinct(term_name), .groups = 'drop')
  
  # Get the clusters that have one unique biological function for each source
  single_clusters <- term.count |>
    dplyr::filter(n_terms == 1) |>
    dplyr::pull(clustr)
  
  # Get the unique biological function names for each source
  single_terms <- unique(clustr.ora.filtr.data[clustr.ora.filtr.data$clustr %in% single_clusters, ]$term_name)
  
  # Initialize an empty vector to store clusters to be merged
  single_terms_clust <- c()
  
  # Loop through unique biological functions
  for (i in single_terms) {
    
    for (source_type in unique(clustr.ora.filtr.data$source)) {
      # Check if more than 1 cluster is associated with the biological function for the current source
      if (length(unique(clustr.ora.filtr.data |> 
                        dplyr::filter(term_name == i & clustr %in% single_clusters & source == source_type) |>
                        dplyr::pull(clustr))) > 1) {
        # Append clusters to be merged to the vector
        single_terms_clust <- append(single_terms_clust, as.character(unique(clustr.ora.filtr.data |>
                                                                               dplyr::filter(term_name == i & clustr %in% single_clusters & source == source_type) |>
                                                                               dplyr::pull(clustr))))
      }
    }
  }
  
  # Create a new data frame to store the merged clusters
  temp.clustr.ora.filtr.data <- clustr.ora.filtr.data
  
  # Set to keep track of fused clusters
  fused_clusters <- c()
  
  # Loop through the single pathways and merge corresponding clusters for each source
  for (term in single_terms) {
    for (source_type in unique(temp.clustr.ora.filtr.data$source)) {
      clusters_to_merge <- temp.clustr.ora.filtr.data |> 
        dplyr::filter(term_name == term & clustr %in% single_terms_clust & source == source_type) |>
        dplyr::pull(clustr)
      
      # Check if any of the clusters have already been fused
      if (any(clusters_to_merge %in% fused_clusters)) {
        cat("Clusters for", source_type, "term", "*", term, "* have already been merged.\n")
        next  # Skip to the next source_type
      }
      
      if (length(unique(clusters_to_merge)) > 1) {
        
        # Find the minimum numeric cluster : fused cluster takes the smallest cluster id 
        new_cluster <- min(clusters_to_merge)
        
        sharing.clusters <- temp.clustr.ora.filtr.data |>
          dplyr::filter(clustr %in% clusters_to_merge)
        
        sharing.clusters$clustr <- new_cluster
        
        # Print the merged clusters, and their shared enrichment
        cat("Merged clusters enriching the", source_type, "term", "*", term, "* :", unique(clusters_to_merge), "\n")
        
        # Concatenate the information from the clusters being merged
        temp.clustr.ora.filtr.data <- dplyr::bind_rows(temp.clustr.ora.filtr.data, sharing.clusters)
        
        # Keep track of fused clusters
        fused_clusters <- union(fused_clusters, clusters_to_merge)
        
        # Identify other clusters with shared biological annotations
        clusters_with_shared_annotation <- temp.clustr.ora.filtr.data |>
          dplyr::filter(term_name %in% sharing.clusters$term_name & term_name %in% single_terms) |>
          dplyr::filter(!clustr %in% clusters_to_merge) |> 
          dplyr::pull(clustr)
        
        if (length(unique(clusters_with_shared_annotation)) > 0) {
          
          # Allow clusters with shared annotation to participate in subsequent merges
          fused_clusters <- setdiff(fused_clusters, clusters_to_merge)
        }
        
        # Keep the vector of cluster ids that have been merged to another, and therefor will be removed as cluster
        clusters_to_merge <- setdiff(clusters_to_merge, min(clusters_to_merge))
        
        # Remove clusters that have merged and are not the final cluster id
        temp.clustr.ora.filtr.data <- temp.clustr.ora.filtr.data |>
          dplyr::filter(!clustr %in% clusters_to_merge)
        
        temp.clustr.ora.filtr.data <- unique(temp.clustr.ora.filtr.data)
        
        temp.clustr.ora.filtr.data <- temp.clustr.ora.filtr.data[order(temp.clustr.ora.filtr.data$clustr),]
        
      }
    }
  }
  
  return(temp.clustr.ora.filtr.data)
  
}
