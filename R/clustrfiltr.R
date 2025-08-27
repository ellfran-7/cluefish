#' Filter clusters based on their gene-set size
#' 
#' @description
#' This function filters clusters based on size, selecting genes that belong to clusters meeting the user-defined size criteria. 
#' 
#' @param getclustrs_data A `dataframe` of type *t* that typically corresponds to the output of the `getclustrs()` function. This input holds at least the columns named `gene_id` and `clustr` respectively holding Ensembl gene and cluster identifiers for the deregulated genes.
#' @param size_filtr The minimum number of genes required for a cluster to be retained (by default: 3). 
#' @return A named `list` holding 2 components, where : 
#'      - `kept` is a dataframe of type *t* similar to the *getclustrs_data* dataframe input with the rows kept after the filter (the genes are part of a cluster equal and over the size limit)
#'      -`removed` is a dataframe of type *t* similar to the *getclustrs_data* dataframe input with the rows removed after the filter (the genes are part of clusters under the size limit)
#'      - `params` is a list of the main parameters used
#' 
#' @export
#'
#' @examples

clustrfiltr <- function(
    getclustrs_data,
    size_filtr = 3
    )
{ 
  # Select necessary columns: Ensembl gene ID and cluster ID, removing duplicate rows in order to pass from transcript per row to gene per row
  dr_g_clustrs <- getclustrs_data |> 
    dplyr::select(gene_id, clustr) |> 
    dplyr::distinct()
  
  # Calculate the distribution of genes within clusters
  cluster_sizes <- table(dr_g_clustrs$clustr)
  
  # Identify clusters with sizes greater than or equal to the specified size filter
  big_clusters <- names(cluster_sizes[cluster_sizes >= size_filtr])
  
  # Print the count of clusters before and after filtering
  message("Total clusters kept: ", length(big_clusters), "/", length(unique(dr_g_clustrs$clustr)))
  
  # Separate clusters into kept and removed based on the size filter to create the list result
  dr_t_clustrs_filtr <- list(kept = getclustrs_data[getclustrs_data$clustr %in% big_clusters, ],
                             removed =  getclustrs_data[!getclustrs_data$clustr %in% big_clusters, ])
  
  # Order the results by cluster number
  dr_t_clustrs_filtr$kept <- dr_t_clustrs_filtr$kept |> 
    dplyr::arrange(clustr)
  
  dr_t_clustrs_filtr$removed <- dr_t_clustrs_filtr$removed |> 
    dplyr::arrange(clustr)
  
  # Reset row names for clarity
  rownames(dr_t_clustrs_filtr$kept) <- NULL
  rownames(dr_t_clustrs_filtr$removed) <- NULL
  
  # Create a list of the results
  clustrfiltr_res <- list(
    
    kept = dr_t_clustrs_filtr$kept,
    
    removed = dr_t_clustrs_filtr$removed,
    
    params = list(
      size_filtr = size_filtr
      )
    
    )
  # Define the result as a structured object with class "clustrfiltres"
  clustrfiltr_res <- structure(clustrfiltr_res, class = "clustrfiltres")
                          
  return(clustrfiltr_res)
}
