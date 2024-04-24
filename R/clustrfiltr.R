#' Gene-set size based filtering of clusters
#' 
#' @description
#' This function removes genes belonging to clusters that are considered "too small" by the user.
#' 
#' @param getclustrs_data A `dataframe` of type *t* that typically corresponds to the output of the `getclustrs()` function. This input holds at least the columns named `ensembl_gene_id` and `clustr` respectively holding Ensembl gene and cluster identifiers for the deregulated genes.
#' @param size_filtr The number of genes in a cluster to consider sufficiently large enough to be a biological complex (by default: 3). 
#' @return A named `list` holding 2 components, where : 
#'      -`kept` is a dataframe of type *t* similar to the *getclustrs_data* dataframe input with the rows kept after the filter (the Ensembl genes are part of a cluster over the size limit)
#'      -`removed` is a dataframe of type *t* similar to the *getclustrs_data* dataframe input with the rows removed after the filter (the Ensembl genes are par tof clusters under the size limit)
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
    dplyr::select(ensembl_gene_id, clustr) |> 
    dplyr::distinct()
  
  # Calculate the distribution of genes within clusters
  cluster_sizes <- table(dr_g_clustrs$clustr)
  
  # Identify clusters with sizes greater than or equal to the specified size filter
  big_clusters <- names(cluster_sizes[cluster_sizes >= size_filtr])
  
  # Print the count of clusters before and after filtering
  cat("Total clusters kept:", length(big_clusters), "/", length(unique(dr_g_clustrs$clustr)), "\n")
  
  # Separate clusters into kept and removed based on the size filter to create the list result
  dr_g_clustrs_filtr <- list(kept = getclustrs_data[getclustrs_data$clustr %in% big_clusters, ],
                             removed =  getclustrs_data[!getclustrs_data$clustr %in% big_clusters, ])
  
  # Reset row names for clarity
  rownames(dr_g_clustrs_filtr$kept) <- NULL
  rownames(dr_g_clustrs_filtr$removed) <- NULL
  
  # Define the result as a structured object with class "clustrfiltres"
  dr_g_clustrs_filtr <- structure(dr_g_clustrs_filtr, class = "clustrfiltres")
                          
  return(dr_g_clustrs_filtr)
}
