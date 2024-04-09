#' Gene-set size based filtering of clusters
#' 
#' @description
#' This function removes genes belonging to clusters that are considered "too small" by the user.
#' 
#' @param clustr.data A table with at least a column named `ensembl_gene_id` and `clustr` 
#' respectively holding Entrez and cluster identifiers for the selected genes.
#' The input should be the output of the `getclustrs` function.
#' @param size.filtr The number of genes in a cluster to be considered sufficiently
#' coherent in a biological perspective (by default: 2).
#' @return A `dataframe` file similar to the input, but only keeping genes part of clusters considered as "large enough" by the user.
#' 
#' @export
#'
#' @examples

clustrfiltr <- function(
    clustr.data,
    size.filtr = 2
)
{ 
  # Loop for removing gene rows part of clusters under the predetermined cluster size limit
  big.clusters <- c() # Empty vector to which the sufficiently large clusters will be appended
  
  for (i in 1:length(unique(clustr.data$clustr)))
  {
    # If the number of gene rows associated with the cluster i is superior to the size filter
    if (length(unique(clustr.data[clustr.data$clustr %in% i, ]$ensembl_gene_id)) > size.filtr)
    {
      # Append the cluster i ID to "big.clusters"
      big.clusters <- append(big.clusters, i)
      
    } else {
      
      # Print the id of filtered-out clusters
      cat("Filtered out cluster", i, "\n")
      
    }
  }
  
  # Print the ratio of clusters before and after the size filter
  cat("Total clusters kept :", length(unique(big.clusters)), "/", length(unique(clustr.data$clustr)))
  
  # Create a dataframe filtered with remaining gene rows part of sufficiently large clusters (> size.filter)
  filtered.clustr <- clustr.data[clustr.data$clustr %in% big.clusters, ]
  
  return(filtered.clustr)
}
