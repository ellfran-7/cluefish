#' Retreive and integrate the clustered Protein-Protein Interaction Network (PPIN) data 
#' 
#' @description
#' This function retrieves and reformats the clustered PPIN data from the StringApp in Cytoscape. 
#'
#' @param getregs_data A dataframe that can correspond to the output of the `getregs()` function. This input holds at least one column named ”ensembl_gene_id” holding Ensembl identifiers for the deregulated genes.
#' @param string_clustr_file The path with the filename of the `node table` .csv file previously downloaded following clustering of the PPIN made using the StringApp in Cytoscape. This node table must hold a `query.term` and `X__mclCluster` column.
#' @return A `dataframe` similar to the *getregs_data* dataframe input with an added column indicating to which cluster belongs each gene (if no cluster is associated : NA)
#' 
#' @export
#'
#' @examples

getclustrs <- function(
    getregs_data, 
    path,
    nodetable_filename
    )
{
  # Read the 'node table' CSV file of the created network from StringApp in Cytoscape
  dr_g_string_clustr <- read.csv(file.path(path, nodetable_filename))
  
  # Reformat the query term column by removing STRING identifiers from gene names (e.g., "\"ENSDARG00000042520\"" -> "ENSDARG00000042520")
  dr_g_string_clustr$query.term <- stringr::str_replace_all(dr_g_string_clustr$query.term, "\"", "")
  
  # Create a dataframe for seamless merge with the 'getregs_data' dataframe
  dr_g_string_clustr <- data.frame(ensembl_gene_id = dr_g_string_clustr$query.term, 
                                   clustr = dr_g_string_clustr$X__mclCluster)
  
  # Create a 'clustr_data' dataframe similar to 'getregs_data' but with an added cluster ID column. This allows us to have dose-response modelling metrics to illustrate the PPIN network
  dr_t_clustr_data <- merge(getregs_data, dr_g_string_clustr, by = "ensembl_gene_id")
  
  # Remove genes not associated with a cluster 
  dr_t_clustr_data <- subset(dr_t_clustr_data, !is.na(clustr)) 
  
  return(dr_t_clustr_data)
}
