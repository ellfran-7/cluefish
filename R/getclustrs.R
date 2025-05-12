#' Retreive and integrate the clustered Protein-Protein Interaction Network (PPIN) data 
#' 
#' @description
#' This function retrieves and reformats the clustered PPIN data from the StringApp in Cytoscape. 
#'
#' @param gene_data A dataframe of type *t* that typically corresponds to the output of the `getids()` or `getregs()` function. This input holds at least one column named "gene_id” holding gene identifiers for the deregulated genes.
#' @param colname_for_merge The identifier column used as query in order to create the PPIN in Cytoscape using the StringApp (e.g. "ensembl_gene_id" or "uniprotsptrembl")
#' @param path Folder where the nodetable exported from Cytoscape is found
#' @param nodetable_filename Filename of the nodetable
#'
#' @return A `dataframe` of type *t* similar to the *gene_data* dataframe input with an added column indicating to which cluster belongs each gene (if no cluster is associated : NA)
#' 
#' @export
#'
#' @examples

getclustrs <- function(
    gene_data,
    colname_for_merge,
    path,
    nodetable_filename
    )
{
  # Read the 'node table' CSV file of the created network from StringApp in Cytoscape
  dr_g_string_clustr <- read.csv(file.path(path, nodetable_filename))
  
  # Check if the following columns are present in the table exported from cytoscape, if not stop the function
  if (any(!c("query.term", "X__mclCluster") %in% colnames(dr_g_string_clustr))) {
    
    stop(
      paste(
      "The exported file from Cytoscape is not the expected node table.",
      "Missing columns:", paste(c("query.term", "X__mclCluster"), collapse = ", "),
      "Check that you exported the node table and not the Network or Edge table.",
      sep = "\n"
      )
    )
    
  }
  
  # Reformat the query term column by removing STRING identifiers from gene names (e.g., "\"ENSDARG00000042520\"" -> "ENSDARG00000042520")
  dr_g_string_clustr$query.term <- stringr::str_replace_all(dr_g_string_clustr$query.term, "\"", "")
  
  # Create a dataframe for seamless merge with the 'getregs_data' dataframe
  dr_g_string_clustr <- data.frame(query_term = dr_g_string_clustr$query.term, 
                                   clustr = dr_g_string_clustr$X__mclCluster)
  
  # Rename the column in the gene_data dataframe to match the 'query_term' column in dr_g_string_clustr
  colnames(dr_g_string_clustr)[1] <- colname_for_merge
  
  # Create a 'clustr_data' dataframe similar to 'getregs_data' but with an added cluster ID column. This allows us to have dose-response modelling metrics to illustrate the PPIN network
  dr_t_clustr_data <- merge(gene_data, dr_g_string_clustr, by = colname_for_merge)
  
  # Remove genes not associated with a cluster 
  dr_t_clustr_data <- subset(dr_t_clustr_data, !is.na(clustr)) 
  
  return(dr_t_clustr_data)
}
