#' Retrieve item identifiers from Ensembl using BioMart
#' 
#' @description
#' This function connects to the Ensembl database, queries a specified species dataset using the biomaRt package, and returns the results of the query. If the query includes the 'external_gene_name' attribute, representing readable gene symbols, the function examines whether any duplicates of this identifier exist across various Ensembl transcript and gene IDs rows. This situation occurs when a gene symbol corresponds to different Ensembl gene IDs and Ensembl transcript IDs. If duplicates are detected, the function modifies the external gene name to address this discrepancy.
#'
#' @param id_query A vector of transcript IDs that typically corresponds to the background transcript list
#' @param species_dataset The name of the species dataset desired on `ensembl.org`.
#' @param id_filter The transcript ID type used in the query (e.g., "ensembl_transcript_id_version")
#' @param id_attribut A vector of new ID types to retrieve from the Ensembl dataset of the specified species (e.g., c("ensembl_gene_id", "external_gene_name"))
#'
#' @return A `dataframe` of type *t* containing the biomaRt query results.
#' 
#' @export
#'
#' @examples

getids <- function(
    id_query, 
    species_dataset, 
    id_filter,
    id_attribut
    )
{ 
  # Connection to the ENSEMBL database and the species dataset 
  ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = species_dataset) 
  
  # Retrieves the specified attributes from the BioMart database given a set of filters and corresponding query values
  query_results <- biomaRt::getBM(
    attributes = id_attribut,
    filters = id_filter,
    values = id_query,
    mart = ensembl
    )
  
  # Check if the "external_gene_name" identifier is requested in the query
  if ("external_gene_name" %in% id_attribut) {
    ## Assign new names to transcripts lacking gene symbols. In such cases, we assign a new name to these transcripts' genes as "unknown" followed by a number based on the order of appearance of these blank names. For instance, if there are 8 occurances of unknown gene names, the 5th gene symbol will be named "unknown.5".
    
    # Assign "Unknown" names to genes lacking a symbol and number them sequentially
    indices <- which(query_results$external_gene_name == "")
    num_duplicates <- sum(query_results$external_gene_name == "")
    new_names <- paste("unknown", seq_len(num_duplicates), sep = ".")
    query_results$external_gene_name[indices] <- new_names
    
    ## Differentiate gene symbols for multiple transcripts associated with the same gene. To achieve this
    ## We append a unique numerical value to gene symbol names based on duplicate positions in Ensembl gene IDs and Ensembl transcript IDs transcript ids. If two transcripts are linked to the same gene ID and symbol "A", their new names will be "A_1.1" and "A_1.2". If they are linked to different gene IDs but the same symbol "A", their new names will be "A_1.1" and "A_2.1" The numbering does not follow a specific order other than alphabetical..
    
    # Order the gene symbols and retrieve the duplicates
    gene_order <- sort(query_results$external_gene_name)
    
    query_results <- query_results[order(gene_order),]
    
    duplicated_genes <- query_results$external_gene_name[duplicated(query_results$external_gene_name)]
    
    # Check if there are duplicated gene symbols in the query results
    if (length(duplicated_genes) > 0) {
      # Iterate over each duplicated gene symbol
      for (i in 1:length(duplicated_genes)) {
        
        # Determine the indices of the rows associated with the duplicated gene symbol
        indices <- which(query_results$external_gene_name == duplicated_genes[i])
        
        # Extract unique Ensembl gene identifiers in the order they appear
        unique_genes <- unique(query_results$ensembl_gene_id[indices])
        
        # Count the appearances of each Ensembl transcript within each Ensembl gene. The ave() function in R to applies the seq_along() function to each group of Ensembl transcripts associated with the same Ensembl gene id.
        transcript_counts <- with(query_results, ave(ensembl_transcript_id_version[indices], ensembl_gene_id[indices], FUN = seq_along))
        
        if (length(unique_genes) > 1) {
          
          # Map each Ensembl gene to a unique numeric identifier based on the order of appearance
          gene_id_map <- as.numeric(factor(query_results$ensembl_gene_id[indices], levels = unique_genes))
          
          # Add ensembl transcript and gene appearance numbering to the gene symbol column 
          query_results$external_gene_name[indices] <- paste0(query_results$external_gene_name[indices], "_g", gene_id_map, "t", transcript_counts)
          
        } else {
          
          # Add ensembl transcript and gene appearance numbering to the gene symbol column 
          query_results$external_gene_name[indices] <- paste0(query_results$external_gene_name[indices], "_t", transcript_counts)
          
          
        }
      }
    }
  }
  

  return(query_results)
}
