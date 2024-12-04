#' Retrieve gene identifiers from Ensembl/Ensembl Genomes
#' 
#' @description
#' This function connects to the Ensembl or Ensembl Metazoa database, queries a specified species dataset using the biomaRt package, and returns the results of the query. If the query includes the 'external_gene_name' attribute, representing readable gene symbols, the function examines whether any duplicates of this identifier exist across various Ensembl transcript and gene IDs rows. This situation occurs when a gene symbol corresponds to different Ensembl gene IDs and Ensembl transcript IDs. If duplicates are detected, the function modifies the external gene name to address this discrepancy.
#'
#' @param id_query A vector of transcript IDs that typically corresponds to the background transcript list
#' @param biomart_db The name of the BioMart database hosted by Ensembl or Ensembl Metazoa. Use `listEnsembl()` to view available Ensembl databases, or `listEnsemblGenomes()` to view available Ensembl Metazoa databases.
#' @param species_dataset The name of the species dataset desired on `ensembl.org`.
#' @param version The Ensembl version to connect to when wanting to connect to an archived Ensembl version
#' @param transcript_id The transcript identifier for the deregulated transcripts used in the query (e.g., "ensembl_transcript_id_version")
#' @param gene_id The gene identifier to be retrieved from the BioMart dataset of the specified species (preferably "ensembl_gene_id"  as it the identifier used in g:profiler)
#' @param gene_name A human-readable gene name identifier to be retrieved from the BioMart dataset of the specified species.
#' @param other_ids One or more additional identifiers or attributes to retrieve from the BioMart dataset of the specified species (e.g., "external_gene_name", "uniprotsptrembl", "string"). Ensure that the retrieved identifier is supported for the organism in the STRING database for subsequent analysis.

#'
#' @return A `dataframe` of type *t* containing the biomaRt query results with a modified "external_gene_name" column if itself present and duplicates exist.
#' 
#' @export
#'
#' @examples

getids <- function(
    id_query, 
    biomart_db,
    species_dataset,
    version = NULL,
    transcript_id,
    gene_id,
    gene_name = NULL,
    other_ids = NULL 
)
{ 
  # Check if the BioMart database is hosted by Ensembl
  if (biomart_db %in% biomaRt::listEnsembl()$biomart) {
    
    # Connection to a BioMart database and the species dataset
    ensembl <- biomaRt::useEnsembl(biomart = biomart_db, dataset = species_dataset, version = version)
    
    # Check if the BioMart database is hosted by Ensembl Genomes
  } else if (biomart_db %in% biomaRt::listEnsemblGenomes()$biomart) {
    
    # Connection to an Ensembl Genome in the BioMart database and the species dataset
    ensembl <- biomaRt::useEnsemblGenomes(biomart = biomart_db, dataset = species_dataset)
    
  } else {
    
    # Message saying that the BioMart or dataset given is not found
    message("The provided BioMart database or dataset is not found.")
  }
  
  # Retrieves the specified attributes from the BioMart database given a set of filters and corresponding query values
  query_results <- biomaRt::getBM(
    attributes = c(transcript_id, gene_id, gene_name, other_ids),
    filters = transcript_id,
    values = id_query,
    mart = ensembl
  )
  
  # Rename the specific columns so that the column adapted to the workflow
  names(query_results)[names(query_results) == transcript_id] <- "transcript_id"
  names(query_results)[names(query_results) == gene_id] <- "gene_id"
  names(query_results)[names(query_results) == gene_name] <- "gene_name"
  
  # Check if the "external_gene_name" identifier is requested in the query under the "gene_name"
  if ("external_gene_name" %in% gene_name) {
    ## Assign new names to transcripts lacking gene symbols. In such cases, we assign a new name to these transcripts' genes as "unknown" followed by a number based on the order of appearance of these blank names. For instance, if there are 8 occurrences of unknown gene names, the 5th gene symbol will be named "unknown.5".
    
    # Assign "Unknown" names to genes lacking a symbol and number them sequentially
    indices <- which(query_results$gene_name == "" | is.na(query_results$gene_name))
    num_unknown <- sum(query_results$gene_name == "" | is.na(query_results$gene_name))
    new_names <- paste("unknown", seq_len(num_unknown), sep = ".")
    query_results$gene_name[indices] <- new_names
    
    ## Differentiate gene symbols for multiple transcripts associated with the same gene. To achieve this
    ## We append a unique numerical value to gene symbol names based on duplicate positions in Ensembl gene IDs and Ensembl transcript IDs transcript ids. If two transcripts are linked to the same gene ID and symbol "A", their new names will be "A_1.1" and "A_1.2". If they are linked to different gene IDs but the same symbol "A", their new names will be "A_1.1" and "A_2.1" The numbering does not follow a specific order other than alphabetical..
    
    # Order the gene symbols and retrieve the duplicates
    gene_order <- sort(query_results$gene_name)
    
    query_results <- query_results[order(gene_order),]
    
    duplicated_genes <- query_results$gene_name[duplicated(query_results$gene_name)]
    
    # Check if there are duplicated gene symbols in the query results
    if (length(duplicated_genes) > 0) {
      # Iterate over each duplicated gene symbol
      for (i in 1:length(duplicated_genes)) {
        
        # Determine the indices of the rows associated with the duplicated gene symbol
        indices <- which(query_results$gene_name == duplicated_genes[i])
        
        # Extract unique Ensembl gene identifiers in the order they appear
        unique_genes <- unique(query_results$gene_id[indices])
        
        # Count the appearances of each Ensembl transcript within each Ensembl gene. The ave() function in R to applies the seq_along() function to each group of Ensembl transcripts associated with the same Ensembl gene id.
        transcript_counts <- with(query_results, ave(transcript_id[indices], gene_id[indices], FUN = seq_along))
        
        if (length(unique_genes) > 1) {
          
          # Map each Ensembl gene to a unique numeric identifier based on the order of appearance
          gene_id_map <- as.numeric(factor(query_results$gene_id[indices], levels = unique_genes))
          
          # Add ensembl transcript and gene appearance numbering to the gene symbol column 
          query_results$gene_name[indices] <- paste0(query_results$gene_name[indices], "_g", gene_id_map, "t", transcript_counts)
          
        } else {
          
          # Add ensembl transcript and gene appearance numbering to the gene symbol column 
          query_results$gene_name[indices] <- paste0(query_results$gene_name[indices], "_t", transcript_counts)
          
          
        }
      }
    }
  }
  
  
  return(query_results)
}
