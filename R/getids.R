#' Retrieve item identifiers from Ensembl using BioMart
#' 
#' @description
#' This function connects to the Ensembl database, queries a specified species dataset using the biomaRt package, and returns the results of the query.
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
  # Connection to the ENSEMVL database and the species dataset 
  ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = species_dataset) 
  
  # Retrieves the specified attributes from the BioMart database given a set of filters and corresponding query values
  query_results <- biomaRt::getBM(id_attribut, 
                            id_filter,
                            id_query,
                            mart = ensembl)

  return(query_results)
}
