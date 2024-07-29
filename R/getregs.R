#' Retrieve the regulatory status of deregulated genes
#' 
#' @description
#' This function retrieves information on the regulatory status of genes using the files downloaded with the `dl_regulation_data()` function.
#' 
#' @param getids_data A dataframe of type *t* that typically corresponds to the output of the `getids()` function. This input holds at least one column named "gene_id" holding gene identifiers for the deregulated genes.  It should be a subset of the output of the `getids()` function with only the significantly deregulated transcripts as rows and transcript/gene identifiers as columns.
#' @param regulator_file The previously downloaded `.txt` Transcription Factor file 
#' from the AnimalTFDB4 database (https://guolab.wchscu.cn/AnimalTFDB4/#/).
#' @param coregulator_file The previously downloaded `.txt` Transcription Co-Factor file from the AnimalTFDB4 database (https://guolab.wchscu.cn/AnimalTFDB4/#/).
#' 
#' @return A `dataframe` of type *t* similar to the *getids_data* input with an added column indicating the (co-)regulator status of each transcript's deregulated gene.
#' 
#' @export
#'
#' @examples

getregs <- function(
    getids_data, 
    regulator_file,
    coregulator_file 
    )
{
  # Read the files containing Transcription Factors and Co-Regulators data
  g_coregs_data <- read.table(coregulator_file, sep = "\t", header = TRUE)
  g_regs_data <- read.table(regulator_file, sep = "\t", header = TRUE)
  
  # Create the vector of the Ensembl genes that are co-regulators
  coregs_items <- as.vector(g_coregs_data$Ensembl) 
  
  # Create the vector of the Ensembl genes that are regulators
  regs_items <- as.vector(g_regs_data$Ensembl) 
  
  # Combine both vectors to include all regulators and co-regulators
  reg_coreg_items <- c(coregs_items, regs_items) 
  
  # Create a logical column 'TF' indicating whether the gene is either a known regulator or co-regulator
  getids_data$TF <- getids_data$gene_id %in% reg_coreg_items
  
  # Return the results
  return(getids_data)
}
