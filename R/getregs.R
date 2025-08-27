#' Retrieve the regulatory status of deregulated genes
#' 
#' @description
#' This function retrieves information on the regulatory status of genes using the files downloaded with the `dl_regulation_data()` function. It firstly checks if the specifed files have beel downloaded, and whether they contain the è `Ensembl`
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
#' \dontrun{
#' # Create example gene data (normally from getids())
#' example_getids_res <- data.frame(
#'   transcript_id = c("ENSDART00000000069.8", "ENSDART00000002164.9", 
#'                     "ENSDART00000001691.8", "ENSDART00000000070.7"),
#'   gene_id = c("ENSDARG00000000068", "ENSDARG00000000069",
#'               "ENSDARG00000001463", "ENSDARG00000008433")
#'   gene_name = c("NHERF1", "DAP", "tdh2", "UNC45B"),
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Add regulatory status (requires downloaded TF/CoTF files)
#' gene_regs <- getregs(
#'   getids_data = example_getids_res,
#'   regulator_file = file.path(tempdir(), "Danio_rerio_TF_example.txt"),
#'   coregulator_file = file.path(tempdir(), "Danio_rerio_Cof_example.txt")
#' )
#' }

getregs <- function(
    getids_data, 
    regulator_file,
    coregulator_file 
    )
{
  # Check if the regulator and co-regulator files exist
  if (!file.exists(regulator_file)) {
    
    stop("The file specified for TF annotations does not exist in this path: ", regulator_file)
    
  }
  
  if (!file.exists(coregulator_file)) {
    
    stop("The file specified for coTF annotations does not exist in this path: ", coregulator_file)
    
  }
  
  # Read the transcription factor and co-factor data files
  g_coregs_data <- read.table(coregulator_file, sep = "\t", header = TRUE)
  g_regs_data <- read.table(regulator_file, sep = "\t", header = TRUE)
  
  # Check if the required columns are present
  if (!"Ensembl" %in% colnames(g_regs_data)) {
    
    stop("The TF file is missing the 'Ensembl' column.")
    
  }
  
  if (!"Ensembl" %in% colnames(g_coregs_data)) {
    
    stop("The coTF file is missing the 'Ensembl' column.")
    
  }
  
  
  # Extract vectors of Ensembl gene identifiers for co-regulators and regulators
  coregs_items <- as.vector(g_coregs_data$Ensembl) 
    regs_items <- as.vector(g_regs_data$Ensembl) 
  
  # Combine both vectors into one for all regulators and co-regulators
  reg_coreg_items <- c(coregs_items, regs_items) 
  
  # Add a logical column 'TF' to indicate if the gene is a regulator or co-regulator
  getids_data$TF <- getids_data$gene_id %in% reg_coreg_items
  
  # Return the updated dataframe
  return(getids_data)
}
