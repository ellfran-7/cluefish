#' Download transcription factors and co-factors data from AnimalTFDB
#'
#' @description
#' This function downloads both transcription factor and cofactor data from the AnimalTFDB4 database. 
#' The downloaded files are saved as .txt files in the specified directory. 
#' If the files already exist locally and `overwrite` is not set to TRUE, the function will not download them again.
#' The files will be stored in the `data/derived-data/` folder, which will be created if necessary.
#' 
#' @param url_tf The URL for the transcription factor data.
#' @param url_cof The URL for the transcription co-factor data.
#' @param path The directory for saving the downloaded files.
#' @param filename_tf The name associated with the downloaded transcription factor file.
#' @param filename_cof The name associated with the downloaded co-factor factor file.
#' @param overwrite If `TRUE`, the function overwrites existing output files; otherwise, it reads the existing file. (default is set to `FALSE`).
#'
#' @return No return value; files are downloaded and stored in the specified directory.
#' @export
#'
#' @examples

dl_regulation_data <- function(
    url_tf,
    url_cof,
    path,
    filename_tf,
    filename_cof,
    overwrite = FALSE
    ) 
  {
  
  # Check if the directory exists, and create it if it does not
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    # `recursive = TRUE` creates intermediate directories as needed
    
    # Message stating that the directory is created
    message("Directory '", path, "' created.")
  }
  
  # Check if files exist locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, filename_tf)) & !overwrite) {
    
    # Message stating that the files already exist, and 'overwrite = TRUE' as to be set in order to replace them
    message("The TF and Cof filenames already exists. Use 'overwrite = TRUE' to replace them")
    
  } else {
    
    ## Create destination folder for data if it doesn't exist
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    
    ## Download TF file
    utils::download.file(url      = paste0(url_tf),
                         destfile = file.path(path, filename_tf),
                         mode     = "wb")
    
    ## Download Cof file
    utils::download.file(url      = paste0(url_cof),
                         destfile = file.path(path, filename_cof),
                         mode     = "wb")
    
  }
  
  ## No return value
  invisible(NULL) 
  
}


