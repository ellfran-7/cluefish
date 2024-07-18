#' Generate a workflow summary table for exploration
#'
#' @description
#' This function creates a concise summary table capturing the key details from the workflow results from the `lonelyfishing()` function. It encompasses all essential information for result exploration, striking a balance by avoiding an overwhelming amount of data that might hinder ease of exploration.
#' 
#' @param lonelyfishing_data The named `list` output of the `lonelyfishing()` function.
#' @param bmdboot_data The DRomics bmdboot `dataframe` results after DRomics::bmdfilter() 
#' @param path Destination folder for the output data results.
#' @param output_filename Output CSV filename.
#' @param overwrite If `TRUE`, the function overwrites existing output files; otherwise, it reads the existing file. (default is set to `FALSE`).
#' 
#' @return No return value; the `.csv` file is downloaded and stored in the specified director. It holds the results of the final step of the workflow with each row being a combination of Ensembl transcript, cluster and biological function annotation. This is to be used as a support for exploration and mechanism discovery.
#' 
#' @export
#'
#' @examples
#' 

results_to_csv <- function(
    lonelyfishing_data,
    bmdboot_data,
    path,
    output_filename,
    overwrite = TRUE
    ) 
{
  
  # Check if the output file already exists locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, output_filename)) && !overwrite) {
    
    message("The summary table csv file already exists. Use 'overwrite = TRUE' to replace it.")
    
  } else {
    
    # Merge lonelyfishing results with bmdboots results after DRomics::bmdfilter()
    dr_t_c_a_workflow_res <- merge(lonelyfishing_data$dr_t_c_a_fishing, bmdboot_data, by.x = "transcript_id", by.y = "id")
    
    
    # Prepare the structure for the summary dataframe
    dr_t_c_a_summary <- data.frame(
      transcript_id = dr_t_c_a_workflow_res$transcript_id,
      gene_id = dr_t_c_a_workflow_res$gene_id,
      newcluster = dr_t_c_a_workflow_res$new_clustr,
      friendliness = dr_t_c_a_workflow_res$friendliness,
      term_name = dr_t_c_a_workflow_res$term_name,
      source = dr_t_c_a_workflow_res$source,
      bmd = as.numeric(dr_t_c_a_workflow_res$BMD.zSD),
      trend = dr_t_c_a_workflow_res$trend
    )
    
    ## Insert columns conditionally if they exist in the input data
    optional_columns <- c("gene_name", "description", "TF")
    for (col in optional_columns) {
      if (col %in% names(dr_t_c_a_workflow_res)) {
        dr_t_c_a_summary[[col]] <- dr_t_c_a_workflow_res[[col]]
      }
    }
    
    # Use relocate to move optional columns after 'gene_id'
    dr_t_c_a_summary <- dr_t_c_a_summary |> 
      dplyr::relocate(dplyr::any_of(optional_columns), .after = gene_id)
    
    # Round BMD.zSD to the tenth for easier reading
    dr_t_c_a_summary$bmd <- round(dr_t_c_a_summary$bmd, 2) 
    
    # Remove repeated rows
    dr_t_c_a_summary <- unique(dr_t_c_a_summary)
    
    # Fast save to csv file (using multiple CPUs)
    data.table::fwrite(x = dr_t_c_a_summary, 
                       file = file.path(path, output_filename), 
                       sep = ";", 
                       dec=  ",",
                       col.names = TRUE, 
                       na = NA) 
    
  }
  
}
