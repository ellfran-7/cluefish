#' Generate a workflow summary table for exploration
#'
#' @description
#' This function creates a concise summary table capturing the key details from the workflow results from the `lonelyfishing()` function. It encompasses all essential information for result exploration, striking a balance by avoiding an overwhelming amount of data that might hinder ease of exploration.
#' 
#' @param lonelyfishing_data The named `file` output of the `lonelyfishing()` function.
#' @param bmdboot_data The DRomics bmdboot dataframe results after DRomics::bmdfilter() 
#' @param path Destination folder for the output data results.
#' @param output_filename Output CSV filename.
#' @param overwrite If `TRUE`, the function overwrites existing output files; otherwise, it reads the existing file. (default is set to `FALSE`).
#' 
#' @return A `.csv` file holding the results of the final step of the workflow with each row being a combination of Ensembl gene, cluster and biological function annotation. This is to be used as a support for exploration and mechanism discovery.
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
    dr_t_a_workflow_res <- merge(lonelyfishing_data$dr_t_c_a_fishing, bmdboot_data, by.x = "ensembl_transcript_id_version", by.y = "id")
    
    # Prepare the structure for the summary dataframe
    dr_t_a_summary <- data.frame(
      ensembl_transcript_id = dr_t_a_workflow_res$ensembl_transcript_id_version,
      ensembl_gene_id = dr_t_a_workflow_res$ensembl_gene_id,
      external_gene_name = dr_t_a_workflow_res$external_gene_name,
      NewCluster = dr_t_a_workflow_res$new_clustr,
      Friendliness = dr_t_a_workflow_res$friendliness,
      Term_name = dr_t_a_workflow_res$term_name,
      Source = dr_t_a_workflow_res$source,
      TF = dr_t_a_workflow_res$TF,
      BMD.zSD = as.numeric(dr_t_a_workflow_res$BMD.zSD),
      Trend = dr_t_a_workflow_res$trend
    )
    
    # Round BMD.zSD to the tenth for easier reading
    dr_t_a_summary$BMD.zSD <- round(dr_t_a_summary$BMD.zSD, 1) 
    
    # Remove repeated rows
    dr_t_a_summary <- unique(dr_t_a_summary)
    
    # Fast save to csv file (using multiple CPUs)
    data.table::fwrite(x = dr_t_a_summary, 
                       file = file.path(path, output_filename), 
                       sep = ";", 
                       dec=  ",",
                       col.names = TRUE, 
                       na = NA) 
    
  }
  
}
