#' Generate a workflow summary table for exploration
#'
#' @description
#' This function creates a concise summary table capturing the key details from the workflow results from the `lonelyfishing()` function. It encompasses all essential information for result exploration, striking a balance by avoiding an overwhelming amount of data that might hinder ease of exploration.
#' 
#' @param lonelyfishing_data The named `list` output of the `lonelyfishing()` function.
#' @param bmdboot_data The DRomics bmdboot `dataframe` results after DRomics::bmdfilter() 
#' @param readable_gene_id Name of the readable gene ID column to add to the summary table as "gene_name" (e.g. "external_gene_name". If none is provided, no additional column will be created (default is set to "NULL")
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
    readable_gene_id = NULL,
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
    dr_t_c_a_workflow_res <- merge(lonelyfishing_data$dr_t_c_a_fishing, bmdboot_data, by.x = "ensembl_transcript_id_version", by.y = "id")
    
    
    # Prepare the structure for the summary dataframe
    dr_t_c_a_summary <- data.frame(
      ensembl_transcript_id = dr_t_c_a_workflow_res$ensembl_transcript_id_version,
      ensembl_gene_id = dr_t_c_a_workflow_res$ensembl_gene_id,
      gene_name = dr_t_c_a_workflow_res$readable_gene_id,
      NewCluster = dr_t_c_a_workflow_res$new_clustr,
      Friendliness = dr_t_c_a_workflow_res$friendliness,
      Term_name = dr_t_c_a_workflow_res$term_name,
      Source = dr_t_c_a_workflow_res$source,
      BMD.zSD = as.numeric(dr_t_c_a_workflow_res$BMD.zSD),
      Trend = dr_t_c_a_workflow_res$trend
    )
    
    # Add readable gene ID column if provided and exists in the data
    if (!is.null(readable_gene_id) && readable_gene_id %in% names(dr_t_c_a_workflow_res)) {
      
      dr_t_c_a_summary$gene_name <- dr_t_c_a_workflow_res[[readable_gene_id]]
      
      dr_t_c_a_summary <- dr_t_c_a_summary[, c("ensembl_transcript_id", "ensembl_gene_id", "gene_name", 
                                               "NewCluster", "Friendliness", "Term_name", "Source", 
                                               "BMD.zSD", "Trend")]
    }
    
    
    # If TF column exists in input data, add it to the summary dataframe
    if ("TF" %in% names(dr_t_c_a_workflow_res)) {
      
      dr_t_c_a_summary$TF <- dr_t_c_a_workflow_res[[tf_column]]
      
    }
    
    # Round BMD.zSD to the tenth for easier reading
    dr_t_c_a_summary$BMD.zSD <- round(dr_t_c_a_summary$BMD.zSD, 1) 
    
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
