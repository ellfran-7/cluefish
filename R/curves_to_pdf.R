#' Generating curvesplots per-cluster to PDF
#' 
#' @description
#' This function generates a PDF file containing a plot of dose-response curves for each cluster of deregulated genes, with each plot labeled with the cluster ID and the number of transcripts in that cluster. The curves are color-coded according to whether the trend is increasing, decreasing, U-shaped, or bell-shaped. The plot axes are labeled with "Dose (µg/L)" and "Signal", and the Y-axis is scaled to be the same across all plots.
#' 
#' @param lonelyfishing_data The named `list` output of the `lonelyfishing()` function.
#' @param bmdboot_data The DRomics bmdboot dataframe results after DRomics::bmdfilter() 
#' @param clustrfusion_data The named `list` output of the `clustrfusion()` function. 
#' @param id_col_for_curves The column giving the identification of each curve (default is "id")
#' @param tested_doses A vector of the tested doses that can be found in the output of the `DRomics::drcfit()`function (unique(f$omicdata$dose))
#' @param annot_order A vector specifying the prioritized order of annotation sources used in the `clustrenrich()` function. This order determines which term or pathway will be used as primary descriptor for each cluster. The relevance of each source can be based on the quantity or quality of information they provide, tailored to the specific case study. For example, for Danio rerio, you might prioritize sources like c("GO:BP", "KEGG", "WP") based on the abundance of information they offer.
#' @param ... Additional arguments passed to the `DRomics::curvesplot()` function
#' @param xunit Unit for the x scale
#' @param xtitle X-axis title 
#' @param ytitle Y-axis title
#' @param colors A vector of colors for different trends (default is set at c("inc" = "#1B9E77", "dec" = "#D95F02", "U" = "#7570B3", "bell" = "#E7298A"))
#' @param path Destination folder for the output data results.
#' @param output_filename Output PDF filename.
#' @param overwrite If `TRUE`, the function overwrites existing output files; otherwise, it reads the existing file. (default is set to `FALSE`).
#'
#' @return No return value; the `.pdf` file is downloaded and stored in the specified director. It holds a curvesplot per slide for each cluster, and those for the Friendly (if created), Lonely and All deregulated gene transcripts.
#' @export
#'
#' @examples
#' 

curves_to_pdf <- function(
    lonelyfishing_data, 
    bmdboot_data,
    clustrfusion_data,
    id_col_for_curves = "id",
    tested_doses,
    annot_order,
    ...,
    xunit,
    xtitle,
    ytitle,
    colors = c("inc" = "#1B9E77", "dec" = "#D95F02", "U" = "#7570B3", "bell" = "#E7298A"),
    path,
    output_filename,
    overwrite = TRUE
    )
{
  
  # Check if the directory exists, and create it if it does not
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    # `recursive = TRUE` creates intermediate directories as needed
    
    # Message stating that the directory is created
    message("Directory '", path, "' created.")
  }
  
  
  # Check if the output file already exists locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, output_filename)) && !overwrite) {
    
    # Message stating that the results file already exists and will be read
    message("The curvesplot pdf file already exists. Use 'overwrite = TRUE' to replace it.")
    
  } else {
    
    # Extract the lonely fishing results 
    dr_t_c_a_fishing4curvesplot <- lonely_fishres$dr_t_c_a_fishing
    
    # Turn the "ensembl_transcript_id" version to "id" for compatibility with the DRomics bmdboot results and for the curvesplot
    names(dr_t_c_a_fishing4curvesplot)[names(dr_t_c_a_fishing4curvesplot) == "transcript_id"] <- "id"
    
    # Merge loenlyfishing results with bmdboots results after DRomics::bmdfilter()
    dr_t_c_a_workflow_res <- merge(dr_t_c_a_fishing4curvesplot, bmdboot_data, by = "id")
    
    ## Subset columns from merged data conditionally if they exist in the input data, removing non-essential columns that can cause redundancy for the dose-response curves
    subset_columns <- c("transcript_id", "gene_id", "gene_name", "description", "TF", 
                        "old_clustr", "friendliness", "term_name", "term_id", "source")
    
    for (col in subset_columns) {
      
      if (col %in% names(dr_t_c_a_workflow_res)) {
        
        dr_t_c_a_workflow_res <- dr_t_c_a_workflow_res[, !names(dr_t_c_a_workflow_res) %in% col, drop = FALSE]
      }
      
    }
    
    # Remove duplicate rows to pass from "t_c_a" to "t".
    dr_t_workflow_res <- unique(dr_t_c_a_workflow_res)
    
    # Add the new_clustr column to the bmdboot results and assign "All" to every row.
    # Assign "All" to represent the entire group of deregulated gene transcripts
    bmdboot_data$new_clustr <- "All" 
    
    # Create a dataframe with all clusters, including the "All" cluster with all deregulated gene transcripts
    dr_t_data4pdf <- rbind(dr_t_workflow_res, bmdboot_data)
    
    # Vector of cluster IDs in increasing order 
    cluster_names <- sort(unique(as.numeric(dr_t_data4pdf$new_clustr)))
    
    # Check if "Friendly" cluster exists in the data, if not, only include "Lonely" and "All" clusters
    if ("Friendly" %in% lonelyfishing_data$dr_g_a_fishing$new_clustr) {
      
      # Adding the "Friendly", "Lonely" and "All" clusters to the end of the cluster ID vector
      cluster_names <- append(cluster_names, c("Friendly", "Lonely", "All"))
      
    } else {
      
      # Adding the "Lonely" and "All" clusters to the end of the cluster ID vector
      cluster_names <- append(cluster_names, c("Lonely", "All"))
      
    }
    
    # Arrange the data by cluster and by source order : GO:BP, KEGG and then WP. GO:BP terms take precedence and are pasted first in the curvesplot title. If there are no enriched GO:BP terms, a KEGG term is pasted. If neither GO:BP nor KEGG terms are present, a WP term is pasted.
    dr_g_a_fusion_ordered <- clustrfusion_data$dr_g_a_fusion |> 
      dplyr::arrange(new_clustr, match(source, annot_order))
    
    # Calculate the first quartile value for each cluster
    groupby <- as.factor(dr_t_data4pdf[, "new_clustr"])
    variable <- dr_t_data4pdf[, "BMD.zSD"]
    firstquartilefun <- function(x) stats::quantile(x, probs = 0.25, na.rm = TRUE)
    dnb <- as.data.frame(table(groupby))
    colnames(dnb) <- c("groupby", "nb_of_items")
    dnb$firstquartile <- tapply(variable, groupby, firstquartilefun)
    dnb$firstquartile <- lapply(dnb$firstquartile, round, 2)
    
    # Avoid scientific notation in x-axis
    options(scipen = 999) 
    
    # Create a blank PDF document to which the function will print each cluster's curvesplot 
    pdf(file = file.path(path, output_filename), width = 5, height = 3, onefile = TRUE)
    
    # Loop through each cluster to create a curvesplot and print to the PDF file
    for (i in 1:length(cluster_names))
    {
      df_clustr <- dr_t_data4pdf[dr_t_data4pdf$new_clustr %in% cluster_names[i],]
      cplot <- DRomics::curvesplot(df_clustr, ...) +
        geom_vline(xintercept = tested_doses,  linetype = 2,  colour = "#302533", 
                   alpha = 0.7,  linewidth = 0.3) +
        labs(title = paste("Cluster", cluster_names[i], ":", dr_g_a_fusion_ordered[dr_g_a_fusion_ordered$new_clustr %in% cluster_names[i],]$term_name, "-", dr_g_a_fusion_ordered[dr_g_a_fusion_ordered$new_clustr %in% cluster_names[i],]$source), 
             subtitle = paste(length(unique(df_clustr$id)), "items", "with 25% <", dnb[dnb$groupby %in% cluster_names[i],]$firstquartile, xunit)) +
        xlab(xtitle) + 
        ylab(ytitle) +
        scale_color_manual(values = colors) +
        theme_light()
      
      print(cplot) 
      
    }
    
    dev.off() # Close the finalized PDF file 
    
  }
  
}
