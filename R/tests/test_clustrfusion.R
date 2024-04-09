clustrfusion2 <- function(
    clustr.ora.filtr.data
)
{
  
  # Associate input to different variable to avoid conflict when manipulation in the for loop
  input.data <- clustr.ora.filtr.data
  
  # Select the term_name column (for the merging process) and the source column (for fusion information) then order the data based on the alphabetical order of the source column. We want to order this column with GO, KEGG then WP. This is reasoned in the README
  term.info <- clustr.ora.filtr.data |> 
    dplyr::select(term_name, source) |> 
    dplyr::arrange(source)
  
  
  # Retrieve the unique terms in alphabetical order
  term.ids <- unique(term.info)
  
  
  # Initialize an empty list to store the information about fused clusters (clusters: annotation)
  fused_clusters_info <- list()
  
  # For loop for merging
  for (i in 1:length(term.ids$term_name)){
    
    # The data associated to the unique term i
    unique.term.data <- input.data[input.data$term_name == term.ids$term_name[i],]
    
    # Remove repetitions
    unique.clustrs <- unique(unique.term.data$clustr)
    
    # If there are multiple clusters associated the the term i
    if (length(unique.clustrs) > 1){
      
      # Merge the clusters by associating the smallest cluster id to the entire column of the subsetted data
      input.data[input.data$clustr %in% unique.clustrs,]$clustr <- min(unique.clustrs)
      
      # Store the information about fused clusters
      fused_clusters_info[[length(fused_clusters_info) + 1]] <- paste(paste(unique.clustrs, collapse = ","), ": ", term.ids$term_name[i], "--", term.ids$source[i])
      
    } else {
      next
    }
  }
  
  # Count the number of clusters before the fusion
  total_clusters_before_fusion <- length(unique(clustr.ora.filtr.data$clustr))
  
  # Count the number of cluster after fusion
  total_clusters_after_fusion <- length(unique(input.data$clustr))
  
  # Calculate the number of clusters that participated in fusion
  clusters_participated_in_fusion <- total_clusters_before_fusion - total_clusters_after_fusion
  
  # Print the ratio of clusters after/before the process so that the user gets an idea of the impact of the process
  cat(total_clusters_after_fusion, "/", total_clusters_before_fusion, "clusters left after the fusion process. \n")
  
  # Remove repeated rows if any
  input.data <- unique(input.data)
  
  # Create the clustr.fusion object
  clustr.fusionres <- list(res = input.data,
                           fusions = fused_clusters_info)
  
  clustr.fusionres <- structure(clustr.fusionres, format = "clustrfusion")
  
  return(clustr.fusionres)
}
