#' Cluster-specific functional enrichment
#' 
#' @description
#' This function performs Over-Representation Analysis (ORA) on each cluster. 
#' An option is added to filter pathways that are enriched by "too few" genes.
#'
#' @param clustr.data The output dataframe from 'clustr.size.filtr' function (with 'ensembl_gene_id' and 'clustr' 
#' columns needed for GO, KEGG, and WP based ORA)
#' @param background The background Ensembl gene list (preferably from the experiment) (output of `get.ids` function)
#' @param background.type Type of background to be used: if non-given = choice between "whole" and "annotated", 
#' if given = choice between "custom" and "custom_annotated".
#' @param databases A vector of data sources to use. Currently, these include GO (GO:BP, GO:MF, GO:CC to select a 
#' particular GO branch), KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP
#' @param org Organism ID defined for the chosen databases: if zebrafish = "Danio rerio"
#' @param pcutoff Adjusted p-value cutoff for Over-Representation analysis (default at 0.05 in gost() function)
#' @param padjmethod P-value adjustment method: one of “gSCS” ,“fdr” and “bonferroni
#' @param exclude.iea Exclude GO electronic annotations (IEA) (default set as "FALSE")
#' @param enrich.size.filtr Option to filter out enriched terms that have too little (min.term.size) and/or too large (max.term.size) gene sets. Recommended for well-annotated organisms.
#' @param min.term.size Minimum gene set size to consider a biological pathway as relevant to the analysis 
#' @param max.term.size Maximum gene set size to consider a biological pathway as relevant to the analysis
#' @param only.highlighted.GO Option to only keep highlighted driver GO terms in the analysis results
#' @param n.genes.filter Minimum number of genes to consider insufficient for an enriched pathway to be kept in the output
#' @param path Destination folder for the output data
#' @param output.filename Output enrichment result filename
#' @param overwrite If TRUE, the file will be downloaded again, and the previous version will be replaced.
#'
#' @return A `.txt` file of the functional enrichment results, similar to the input (`clustr.data`) structure.
#' 
#' @export
#'
#' @examples

clustrenrich <- function(
    clustr.data,
    background, 
    background.type = "custom_annotated", 
    databases, 
    org, 
    pcutoff = 0.05,  
    padjmethod = "fdr",
    exclude.iea = FALSE, 
    enrich.size.filtr = TRUE, 
    min.term.size = 5,
    max.term.size = 500, 
    only.highlighted.GO = TRUE,
    n.genes.filter = 2, 
    path, 
    output.filename, 
    overwrite = FALSE 
)

{
  
  # Check if the output file already exists locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, output.filename)) && !overwrite) {
    
    message("The enrichment results file already exists. Reading this file. Use 'overwrite = TRUE' to replace it.")
    
    # Read the existing file and return it 
    read.table(file.path(path, output.filename))
    
  } else {
    
    # Select the two columns of interest : query genes and cluster ids
    # Then order this by cluster number (1, 2, 3 ... n)
    modif.cluster.data <- clustr.data |> 
      dplyr::select(ensembl_gene_id, clustr) |> 
      dplyr::group_by(clustr) 
    
    # Set the “clustr” column values to factors for appropriate list names
    modif.cluster.data$clustr <- factor(modif.cluster.data$clustr, 
                                        levels = unique(modif.cluster.data$clustr))
    
    # Turn the “ensembl_gene_ids” column values to characters
    modif.cluster.data$ensembl_gene_id <- as.character(modif.cluster.data$ensembl_gene_id)
    
    # Create the list by dividing the gene set data (“ensembl_gene_id”) into the groups (“clustr”). This is reassembled in the form of a list of vectors containing the values for the groups.
    cluster.list <- split(modif.cluster.data$ensembl_gene_id, 
                          modif.cluster.data$clustr)
    
    
    # Perform Over-Representation Analysis (ORA) using gprofiler2::gost
    multi.gostres <- gprofiler2::gost(
      query = cluster.list, 
      organism = org, 
      ordered_query = FALSE, 
      multi_query = FALSE, 
      significant = TRUE, 
      exclude_iea = exclude.iea, 
      measure_underrepresentation = FALSE, 
      evcodes = TRUE, 
      user_threshold = 0.05, 
      correction_method = padjmethod, 
      domain_scope = background.type,
      custom_bg = background, 
      numeric_ns = "", 
      sources = databases, 
      as_short_link = FALSE, 
      highlight = TRUE 
    ) 
    
    
    # Format the results in order for more conveniant exploration
    multi.gostres$result <- multi.gostres$result |> 
      dplyr::arrange(query)
    
    # Remove non-highlighted GO terms from further analysis
    if (only.highlighted.GO == TRUE) {
      
      multi.gostres$result <- multi.gostres$result |> 
        dplyr::filter(
          ((grepl("GO", source) & highlighted == TRUE) | (!(grepl("GO", source))))
        )
    }
    
    # Remove too little and too large biological function gene sets from further analysis  
    multi.gostres$result <- multi.gostres$result |> 
      dplyr::filter(
        min.term.size < term_size & term_size < max.term.size
      )
    
    # Transform the dataframe from "cluster per row" to "gene per row" :
    multi.gostres.by.gene <- multi.gostres$result |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::select(intersection, query, term_name, term_id, source) |> 
      dplyr::rename(ensembl_gene_id = intersection, 
                    clustr = query)
    
    if (enrich.size.filtr == TRUE) {
      
      # Remove non-highlighted GO terms and keep the biological functions that are sufficiently enriched by enough genes
      clustr.term.count <- multi.gostres.by.gene  |> 
        dplyr::group_by(clustr, term_name) |> 
        dplyr::summarise(n_genes = length(unique(ensembl_gene_id)), .groups = "drop") |> 
        dplyr::ungroup()  
      
      # Remove biological functions that are enriched by less than 3 genes
      clustr.term.filtr <- subset(clustr.term.count, n_genes > n.genes.filter)
      
      clustr.kept.term <- merge(multi.gostres.by.gene, clustr.term.filtr, by = c("clustr", "term_name"))
      
      clustr.kept.term$n_genes <- NULL
      clustr.kept.term <- unique(clustr.kept.term)
      clustr.kept.term <- clustr.kept.term[, c(3, 1, 2, 4, 5)]
      
      # Print the number of enriching and non-enriching clusters
      cat(length(unique(clustr.kept.term$clustr)), "/",  length(unique(clustr.data$clustr)), "clusters participating in enrichment.", "\n")
      
      # Print the number of terms removed because of filter
      cat(length(unique(clustr.kept.term$term_name)), "/",  length(unique(multi.gostres.by.gene$term_name)), "enriched terms kept after term filter.", "\n")
      
      clustr.kept.term$clustr <- as.numeric(clustr.kept.term$clustr)
      clustr.kept.term <- clustr.kept.term[order(clustr.kept.term$clustr), ]
      
      
      write.table(clustr.kept.term, file.path(path, output.filename))
      return(clustr.kept.term)
      
    } else { 
      
      multi.gostres.by.gene$clustr <- as.numeric(multi.gostres.by.gene$clustr)
      multi.gostres.by.gene <- multi.gostres.by.gene[order(multi.gostres.by.gene$clustr), ]
      
      write.table(multi.gostres.by.gene, file.path(path, output.filename))
      return(multi.gostres.by.gene)
      
    }
    
  }
  
}
