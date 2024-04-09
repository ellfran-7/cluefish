#' Fishing lonely genes into clusters based mapping of annotation to enrichment
#'
#' @description
#' This function identifies and integrates lonely genes (genes not found in any cluster) into existing clusters. The integration is based on shared gene pathway annotations and cluster pathway enrichments. The lonely genes are identified from the input gene data (`responsiv.data`) and subsequently associated with clusters based on shared biological annotations from pathways.
#' 
#' @param responsiv.data A data frame with at least two columns: `ensembl_gene_id` and `term_name`, containing Entrez gene identifiers and pathway annotations for selected genes. This data frame is typically the output of the `get_reg_annot` function.
#' @param responsiv.annot.genes A data frame representing the output of the 'annot.genes' function, encompassing gene annotations from various databases (GO, KEGG, WP).
#' @param filtr.clustr.data A data frame with at least two columns: `ensembl_gene_id` and `clustr`, containing Entrez gene identifiers and cluster identifiers for selected genes. Usually, the output of the `clustr_size_filtr` function.
#' @param clustr.ora.filtr.data A data frame with at least columns named `ensembl_gene_id` and `go_term`, `kegg_pathway` and/or `wiki_pathway`. The input should be the output of the `clustr_enrich` function.
#' @param clustr.fusion.data A data frame with at least two columns: `ensembl_gene_id` and `clustr`, containing Ensembl gene identifiers and cluster identifiers for selected genes. Typically, this data frame is the output of the `clustr_fusion` function.
#' @param overfriendly.limit The maximum number of clusters a gene can be part of to be considered "Friendly". Genes exceeding this limit are assigned to a separate "Friendly" cluster.
#' @param path The destination folder for storing intermediate data and the output file.
#' @param output.filename The filename for the output file containing the expanded clusters.
#' @param overwrite If `TRUE`, the function overwrites existing output files; otherwise, it reads the existing file. Default is `FALSE`.
#' 
#' @return A `.txt` file containing the resulting expanded clusters with detailed information for each gene, including Ensembl transcript ID, Ensembl gene ID, Entrez gene ID, Gene Symbol, Cluster ID, Biological function annotations, and (co)regulator information.
#' 
#' @export
#'
#' @examples
#' 

lonelyfishing <- function(
    responsiv.data,
    responsiv.annot.genes, 
    clustr.fusion.data, 
    overfriendly.limit = 0, 
    path, 
    output.filename,
    overwrite = FALSE 
)
  
{
  
  # Check if the output file already exists locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, output.filename)) && !overwrite) {
    
    message("The lonely fishing results file already exists. Reading this file. Use 'overwrite = TRUE' to replace it.")
    
    # Read the existing file and return it 
    read.table(file.path(path, output.filename))
    
  } else {
    
    # If the output file does not exist or overwrite is set to TRUE, proceed with the fishing
    
    # Step 1: Filtering lonely genes and small clusters
    responsiv.data.no.clustr.genes <- responsiv.data[!(responsiv.data$ensembl_gene_id %in% clustr.fusion.data$ensembl_gene_id),] 
    
    
    # Step 2: Filtering annotations based on enriched biological functions
    
    ## In the fishing process, we want to consider only terms and pathways that are enriched.
    ## I want to keep rows where at least one of the values in go_term, kegg_pathway, or wiki_pathway is present in clustr.fusion.data$term_name. If any of these columns has NA, you still want to keep the row if one of the other two columns has a valid value
    
    responsiv.annot.genes.filtered <- responsiv.annot.genes |> 
      dplyr::filter(
        (
          (
            is.na(go_term) & kegg_pathway %in% responsiv.annot.genes$kegg_pathway & wiki_pathway %in% responsiv.annot.genes$wiki_pathway | 
              is.na(kegg_pathway) & go_term %in% clustr.fusion.data$term_name & wiki_pathway %in% responsiv.annot.genes$wiki_pathway |
              is.na(wiki_pathway) & kegg_pathway %in% responsiv.annot.genes$kegg_pathway & go_term %in% clustr.fusion.data$term_name
          ) |
            (
              is.na(go_term) & is.na(kegg_pathway) & wiki_pathway %in% responsiv.annot.genes$wiki_pathway | 
                is.na(kegg_pathway) & is.na(wiki_pathway) & go_term %in% clustr.fusion.data$term_name |
                is.na(wiki_pathway) & is.na(go_term) & kegg_pathway %in% responsiv.annot.genes$kegg_pathway
            ) |
            (
              go_term %in% clustr.fusion.data$term_name & kegg_pathway %in% responsiv.annot.genes$kegg_pathway &  wiki_pathway %in% responsiv.annot.genes$wiki_pathway
            )
        )
      )
    
    # Remove unecessary columns
    responsiv.annot.genes.filtered <- responsiv.annot.genes.filtered |> 
      dplyr::select(-c(entrezgene_id,
                       ensembl_transcript_id_version,
                       external_gene_name))
    
    # Step 3: Merge lonely genes with filtered biological function annotations
    
    responsiv.data.no.clustr.genes <- merge(responsiv.data.no.clustr.genes, responsiv.annot.genes.filtered, by = "ensembl_gene_id")
    
    # Pivot the data to long format, filter and arrange
    responsiv.data.no.clustr.genes <- responsiv.data.no.clustr.genes |>
      tidyr::pivot_longer(cols = c(go_term, kegg_pathway, wiki_pathway), names_to = "term_type", values_to = "term_name") |>
      dplyr::filter(!is.na(term_name) & term_name %in% clustr.fusion.data$term_name) |>
      dplyr::arrange(ensembl_gene_id)
    
    # Create a dataframe with relevant columns
    responsiv.data.no.clustr.genes <- data.frame(entrezgene_id = responsiv.data.no.clustr.genes$entrezgene_id,
                                                 ensembl_gene_id = responsiv.data.no.clustr.genes$ensembl_gene_id,
                                                 term_name = responsiv.data.no.clustr.genes$term_name)
    
    responsiv.data.no.clustr.genes <- unique(responsiv.data.no.clustr.genes)
    
    ## This is the dataframe of lonely genes (= rows) with two gene ID columns and 
    ## the associated GO/KEGG/WP column. These are the genes that are susceptible
    ## to be fished into clusters.
    
    ## Now expanding the clusters by fishing lonely genes sharing the 
    ## same driver GO, KEGG and/or WP annotations as a cluster !
    
    # Initialize empty vectors
    n_genes <- c() # to store genes sharing the cluster term annotation
    clustr_id <- c() # to store the expanded cluster(s) to which the genes are fished 
    
    # Progress bar for fishing process
    pb <- progress::progress_bar$new(
      format = "  fishing [:bar] :percent elapsed: :elapsed",
      clear = FALSE,
      total = length(clustr.fusion.data$term_name)
    )
    
    # Loop through each biological annotation
    for (i in 1:length(clustr.fusion.data$term_name))
      
    {
      # Loop through each gene within the biological annotation
      for (j in 1:dim(responsiv.data.no.clustr.genes)[1])
        
      {
        # Check if gene shares a biological annotation with a cluster
        if (responsiv.data.no.clustr.genes$term_name[j] == clustr.fusion.data$term_name[i])
          
        {
          clustr_id <- append(clustr_id, clustr.fusion.data$clustr[i]) 
          n_genes <- append(n_genes, responsiv.data.no.clustr.genes[j,]$ensembl_gene_id) 
        }
      }
      
      # Update progress bar
      pb$tick()
      
    }
    
    # Finish progress bar
    pb$terminate()
    
    # Create the dataframe of fished genes for rbind 
    added.fish.genes.4rbind <- data.frame(clustr = clustr_id,
                                          ensembl_gene_id = n_genes)
    
    # Removed repeated rows
    added.fish.genes.4rbind <- unique(added.fish.genes.4rbind)
    
    
    # Original genes without clusters
    no.clustr.data <- responsiv.data |> 
      dplyr::filter((!ensembl_gene_id %in% clustr.fusion.data$ensembl_gene_id) &
                      ensembl_gene_id %in% responsiv.annot.genes$ensembl_gene_id)
    
    # Print information about lonely genes fished 
    cat(length(unique(added.fish.genes.4rbind$ensembl_gene_id)), "/", length(unique(no.clustr.data$ensembl_gene_id)), "lonely annotated genes fished !", "\n")
    
    cat(length(unique(no.clustr.data[!no.clustr.data$ensembl_gene_id %in% added.fish.genes.4rbind$ensembl_gene_id,]$ensembl_gene_id)), "remain lonely !")
    
    ## This is a dataframe of clusters after size filtering, biological annotation filtering, cluster fusion and lonely gene fishing. This only contains genes that are annotated, meaning that we need to retrieve the non-annotated genes from the original string clusters. 
    
    ## First we need to fuse the original string clusters based on the same criteria as we created the 'clustr.fusion' function output. This means that we proceed with cluster fusion of cluster with shared single enriching terms, all the while merging origin string clusters based on the results. So here we apply the same procedure for fusion.
    
    
    # Group by cluster and count the number of unique terms for each cluster
    term.count <- clustr.ora.filtr.data |>
      dplyr::group_by(clustr, source) |>
      dplyr::summarize(n_term = dplyr::n_distinct(term_name), .groups = 'drop')
    
    # Get the cluster names that only have one unique biological annotation
    single_clusters <- term.count |>
      dplyr::filter(n_term == 1) |>
      dplyr::pull(clustr)
    
    # Get unique terms from single clusters
    single_terms <- unique(clustr.ora.filtr.data[clustr.ora.filtr.data$clustr %in% single_clusters,]$term_name)
    
    # Initialize vector to store clusters to be merged
    single_terms_clust <- c()
    
    # Loop through single terms
    for (i in single_terms) {
      
      # Loop through each source type (e.g., GO:BP, KEGG)
      for (source_type in unique(clustr.ora.filtr.data$source)) {
        
        # Check if clusters with shared annotation should be merged
        if (length(unique(clustr.ora.filtr.data |> 
                          dplyr::filter(term_name == i & clustr %in% single_clusters & source == source_type) |>
                          dplyr::pull(clustr))) > 1) {
          
          # Append clusters to be merged
          single_terms_clust <- append(single_terms_clust, as.character(unique(clustr.ora.filtr.data |>
                                                                                 dplyr::filter(term_name == i & clustr %in% single_clusters & source == source_type) |>
                                                                                 dplyr::pull(clustr))))
        }
      }
    }
    
    # Create a temporary data frame to store merged enriching clusters
    temp.clustr.ora.filtr.data <- clustr.ora.filtr.data
    
    # Create a temporary data frame to store merged original string clusters
    temp.filtr.clustr.data <- filtr.clustr.data 
    
    
    # Initialize vector to track fused clusters
    fused_clusters <- c()
    
    # Loop through the single terms and merge corresponding clusters for each source
    for (term in single_terms) {
      for (source_type in unique(temp.clustr.ora.filtr.data$source)) {
        clusters_to_merge <- temp.clustr.ora.filtr.data |> 
          dplyr::filter(term_name == term & clustr %in% single_terms_clust & source == source_type) |>
          dplyr::pull(clustr)
        
        # Check if any of the clusters have already been fused
        if (any(clusters_to_merge %in% fused_clusters)) {
          next  # Skip to the next source_type
        }
        
        # Check if clusters should be merged
        if (length(unique(clusters_to_merge)) > 1) 
        {  
          
          # Find the minimum numeric cluster
          new_cluster <- min(clusters_to_merge)
          
          # Extract clusters sharing a single biological annotatio
          sharing.clusters <- temp.clustr.ora.filtr.data |>
            dplyr::filter(clustr %in% clusters_to_merge)
          
          # Replace cluster numbers with the new fused cluster number
          sharing.clusters$clustr <- new_cluster
          
          # Concatenate information from the enriching clusters being merged
          temp.clustr.ora.filtr.data <- dplyr::bind_rows(temp.clustr.ora.filtr.data, sharing.clusters)
          
          # Concatenate information from the original string clusters being merged
          temp.filtr.clustr.data <- dplyr::bind_rows(temp.filtr.clustr.data, sharing.clusters)
          
          # Keep track of fused clusters
          fused_clusters <- union(fused_clusters, clusters_to_merge)
          
          # Identify other clusters with shared biological annotations
          clusters_with_shared_annotation <- temp.clustr.ora.filtr.data |>
            dplyr::filter(term_name %in% sharing.clusters$term_name & term_name %in% single_terms) |>
            dplyr::filter(!clustr %in% clusters_to_merge) |> 
            dplyr::pull(clustr)
          
          if (length(unique(clusters_with_shared_annotation)) > 0) {
            
            # Allow clusters with shared annotation to participate in subsequent merges
            fused_clusters <- setdiff(fused_clusters, clusters_to_merge)
          }
          
          # Keep the vector of cluster ids that have been merged to another
          clusters_to_merge <- setdiff(clusters_to_merge, min(clusters_to_merge))
          
          # Remove clusters from temporary enrichres df res that have merged and are not the final cluster id
          temp.clustr.ora.filtr.data <- temp.clustr.ora.filtr.data |>
            dplyr::filter(!clustr %in% clusters_to_merge)
          
          temp.clustr.ora.filtr.data <- unique(temp.clustr.ora.filtr.data)
          
          # Remove the same clusters from the temporary string cluster df that have merged and are not the final cluster id
          temp.filtr.clustr.data <- temp.filtr.clustr.data |>
            dplyr::filter(!clustr %in% clusters_to_merge)
          
          temp.filtr.clustr.data <- unique(temp.filtr.clustr.data)
          
          
        } else {
          next
        }
      }
    }
    
    # Filter and select relevant columns
    temp.filtr.clustr.data <- temp.filtr.clustr.data[temp.filtr.clustr.data$clustr %in% clustr.fusion.data$clustr, ] 
    
    # Create data frame for rbind
    temp.filtr.clustr.4rbind <- data.frame(ensembl_gene_id = temp.filtr.clustr.data$ensembl_gene_id,
                                           clustr = temp.filtr.clustr.data$clustr)
    
    ## These are the clusters from STRING after size filter and cluster fusion (fusion based on clusters sharing single biological enrichment per source). This dataframe contains only the original gene clusters without the annotation steps etc.This allows us to complete the final expanded clusters.
    
    # We also need to cross the info from the fusion step based, as there is info that is lost when only fusing the naked string clusters 
    
    # Cross information from the fusion step
    temp.clustr.fusion.4.rbind <- data.frame(ensembl_gene_id = temp.clustr.ora.filtr.data$ensembl_gene_id,
                                             clustr = temp.clustr.ora.filtr.data$clustr)
    
    bigger.clustr.data <- rbind(added.fish.genes.4rbind, temp.filtr.clustr.4rbind, temp.clustr.fusion.4.rbind) 
    
    # Drop duplicates
    bigger.clustr.data <- unique(bigger.clustr.data)
    
    # This dataframe contains all kept clusters with all genes types of genes inside : 
    # - genes in original cluster but not annotated
    # - genes in original cluster enriching a kept biological annotation
    # - genes from a different original cluster but fusioned together 
    # - lonely genes not originally in cluster but fishing based shared biological annotation
    
    ## But what are the left genes that haven't been associated with a cluster?
    
    # Select the remaining lonely data
    lonely.data <- responsiv.data[!(responsiv.data$ensembl_gene_id %in% bigger.clustr.data$ensembl_gene_id), ]
    
    # Prepare data for rbind
    lonely.data.4rbind <- data.frame(ensembl_gene_id = lonely.data$ensembl_gene_id,
                                     clustr = "Lonely")
    
    exp.clustr.data <- rbind(bigger.clustr.data, lonely.data.4rbind)
    
    # Select only relevant columns for merge, to gain TF information
    responsiv.annot.genes.filtered <- responsiv.annot.genes.filtered |> 
      dplyr::select(ensembl_gene_id,
                    go_term,
                    kegg_pathway,
                    wiki_pathway)
    
    responsiv.annot.genes.filtered <- merge(responsiv.data, responsiv.annot.genes.filtered, by = "ensembl_gene_id", all = TRUE)
    
    # Select only relevant columns for another merge
    responsiv.annot.genes.filtered <- responsiv.annot.genes.filtered |> 
      dplyr::select(ensembl_gene_id,
                    go_term,
                    kegg_pathway,
                    wiki_pathway,
                    TF)
    
    exp.clustr.data <- merge(exp.clustr.data, responsiv.annot.genes.filtered, by = "ensembl_gene_id", all = TRUE)
    
    # We now have the expanded clusters dataframe, similar to 'bigger.clustr.data'
    # but with the added lonely genes associated to the "Lonely" cluster
    
    ### But what about gene that are part of many clusters ?
    
    # When fishing lonely genes, we retrieve alot of information, especially when we are taking into account multiple annotations from multiple sources (GO, KEGG and WP). This means that we will retrieve a significant amount of genes scattered around maybe 3, 4 or more clusters. These genes can be named "Friendly" genesas they are part of a multitude of clusters, in the contrary to the "Lonely" genes. We can therefore associate these friendly genes to a "Friendly" cluster based a limit upon how many clusters a gene has to be part of to be considered friendly or not. But if the gene is a transcription (co-)factor, then should be kept in the cluster as they are interesting to evaluate
    
    exp.clustr.data <- exp.clustr.data |> 
      dplyr::group_by(ensembl_gene_id) |> 
      dplyr::mutate(clustr_count = dplyr::n_distinct(clustr)) |> 
      dplyr::mutate(clustr = dplyr::if_else(clustr_count > overfriendly.limit & TF == FALSE, "Friendly", as.character(clustr)))
    
    
    
    # We now need to retrieve the dose-response modeling info from
    # the 'dromics.workflow' ! 
    exp.clustr.data <- merge(exp.clustr.data[!is.na(exp.clustr.data$ensembl_gene_id), ],
                             responsiv.data[!is.na(responsiv.data$ensembl_gene_id), ], by = c("ensembl_gene_id", "TF"))
    # merging the output of 'get.ids' with the expanded cluster dataframe for 
    # adding the different gene IDs and the (co)regulator information but making
    # sure that no item with ensembl_gene_id = "NA" merges (!is.na() function)
    
    colnames(exp.clustr.data)[colnames(exp.clustr.data) == "ensembl_transcript_id_version"] <- "id"
    # modifying the column id so that DRomics plots recognize the column correctly
    # (the plots in this package take "id" as the column for plotting)
    
    # We nearly have all the info, but we are missing the items that are not annotated with the Ensembl identification, therefor they have not been fished. We need to add these genes to hte results : 
    
    if (dim(responsiv.data[!(responsiv.data$ensembl_transcript_id_version %in% exp.clustr.data$id),])[1] != 0) {
      
      lonely.data.no.ensembl <- responsiv.data[!(responsiv.data$ensembl_transcript_id_version %in% exp.clustr.data$id),]
      
      lonely.data.no.ensembl.4rbind <- data.frame(entrezgene_id = NA,
                                                  ensembl_gene_id = lonely.data.no.ensembl$ensembl_gene_id,
                                                  clustr = "Lonely",
                                                  go_term = NA,
                                                  kegg_pathway = NA,
                                                  wiki_pathway = NA,
                                                  clustr_count = 1,
                                                  id = lonely.data.no.ensembl$ensembl_transcript_id_version,
                                                  external_gene_name = lonely.data.no.ensembl$external_gene_name,
                                                  TF = lonely.data.no.ensembl$TF)
      
      exp.clustr.data <- rbind(exp.clustr.data, lonely.data.no.ensembl.4rbind)
      
      exp.clustr.data <- unique(exp.clustr.data) # remove repeated rows
      
    }
    
    exp.clustr.data <- exp.clustr.data |> 
      dplyr::select(id, ensembl_gene_id, entrezgene_id, external_gene_name, TF, clustr, clustr_count, everything()) |> 
      dplyr::arrange(clustr)
    
    write.table(exp.clustr.data, file.path(path, output.filename))
    
    return(exp.clustr.data)
    
    
  }
}