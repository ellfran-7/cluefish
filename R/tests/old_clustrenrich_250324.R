#' Cluster-specific functional enrichment
#' 
#' @description
#' Firstly, this function performs Over-Representation Analysis (ORA) on each cluster. 
#' An option is added to filter pathways that are enriched by "too few" genes. Secondly, it adds the data for genes in a string cluster but not participating in the enrichment. Lastly, it retrieves up-to-date gene sets for GO, KEGG, or WikiPathways (WP) from within the g:profiler database. 
#'
#' @param clustr.data The $generes dataframe from output of the 'clustrfiltr()' function with at least  'ensembl_gene_id' and 'clustr' columns needed for GO, KEGG, and WP based ORA.
#' @param responsiv.genes The dataframe output of the 'getregs()' function with at least the 'ensembl_gene_id' column for ORA.
#' @param background The background Ensembl gene list (preferably from the experiment) found in the output of DRomics::drcfit().
#' @param background.type Type of background to be used: if non-given = choice between "whole" and "annotated", 
#' if given = choice between "custom" and "custom_annotated".
#' @param databases A vector of data sources to use. Currently, these include GO (GO:BP, GO:MF, GO:CC to select a 
#' particular GO branch), KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP
#' @param org Organism ID defined for the chosen databases: if zebrafish = "Danio rerio"
#' @param pcutoff Adjusted p-value cutoff for Over-Representation analysis (default at 0.05 in gost() function)
#' @param padjmethod P-value adjustment method: one of “gSCS” ,“fdr” and “bonferroni
#' @param exclude.iea Exclude GO electronic annotations (IEA) (default set as "FALSE")
#' @param enrich.size.filtr Option to filter out enriched terms that have too little (min.term.size) and/or too large (max.term.size) gene sets. Recommended for well-annotated organisms (default = TRUE).
#' @param min.term.size Minimum gene set size to consider a biological pathway as relevant to the analysis 
#' @param max.term.size Maximum gene set size to consider a biological pathway as relevant to the analysis
#' @param only.highlighted.GO Option to only keep highlighted driver GO terms in the analysis results (default = TRUE)
#' @param n.genes.filter Minimum number of genes to consider sufficient for an enriched pathway to be kept in the output
#' @param path Destination folder for the output data
#' @param output.filename Output enrichment result filename
#' @param overwrite If TRUE, the file will be downloaded again, and the previous version will be replaced.
#'
#' @return A named list file holding 4 components : 
#' - "generes": a dataframe of the results with each row being a gene ("ensembl_gene_id") associated to a specific biological function ("term_name")
#' - "gostres": a named list where 'result' contains the data frame with enrichment analysis results, and 'meta' contains metadata necessary for creating a Manhattan plot. This is the original output of a gprofiler2::gost()
#' - "geneannots": a dataframe holding all the annotations found from the g:profiler database for all the deregulated genes
#' - "simplifylog": a dataframe holding the number of biological functions enriched per cluster, before and after each filtering step for each source.
#' 
#' @export
#'
#' @examples

clustrenrich <- function(
    clustr.data,
    responsiv.genes,
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
    
    # Select the two columns of interest : query genes ("ensembl_gene_id") and cluster ids ("clustr")
    modif.cluster.data <- clustr.data |> 
      dplyr::select(ensembl_gene_id, clustr) |> 
      dplyr::group_by(clustr) 
    
    # Set the “clustr” column values to factors for appropriate list names
    modif.cluster.data$clustr <- factor(modif.cluster.data$clustr, 
                                        levels = unique(modif.cluster.data$clustr))
    
    # Turn the “ensembl_gene_ids” column values to characters
    modif.cluster.data$ensembl_gene_id <- as.character(modif.cluster.data$ensembl_gene_id)
    
    # Create the list by dividing the gene set data (“ensembl_gene_id”) into the groups (“clustr”). This is reassembled in the form of a list of vectors containing the values for the groups. It will be the input of the gprofiler2::gost() function and this format allows us to perform a gost() on each cluster
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
      user_threshold = pcutoff, 
      correction_method = padjmethod, 
      domain_scope = background.type,
      custom_bg = background, 
      numeric_ns = "", 
      sources = databases, 
      as_short_link = FALSE, 
      highlight = TRUE 
    ) 
    
    # Format the results in order for more convenient exploration
    multi.gostres$result <- multi.gostres$result |> 
      dplyr::mutate(query = as.numeric(query)) |> 
      dplyr::arrange(query)
    
    
    # ---------------
    # using the "multi.gostres$result" we create one of the components the output of clustrenrich(), which is the $simplifylog. It holds the number of biological functions enriched per cluster, before and after each filtering step for each source. In this first block, we create the base dataframe with no filters performed yet.
    
    # Create a dataframe holding the number of terms enriched by clusters before and after each condition. This will be updated after each sub-step of this function.
    all.enriched.terms <- multi.gostres$result |> 
      dplyr::select(query, term_name, source)
    
    # Group the data by query and source, count the number of distinct term names per combination of query and source, turn the source column categories into columns associated to their proper count by cluster, rename the columns for coherence and turn the tibble into a dataframe.
    term.count.per.clustr <- all.enriched.terms |>
      dplyr::group_by(query, source) |>
      dplyr::summarize(count = dplyr::n_distinct(term_name), .groups = "drop") |>
      tidyr::pivot_wider(names_from = source, values_from = count, values_fill = 0) |> 
      dplyr::rename(clustr = query, all.GO = 'GO:BP', all.KEGG = KEGG, all.WP = WP) |> 
      as.data.frame()
    # ---------------
    
    
    # Condition whether to remove non-highlighted GO terms from further analysis
    if (only.highlighted.GO == TRUE) {
      
      multi.gostres$result <- multi.gostres$result |> 
        dplyr::filter(
          ((grepl("GO", source) & highlighted == TRUE) | (!(grepl("GO", source))))
        )
      
      # ----------------
      # When we filter out non-highlighted GO terms, we want to know how many are left, and add it to our "term.count.per.clustr" dataframe
      
      # select only the three columns needed from the results
      high.enriched.go <- multi.gostres$result |> 
        dplyr::select(query, term_name, source)
      
      # Select only terms associated to the GO:BP source, remove the source column, group the data by the query column, summarise the groups with a new column "driver.GO" thats counts the number of unique occurences of term names within each query group, rename the query column to "clustr" and turn the tibble into a dataframe.
      high.go.count.per.clustr <- high.enriched.go |>
        dplyr::filter(source == "GO:BP") |> 
        dplyr::select(-source) |> 
        dplyr::group_by(query) |> 
        dplyr::summarize(driver.GO = dplyr::n_distinct(term_name), .groups = "drop") |> 
        dplyr::rename(clustr = query) |> 
        as.data.frame()
      
      # Crush the previous "term.count.per.clustr" with the newly merged dataframe also holding driver.GO.term occurence counting
      term.count.per.clustr <- merge(term.count.per.clustr, high.go.count.per.clustr, by = "clustr", all = TRUE)
      # ----------------
      
    }
    
    # Depending on sizes chosen, remove too little and too large biological function gene sets from further analysis  
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
    
    # Condition whether to remove biological functions that are not insufficiently enriched
    if (enrich.size.filtr == TRUE) {
      
      # Remove non-highlighted GO terms and keep the biological functions that are sufficiently enriched by enough genes
      clustr.term.count <- multi.gostres.by.gene  |> 
        dplyr::group_by(clustr, term_name) |> 
        dplyr::summarise(n_genes = length(unique(ensembl_gene_id)), .groups = "drop") |> 
        dplyr::ungroup()  
      
      # Remove biological functions that are enriched by less than 3 genes
      clustr.term.filtr <- subset(clustr.term.count, n_genes >= n.genes.filter)
      
      # Combine this with the gprofiler2::gost() results
      clustr.kept.term <- merge(multi.gostres.by.gene, clustr.term.filtr, by = c("clustr", "term_name"))
      
      # Remove the "n_genes" column
      clustr.kept.term$n_genes <- NULL
      
      # Remove repeated rows
      clustr.kept.term <- unique(clustr.kept.term)
      
      # Order and keep the selected columns
      clustr.kept.term <- clustr.kept.term |> 
        dplyr::select(ensembl_gene_id, clustr, term_name, term_id, source)
      
      # Print the number of enriching and non-enriching clusters
      cat(length(unique(clustr.kept.term$clustr)), "/",  length(unique(clustr.data$clustr)), "clusters participating in enrichment.", "\n")
      
      # Print the number of terms removed because of filter
      cat(length(unique(clustr.kept.term$term_name)), "/",  length(unique(multi.gostres.by.gene$term_name)), "enriched terms kept after term filter.", "\n")
      
      # Transform the clustr column to numeric to order the data by clustr
      clustr.kept.term$clustr <- as.numeric(clustr.kept.term$clustr)
      clustr.kept.term <- clustr.kept.term[order(clustr.kept.term$clustr), ]
      
      multi.gostres.by.gene <- clustr.kept.term
      
    } else { 
      
      # Transform the clustr column to numeric to order the data by clustr
      multi.gostres.by.gene$clustr <- as.numeric(multi.gostres.by.gene$clustr)
      multi.gostres.by.gene <- multi.gostres.by.gene[order(multi.gostres.by.gene$clustr), ]
      
    }
    
    # --------------
    # After filtering the gene set sizes and the enrichment sizes of biological functions, we also want to get the number of occurences for each cluster and source combination.
    
    # Select only the columns of interest
    after.filtr.enriched.terms <- multi.gostres.by.gene |> 
      dplyr::select(clustr, term_name, source)
    
    # Here we group the data by clustr and source, count the number of distinct term names per combination of clustr and source, turn the source column categories into columns associated to their proper count by cluster, rename the columns for coherence, remove the columns source and term name and turn the tibble into a dataframe.
    kept.term.count.per.clustr <- after.filtr.enriched.terms |>
      dplyr::group_by(clustr, source) |>
      dplyr::summarize(count = dplyr::n_distinct(term_name), .groups = "drop") |>
      tidyr::pivot_wider(names_from = source, values_from = count, values_fill = 0) |> 
      dplyr::rename(kept.driver.GO = 'GO:BP', kept.KEGG = KEGG, kept.WP = WP) |> 
      as.data.frame()
    
    # Here we merge crush the previous term.count.per.clustr with the newly merged dataframe also holding driver.GO.term occurence counting
    term.count.per.clustr <- merge(term.count.per.clustr, kept.term.count.per.clustr, by = 'clustr', all = TRUE)
    
    # Rearranging the rows in source then step within the function order
    term.count.per.clustr <- term.count.per.clustr |> 
      dplyr::select(clustr, all.GO, driver.GO, kept.driver.GO, all.KEGG, kept.KEGG, all.WP, kept.WP)
    
    # Replace all NAs with 0 (NAs appear because of the filtering steps, removing clusters therefor they were originally enriched terms but now they don't)
    term.count.per.clustr <- term.count.per.clustr |> 
      dplyr::mutate_all(~replace(., is.na(.), 0))
    # --------------
    
    
    # As we have performed enrichment, we only have the genes and data participating in the enrichment of terms. Therefor all genes that are in a string cluster but not participating in any enrichment are not in the direct results. They need to be put back into the results.
    
    # Retrieve the data for all ensembl gene ids not participating in any enrichment
    no.enrich.data <- clustr.data[!(clustr.data$ensembl_gene_id %in% multi.gostres.by.gene$ensembl_gene_id),]
    
    # Construct a df for the rbind with the gost enrichment output
    no.enrich.data.4.rbind <- data.frame(ensembl_gene_id = no.enrich.data$ensembl_gene_id,
                                         clustr = no.enrich.data$clustr,
                                         term_name = NA,
                                         term_id = NA,
                                         source = NA)
    
    # Combine both results : enrichment and non-enrichment genes
    full.data <- rbind(multi.gostres.by.gene, no.enrich.data.4.rbind)
    
    # Remove repeated rows
    full.data <- unique(full.data)
    
    # Order the data by clustr
    full.data <- full.data[order(full.data$clustr),]
    
    
    # STEP 7 bis - Use the gprofiler2::gost() function to retrieve annotations for the deregulated genes
    
    # Retrieve GO:BP, KEGG and WP annotations using the gprofiler::gost() function :
    annot.gostres <- gprofiler2::gost(
      query = responsiv.genes, 
      organism = org, 
      ordered_query = FALSE, 
      multi_query = FALSE, 
      significant = FALSE, 
      exclude_iea = FALSE, 
      measure_underrepresentation = FALSE, 
      evcodes = TRUE,
      sources = databases, 
      as_short_link = FALSE, 
      highlight = FALSE 
    ) 
    
    # Transform the dataframe from "cluster per row" to "gene per row" :
    annot.genes <- annot.gostres$result |> 
      tidyr::separate_rows(intersection, sep = ",") |> 
      dplyr::select(intersection, term_name, term_id, source) |> 
      dplyr::rename(ensembl_gene_id = intersection)
    
    # Remove repeated rows
    annot.genes <- unique(annot.genes)
    
    # Now create the resulting list combining the enrichment results and the annotations
    clustr.enrichres <- list(dr_g_a_enrich = full.data,
                             gostres = multi.gostres,
                             dr_g_a_whole = as.data.frame(annot.genes),
                             c_simplifylog = term.count.per.clustr)
    
    # Name the structure of the result
    clustr.enrichres <- structure(clustr.enrichres, class = "clustrenrichres")
    
    # Save the object
    saveRDS(clustr.enrichres, file.path(path, output.filename))
    
    return(clustr.enrichres)
    
  }
  
}
