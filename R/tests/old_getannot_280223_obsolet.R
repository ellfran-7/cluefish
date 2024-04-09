#' Retrieve biological and functional annotations
#'
#' @description
#' This function retrieves up-to-date gene sets for GO, KEGG, or WikiPathways (WP) database. 
#' The WP file needs to be in the same folder as the saved output of this function.
#' 
#' @param query.ids A table with at least one column named `entrezgene_id` holding  
#' Entrez identifiers for the sleected genes. 
#' The input should be the output of the `get_ids` function, with item rows and identifier columns.
#' @param database The database from which annotations are retreived ("GO", "KEGG" or "WP").
#' @param org.go The organism dataset namespace for biomaRt.
#' @param select.go.category is an option to select one or multiple of the 3 GO term categories ("biological_process", "molecular function" and "cellular_compartiment")
#' @param go.category When `select.go.category == TRUE`, choose the GO term category(ies) to be kept.
#' @param org.kegg The organism ID reference in the KEGG annotation database (e.g., "dre" for Danio Rerio)
#' @param org.wp The organism ID reference in the WP annotation database (e.g., "Danio rerio")
#' @param path The destination folder for wp.data and the output data.
#' @param wp.filename The GMT species file name downloaded from the WP database 
#' (https://data.wikipathways.org/current/gmt/) (only use if `database = "WP"`).
#' @param output.filename The output data filename
#' @param overwrite A `logical` value. If `TRUE`, the file will be downloaded again and the previous version will be replaced.
#'
#' @return A `.txt` file with GO, KEGG and/or Wikipathway annotations for the queried genes.
#' 
#' @export
#'
#' @examples
#' 

getannot <- function(
    query.ids, 
    database, 
    org.go, 
    select.go.category = FALSE, 
    go.category,
    org.kegg, 
    org.wp,
    path, 
    wp.filename, 
    output.filename, 
    overwrite = FALSE 
)
{
  # Check if the output file already exists locally and overwrite is not set to TRUE
  if (file.exists(file.path(path, output.filename)) && !overwrite) {
    
    message("The responsive annotation file already exists. Reading this file. Use 'overwrite = TRUE' to replace it.")
    
    # Read the existing file and return it
    read.table(file.path(path, output.filename))
    
  } else {
    
    # Create a "annot.results" dataframe to store annotation results. As we are interesting in gene annotations relative to the experiment, we will use the background gene id list as the mapping dataframe.
    annot.results <- query.ids
    
    # Retrieve data from the GO database if specified
    if ("GO" %in% database){ 
      
      ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = org.go) 
      
      
      # Retrieve up-to-date GO annotations gene sets using "biomaRt"
      golist <- biomaRt::getBM(mart = ensembl,
                               attributes = c("name_1006", 
                                              "namespace_1003", 
                                              "entrezgene_id"),
                               filters = "ensembl_transcript_id_version",
                               values = query.ids$ensembl_transcript_id_version)
      
      # Remove non-annotated genes from the background
      golist <- golist[golist$name_1006 != "", ] 
      
      # Select the GO terms to be kept, based on the chosen GO category selected
      if (select.go.category == TRUE) { 
        golist <- golist |> 
          dplyr::filter(namespace_1003 %in% go.category) 
      }
      
      # Rename columns accordingly
      colnames(golist) <- c("go_term", 
                            "go_category",
                            "entrezgene_id") 
      
      annot.results <- merge(annot.results, golist, by = "entrezgene_id", all.x = TRUE)
      
    }
    
    # Retrieve data from the KEGG database if specified
    if ("KEGG" %in% database){
      
      pathways <- KEGGREST::keggList("pathway", org.kegg)
      
      # Order the list alphabetically
      pathways <- pathways[sort(names(pathways))] 
      
      # Retrieve KEGG pathway to gene links
      pathways2gene <- KEGGREST::keggLink(org.kegg, "pathway")
      
      # Remove "path:" prefix
      names(pathways2gene) <- sub("^path:", "", names(pathways2gene)) 
      
      # Function to return all genes associated with the corresponding KEGG pathways while removing the org.kegg prefix
      .getGenes <- function(pathway)
      { 
        genes <- pathways2gene[names(pathways2gene) == pathway]
        genes <- sub(paste0("^", org.kegg, ":"), "", genes)
        genes <- unname(sort(genes))
        return(genes)
      }
      
      # Apply the function to all the names in the pathways list
      gs <- lapply(names(pathways), .getGenes) 
      
      # Build first gmt column: the ID (format: <pathway.nr>_<pathway.title>)
      .makeGSNames <- function(ids, titles)
      {
        ids <- sub("^path:", "", ids)
        titles <- vapply(titles, 
                         function(title) unlist(strsplit(title, " - "))[1], 
                         character(1))
        titles <- sub("^ +", "", titles)
        titles <- sub(" +$", "", titles)
        titles <- gsub(" ", "_", titles)
        ids <- paste(ids, titles, sep="_")
        return(ids)
      }
      
      titles <- vapply(pathways, 
                       function(title) unlist(strsplit(title, " - "))[1], 
                       character(1))
      
      # Rename the gene set list with the new titles
      names(gs) <- titles
      
      # Turn the list of gene sets into a dataframe with col1: associated KEGG paths and col2: all KEGG annotated gene IDs for the organism. Rows represent the genes.
      kegglist <- data.frame(
        pathway = rep(names(gs), sapply(gs, length)),
        gene = unlist(gs),
        stringsAsFactors = FALSE
      ) 
      
      # Rename the columns accordingly
      colnames(kegglist) <- c("kegg_pathway", "entrezgene_id") 
      
      annot.results <- merge(annot.results, kegglist, by = "entrezgene_id", all.x = TRUE)
      
    } 
    
    # Retrieve data from the WP database if specified
    if ("WP" %in% database) { 
      # After downloading the GMT file from the WikiPathways db for the specific organism model The 'readPathwayGMT' function in the 'rWikiPathways' package is used to convert the GMT file to a dataframe of pathway-gene associations
      wikidata <- rWikiPathways::readPathwayGMT(file.path(path, wp.filename))
      wikilist <- data.frame(wiki_pathway = wikidata$name, 
                             entrezgene_id = wikidata$gene)
      
      annot.results <- merge(annot.results, wikilist, by = "entrezgene_id", all.x = TRUE)
      
    }
    
    # Write the annotation results to a file
    write.table(annot.results, file.path(path, output.filename))
    
    return(annot.results)
    
  }
  
}




# ========================================================================





# STEP 3 - BIOLOGICAL FUNCTION ANNOTATION RETRIEVAL 
# -----------------------------------------------------------------------------
#
# This script retrieves biological function annotations using the 'getannot' function.
#
# The function used in the script has been developed for this project
# and can be found in the folder R/, under the name '2-getannot.R'.

# We only need the annotations for the selected transcripts/genes from the DRomics workflow :
dr_t_ids <- bg_t_ids[bg_t_ids$ensembl_transcript_id_version %in% b$res$id,]

responsiv_bio_annotres <-  getannot(
  query.ids = dr_t_ids,
  database = c("GO", "KEGG", "WP"),
  org.go = "drerio_gene_ensembl",
  select.go.category = TRUE,
  go.category = "biological_process",
  org.kegg = "dre",
  path = "data/derived-data/",
  wp.filename = "wikipathways-20240210-gmt-Danio_rerio.gmt",
  output.filename = "responsiv_bio_annotres.txt",
  overwrite = TRUE
)



