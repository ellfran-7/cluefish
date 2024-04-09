responsiv.data = dr_t_regs
responsiv.annot.genes = clustr_enrichres$geneannots
clustr.fusion.data = clustr_fusionres$generes
friendly.limit = 0
path = "data/derived-data"
output.filename = "lonely_fishres"




# If the output file does not exist or overwrite is set to TRUE, proceed with the fishing

# Step 1: Filtering lonely genes and small clusters
responsiv.data.no.clustr.genes <- responsiv.data[!(responsiv.data$ensembl_gene_id %in% clustr.fusion.data$ensembl_gene_id),] 

# Step 2: Filtering annotations based on enriched biological functions

## In the fishing process, we want to consider only terms and pathways that are enriched.
## I want to keep rows where at least one of the values in go_term, kegg_pathway, or wiki_pathway is present in clustr.fusion.data$term_name. If any of these columns has NA, you still want to keep the row if one of the other two columns has a valid value
responsiv.annot.genes.filtered <- responsiv.annot.genes |> 
  dplyr::filter(term_name %in% clustr.fusion.data$term_name)

# Step 3: Merge lonely genes with filtered biological function annotations
responsiv.data.no.clustr.genes <- merge(responsiv.data.no.clustr.genes, responsiv.annot.genes.filtered, by = "ensembl_gene_id")


# Select columns of interest in fishing
responsiv.data.no.clustr.genes <- responsiv.data.no.clustr.genes |> 
  dplyr::select(ensembl_gene_id, term_name)

# Remove repetitions
responsiv.data.no.clustr.genes <- unique(responsiv.data.no.clustr.genes)

## This is the dataframe of lonely genes (= rows) with an Ensembl gene id column and 
## the associated  term_name column column. These are the genes that are susceptible
## to be fished into clusters.

## Now expanding the clusters by fishing lonely genes sharing the 
## same driver GO, KEGG and/or WP annotations as a cluster !

modified.clustr.fusion.data <- clustr.fusion.data |> 
  dplyr::select(-ensembl_gene_id)

modified.clustr.fusion.data <- unique(modified.clustr.fusion.data)

lonely.merged.data <- merge(responsiv.data.no.clustr.genes, modified.clustr.fusion.data, by = "term_name")


lonely.fish.result <- rbind(clustr.fusion.data, lonely.merged.data)


# Original genes without clusters
no.clustr.data <- responsiv.data |> 
  dplyr::filter((!ensembl_gene_id %in% clustr.fusion.data$ensembl_gene_id) &
                  ensembl_gene_id %in% responsiv.annot.genes$ensembl_gene_id)

# Print information about lonely genes fished 
cat(length(unique(lonely.merged.data$ensembl_gene_id)), "/", length(unique(no.clustr.data$ensembl_gene_id)), "lonely annotated genes fished !", "\n")

cat(length(unique(no.clustr.data[!no.clustr.data$ensembl_gene_id %in% lonely.merged.data$ensembl_gene_id,]$ensembl_gene_id)), "remain lonely !")

# This dataframe contains all kept clusters with all genes types of genes inside : 
# - genes in original cluster but not annotated
# - genes in original cluster enriching a kept biological annotation
# - genes from a different original cluster but fusioned together 
# - lonely genes not originally in cluster but fishing based shared biological annotation

## But what are the left genes that haven't been associated with a cluster?

# Select the remaining lonely data
lonely.data <- responsiv.data[!(responsiv.data$ensembl_gene_id %in% lonely.fish.result$ensembl_gene_id), ]

# Prepare data for rbind
lonely.data.4rbind <- data.frame(ensembl_gene_id = lonely.data$ensembl_gene_id,
                                 old_clustr = "Lonely",
                                 new_clustr = "Lonely")

lonely.data.4rbind <- merge(lonely.data.4rbind, responsiv.annot.genes, by = "ensembl_gene_id")

exp.clustr.data <- rbind(lonely.fish.result, lonely.data.4rbind)

exp.clustr.data <- exp.clustr.data |> 
  dplyr::select(ensembl_gene_id, old_clustr, new_clustr, 
                term_name, term_id, source) |>
  dplyr::arrange(new_clustr, term_name)

exp.clustr.data <- unique(exp.clustr.data)


# We now have the expanded clusters dataframe, similar to 'lonely.fish.result'
# but with the added lonely genes associated to the "Lonely" cluster

### But what about gene that are part of many clusters ?

# When fishing lonely genes, we retrieve alot of information, especially when we are taking into account multiple annotations from multiple sources (GO, KEGG and WP). This means that we will retrieve a significant amount of genes scattered around maybe 3, 4 or more clusters. These genes can be named "Friendly" genesas they are part of a multitude of clusters, in the contrary to the "Lonely" genes. We can therefore associate these friendly genes to a "Friendly" cluster based a limit upon how many clusters a gene has to be part of to be considered friendly or not. But if the gene is a transcription (co-)factor, then should be kept in the cluster as they are interesting to evaluate

if (friendly.limit != 0) {
  
  exp.clustr.data <- exp.clustr.data |> 
    dplyr::group_by(ensembl_gene_id) |> 
    dplyr::mutate(friendliness = dplyr::n_distinct(new_clustr)) |> 
    dplyr::mutate(new_clustr = dplyr::if_else(friendliness > overfriendly.limit & TF == FALSE, "Friendly", as.character(new_clustr))) |> 
    dplyr::ungroup()

  } else {
  
  exp.clustr.data <- exp.clustr.data |> 
    dplyr::group_by(ensembl_gene_id) |> 
    dplyr::mutate(friendliness = dplyr::n_distinct(new_clustr)) |> 
    dplyr::ungroup()
}
  

# We nearly have all the info, but we are missing the items that are not annotated with the Ensembl identification, therefor they have not been fished. We need to add these genes to hte results : 

if (dim(responsiv.data[!(responsiv.data$ensembl_gene_id %in% exp.clustr.data$ensembl_gene_id),])[1] != 0) {
  
  lonely.data.no.ensembl <- responsiv.data[!(responsiv.data$ensembl_gene_id %in% exp.clustr.data$ensembl_gene_id),]
  
  lonely.data.no.ensembl.4rbind <- data.frame(ensembl_gene_id = lonely.data.no.ensembl$ensembl_gene_id,
                                              old_clustr = "Lonely",
                                              new_clustr = "Lonely",
                                              term_name = NA,
                                              term_id = NA,
                                              source = NA,
                                              friendliness = 1)
  
  exp.clustr.data <- rbind(exp.clustr.data, lonely.data.no.ensembl.4rbind)
  
  exp.clustr.data <- unique(exp.clustr.data) # remove repeated rows
  
}

# We can finally add all the remaining information : other gene IDs/transcript ID and TF status
exp.clustr.data <- merge(exp.clustr.data, responsiv.data, by = "ensembl_gene_id", all = TRUE)

exp.clustr.data <- exp.clustr.data |> 
  dplyr::select(ensembl_transcript_id_version, ensembl_gene_id, entrezgene_id, external_gene_name, TF, old_clustr, new_clustr, friendliness, everything()) |> 
  dplyr::arrange(new_clustr, source, term_name)

saveRDS(exp.clustr.data, file.path(path, output.filename))

return(exp.clustr.data)










# --------------------- FULL DPLYR VERSION (nearly)

# Step 1: Filtering lonely genes and small clusters
responsiv.data.no.clustr.genes <- responsiv.data |> 
  dplyr::anti_join(clustr.fusion.data, by = "ensembl_gene_id")

# Step 2: Filtering annotations based on enriched biological functions
responsiv.annot.genes.filtered <- responsiv.annot.genes |> 
  dplyr::semi_join(clustr.fusion.data, by = c("term_name" = "term_name"))

# Step 3: Merge lonely genes with filtered biological function annotations
responsiv.data.no.clustr.genes <- dplyr::inner_join(responsiv.data.no.clustr.genes, responsiv.annot.genes.filtered, by = "ensembl_gene_id") |> 
  dplyr::distinct()

# Expand clusters by fishing lonely genes sharing the same driver GO, KEGG, and/or WP annotations as a cluster
modified.clustr.fusion.data <- clustr.fusion.data |> 
  dplyr::select(-c(ensembl_gene_id, term_id, source)) |> 
  dplyr::distinct()
lonely.merged.data <- dplyr::inner_join(responsiv.data.no.clustr.genes, modified.clustr.fusion.data, by = "term_name") 
lonely.fish.result <- dplyr::bind_rows(clustr.fusion.data, lonely.merged.data)



# Retrieve the DRomics BMD computed transcript genes not associated to clusters
no.clustr.data <- responsiv.data |> 
  dplyr::filter((!ensembl_gene_id %in% clustr.fusion.data$ensembl_gene_id) &
                  ensembl_gene_id %in% responsiv.annot.genes$ensembl_gene_id)

# Print information about lonely genes fished 
cat(length(unique(lonely.merged.data$ensembl_gene_id)), "/", length(unique(no.clustr.data$ensembl_gene_id)), "lonely annotated genes fished !", "\n")

cat(length(unique(no.clustr.data[!no.clustr.data$ensembl_gene_id %in% lonely.merged.data$ensembl_gene_id,]$ensembl_gene_id)), "remain lonely !")

# This dataframe contains all kept clusters with all genes types of genes inside : 
# - genes in original cluster but not annotated
# - genes in original cluster enriching a kept biological annotation
# - genes from a different original cluster but fusioned together 
# - lonely genes not originally in cluster but fishing based shared biological annotation

## But what are the left genes that haven't been associated with a cluster?

# Select the remaining lonely data
lonely.data <- responsiv.data |> 
  dplyr::anti_join(lonely.fish.result, by = "ensembl_gene_id")
lonely.data.4rbind <- data.frame(ensembl_gene_id = lonely.data$ensembl_gene_id,
                                 old_clustr = "Lonely",
                                 new_clustr = "Lonely")

lonely.data.4rbind <- dplyr::inner_join(lonely.data.4rbind, responsiv.annot.genes, by = "ensembl_gene_id", relationship = "many-to-many")

exp.clustr.data <- dplyr::bind_rows(lonely.fish.result, lonely.data.4rbind)

exp.clustr.data <- exp.clustr.data |> 
  dplyr::select(ensembl_gene_id, old_clustr, new_clustr, 
                term_name, term_id, source) |>
  dplyr::arrange(new_clustr, term_name) |> 
  dplyr::distinct()

### But what about gene that are part of many clusters ?
# Handling "Friendly genes

if (friendly.limit != 0) {
  
  exp.clustr.data <- exp.clustr.data |> 
    dplyr::group_by(ensembl_gene_id) |> 
    dplyr::mutate(friendliness = dplyr::n_distinct(new_clustr)) |> 
    dplyr::mutate(new_clustr = dplyr::if_else(friendliness > overfriendly.limit & TF == FALSE, "Friendly", as.character(new_clustr))) |> 
    dplyr::ungroup()
  
} else {
  
  exp.clustr.data <- exp.clustr.data |> 
    dplyr::group_by(ensembl_gene_id) |> 
    dplyr::mutate(friendliness = dplyr::n_distinct(new_clustr)) |> 
    dplyr::ungroup()
}


# Include remaining lonely genes without Ensembl IDs

lonely.data.no.ensembl <- responsiv.data[!(responsiv.data$ensembl_gene_id %in% exp.clustr.data$ensembl_gene_id),]

lonely.data.no.ensembl <- dplyr::anti_join(responsiv.data, exp.clustr.data, by = "ensembl_gene_id")

lonely.data.no.ensembl.4rbind <- data.frame(ensembl_gene_id = lonely.data.no.ensembl$ensembl_gene_id,
                                            old_clustr = "Lonely",
                                            new_clustr = "Lonely",
                                            term_name = NA,
                                            term_id = NA,
                                            source = NA,
                                            friendliness = 1)

exp.clustr.data <- dplyr::bind_rows(exp.clustr.data, lonely.data.no.ensembl.4rbind) |> 
  dplyr::distinct()



# Merge with the original data
exp.clustr.data <- dplyr::full_join(exp.clustr.data, responsiv.data, by = "ensembl_gene_id") |>
  dplyr::select(ensembl_transcript_id_version, ensembl_gene_id, entrezgene_id, external_gene_name, TF, old_clustr, new_clustr, friendliness, everything()) |> 
  dplyr::arrange(new_clustr, source, term_name) |> 
  dplyr::distinct()
