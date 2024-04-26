
require(ggplot2)
test <- DRomics::bmdplot(b_fishres_no_redund[b_fishres_no_redund$new_clustr %in% 1, ], 
        add.CI = TRUE,
        colorby = "trend",
        shapeby = "TF",
        point.size = 3,
        point.alpha = 1,
        line.size = 1,
        line.alpha = 0.7,
        add.label = TRUE,
        BMD_log_transfo = TRUE) +
  xlab("BMD.zSD (μg/L)") + 
  scale_colour_manual(values = c("inc" = "#53e078", "dec" = "#fda547", 
                                "U" = "#9590c4", "bell" = "#fc67a2")) +
  theme_light()


# TEST
# As a gene can be associated to multiple transcripts, we want to be able to differentiate these gene symbols for visualisations.

test <- data.frame(transcript = c("a", "b", "c"),
                   gene = c("A", "B", "A"))
                   
gene_order <- sort(test$gene)

test <- test[order(gene_order),]

duplicated_genes <- test$gene[duplicated(test$gene)]

if (length(duplicated_genes) > 0) {
  for (i in 1:length(duplicated_genes)) {
    indices <- which(test$gene == duplicated_genes[i])
    num_duplicates <- sum(test$gene == duplicated_genes[i])
    new_names <- paste(duplicated_genes[i], seq_len(num_duplicates), sep = ".")
    test$gene[indices] <- new_names
  }
}

# Load DRomics bmdboot object
b <- readRDS(file = "data/raw-data/bootres_zebrafish_phtalate_UF_seed3_5000iter.rds")

# We filter the bmdboot result by selecting only transcripts with a defined confidence interval around the BMD
BMDres_definedCI <- DRomics::bmdfilter(b$res, BMDfilter = "definedCI")

fishres <- readRDS("outputs/lonely_fishres_2024_04_22.rds")

b_fishres <- merge(fishres$dr_t_c_a_fishing, BMDres_definedCI, by.x = "ensembl_transcript_id_version", by.y = "id")

# Some transcripts lack a gene symbol (a readable gene name), which complicates visualizations. In such cases, we assign a new name to these transcripts' genes as "Unknown" followed by a position number based on the number of blank names present. For instance, if there are 8 instances of unknown gene names, the 5th gene will be named "Unknown.5".

indices <- which(b_fishres$external_gene_name == "")
num_duplicates <- sum(b_fishres$external_gene_name == "")
new_names <- paste("Unknown", seq_len(num_duplicates), sep = ".")
b_fishres$external_gene_name[indices] <- new_names

# Since a gene can be associated with multiple transcripts, it's essential to distinguish these gene symbols for visualizations. To achieve this, we append a unique numerical value to the gene names corresponding to their duplicate positions. For example, if two transcripts are associated with the gene "A", the new names for each transcript's gene will be "A.1" and "A.2". The numbering does not follow a specific order.

gene_order <- sort(b_fishres$external_gene_name)

b_fishres <- b_fishres[order(gene_order),]

duplicated_genes <- b_fishres$external_gene_name[duplicated(b_fishres$external_gene_name)]

duplicated_genes <- duplicated_genes[duplicated_genes != ""]

unknown_count = 0

if (length(duplicated_genes) > 0) {
  for (i in 1:length(duplicated_genes)) {
    indices <- which(b_fishres$external_gene_name == duplicated_genes[i])
    num_duplicates <- sum(b_fishres$external_gene_name == duplicated_genes[i])
    new_names <- paste(duplicated_genes[i], seq_len(num_duplicates), sep = ".")
    b_fishres$external_gene_name[indices] <- new_names
  }
}


# Connection to the ENSEMVL database and the species dataset 
ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl") 

# Retrieves the specified attributes from the BioMart database given a set of filters and corresponding query values
query_results <- biomaRt::getBM(attributes = c("ensembl_transcript_id_version",
                                  "ensembl_gene_id", 
                                  "external_gene_name"), 
                                filters = "ensembl_transcript_id_version",
                                values = f$omicdata$item,
                                mart = ensembl)

# Assign new names to transcripts lacking gene symbols. In such cases, we assign a new name to these transcripts' genes as "unknown" followed by a position number based on the number of blank names present. For instance, if there are 8 instances of unknown gene names, the 5th gene will be named "Unknown.5".

indices <- which(query_results$external_gene_name == "")
num_duplicates <- sum(query_results$external_gene_name == "")
new_names <- paste("unknown", seq_len(num_duplicates), sep = ".")
query_results$external_gene_name[indices] <- new_names

# Differentiate gene symbols for multiple transcripts associated with the same gene. To achieve this
# We append a unique numerical value to gene symbol names based on duplicate positions in Ensembl gene IDs and Ensembl transcript IDs transcript ids. If two transcripts are linked to the same gene ID and symbol "A", their new names will be "A_1.1" and "A_1.2". If they are linked to different gene IDs but the same symbol "A", their new names will be "A_1.1" and "A_2.1" The numbering does not follow a specific order other than alphabetical..

# Append unique numerical values to gene symbols associated with duplicate transcripts
gene_order <- sort(query_results$external_gene_name)

query_results <- query_results[order(gene_order),]

duplicated_genes <- query_results$external_gene_name[duplicated(query_results$external_gene_name)]

duplicated_genes <- duplicated_genes[duplicated_genes != ""]

unknown_count = 0

if (length(duplicated_genes) > 0) {
  for (i in 1:length(duplicated_genes)) {
    indices <- which(query_results$external_gene_name == duplicated_genes[i])
    num_duplicates <- sum(query_results$external_gene_name == duplicated_genes[i])
    new_names <- paste(duplicated_genes[i], seq_len(num_duplicates), sep = ".")
    query_results$external_gene_name[indices] <- new_names
  }
}


test1 <- data.frame(gene = c("A", "B", "A", "B"),
                    transcript = c("a", "b", "c", "d"),
                    symbol = c("GeneX_1.1", "Gene_2.1", "Gene_1.2", "Gene_2.2"))

test2 <- data.frame(gene = c("A", "B", "A", "B"),
                    transcript = c("a", "b", "c", "d"),
                    symbol = c("GeneX", "GeneY", "GeneX", "GeneY"))

# Some transcripts lack a gene symbol (a readable gene name), which complicates visualizations. In such cases, we assign a new name to these transcripts' genes as "Unknown" followed by a position number based on the number of blank names present. For instance, if there are 8 instances of unknown gene names, the 5th gene will be named "Unknown.5".

indices <- which(b_lonely_fishres$external_gene_name == "")
num_duplicates <- sum(b_lonely_fishres$external_gene_name == "")
new_names <- paste("Unknown", seq_len(num_duplicates), sep = ".")
b_lonely_fishres$external_gene_name[indices] <- new_names

# Since a gene can be associated with multiple transcripts, it's essential to distinguish these gene symbols for visualizations. To achieve this, we append a unique numerical value to the gene names corresponding to their duplicate positions. For example, if two transcripts are associated with the gene "A", the new names for each transcript's gene will be "A.1" and "A.2". The numbering does not follow a specific order.

symbol_order <- sort(test2$symbol)
gene_order <- sort(test2$gene)
test2 <- test2[order(symbol_order, gene_order),]

test2_gene_level <- test2 |> 
  dplyr::select(gene, symbol)

duplicated_genes <- unique(test2_gene_level$symbol[duplicated(test2_gene_level$symbol)])


if (length(duplicated_genes) > 0)
  
indices <- which(test2$symbol == unique(duplicated_genes[1]))
num_duplicates <- sum(test2 == duplicated_genes[1])
new_names <- paste(duplicated_genes[1], seq_len(num_duplicates), sep = "_")
test2$symbol[indices] <- new_names


test2_transcript_level <- test2 |> 
  dplyr::select(transcript, symbol)

duplicated_genes <- unique(test2_transcript_level$symbol[duplicated(test2_transcript_level$symbol)])


if (length(duplicated_genes) > 0)
  
indices <- which(test2$symbol == unique(duplicated_genes[1]))
num_duplicates <- sum(test2 == duplicated_genes[1])
new_names <- paste(duplicated_genes[1], seq_len(num_duplicates), sep = "_")
test2$symbol[indices] <- new_names











# Sample dataframe
df <- data.frame(
  gene = c("A", "A", "A", "B"),
  transcript = c("a", "b", "c", "d"),
  symbol = c("GeneX", "GeneX", "GeneX", "GeneY")
)

# Function to add numbering
add_numbering <- function(df) {
  # Get unique gene identifiers in the order they appear
  unique_genes <- unique(df$gene)
  
  # Map each gene to a unique numeric identifier based on order of appearance
  gene_id_map <- as.numeric(factor(df$gene, levels = unique_genes))
  
  # Create a unique identifier for each transcript within each gene
  transcript_counts <- with(df, ave(transcript, gene, FUN = seq_along))
  transcript_duplicates <- with(df, ave(transcript, gene, FUN = seq_along))
  
  # Add numbering to the symbol column
  df$symbol <- paste(df$symbol, gene_id_map, transcript_duplicates, sep = "_")
  
  return(df)
}

# Apply the function
df <- add_numbering(df)

# Output the modified dataframe
print(df)








































# Connection to the ENSEMVL database and the species dataset 
ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl") 

# Retrieves the specified attributes from the BioMart database given a set of filters and corresponding query values
query_results <- biomaRt::getBM(attributes = c("ensembl_transcript_id_version",
                                               "ensembl_gene_id", 
                                               "external_gene_name"), 
                                filters = "ensembl_transcript_id_version",
                                values = f$omicdata$item,
                                mart = ensembl)

## Assign new names to transcripts lacking gene symbols. In such cases, we assign a new name to these transcripts' genes as "unknown" followed by a number based on the order of appearance of these blank names. For instance, if there are 8 occurances of unknown gene names, the 5th gene symbol will be named "unknown.5".

# Assign "Unknown" names to genes lacking a symbol and number them sequentially
indices <- which(query_results$external_gene_name == "")
num_duplicates <- sum(query_results$external_gene_name == "")
new_names <- paste("unknown", seq_len(num_duplicates), sep = ".")
query_results$external_gene_name[indices] <- new_names

## Differentiate gene symbols for multiple transcripts associated with the same gene. To achieve this
## We append a unique numerical value to gene symbol names based on duplicate positions in Ensembl gene IDs and Ensembl transcript IDs transcript ids. If two transcripts are linked to the same gene ID and symbol "A", their new names will be "A_1.1" and "A_1.2". If they are linked to different gene IDs but the same symbol "A", their new names will be "A_1.1" and "A_2.1" The numbering does not follow a specific order other than alphabetical..

# Order the gene symbols and retrieve the duplicates
gene_order <- sort(query_results$external_gene_name)

query_results <- query_results[order(gene_order),]

duplicated_genes <- query_results$external_gene_name[duplicated(query_results$external_gene_name)]

# Check if there are duplicated gene symbols in the query results
if (length(duplicated_genes) > 0) {
  # Iterate over each duplicated gene symbol
  for (i in 1:length(duplicated_genes)) {
    
    # Determine the indices of the rows associated with the duplicated gene symbol
    indices <- which(query_results$external_gene_name == duplicated_genes[i])
    
    # Extract unique Ensembl gene identifiers in the order they appear
    unique_genes <- unique(query_results$ensembl_gene_id[indices])
    
    # Map each Ensembl gene to a unique numeric identifier based on the order of appearance
    gene_id_map <- as.numeric(factor(query_results$ensembl_gene_id[indices], levels = unique_genes))
    
    # Count the appearances of each Ensembl transcript within each Ensembl gene
    transcript_counts <- with(query_results, ave(ensembl_transcript_id_version[indices], ensembl_gene_id[indices], FUN = seq_along))
    
    # Add ensembl transcript and gene appearance numbering to the gene symbol column 
    query_results$external_gene_name[indices] <- paste0(query_results$external_gene_name[indices], "_", gene_id_map, ".", transcript_counts)
    
  }
}

