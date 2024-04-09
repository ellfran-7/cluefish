# SUMMARY TABLES


b <- readRDS("data/derived-data/zebra_boots_05")
lonely_fishres <- read.table("data/derived-data/lonely_fishres.txt")


b_fishres <- merge(lonely_fishres, b$res, by = "id")

pergene.res.summary <- data.frame(
  ensembl_transcript_id = b_fishres$id,
  ensembl_gene_id = b_fishres$ensembl_gene_id,
  entrezgene_id = b_fishres$entrezgene_id,
  external_gene_name = b_fishres$external_gene_name,
  Cluster = b_fishres$clustr,
  Overfriendliness = b_fishres$clustr_count,
  GO.term = b_fishres$go_term,
  KEGG.pathway = b_fishres$kegg_pathway,
  Wikipathway = b_fishres$wiki_pathway,
  TF = b_fishres$TF,
  BMD.zSD = as.numeric(b_fishres$BMD.zSD),
  Trend = b_fishres$trend
  )


pergene.res.summary$BMD.zSD <- round(pergene.res.summary$BMD.zSD, 1) 
# easier reading of BMD 

pergene.res.summary <- unique(pergene.res.summary)

# fast save to csv file (using multiple CPUs) for easier interpretation analysis
data.table::fwrite(pergene.res.summary, "outputs/lonely_fishres_summary.csv", sep = ";", 
       col.names = TRUE, na = NA) # the string to use for missing values in the data 



# CURVESPLOT FOR ALL TO PDF

b_fishres_no_redund <- subset(b_fishres, select = -c(entrezgene_id,
                                                     ensembl_gene_id,
                                                     external_gene_name,
                                                     clustr_count,
                                                     go_term,
                                                     kegg_pathway,
                                                     wiki_pathway,
                                                     TF))

# these are the column numbers that have been removed !

b_fishres_no_redund <- unique(b_fishres_no_redund)




curves2pdf <- function(b.exp.clustr.data, # output of first or second 'lonely.fishing'
                       responsiv.data, # output of 'regs.in.data'
                       tested.doses, # vector of tested doses
                       log_transfo = TRUE, # If TRUE a log transformation of the
                       # dose is used in the plot. This option 
                       # needs a definition of a strictly positive
                       # value of xmin in input.
                       curvesplot.pdf.filename)
  # filename of curvesplots for each cluster
{
  
  # formatting the responsiv.data to be rbinded to b.exp.clustr.data (this lets use 
  # add a curvesplot of ALL data)
  
  responsiv.data$clustr <- "All" # clustr id for All 
  
  # now creating a df with all clusters, and an added cluster "All" with all 
  # responsive genes
  df.4pdf <- rbind(b.exp.clustr.data, responsiv.data)
  
  
  # Order the data by cluster
  df.4pdf <- df.4pdf[order(suppressWarnings(as.numeric(gsub("Cluster_", "", df.4pdf$clustr)))), ]
  
  # Remove the temporary column used for sorting
  df.4pdf$cluster_numeric <- NULL
  
  # Remove grouping
  df.4pdf <- as.data.frame(df.4pdf)
  
  cluster.names <- unique(df.4pdf$clustr)
  # adding the 2 character clusters to the end (because the as.numeric deleted
  # these ones and we want them at the end)
  cluster.names <- append(cluster.names, c("Lonely", "All"))
  
  
  # create a blank pdf file document to which we shall print each curvesplot 
  pdf(curvesplot.pdf.filename, width = 5, height = 3, onefile = TRUE)
  
  for (i in cluster.names) 
  {
    df.clust <- df.4pdf[df.4pdf$clustr %in% i,]
    cplot <- DRomics::curvesplot(df.clust, addBMD = TRUE, scaling = TRUE, colorby = "trend",
                        npoints= 100, free.y.scales = FALSE,
                        xmin = 0.01, xmax = 100, dose_log_transfo = log_transfo, line.size = 1, line.alpha = 0.4, point.size = 2, point.alpha = 0.4) +
      geom_vline(xintercept = tested.doses, linetype = 2) +
      labs(title = i, subtitle = paste("with", length(unique(df.clust$id)), "items")) +
      xlab("Dose (µg/L)") + ylab("Signal") +
      scale_colour_manual(values = c("inc" = "#1B9E77", "dec" = "#D95F02", "U" = "#7570B3", "bell" = "#E7298A")) +
      theme_light()
    
    print(cplot) # print to pdf file
    
  }
  
  dev.off() # close the finalise pdf file 
  
}

f <- readRDS("data/derived-data/zebra_fitted_05")
tested.doses <- unique(f$omicdata$dose)

require(ggplot2)
curves2pdf(b.exp.clustr.data = b_fishres_no_redund, 
           responsiv.data = b$res, 
           tested.doses = unique(f$omicdata$dose),
           curvesplot.pdf.filename = "outputs/lonely_fishres_curvesplots.pdf")


## modded curvesplot

curves2pdfmod(b.exp.clustr.data = b_fishres_no_redund, 
              responsiv.data = b$res, 
              clustr.fusion.data = clustr_fusionres,
              tested.doses = unique(f$omicdata$dose),
              curvesplot.pdf.filename = "outputs/lonely_fishres_curvesplots_modded.pdf")




sensitivityplot2(b_fishres_no_redund[!b_fishres_no_redund$clustr %in% c("Friendly", "Lonely"), ], 
                 group = "clustr",
                 BMDsummary = "median.and.IQR",
                 BMD_log_transfo = TRUE) +
  geom_point(color = "#302533") +
  guides(size = guide_legend(title = "nb of transcripts")) +
  ylab("BMD 25th quantiles (μg/L)") +
  xlab("Expanded clusters") +
  theme_light()



DRomics::bmdplot(b_fishres_no_redund[b_fishres_no_redund$clustr %in% "Lonely", ], 
        add.CI = TRUE,
        colorby = "trend",
        BMD_log_transfo = TRUE) +
  xlab("BMD.zSD (μg/L)") + 
  scale_colour_manual(values = c("inc" = "#53e078", "dec" = "#fda547", 
                                "U" = "#9590c4", "bell" = "#fc67a2")) +
  theme_light()



