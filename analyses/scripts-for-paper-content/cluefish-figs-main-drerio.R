#> ==================================================================
#> Script to generate the main figures in Franklin et al. (submitted)
#> ==================================================================
#> Case: drerio - RNA-seq ****
#> Run : 12/04/24 ****



#> Load and set everything necessary for the script (e.g., functions, variables)
#> -----------------------------------------------------------------------------

#> Load the project functions --------------------

devtools::load_all(here::here())

#> Load packages --------------------

require(ggplot2)
require(patchwork)
require(here)

#> Load fonts --------------------

require(extrafont)
# First, import the fonts installed on the system, if not already done
# font_import()

# Second, register the fonts with the PDF output device (only necessary in session where you ran font_import())
loadfonts()



#> Set the file_date variable --------------------
file_date = "2024-12-04"



#> Choosing Colors --------------------

# Diverging and Sequential Palette for General Plots ---

# Colorblind-safe palettes were created using:
#   - Chroma.js by Gregor Aisch: [https://gka.github.io/palettes/#/9|s|00429d,96ffea,ffffe0|ffffe0,ff005e,93003a|1|1]
#   - Viz Palette by Elijah Mekks and Susie Lu: [https://projects.susielu.com/viz-palette?colors=[%22#4cb0af%22,%22#e4b5ff%22]&backgroundColor=%22white%22&fontColor=%22black%22&mode=%22normal%22]

# 3 colours for GO, KEGG and WP (following Paul Tol's colours palettes):

# Choice 1:
paul_tol_pal <- c('#77aadd', '#44bb99', '#eedd88')

# Choice 2:
paul_tol_pal <- c('#88CCEE', '#CC6677', '#eedd88')


# Diverging palette ranging from dark teal to dark purple/orchid:
cluefish_pal_div9 <- c('#00393a', '#086969', '#2b9c9b', '#6cd0ce', '#f5f5f5', '#dcadf5', '#a67abf', '#724a8c', '#421d5b')

# Sequential  grey palette:
grey_pal_seq9 <- c('#1f1f1f', '#353535', '#4d4d4d', '#666666', '#808080', '#9b9b9b', '#b7b7b7', '#d4d4d4', '#f2f2f2')

# Qualitative Color Scheme for Curves Plots ---

# Colorblind-safe qualitative colors inspired by Paul Tol's quantitative color schemes:
#   - Reference: [https://personal.sron.nl/~pault/#sec:qualitative]
curves_trend_cols <- c("inc" = "#009988", "dec" = "#EE7733", "U" = "#0077BB", "bell" = "#EE3377")




#> Load the necessary data --------------------

## DRomics pipeline - Transcriptomics data ----

# Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = here::here("data", "derived-data", "fitres_zebrafish_phtalate.rds"))

# Load DRomics "bmdboot" results filtered with only transcripts with a defined confidence interval around the BMD
b_definedCI <- readRDS(file = here::here("data", "derived-data", "bootres_zebrafish_phtalate_UF_seed3_5000iter_definedCI.rds"))



## Cluefish workflow ----

# Load the clustrenrich() result
clustr_enrichres <- readRDS(here::here("outputs", file_date, paste0("clustr_enrichres_", file_date, ".rds")))

# Retrieve the clustrfusion() result
clustr_fusionres <- clustrfusion(
  clustrenrich_data = clustr_enrichres
)

# Load the lonelyfishing() result
lonelyfishing_data <- readRDS(here::here("outputs", file_date, paste0("lonely_fishres_", file_date, ".rds")))

# Modify the "ensembl_transcript_id_version" column to "id" for easier merge and integration to DRomics visualisations
names(lonelyfishing_data$dr_t_c_a_fishing)[names(lonelyfishing_data$dr_t_c_a_fishing) == "transcript_id"] <- "id"

# Combine DRomics info and the workflow info
b_lonely_fishres <- merge(lonelyfishing_data$dr_t_c_a_fishing, b_definedCI,  by = "id")

# Load the lonely cluster analysis simplenrich() result
lonely_cluster_analysis_res <- readRDS(here::here("outputs", file_date, paste0("lonely_clustr_analysis_res_", file_date, ".rds")))

## Standard workflow ----

# Load the standard workflow results
stand_res <- readRDS(here::here("outputs", file_date, paste0("standard_pipeline_res_", file_date, ".rds")))

# ------------------------------------------------------------------------------



#> -----------------------------------------------------------------------------
#> **Main Figure 2.**
#> Title : Results of the Cluefish workflow applied to the zebrafish dataset. 
#> -----------------------------------------------------------------------------

#> In this section, we will generate and combine Venn diagrams, sensitivity plots, and curve plots into a single composite plot. 
#> While the sensitivity and curve plots can be seamlessly integrated using the <patchwork> package, the Venn diagrams are saved as .png files and must be inserted separately. Therefore, during the layout formatting with <patchwork>, space will be allocated specifically for placing the Venn diagrams.
#> -----------------------------------------------------------------------------



#> Generate Venn Diagrams of enriched term content between both workflow =======

#> Create functions to streamline Venn diagram plots ----------

# Venn diagrams ----

venns4paper <- function(
  data1,
  data2,
  data3,
  variable,
  source_id,
  ...,
  category.names = c("" , ""),
  filename = "venn-default-output.png",
  output = TRUE
  )

{

  # Enriched terms kept at the end of the standard approach
  data1_terms <- unique((data1 |> 
                           dplyr::filter(source == source_id))[,variable])
  
  # Enriched terms kept at the end of cluefish
  data2_terms <- unique((data2 |> 
                           dplyr::filter(source == source_id))[,variable])
  
  # Enriched terms kept at the end of the lonely cluster analysis
  data3_terms <- unique((data3 |> 
                           dplyr::filter(source == source_id))[,variable])
  
  # Combine enriched terms from cluefish and the lonely cluster analysis
  data2_and_data3_terms <- unique(c(data2_terms, data3_terms))
  
  # Generate the venn diagram
  VennDiagram::venn.diagram(
    x = list(
      data1_terms, 
      data2_and_data3_terms),
    category.names = category.names,
    filename = filename,
    output = output,
    ...
    )
  
}


## Select the data ----

# the filtered standard dataframe (simplenrich() [highlighted, gene set size limits, ngenes enrich limit])
stand_results <- stand_res$filtered$dr_a

# the final Cluefish dataframe (clustrenrich() [highlighted, gene set size limits, ngenes enrich limit], clustrfusion(), lonelyfishing())
# !!!!!!!!!! This only what is enriched and kept after the various steps !!!!!!!!!!!
cluefish_results <- b_lonely_fishres |> 
  dplyr::filter(term_name %in% clustr_fusionres$dr_g_a_fusion$term_name)

# the lonely cluster analysis dataframe 
lonelycluster_analysis_results <- lonely_cluster_analysis_res$filtered$dr_a

# Create the source vector
source_vector <- unique(na.omit(cluefish_results$source))


## Generating the figures for each source ----

# For each source in the source vector, run the venns4paper() function

# GO:BP
venns4paper(
  data1 = stand_results,
  data2 = cluefish_results,
  data3 = lonelycluster_analysis_results,
  variable = "term_name",
  source_id = source_vector[1],
  filename = here::here("figures", "for-combo", paste0("fig-venndiagramm-1-for-layout-", gsub(":", "", source_vector[1]), "-drerio.png")),
  height = 5, 
  width = 5, 
  units= "cm",
  resolution = 600,
  lwd = c(1, 1),
  lty = c(2, 1),
  fill = c(paul_tol_pal[1], paul_tol_pal[1]), 
  col = c(grey_pal_seq9[2], grey_pal_seq9[2]), 
  cex = .9,
  fontface = "bold",
  fontfamily = "Arial",
  label.col = c(grey_pal_seq9[2], grey_pal_seq9[2], grey_pal_seq9[2]),
  cat.cex = .2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.col = c(grey_pal_seq9[3], cluefish_pal_div9[2]),
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "Arial"
)

# KEGG
venns4paper(
  data1 = stand_results,
  data2 = cluefish_results,
  data3 = lonelycluster_analysis_results,
  variable = "term_name",
  source_id = source_vector[2],
  filename = here::here("figures", "for-combo", paste0("fig-venndiagramm-2-for-layout-", gsub(":", "", source_vector[2]), "-drerio.png")),
  height = 4.5, 
  width = 4.5, 
  units= "cm",
  resolution = 900,
  lwd = c(1, 1),
  lty = c(2, 1),
  fill = c(paul_tol_pal[2], paul_tol_pal[2]), 
  col = c(grey_pal_seq9[2], grey_pal_seq9[2]), 
  cex = .9,
  fontface = "bold",
  fontfamily = "Arial",
  label.col = c("white", grey_pal_seq9[2], grey_pal_seq9[2]),
  cat.cex = .2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.col = c(grey_pal_seq9[3], cluefish_pal_div9[2]),
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "Arial"
)

# WP
venns4paper(
  data1 = stand_results,
  data2 = cluefish_results,
  data3 = lonelycluster_analysis_results,
  variable = "term_name",
  source_id = source_vector[3],
  filename = here::here("figures", "for-combo", paste0("fig-venndiagramm-3-for-layout-", gsub(":", "", source_vector[3]), "-drerio.png")),
  height = 2.5, 
  width = 2.5, 
  units= "cm",
  resolution = 900,
  lwd = c(1, 1),
  lty = c(2, 1),
  fill = c(paul_tol_pal[3], paul_tol_pal[3]), 
  col = c(grey_pal_seq9[2], grey_pal_seq9[2]), 
  cex = .9,
  fontface = "bold",
  fontfamily = "Arial",
  label.col = c("white", grey_pal_seq9[2], grey_pal_seq9[2]),
  cat.cex = .2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.col = c(grey_pal_seq9[3], cluefish_pal_div9[2]),
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "Arial"
)

#> ==============




#> Generate the sensitivity plot of clusters ===================================


## Select and Prepare Data ----

# Generate the sensitivity plot using the cleaned and filtered data
sp_cl <- DRomics::sensitivityplot(
  b_lonely_fishres_no_redund_selected, 
  group = "new_clustr",
  BMDsummary = "first.quartile",
  BMD_log_transfo = TRUE,
  line.size = 0.7, 
  line.alpha = 0.3
) +
  xlab("Clusters") + 
  ylab("BMD 25th quantiles (μg/L) (in log scale)") +
  labs(size = "n° of transcripts") +
  geom_point(color = "#333333") +
  theme_light() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.15),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = grey_pal_seq9[6]))

sp_cl

#> ==============




#> Generate the curvesplots of the four most sensitivity clusters ==============


#> Generate Curvesplots of Interest ----------

## Create Curvesplots for Specific Terms ----


# Filter data for clustr "43" and generate the curves plot -
clustr_choice <- "43"

term_df4drplots <- b_lonely_fishres |> 
  dplyr::filter(new_clustr == clustr_choice)

cp1 <- DRomics::curvesplot(
  term_df4drplots, 
  addBMD = TRUE, 
  scaling = TRUE, 
  colorby = "trend",
  npoints= 100, 
  free.y.scales = FALSE,
  xmin = 0.1, 
  xmax = 100,
  dose_log_transfo = TRUE,
  line.size = 0.5, 
  line.alpha = 0.6,
  point.size = 2, 
  point.alpha = 0.6
) +
  geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.7) +
  labs(x = NULL) +
  xlab("Dose (μg/L) (in log scale)") +
  ylab("Signal (scaled)") +
  scale_color_manual(values = c("inc" = "#009988", "dec" = "#EE7733", "U" = "#0077BB", "bell" = "#EE3377")) +
  theme_light() +
  annotate(geom = "text", 
           x = 0.4,
           y = 0.7, 
           size  = 3,
           label = paste0("Cluster ", clustr_choice))


# Filter data for Cluster #37 and generate the curves plot -
clustr_choice <- "37"

term_df4drplots <- b_lonely_fishres |> 
  dplyr::filter(new_clustr == clustr_choice)

cp2 <- DRomics::curvesplot(
  term_df4drplots, 
  addBMD = TRUE, 
  scaling = TRUE, 
  colorby = "trend",
  npoints= 100, 
  free.y.scales = FALSE,
  xmin = 0.1, 
  xmax = 100,
  dose_log_transfo = TRUE,
  line.size = 0.5, 
  line.alpha = 0.6,
  point.size = 2, 
  point.alpha = 0.6
) +
  geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.7) +
  xlab("Dose (μg/L) (in log scale)") +
  ylab("Signal (scaled)") +
  # labs(x = NULL, y = NULL) +
  guides(color = "none") +  # Hide legend
  scale_color_manual(values = c("inc" = "#009988", "dec" = "#EE7733", "U" = "#0077BB", "bell" = "#EE3377")) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_light() +
  annotate(geom = "text", 
           x = 0.4,
           y = 0.7, 
           size  = 3,
           label = paste0("Cluster ", clustr_choice))


# Filter data for clustr "10" and generate the curves plot -
clustr_choice <- "10"

term_df4drplots <- b_lonely_fishres |> 
  dplyr::filter(new_clustr == clustr_choice)

cp3 <- DRomics::curvesplot(
  term_df4drplots, 
  addBMD = TRUE, 
  scaling = TRUE, 
  colorby = "trend",
  npoints= 100, 
  free.y.scales = FALSE,
  xmin = 0.1, 
  xmax = 100,
  dose_log_transfo = TRUE,
  line.size = 0.5, 
  line.alpha = 0.6,
  point.size = 2, 
  point.alpha = 0.6
) +
  geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.7) +
  xlab("Dose (μg/L) (in log scale)") + 
  ylab("Signal (scaled)") +
  # labs(y = NULL) +
  guides(color = "none") +  # Hide legend
  scale_color_manual(values = c("inc" = "#009988", "dec" = "#EE7733", "U" = "#0077BB", "bell" = "#EE3377")) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_light() +
  annotate(geom = "text", 
           x = 0.4,
           y = 0.7, 
           size  = 3,
           label = paste0("Cluster ", clustr_choice))


# Filter data for clustr '15' and generate the curves plot -
clustr_choice <- "15"

term_df4drplots <- b_lonely_fishres |> 
  dplyr::filter(new_clustr == clustr_choice)

cp4 <- DRomics::curvesplot(
  term_df4drplots, 
  addBMD = TRUE, 
  scaling = TRUE, 
  colorby = "trend",
  npoints= 100, 
  free.y.scales = FALSE,
  xmin = 0.1, 
  xmax = 100,
  dose_log_transfo = TRUE,
  line.size = 0.5, 
  line.alpha = 0.6,
  point.size = 2, 
  point.alpha = 0.6
) +
  geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.7) +
  xlab("Dose (μg/L) (in log scale)") + 
  ylab("Signal (scaled)") +
  # labs(y = NULL) +
  guides(color = "none") +  # Hide legend
  scale_color_manual(values = c("inc" = "#009988", "dec" = "#EE7733", "U" = "#0077BB", "bell" = "#EE3377")) +
  theme_light() +
  annotate(geom = "text", 
           x = 0.4,
           y = 0.7, 
           size = 3,
           label = paste0("Cluster ", clustr_choice))


#> ==============




#> Combine the sensitvity and curves plots with <patchwork> ====================


## Arrange the Curvesplots into a single figure ----

# Combine curves plots into a single layout
cp_plots <- cp1 + cp2 + cp3 + cp4 +
  patchwork::plot_layout(axes = "collect",
                         guides = "collect") +  # Collect guides (legends) into one place
  # patchwork::guide_area() + 
  patchwork::plot_annotation(tag_levels = "A") &  # Add tags (A, B, C) to the plots
  theme(
    text = element_text(family = "Segoe UI"),
    axis.title.x = element_text(size = 9, margin = margin(t = 5)),   # Adjust x-axis label size
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),   # Adjust y-axis label size
    axis.text.x = element_text(size = 8),      # Adjust x-axis tick labels size
    axis.text.y = element_text(size = 8),      # Adjust y-axis tick labels size
    legend.title = element_text(size = 10, hjust = 0.5),   # Adjust legend title size
    legend.text = element_text(size = 9),       # Adjust legend text size
    legend.direction = "vertical",
    legend.box.just = "center"
  ) 

cp_plots

#> Combine the all plots ----

# Create the layout design
layout <- "AAAAA#######
           AAAAA#######
           AAAAABBBBBBB
           AAAAABBBBBBB
           AAAAABBBBBBB
           AAAAABBBBBBB"

sp_cp_plots <- sp_cl + cp_plots + 
  patchwork::plot_annotation(tag_levels = list(c("A", "C", "D", "E", "F"))) +  # Add tags (A, B, C) to the plots
  patchwork::plot_layout(design = layout) &  
  theme(
    text = element_text(family = "Segoe UI"),
    axis.title.x = element_text(size = 10, margin = margin(t = 5)),   # Adjust x-axis label size
    axis.title.y = element_text(size = 10, margin = margin(r = 5)),   # Adjust y-axis label size
    axis.text.x = element_text(size = 8),      # Adjust x-axis tick labels size
    axis.text.y = element_text(size = 8),      # Adjust y-axis tick labels size
    legend.title = element_text(size = 9, hjust = 0.5),   # Adjust legend title size
    legend.text = element_text(size = 8),       # Adjust legend text size
    plot.tag = element_text(size = 14),
    legend.direction = "vertical",
    legend.box.just = "center"
  )

sp_cp_plots

# Path and Filename for saving the figure
path <- here::here("figures", "for-combo")
filename <- paste0("fig-combined-drerio", Sys.Date())

# Save the combined plot to a PowerPoint document as an editable vector layer
require(rvg)
require(officer)

# Convert the plot object to a class dml (for editable graphics in PowerPoint)
p_cp <- rvg::dml(ggobj = sp_cp_plots)

# Initialize PowerPoint presentation and add a slide with the plot
my_figs <- officer::read_pptx() |> 
  officer::add_slide() |> 
  officer::ph_with(p_cp, ph_location(left = 0, top = 0, height = 7.48, width = 8.66)) |> 
  base::print(target = here::here(path, paste0(filename, ".pptx")))

# Note: The figure was exported to a .pptx format to insert the Venn diagrams and adjust the layout to better align with the paper's guidelines.

#> ==============


# ------------------------------------------------------------------------------






#> -----------------------------------------------------------------------------
#> **Main Figure 3.**
#> Title : Dose-response curves and computed benchmark dose (BMD) points for deregulated zebrafish transcripts and individual endpoints.
#> -----------------------------------------------------------------------------

#> In this section, we will generate and combine dose-response curves plots of :
#>  - Hox gene transcripts
#>  - Transcripts associated with eye development and function, with transcripts’ genes from the crystallin family highlighted in pink and other eye-related genes in light grey
#>  - Body length measurements
#>  - Eye surface area measurements
#> -----------------------------------------------------------------------------


#> Curvesplot of eye surface measurements ======================================

# Filter apical data for eye area and generate the curves plot
eyearea_df4drplots <- b_api$res |> 
  dplyr::filter(id == "eye_area") |> 
  dplyr::mutate(color = "#555555") # Adding a uniform colour column

(
  cp_eyearea <- DRomics::curvesplot(
    eyearea_df4drplots, 
    addBMD = TRUE, 
    scaling = FALSE, 
    colorby = "color",
    npoints= 100, 
    free.y.scales = TRUE,
    xmin = 0.1, 
    xmax = 100,
    dose_log_transfo = TRUE,
    line.size = 0.5, 
    line.alpha = 0.6,
    point.size = 2,  
    point.alpha = 0.6
  ) +
    geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.7) +
    xlab("Dose (μg/L) (in log scale)") + 
    ylab(expression(paste("Eye area (", mm^{2}, ")"))) +
    guides(color = "none") +
    scale_color_manual(values = eyearea_df4drplots$color[1]) + 
    theme_light() 
)

#> ==============



#> Curvesplot of eye surface measurements ======================================

# Filter apical data for body length and generate the curves plot
bodylen_df4drplots <- b_api$res |> 
  dplyr::filter(id == "body_length") |> 
  dplyr::mutate(color = "#555555") # this is just to create a uniform column 

(
  cp_bodylen <- DRomics::curvesplot(
    bodylen_df4drplots, 
    addBMD = TRUE, 
    scaling = FALSE, 
    colorby = "color",
    npoints= 100, 
    free.y.scales = FALSE,
    xmin = 0.1, 
    xmax = 100,
    dose_log_transfo = TRUE,
    line.size = 0.5, 
    line.alpha = 0.6,
    point.size = 2,  
    point.alpha = 0.6
  ) +
    geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.7) +
    xlab("Dose (μg/L) (in log scale)") + 
    ylab("Body length (mm)") +
    guides(color = "none") +
    scale_color_manual(values = bodylen_df4drplots$color[1]) + 
    theme_light()
)

#> ==============




#> Curvesplot of deregulated Hox gene transcripts ==============================

# Filter transcriptomic data for Hox genes and generate the curves plot
hoxgenes_df4drplots <- b_lonely_fishres |> 
  dplyr::filter(grepl("hox", gene_name)) |> 
  dplyr::mutate(color = "#555555") # this is just to create a uniform column 

(
  cp_hoxgenes <- DRomics::curvesplot(
    hoxgenes_df4drplots, 
    addBMD = TRUE, 
    scaling = TRUE, 
    colorby = "color", 
    npoints= 100, 
    free.y.scales = FALSE,
    xmin = 0.1, 
    xmax = 100,
    dose_log_transfo = TRUE,
    line.size = 0.5, 
    line.alpha = 0.6,
    point.size = 2, 
    point.alpha = 0.6
  ) +
    geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.7) +
    xlab("Dose (μg/L) (in log scale)") + 
    ylab("Signal (scaled)") +
    guides(color = "none") +
    scale_color_manual(values = hoxgenes_df4drplots$color[1]) + 
    scale_y_continuous(limits = c(-1, 1)) +
    theme_light() +
    annotate(geom = "text", 
             x = 0.5,
             y = -0.6, 
             size  = 2.5,
             label = paste0("Hox genes"))
  
)

#> ==============



#> Curvesplot of deregulated eye development related gene transcripts ==========

# Filter transcriptomic data for eye development related genes and generate the curves plot

# Get all the annotations of all the deregulated genes
whole <- clustr_enrichres$dr_g_a_whole

# GO term list to visualize
eye_terms <- c("lens development in camera-type eye",
               "camera-type eye development",
               "eye development",
               "visual system development",
               "visual perception",
               "sensory perception of light stimulus",
               "sensory system development",
               "sensory organ development")

# Get all the genes associated to these functions
eye_genes <- whole |> 
  dplyr::filter(term_name %in% eye_terms)

# OR pattern matching with of key words

eye_genes <- whole |> 
  dplyr::filter(grepl("lens", term_name) |
                  grepl("eye", term_name) |
                  grepl("visual", term_name) |
                  # grepl("light", term_name) |
                  grepl("retina", term_name))

# Retrieve the gene names for these ensembl genes
eye_gene_ids <- merge(eye_genes, bg_t_ids, by = "gene_id")

# Modify the "ensembl_transcript_id_version" column to "id" for easier merge and integration to DRomics visualisations
names(eye_gene_ids)[names(eye_gene_ids) == "transcript_id"] <- "id"

# Combine DRomics info and the workflow info
b_eyegenes <- merge(eye_gene_ids, b_definedCI,  by = "id")

#  Add a logical column if the gene is a crystallin or not
crystalgenes_df4drplots <- b_eyegenes |> 
  dplyr::mutate(gene_family = ifelse(grepl("cry", gene_name), "cristallin", "other eye-related processes"))

(
  cp_crystaleyegenes <- DRomics::curvesplot(
    crystalgenes_df4drplots, 
    addBMD = TRUE, 
    scaling = TRUE, 
    colorby = "gene_family",
    npoints= 100, 
    free.y.scales = FALSE,
    xmin = 0.1, 
    xmax = 100,
    dose_log_transfo = TRUE,
    line.size = 0.5, 
    line.alpha = 0.6,
    point.size = 2, 
    point.alpha = 0.6
  ) +
    geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.7) +
    xlab("Dose (μg/L) (in log scale)") + 
    ylab("Signal (scaled)") +
    labs(color = "Gene family") +
    scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = c("cristallin" = "#f6618c", "other eye-related processes" = "#aabcc2")) +
    theme_light()
  
)

#> ==============




#> Combine all previous curvesplots into one singular plot =====================

cp_plots <- (cp_hoxgenes + cp_crystaleyegenes) / (cp_bodylen + cp_eyearea) + 
  patchwork::plot_annotation(tag_levels = "A") +  # Add tags (A, B, C) to the plots
  patchwork::plot_layout(width = 1) &  # Collect guides (legends) into one place
  theme(
    text = element_text(family = "Segoe UI"),
    axis.title.x = element_text(size = 8, margin = margin(t = 5)),   # Adjust x-axis label size
    axis.title.y = element_text(size = 8, margin = margin(r = 5)),   # Adjust y-axis label size
    axis.text.x = element_text(size = 7),      # Adjust x-axis tick labels size
    axis.text.y = element_text(size = 7),      # Adjust y-axis tick labels size
    legend.title = element_text(size = 9, hjust = 0.5),   # Adjust legend title size
    legend.text = element_text(size = 7),       # Adjust legend text size
    legend.direction = "vertical",
    legend.box.just = "center"
  )

cp_plots

# Path and filename for saving the figure
path <- here::here("figures", "for-supp")
filename <- paste0("supp-fig-curves-plots-drerio", Sys.Date())

# Save the combined plot as a high-quality PNG
ggsave(
  filename = here::here(path, paste0(filename, ".png")), 
  plot = cp_plots, 
  dpi = 600, 
  height = 12, 
  width = 16, 
  units = "cm"
)

#> ==============