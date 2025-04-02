#> ===========================================================================
#> Script to generate the supplementary figures in Franklin et al. (submitted)
#> ===========================================================================
#> Case:  rnorvegicus - TempO-seq ****
#> Run : 01/04/25 ****



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
file_date = "2025-04-02"



#> Create directory for saving supp. data (if it doesn't already exist) --------
dir_path <- "figures/for-supp"

if (!dir.exists(dir_path)) { # Check if the directory path exists
  
  dir.create(dir_path) # Create it if not
  
}


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
f <- readRDS(file = here::here("data", "derived-data", "fitres_rat_liver_pfoa.rds"))

# The id names are homemade, sciome generated, combining capital gene names and other numbers (of unknown origin)
f_id_mod <- sub("_.*", "", f$omicdata$item) 

# Load DRomics "bmdboot" results filtered with only transcripts with a defined confidence interval around the BMD
b_definedCI <- readRDS(file = here::here("data", "derived-data","bootres_rat_liver_pfoa_seed1234_5000iter_definedCI.rds"))

# Extract and clean gene identifiers from 'b_definedCI'
b_definedCI_mod <- b_definedCI |> 
  dplyr::mutate(id_mod = sub("_.*", "", id)) |> 
  dplyr::select(id, id_mod)

# Normalise all gene IDs to lowercase, except for "LOC" identifiers
b_definedCI_mod <- b_definedCI_onlyids |> 
  dplyr::mutate(
    id_mod = dplyr::if_else(grepl("LOC", id_mod), id_mod, tolower(id_mod))
  )





## Cluefish workflow ----

# Load the getids() result
bg_t_ids <- read.table(paste0("outputs/", file_date, "/bg_t_ids_", file_date, ".txt"))

# Load the clustrenrich() result
clustr_enrichres <- readRDS(here::here("outputs", file_date, paste0("clustr_enrichres_", file_date, ".rds")))

# Retrieve the clustrfusion() result
clustr_fusionres <- clustrfusion(
  clustrenrich_data = clustr_enrichres
)

# Load the lonelyfishing() result
lonelyfishing_data <- readRDS(here::here("outputs", file_date, paste0("lonely_fishres_", file_date, ".rds")))

# Combine DRomics info and the workflow info
b_lonely_fishres <- merge(lonelyfishing_data$dr_t_c_a_fishing, b_definedCI,  by = "id")

# Load the lonely cluster analysis simplenrich() result
lonely_cluster_analysis_res <- readRDS(here::here("outputs", file_date, paste0("lonely_clustr_analysis_res_", file_date, ".rds")))

## Standard workflow ----

# Load the standard workflow results
stand_res <- readRDS(here::here("outputs", file_date, paste0("standard_pipeline_res_", file_date, ".rds")))

# ------------------------------------------------------------------------------



#> -----------------------------------------------------------------------------
#> **Supp. Figure ???.**
#> Title : Principal Component Analysis (PCA) plots showing the data before (A) and after (B) removal of one control sample (replicate) 
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
  filename = here::here("figures", "for-combo", paste0("fig-venndiagramm-1-for-layout-", gsub(":", "", source_vector[1]), "-rnorvegicus.png")),
  height = 2.5, 
  width = 2.5, 
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
  filename = here::here("figures", "for-combo", paste0("fig-venndiagramm-2-for-layout-", gsub(":", "", source_vector[2]), "-rnorvegicus.png")),
  height = 5, #4.5 
  width = 5, #4.5
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
  filename = here::here("figures", "for-combo", paste0("fig-venndiagramm-3-for-layout-", gsub(":", "", source_vector[3]), "-rnorvegicus.png")),
  height = 2.2, 
  width = 2.2, 
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

# Remove unnecessary columns that can cause redundancy in the plot
# Only retain essential columns for the plot generation
b_lonely_fishres_no_redund <- b_lonely_fishres |> 
  dplyr::select(-c(gene_id, TF, old_clustr, friendliness, term_name, term_id, source)) |> 
  dplyr::distinct()

# Exclude the "Lonely" cluster as it is disproportionately large compared to other clusters
# This exclusion helps in visualizing the relative sizes of other clusters more effectively
b_lonely_fishres_no_redund_selected <- b_lonely_fishres_no_redund[!b_lonely_fishres_no_redund$new_clustr %in% "Lonely",]

## Generate the sensitivity plot using the cleaned and filtered data ----
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
    legend.background = element_blank(),
    legend.box.background = element_blank())

sp_cl


# Create an empty plot to serve as a placeholder for your Venn diagrams ----
empty_plot <- ggplot() + 
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Combine the plots using patchwork ----
sp_space_plots <- sp_cl / empty_plot +
  patchwork::plot_layout(heights = c(1, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(family = "Segoe UI"),
    axis.title.x = element_text(size = 10, margin = margin(t = 5)),
    axis.title.y = element_text(size = 10, margin = margin(r = 5)),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(size = 9, hjust = 0.5),
    legend.text = element_text(size = 8),
    plot.tag = element_text(size = 14),
    legend.direction = "vertical",
    legend.box.just = "center"
  )

# Display the combined plot
sp_space_plots


# Path and Filename for saving the figure
path <- here::here("figures", "for-combo")
filename <- paste0("fig-combined-rnorvegicus-", Sys.Date())

# Save the combined plot to a PowerPoint document as an editable vector layer
require(rvg)
require(officer)

# Convert the plot object to a class dml (for editable graphics in PowerPoint)
p_sp <- rvg::dml(ggobj = sp_space_plots)

# Initialize PowerPoint presentation and add a slide with the plot
my_figs <- officer::read_pptx() |> 
  officer::add_slide() |> 
  officer::ph_with(p_sp, ph_location(left = 0, top = 0, height = 7.48, width = 8.66)) |> 
  base::print(target = here::here(path, paste0(filename, ".pptx")))

# Note: The figure was exported to a .pptx format to insert the Venn diagrams and adjust the layout to better align with the paper's guidelines.

#> ==============

