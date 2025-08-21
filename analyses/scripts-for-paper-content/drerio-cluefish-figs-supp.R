#> ===========================================================================
#> Script to generate the supplementary figures in Franklin et al. (2025)
#> ===========================================================================
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

## DRomics pipeline - Anchoring data ----

# Load DRomics "bmdboot" results of apical data
b_api <- readRDS(file = here::here("data", "derived-data", "bootres_apical_zebrafish_phtalate_UF_seed3_5000iter.rds"))



## DRomics pipeline - Transcriptomics data ----

# Load DRomics drcfit object (which holds the background transcript list) 
f <- readRDS(file = here::here("data", "derived-data", "fitres_zebrafish_phtalate.rds"))

# Load DRomics "bmdboot" results filtered with only transcripts with a defined confidence interval around the BMD
b_definedCI <- readRDS(file = here::here("data", "derived-data","bootres_zebrafish_phtalate_UF_seed3_5000iter_definedCI.rds"))



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
#> **Supp. Figure 1.**
#> Title : Principal Component Analysis (PCA) plots showing the data before (A) and after (B) removal of one control sample (replicate) 
#> -----------------------------------------------------------------------------

# Set seed in order to reproduce results  ---------

set.seed(1234) # Fixing the seed to reproduce the results


# Original RNAseq data ---------

# Load the derived data used as input in the dromics pipeline
zebra_dbp_df <- read.table("data/derived-data/zebra_dbp_df.txt")

# Re-run RNAseqdata to reproduce the PCA plot
(og_o <- DRomics::RNAseqdata(zebra_dbp_df, # or use the file directory
                          transfo.method = 'rlog', 
                          round.counts = TRUE)) 

og_pca <- DRomics::PCAdataplot(og_o) + ggplot2::theme_bw()



# RNAseq data with a control sample removed ---------

# Remove the 0B replicate (the second column) from the dataframe.
zebra_dbp_df_rm0 <- zebra_dbp_df[, -2]

# Reset column names to V1, V2, V3, etc. after removing the column.
colnames(zebra_dbp_df_rm0) <- paste0("V", seq_along(colnames(zebra_dbp_df_rm0)))

# Re-run RNAseqdata to reproduce the PCA plot
(o <- DRomics::RNAseqdata(zebra_dbp_df_rm0,  
                          transfo.method = 'rlog', 
                          round.counts = TRUE)
)

pca <- DRomics::PCAdataplot(o) + ggplot2::theme_bw()



# Create a combined plot using <patchwork> ---------

require(patchwork)
require(ggplot2)

pca_plots <- og_pca + pca +
  patchwork::guide_area() + 
  patchwork::plot_annotation(tag_levels = "A") +  # Add tags (A, B, C) to the plots
  patchwork::plot_layout(guides = "collect") &  # Collect guides (legends) into one place
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

pca_plots


# Path and Filename for saving the figure
path <- here::here("figures", "for-supp")
filename <- paste0("supp-fig-pca-before-after-drerio", Sys.Date())

# Save the combined plot as a high-quality PNG
ggsave(
  filename = here::here(path, paste0(filename, ".png")), 
  plot = pca_plots, 
  dpi = 600, 
  height = 8, 
  width = 20, 
  units = "cm"
)

# ------------------------------------------------------------------------------





#> -----------------------------------------------------------------------------
#> **Supp. Figure 2.**
#> Title : Dose-response fitted curve and computed BMD point for the rxraa gene transcript (ENSDART00000080481.6)
#> -----------------------------------------------------------------------------

# Filter the significantly deregulated rxraa transcript
t_rxra_data <- b_definedCI |> 
  dplyr::filter(id == "ENSDART00000080481.6") |> 
  dplyr::mutate(color = "#555555") # This is a color column just so that we can produce a curvesplot() by color

(
  cp_rxraa <- DRomics::curvesplot(
    t_rxra_data, 
    addBMD = TRUE, 
    scaling = TRUE, 
    colorby = "color", 
    npoints= 100, 
    free.y.scales = FALSE,
    xmin = 0.1, 
    xmax = 100,
    dose_log_transfo = TRUE,
    line.size = 1, 
    line.alpha = 0.6,
    point.size = 4, 
    point.alpha = 0.6
  ) +
    geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.7) +
    xlab("Dose (μg/L) (in log scale)") + 
    ylab("Signal (scaled)") +
    guides(color = "none") +
    scale_color_manual(values = t_rxra_data$color[1]) + 
    scale_y_continuous(limits = c(-1, 1)) +
    theme_light()
)

# Path and Filename for saving the figure
path <- here::here("figures", "for-supp")
filename <- paste0("supp-fig-rxraa-drerio", Sys.Date())

ggsave(
  filename = here::here(path, paste0(filename, ".png")), 
  plot = cp_rxraa, 
  dpi = 600, 
  height = 12, 
  width = 12, 
  units = "cm"
)

#> -----------------------------------------------------------------------------





#> -----------------------------------------------------------------------------
#> **Supp. Figure 3.**
#> Title : Plots of the raw dose-response data for the transcripts of three Cyp26 genes: (A) cyp26a1 (ENSDART00000041728.7), (B) cyp26b1 (ENSDART00000110347.3), and (C) cyp26c1 (ENSDART00000077809.5).
#> -----------------------------------------------------------------------------

# Filter the Cyp26a transcript - cyp26a1
cyp26a <- bg_t_ids[grepl("cyp26a", bg_t_ids$gene_name),]$transcript_id

g_cyp26a_data <- DRomics::targetplot(cyp26a, f, dose_log_transfo = FALSE) +
  theme_bw() +
  geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.5) +
  ylab("Signal (scaled)") +
  labs(x = NULL)

# Filter the Cyp26b transcript - cyp26b1
cyp26b <- bg_t_ids[grepl("cyp26b", bg_t_ids$gene_name),]$transcript_id

g_cyp26b_data <- DRomics::targetplot(cyp26b, f, dose_log_transfo = FALSE) +
  theme_bw() +
  geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.5) +
  xlab("Dose (μg/L) (in log scale)") +
  labs(y = NULL)

# Filter the Cyp26c transcript - cyp26c1
cyp26c <- bg_t_ids[grepl("cyp26c", bg_t_ids$gene_name),]$transcript_id

g_cyp26c_data <- DRomics::targetplot(cyp26c, f, dose_log_transfo = FALSE) +
  theme_bw() +
  geom_vline(xintercept = unique(f$omicdata$dose), linetype = 3, colour = "#999999", alpha = 0.7, linewidth = 0.5) +
  xlab("Dose (μg/L) (in log scale)") +
  ylab("Signal (scaled)")

# Create a combined plot foolowing a specifc layout using <patchwork> ---------
layout <- '
AB
C#'

cyp26_plots <- g_cyp26a_data + g_cyp26b_data + g_cyp26c_data +
  patchwork::plot_layout(design = layout) &  
  patchwork::plot_annotation(tag_levels = "A") &  # Add tags (A, B, C) to the plots
  theme(
    text = element_text(family = "Segoe UI"),
    axis.title.x = element_text(size = 9, margin = margin(t = 5)),   # Adjust x-axis label size
    axis.title.y = element_text(size = 9, margin = margin(r = 5)),   # Adjust y-axis label size
    axis.text.x = element_text(size = 8),      # Adjust x-axis tick labels size
    axis.text.y = element_text(size = 8),      # Adjust y-axis tick labels size
    legend.title = element_text(size = 10, hjust = 0.5),   # Adjust legend title size
    legend.text = element_text(size = 8),       # Adjust legend text size
    legend.direction = "vertical",
    legend.box.just = "center"
  )

cyp26_plots

# Path and Filename for saving the figure
path <- here::here("figures", "for-supp")
filename <- paste0("supp-fig-cyp26abc-plots-drerio", Sys.Date())

# Save the combined plot as a high-quality PNG
ggsave(
  filename = here::here(path, paste0(filename, ".png")), 
  plot = cyp26_plots, 
  dpi = 600, 
  height = 14, 
  width = 14, 
  units = "cm"
)

#> -----------------------------------------------------------------------------

