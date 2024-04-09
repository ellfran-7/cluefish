# STEP 12 -  GENERATING CURVESPLOTS FOR EAHC CLUSTER TO PDF
# -----------------------------------------------------------------------------
#
# This script generates output PDF file containing a plot of dose-response curves for each cluster of genes, with each plot labeled with the cluster ID and the number of transcripts in that cluster. The curves are color-coded according to whether the trend is increasing, decreasing, U-shaped, or bell-shaped. The plot axes are labeled with "Dose (µg/L)" and "Signal", and the y-axis is scaled to be the same across all plots.
#
# The function used in the script has been developed for this project and can be found in the folder R/, under the name '9-curves2pdf.R'.
#

curves2pdf(
    lonely.fish.data = lonely_fishres$dr_g_a_fishing,
    dromics.data = BMDres_definedCI, 
    clustr.fusion.data = clustr_fusionres$dr_g_a_fusion,
    tested.doses = unique(f$omicdata$dose), 
    addBMD = TRUE,
    scaling = TRUE,
    npoints= 100,
    free.y.scales = FALSE,
    xmin = 0.01, 
    xmax = 100, 
    dose.log.transfo = TRUE, 
    line.size = 0.7, 
    line.alpha = 0.4, 
    point.size = 2, 
    point.alpha = 0.4,
    unit = "µg/L",
    xtitle = "Dose (µg/L)",
    ytitle = "Signal",
    colors = c("inc" = "#1B9E77", "dec" = "#D95F02", "U" = "#7570B3", "bell" = "#E7298A"),
    curvesplot.pdf.filename = "outputs/lonely_fishres_curvesplots.pdf")
