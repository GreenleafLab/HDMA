#!/share/software/user/open/R/4.1.2/bin/Rscript

# Call marker peaks  --------------------------------------------------

# imports -----------------------------------------------------------------------------------------
#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(here)
})

# Get additional functions, etc.:
scriptPath <- here("code/utils")
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# user inputs ---------------------------------------------------------------------------------------
# set working directory
wd <- here("output/01-preprocessing/03")

# set working directory
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
pointSize <- 0.5
addArchRGenome("hg38")
addArchRThreads(8)


# Call marker peaks ----------------------------------------------------------------
proj <- loadArchRProject(paste0(wd, "/allSamples"))

# Identify Marker Peaks while controling for TSS and Depth Biases
# consider adding a for loop to do pairwise comparisons and save everything 

# # test code to check AG_8 (29 cells) and AG_9 (25 cells), 
# #           which led to the error "Cells in foreground and background are 0"
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "RNA_Clusters",
  useGroups = c("TR_11","AG_8", "AG_9"),
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 1000,
  testMethod = "wilcoxon"
)

markerPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "RNA_Clusters",
  useGroups = c("TR_11"),
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 1000,
  testMethod = "wilcoxon"
)


# exclude small clusters with <50 cells from marker peak calling
excl <- table(proj$RNA_Clusters)[table(proj$RNA_Clusters) < 50] %>% names

markerPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  useGroups = unique(proj$RNA_Clusters)[!unique(proj$RNA_Clusters) %in% excl],
  groupBy = "RNA_Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 1000,
  testMethod = "wilcoxon"
)

saveRDS(markerPeaks, paste0(wd, "/markerPeaks_Clusters.rds"))

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
saveRDS(markerList, paste0(wd, "/markerPeaks_Clusters_FDR0.1_Log2FC0.5.rds"))

#Visualize Markers as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  binaryClusterRows = TRUE,
  clusterCols = FALSE,
  transpose = FALSE
)
draw(heatmapPeaks, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapPeaks, name="Peak-Marker-Heatmap-Clusters", width=10, height=15, ArchRProj=proj, addDOC=FALSE)
