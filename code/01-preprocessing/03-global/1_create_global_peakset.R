#!/share/software/user/open/R/4.1.2/bin/Rscript

# Call peaks on all clusters across organs --------------------------------------------------

# imports -----------------------------------------------------------------------------------------
#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(ggrepel)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(here)
})

scriptPath <- here("code/utils")
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/GO_wrappers.R"))

# user inputs ---------------------------------------------------------------------------------------
# set working directory
wd <- here("output/01-preprocessing/03")
input.dir <- here("output/01-preprocessing/02/")
meta.dir <- here("output/01-preprocessing/02/shared/meta/")

# set working directory
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
pointSize <- 0.5
addArchRGenome("hg38")
addArchRThreads(8)

# Collect all ArrowFiles ---------------------------------------------------------------------------------------
# Manually define which organs to keep
toKeep <- read_tsv(here("code/02-global_analysis/01-organs_keep.tsv"))
organs <- toKeep$organs

# Grab all arrowfile paths
organ.dirs <- sapply(organs, function(organ) {paste0(input.dir, "/", organ, "/atac_preprocess_output/ATAC_obj_clustered_final/ArrowFiles")})
ArrowFiles <- lapply(organ.dirs, function(dir) {list.files(dir, full.names = T)})
ArrowFiles <- ArrowFiles %>% unname() %>% unlist()

# Create ArchRProject ---------------------------------------------------------------------------------------
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = paste0(wd, "/allSamples")
)

saveArchRProject(proj)
proj <- loadArchRProject(paste0(wd, "/allSamples"))

# Transfer annotations ----------------------------------------------------------------

annot.list <- list()
for (i in seq_along(toKeep$code)){
  organ <- toKeep$code[i]
  meta <- readr::read_tsv(paste0(meta.dir, "/", organ, "_meta.txt"))
  df <- data.frame(CellNames=meta$cb, RNA_Clusters=paste0(organ, "_", meta$L1_clusterID))
  annot.list[[organ]] <- df
}
annot <- bind_rows(annot.list)

saveRDS(annot, paste0(wd, "/allSamples_RNA_Clusters.rds"))
annot <- readRDS(paste0(wd, "/allSamples_RNA_Clusters.rds"))

# Add RNA_Clusters to the all sample project
proj <- addCellColData(
  ArchRProj = proj,
  name = "RNA_Clusters",
  cells = annot$CellNames,
  data = paste0(annot$RNA_Clusters),
  force = TRUE
)

# Call Peaks ----------------------------------------------------------------

# Create Group Coverage Files that can be used for downstream analysis
proj <- addGroupCoverages(
  ArchRProj=proj,
  groupBy= "RNA_Clusters",
  minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40, BOR = 50)
)
saveArchRProject(proj)

# Find Path to Macs2 binary
pathToMacs2 <- findMacs2()

# Call Reproducible Peaks w/ Macs2
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "RNA_Clusters",
  peaksPerCell = 500, # The upper limit of the number of peaks that can be identified per cell-grouping in groupBy. (Default = 500)
  pathToMacs2 = pathToMacs2,
  force = TRUE
)
saveArchRProject(proj)

# Add Peak Matrix
proj <- addPeakMatrix(ArchRProj = proj, force = TRUE)

# Get peaks
peaks <- getPeakSet(proj)

# Save peak gRanges object and ArchR project
saveRDS(peaks, paste0(wd, "/peaks_gr.rds"))
saveArchRProject(proj)

# Save per cluster pseudobulked peak matrix
bulkpeak <- getGroupSE(ArchRProj = proj,
                      useMatrix = "PeakMatrix",
                      groupBy = "RNA_Clusters",
                      divideN = F,
                      scaleTo = NULL)
saveRDS(bulkpeak, "allorgan_peakmatrix_pseudobulked.rds")

# Calculate chromVAR deviations ----------------------------------------------------------------
proj <- loadArchRProject(paste0(wd, "/allSamples"))

cisbp2021 <- readRDS(here("data/external/Kartha2022_cisbp/cisBP_human_pfms_2021.rds")) # this is from Kartha et al 2022 Cell Genomics
cisbp_pwms <- TFBSTools::toPWM(cisbp2021)
proj <- addMotifAnnotations(ArchRProj = proj, annoName="Motif", motifPWMs = cisbp_pwms, force = TRUE)

# Add background peaks
proj <- addBgdPeaks(proj, force = TRUE)
saveArchRProject(proj)
#(WARNING: 2h+ use more cores and memory, if it fails go to 1.0.2 release of ArchR)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(proj)

# get deviation matrix and save as rds
message("fetching chromVAR deviations")
dev <- getMatrixFromProject(proj, useMatrix="MotifMatrix", useSeqnames="z")
saveRDS(dev, paste0(wd, "/all_chromvar_z.rds"))

# Count Fragments Per Cluster ----------------------------------------------------------------
proj <- loadArchRProject(paste0(wd, "/allSamples"))

cellTypes <- unique(proj$RNA_Clusters)
min_reads = 5e6

# Count number of total fragments per cluster
readsPerCluster <- proj@cellColData %>% 
  as.data.frame() %>% 
  dplyr::group_by(RNA_Clusters) %>% 
  summarise(total_reads = sum(nFrags))

# reorder and save tsv
readsPerCluster <- readsPerCluster %>% arrange(desc(total_reads))
write_tsv(readsPerCluster, paste0(wd, "/fragmentsPerCluster.tsv"))

# Plot fragments per cluster
readsPerCluster <- read_tsv(paste0(wd, "/fragmentsPerCluster.tsv")) %>% as.data.frame()
nPassing <- length(which(readsPerCluster$total_reads > min_reads))

pdf(paste0(wd, "/fragmentsPerCluster.pdf"), w = 25, h = 4)
ggplot(readsPerCluster, 
       aes(
         x = reorder(RNA_Clusters, -total_reads), 
         y = total_reads)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = min_reads) + 
  theme_ArchR() +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(1, 10e10)
  ) +
  xlab("") +
  ylab("Fragments per cluster") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_text(size = 14)
  ) +
  ggtitle(paste0("Number of clusters with >", min_reads, " reads per cluster = ", nPassing, "/", length(cellTypes)))
dev.off()

# Export Pseudobulked GeneScoreMatrix ----------------------------------------------------------------
bulkGS <- getGroupSE(proj, useMatrix = "GeneScoreMatrix", groupBy = "RNA_Clusters")
saveRDS(bulkGS, paste0(wd, "/GeneScoreByRNAClustersPseodobulkedMatrix.rds"))

