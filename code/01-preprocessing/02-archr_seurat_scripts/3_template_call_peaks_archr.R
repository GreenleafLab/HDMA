#!/share/software/user/open/R/4.1.2/bin/Rscript

# Call peaks on clusters from overall clustering --------------------------------------------------

# imports -----------------------------------------------------------------------------------------
#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(Seurat)
  library(tidyr)
  library(ggrepel)
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
organ_code <- "##"
organ_batch <- "##" #e.g. 11
wd <- here("output/01-preprocessing/02/##ORGAN##/atac_preprocess_output")
outdir <- paste0(wd, "/ATAC_obj_clustered_peaks_final")
rna_wd <- here("output/01-preprocessing/02/##ORGAN##/rna_preprocess_output")
rna_gtf <- here("data/reference/gtf/gencode.v41.annotation.BPfiltered.gtf")

meta.path <- here("output/01-preprocessing/02/shared/SampleMetadata.csv")
cutoffs.path <- here("output/01-preprocessing/02/shared/sample_filter_cutoffs_metadata.csv")

# post RNA preprocess final cell list (not sample meta data)
dir_meta <- paste0(here("output/01-preprocessing/02/shared/meta/"), organ_code, "_meta.txt")

# set working directory
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
pointSize <- 0.5
addArchRGenome("hg38")

# load meta data
meta <- read.csv(meta.path)
rownames(meta) <- meta$SampleName
samp.names <- meta$SampleName[meta$Batch==organ_batch]
samp.sex <- meta[samp.names, "Sex"] # this will be empty until we manually fill it after looking at rna
samp.age <- paste0("PCW",meta[samp.names,"PCW"])
names(samp.sex) <- samp.names
names(samp.age) <- samp.names

# main ----------------------------------------------------------------------------------------------------
# once RNA preprocess is complete, run this script

## Filter using combined (RNA + ATAC) filtered cells ------------------------------------------------------
atac_proj <- loadArchRProject(paste0(wd, "/filtered_r1"))

# Read round 1 white list that's been filtered on both ATAC and RNA
r1_meta <- readr::read_tsv(dir_meta)

cellnames <- atac_proj$cellNames[which(getCellNames(atac_proj) %in% r1_meta$cb)]

atac_proj <- atac_proj[cellnames, ]

# Add age 
atac_proj$age <- samp.age[atac_proj$Sample] %>% unlist()

## Add cluster labels --------------------------------------------------------------------------
meta_rna <- r1_meta %>% as.data.frame() %>% dplyr::select(L1_clusterID)
meta_rna$CellNames <- r1_meta$cb

meta_atac <- atac_proj@cellColData %>% as.data.frame()
meta_atac$CellNames <- rownames(meta_atac)

meta <- inner_join(meta_atac, meta_rna, by = "CellNames")

# Drop any additional cells that needs to be dropped
atac_proj <- atac_proj[meta$CellNames,]

# Transfer cluster labels
atac_proj$RNA_Clusters <- paste0("c", meta$L1_clusterID)

saveArchRProject(atac_proj, outputDirectory = paste0(wd, "/ATAC_obj_clustered_final"), dropCells = TRUE)

## Add gene expression matrix and PCA to ArchRProject -------------------------------------------------------------------
atac_proj <- loadArchRProject(paste0(wd, "/ATAC_obj_clustered_final"))
rna_proj <- paste0(rna_wd, "/RNA_obj_clustered_final.rds") %>% readRDS
gtf.file <- rna_gtf

atac_proj <- addSeuratGeneExpressionMatrix(atac_proj, gtf.file, rna_proj)
atac_proj <- addSeuratPCAtoArchRProj(rna_proj, reduction = "pca", atac_proj, convertCellName = F)

saveArchRProject(atac_proj, outputDirectory = paste0(wd, "/ATAC_obj_clustered_final"), dropCells = TRUE)

# replot ATAC TSS vs nFrag QC plot on final cells 
plotTSSnFrag(atac_proj, cutoffs.path)

# replot the QC violin plots on final cells
plotList <- list()
plotList[[1]] <- plotGroups(ArchRProj = atac_proj, groupBy = "Sample", colorBy = "cellColData", 
                            name = "log10(nFrags)", plotAs = "violin", addBoxPlot = TRUE, alpha = 0.4)
plotList[[2]] <- plotGroups(ArchRProj = atac_proj, groupBy = "Sample", colorBy = "cellColData", 
                            name = "TSSEnrichment", plotAs = "violin", addBoxPlot = TRUE, alpha = 0.4)
plotPDF(plotList = plotList, name = "Violin-TSS-log10(nFrag)-filtered", width = 4, height = 4,  ArchRProj = atac_proj, addDOC = FALSE)

## Copy ArchR project to new directory for adding peak information -----------------------------------------------
#atac_proj <- loadArchRProject(paste0(wd, "/ATAC_obj_clustered_final"))

proj <- saveArchRProject(
  ArchRProj = atac_proj,
  outputDirectory = outdir, load=T)

## ATAC dim reduction with iterative LSI -------------------------------------------------------------------
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", 
                        iterations = 2, varFeatures = 25000, dimsToUse = 1:30)

# calculate correlation of LSI dims with TSS enrichment
cor_LSI_TSS <- lapply(1:dim(proj@reducedDims$IterativeLSI$matSVD)[2], function(n){
  cor.test(proj@reducedDims$IterativeLSI$matSVD[,n],proj$TSSEnrichment,method="spearman")$estimate
}) %>% unlist %>% abs
names(cor_LSI_TSS) <- 1:length(cor_LSI_TSS)
p <- ggplot() + geom_point(aes(x=1:length(cor_LSI_TSS),y=cor_LSI_TSS)) +
  geom_line(aes(x=1:length(cor_LSI_TSS),y=cor_LSI_TSS)) + theme_classic() +
  xlab("LSI dimension") + ylab("|Spearman rho|")
p
plotPDF(p, name="Correlation_LSI_TSSEnrichment", ArchRProj = proj, width=5, height=5, addDOC = F)

cor_LSI_depth <- proj@reducedDims$IterativeLSI$corToDepth$scaled
p <- ggplot() + geom_point(aes(x=1:length(cor_LSI_depth),y=cor_LSI_depth)) +
  geom_line(aes(x=1:length(cor_LSI_depth),y=cor_LSI_depth)) + theme_classic() +
  xlab("LSI dimension") + ylab("correlation")
plotPDF(p, name="Correlation_LSI_depth", ArchRProj = proj, width=5, height=5, addDOC = F)

# cluster
proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat",
  name = "ATAC_Clusters", resolution = 0.8)

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", name = "UMAP", 
                nNeighbors = 30, minDist = 0.5, metric = "cosine", force=T)

# add some meta
# Read round 1 white list that's been filtered on both ATAC and RNA
r1_meta <- readr::read_tsv(dir_meta) %>% as.data.frame
rownames(r1_meta) <- r1_meta$cb
proj$RNA_NamedCluster_L1 <- r1_meta[proj$cellNames,"L1_clusterName"]
proj$RNA_NamedCluster_L1b <- r1_meta[proj$cellNames,"L2_clusterID"]
proj$RNA_NamedCluster_L2 <- r1_meta[proj$cellNames,"L2_clusterName"]

# 
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "ATAC_Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "RNA_Clusters", embedding = "UMAP") # this is RNA cluster
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "RNA_NamedCluster_L1", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "RNA_NamedCluster_L1b", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "RNA_NamedCluster_L2", embedding = "UMAP")

plotPDF(p1,p2,p3,p4,p5, name = "Plot-UMAP-ATAC-RNA-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(proj)

## Call Peaks ---------------------------------------------------------------------------------------------------
#proj <- loadArchRProject(outdir)

# Create Group Coverage Files that can be used for downstream analysis
proj <- addGroupCoverages(
  ArchRProj=proj, 
  groupBy= "RNA_Clusters", 
  #minCells = 50, # The minimum number of cells required in a given cell group to permit insertion coverage file generation. (default = 40, BOR = 50)
  force=TRUE
  )

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

# Add Peak Matrix
proj <- addPeakMatrix(ArchRProj = proj, force = TRUE)

# Save project
saveArchRProject(proj)


## Identifying Marker Peaks -------------------------------------------------------------------------------------

# Identify Marker Peaks while controling for TSS and Depth Biases
markerPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "RNA_Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    maxCells = 1000,
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

# p <- plotBrowserTrack(
#   ArchRProj = proj,
#   groupBy = "RNA_Clusters",
#   geneSymbol = c("MYOCD"),
#   upstream = 50000,
#   downstream = 50000
# )
# grid::grid.draw(p$MYOCD)

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

saveRDS(markerPeaks, paste0(outdir, "/markerPeaks_Clusters.rds"))

## Motif annnotations and chromVAR deviation -------------------------------------------------------------------

# Add motif annotations based on Buenrostro lab's cleaned up cisbp pwms
#proj <- addMotifAnnotations(ArchRProj = proj, name = "Motif", motifSet = "cisbp", force = TRUE)
cisbp2021 <- readRDS(here("data/external/Kartha2022_cisbp/cisBP_human_pfms_2021.rds")) # this is from Kartha et al 2022 Cell Genomics
cisbp_pwms <- TFBSTools::toPWM(cisbp2021)
proj <- addMotifAnnotations(ArchRProj = proj, annoName="Motif", motifPWMs = cisbp_pwms, force = TRUE)

# Add background peaks
proj <- addBgdPeaks(proj, force = TRUE)

#(WARNING: 2h+ use more cores and memory, if it fails go to 1.0.2 release of ArchR)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(proj)

## Plot Motif Enrichment in Marker Peaks ------------------------------------------------------------------------

#Identify Motif Enrichments in marker peaks
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Subset to clusters that have at least some enrichment
log10pCut <- 10

#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
    enrichMotifs, 
    n=5, # tunable param default is 10
    transpose=FALSE, 
    cutOff=log10pCut
)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot")
plotPDF(heatmapEM, name="Motifs-Enriched-Heatmap-Clusters", width=8, height=12, ArchRProj=proj, addDOC=FALSE)

saveArchRProject(proj)


## Plot motif deviations across all the clusters ------------------------------------------------------------

# Plot motifs that are the most variable across the dataset
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(proj)


## Calculate coaccessibility and peak-to-gene linkages ----------------------------------------------------------
# should this be combined LSI?

# Calculate coaccessibility
proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "IterativeLSI"
)

# Calculate peak-to-gene links
proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  useMatrix = "GeneExpressionMatrix",
  #reducedDims = "IterativeLSI",
  reducedDims = "pca" # new: use RNA PCA dim reduct for peak 2 gene links
)

saveArchRProject(proj)

# filter out RNA genes not present in the ArchR annotation (mostly antisense and lncRNA)
p2gtmp <- metadata(proj@peakSet)$Peak2GeneLinks

rna_genes <- metadata(p2gtmp)$geneSet
atac_genes <- getGenes()
genespf <- which(rna_genes$name %in% atac_genes$symbol)
p2gpf <- p2gtmp[p2gtmp$idxRNA %in% genespf,]
metadata(proj@peakSet)$Peak2GeneLinks <- p2gpf


# plot peak 2 gene
p1 <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "RNA_Clusters", k=length(unique(proj$RNA_Clusters)), seed=1)
p2 <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "RNA_Clusters", k=2*length(unique(proj$RNA_Clusters)), seed=1)
p3 <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "RNA_NamedCluster_L2", k=length(unique(proj$RNA_Clusters)), seed=1)

plotPDF(p1,p2,p3, name = "Peak2Gene-Heatmap-filteredgenes", width=12, height=12, ArchRProj=proj, addDOC=F)

## GO enrichment of peak2gene kmeans groups ----------------------------------------------------
p2gMat <- plotPeak2GeneHeatmap(proj, groupBy="RNA_Clusters", returnMatrices=TRUE, k=length(unique(proj$RNA_Clusters)), seed=1)

# Get association of peaks to clusters
kclust_df <- data.frame(kclust=p2gMat$ATAC$kmeansId, peakName=p2gMat$Peak2GeneLinks$peak, gene=p2gMat$Peak2GeneLinks$gene)

# Fix peakname
kclust_df$peakName <- sapply(kclust_df$peakName, function(x){strsplit(x, ":|-")[[1]] %>% paste(.,collapse="_")})

kclust <- unique(kclust_df$kclust) %>% sort()
all_genes <- kclust_df$gene %>% unique() %>% sort()

# Save table of top linked genes per kclust
nGOgenes <- 200
nclust <- length(kclust)
topKclustGenes <- lapply(kclust, function(k){
  kclust_df[kclust_df$kclust == k,]$gene %>% getFreqs() %>% head(nGOgenes) %>% names()
}) %>% do.call(cbind,.)

plotDir <- paste0(getOutputDirectory(proj), "/Plots/Peak2Gene/")
dir.create(plotDir, recursive=T, showWarnings=F)
outfile <- paste0(plotDir, sprintf("/top%s_genes_kclust_k%s_filteredgenes.tsv",nGOgenes, nclust))
write.table(topKclustGenes, file=outfile, quote=FALSE, sep='\t', row.names = FALSE, col.names=TRUE)

# plot highly regulated genes
hrg <- kclust_df$gene %>% getFreqs() %>% as.data.frame
colnames(hrg) <- "NumberLinkedPeaks"
hrg$gene <- rownames(hrg)
hrg$order <- 1:dim(hrg)[1]
p <- ggplot(hrg, aes(x=order, y=NumberLinkedPeaks, label=gene)) + geom_point() + xlab("Ranked gene list") +
  ylab("Number of Linked Peaks") + geom_label_repel(data=hrg[1:20,], max.overlaps=10) +
  theme_classic()
ggsave(paste0(plotDir,"/NumLinkedPeaks_vs_RankedGenes.pdf"), width=5, height=5)

GOresults <- lapply(kclust, function(k){
  message(sprintf("Running GO enrichments on k cluster %s...", k))
  clust_genes <- topKclustGenes[,k]
  upGO <- rbind(
    calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="BP") 
    #calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="MF")
    #calcTopGo(all_genes, interestingGenes=upGenes, nodeSize=5, ontology="CC")
  )
  upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
})

names(GOresults) <- paste0("cluster_", kclust)

# Plots of GO term enrichments:
pdf(paste0(plotDir, sprintf("/kclust_GO_15termsBPonlyBarLim_k%s_filteredgenes.pdf", nclust)), width=10, height=8)
for(name in names(GOresults)){
  goRes <- GOresults[[name]]
  if(nrow(goRes)>1){
    print(topGObarPlot(goRes, cmap = cmaps_BOR$comet, 
                       nterms=15, border_color="black", 
                       barwidth=0.85, title=name, barLimits=c(0, 5)))
  }
}
dev.off()

# browser tracks
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "RNA_NamedCluster_L2", 
  geneSymbol = hrg[1:10, "gene"], 
  upstream = 100000,
  downstream = 100000,
  loops = getPeak2GeneLinks(proj)
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-filteredgenes.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)



## ATAC dim reduct with peak matrix -------------------------------------------------------------------
proj <- loadArchRProject(outdir)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "PeakMatrix", name = "IterativeLSI_Peak", 
                        iterations = 2, varFeatures = 25000, dimsToUse = 1:30)

# cluster
proj <- addClusters(input = proj, reducedDims = "IterativeLSI_Peak", method = "Seurat",
                    name = "ATAC_Clusters_Peak", resolution = 0.8)

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI_Peak", name = "UMAP_Peak", 
                nNeighbors = 30, minDist = 0.5, metric = "cosine", force=T)


p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "ATAC_Clusters", embedding = "UMAP_Peak")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "RNA_Clusters", embedding = "UMAP_Peak")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "RNA_NamedCluster_L1", embedding = "UMAP_Peak")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "RNA_NamedCluster_L1b", embedding = "UMAP_Peak")
p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "RNA_NamedCluster_L2", embedding = "UMAP_Peak")

plotPDF(p1,p2,p3,p4,p5, name = "Plot-UMAP-ATAC-RNA-Clusters-PeakMatrix.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(proj)
