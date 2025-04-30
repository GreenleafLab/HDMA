#!/share/software/user/open/R/4.1.2/bin/Rscript

# Re-add gene expression matrix to ATAC objects and redo peak2gene links --------------------------------------------------
# use decontX assay instead of RNA assay from Seurat object

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
args <- commandArgs(trailingOnly = T)
organ_name <- args[1]

#organ_name <- "Adrenal" # test code 

wd <- sprintf(here("output/01-preprocessing/02/%s/atac_preprocess_output"), organ_name)
outdir <- paste0(wd, "/ATAC_obj_clustered_peaks_final")
rna_wd <- sprintf(here("output/01-preprocessing/02/%s/rna_preprocess_output"), organ_name)
rna_gtf <- here("data/reference/gtf/gencode.v41.annotation.BPfiltered.gtf")

# set working directory
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
pointSize <- 0.5
addArchRGenome("hg38")
addArchRThreads(1)

# main ----------------------------------------------------------------------------------------------------

## Add gene expression matrix to ArchRProject (no peak matrix) -------------------------------------------------------------------
atac_proj <- loadArchRProject(paste0(wd, "/ATAC_obj_clustered_final"))
atac_proj <- saveArchRProject(atac_proj, paste0(wd, "/ATAC_obj_clustered_final_decontx"), load=T)

rna_proj <- paste0(rna_wd, "/RNA_obj_clustered_final.rds") %>% readRDS
gtf.file <- rna_gtf

atac_proj <- addSeuratGeneExpressionMatrix(atac_proj, gtf.file, rna_proj, assay="decontX")

saveArchRProject(atac_proj, dropCells = TRUE)
rm(atac_proj)

## Add gene expression matrix to ArchRProject final with peaks called -----------------------------------------------
proj <- loadArchRProject(outdir)
proj <- saveArchRProject(proj, paste0(wd, "/ATAC_obj_clustered_peaks_final_decontx"), load=T)

proj <- addSeuratGeneExpressionMatrix(proj, gtf.file, rna_proj, assay="decontX")

saveArchRProject(proj, dropCells = TRUE)

## Calculate peak-to-gene linkages ----------------------------------------------------------

# fixing pca cell name order
proj@reducedDims$pca$matDR <- proj@reducedDims$pca$matDR[getCellNames(proj),]

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


# browser tracks
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "RNA_NamedCluster_L2", 
  geneSymbol = hrg[1:10, "gene"], 
  upstream = 250000,
  downstream = 250000,
  loops = getPeak2GeneLinks(proj)
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-filteredgenes-250k.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

# with gene expression
# remove genes with hyphen in the name due to parsing error 
candidates <- hrg[1:10, "gene"]
candidates <- candidates[!grepl("-", candidates)]

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "RNA_NamedCluster_L2", 
  geneSymbol = candidates, 
  upstream = 250000,
  downstream = 250000,
  loops = getPeak2GeneLinks(proj),
  useMatrix = "GeneExpressionMatrix"
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-filteredgenes-250k-gex.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 10, height = 5)

