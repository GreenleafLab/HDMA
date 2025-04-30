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

# Call Peak2Gene Linkages ----------------------------------------------------------------
#proj <- loadArchRProject(paste0(wd, "/allSamples")) # use raw RNA counts
proj <- loadArchRProject(paste0(wd, "/allSamples_decontx")) # or use decontxed RNA counts

# Add Iterative LSI
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "PeakMatrix", name = "IterativeLSI_Peak")

# Calculate coaccessibility
proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "IterativeLSI_Peak"
)

saveArchRProject(proj)

proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  useMatrix = "GeneExpressionMatrix",
  reducedDims = "IterativeLSI_Peak",
  k = 200
)

# try different k and overlap cutoffs
# try subsetting per organ and take the union of p2g links

# Errors w/ --> now solved by changing the RNA matrices to have a common set of features
#3: stop("Error not all FeatureDF for assay is the same!")
#2: .getFeatureDF(ArrowFiles, useMatrix, threads = threads)
saveArchRProject(proj)

# filter out RNA genes not present in the ArchR annotation (mostly antisense and lncRNA)
p2gtmp <- metadata(proj@peakSet)$Peak2GeneLinks

rna_genes <- metadata(p2gtmp)$geneSet
atac_genes <- getGenes()
genespf <- which(rna_genes$name %in% atac_genes$symbol)
p2gpf <- p2gtmp[p2gtmp$idxRNA %in% genespf,]
metadata(proj@peakSet)$Peak2GeneLinks <- p2gpf


# add some annotation columns for plotting
annot.list <- list()
for (i in seq_along(toKeep$code)){
  organ <- toKeep$code[i]
  meta <- readr::read_tsv(paste0(meta.dir, "/", organ, "_meta.txt"))
  df <- data.frame(CellNames=meta$cb, RNA_Clusters=paste0(organ, "_", meta$L1_clusterID), 
                   RNA_NamedCluster_L2=paste0(organ, "_", meta$L1_clusterID, "_", meta$L2_clusterName))
  annot.list[[organ]] <- df
}
annot <- bind_rows(annot.list)

# Add RNA_Clusters to the all sample project
proj <- addCellColData(
  ArchRProj = proj,
  name = "RNA_NamedCluster_L2",
  cells = annot$CellNames,
  data = paste0(annot$RNA_NamedCluster_L2),
  force = TRUE
)

proj$organ_code <- lapply(proj$RNA_Clusters, function(n){strsplit(n, split="_")[[1]][1]}) %>% unlist 

# plot peak 2 gene
p1 <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "RNA_Clusters", k=10, seed=1)
p2 <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "RNA_Clusters", k=20, seed=1)
p3 <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "RNA_NamedCluster_L2", k=10, seed=1)
p4 <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "organ_code", k=10, seed=1)

#plotDir <- paste0(getOutputDirectory(proj), "/Plots/Peak2Gene/")
plotDir <- paste0(getOutputDirectory(proj), "/Plots/Peak2Gene-k200-rerun/")
dir.create(plotDir, recursive=T, showWarnings=F)

plotPDF(p1,p2,p3, name = file.path("Peak2Gene-k200-rerun", "Peak2Gene-Heatmap-filteredgenes"), width=12, height=12, ArchRProj=proj, addDOC=F)
plotPDF(p4, name = file.path("Peak2Gene-k200-rerun", "Peak2Gene-Heatmap-filteredgenes-byorgan"), width=12, height=12, ArchRProj=proj, addDOC=F)

# GO enrichment of peak2gene kmeans groups ----------------------------------------------------
p2gMat <- plotPeak2GeneHeatmap(proj, groupBy="organ_code", returnMatrices=TRUE, k=10, seed=1)

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
write.table(hrg, file.path(plotDir, sprintf("/top%s_hrg_kclust_k%s_filteredgenes.tsv", nGOgenes, nclust)), quote=FALSE, sep="\t", row.names=FALSE)

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
  groupBy = "organ_code", 
  geneSymbol = hrg[1:20, "gene"], 
  upstream = 100000,
  downstream = 100000,
  loops = getPeak2GeneLinks(proj)
)

plotPDF(plotList = p, 
        name = file.path(plotDir, "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-filteredgenes.pdf"), 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)



