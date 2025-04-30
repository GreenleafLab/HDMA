suppressPackageStartupMessages({
  library(Seurat) 
  library(SeuratObject)
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(ggthemes)
  library(BPCells)
  library(celda) 
  library(future)
  library(here)
})
##########################################
# Plot average expression of genes for each chromosome
##########################################
output.dir <- here::here("output/01-preprocessing/02/Brain/rna_preprocess_output/")
setwd(output.dir)

plotDir <- paste0(output.dir, "plots/clustered/SCTdecontX_v4")
obj <- readRDS("RNA_obj_clustered_SCTdecontX_v4.rds")

# log normalize RNA counts
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)

# Get expressed genes
rawCounts <- GetAssayData(object = obj, assay="RNA", slot = "counts")
minUMIs <- 2 # Genes with at least 2 UMIs
minCells <- 5 # Genes expressed in at least 5 cells
expressedGenes <- rownames(rawCounts[rowSums(rawCounts > minUMIs) > minCells,])

require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Hs.eg.db)

# Get ENTREZID for each gene
expressedGenesEntrezID <- select(
  org.Hs.eg.db,
  keys = expressedGenes,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL")

expressedGenesEntrezID <- drop_na(expressedGenesEntrezID)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
geneGR <- GenomicFeatures::genes(txdb)
std_chrm <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
              "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
              "chr21", "chr22", "chrX", "chrY")

geneGR <- geneGR[seqnames(geneGR) %in% std_chrm]

chrByGeneID <- data.frame(chr = seqnames(geneGR), ENTREZID = geneGR$gene_id)
chrByGeneName <- inner_join(chrByGeneID, expressedGenesEntrezID, by = "ENTREZID")

chrByGenes <- list()
chrByGenes <- lapply(std_chrm, function(x) {
  genes <- chrByGeneName %>% filter(chr == x) %>% dplyr::select(SYMBOL) %>% unlist() %>% unname()
})
names(chrByGenes) <- std_chrm

# For each chromosome, append expression for control and trisomy18 for each gene on that chromosome
chrByExpr <- data.frame()
for(ix in seq_along(std_chrm)) {
  features <- chrByGenes[ix] %>% unname() %>% unlist()
  #avgExpr <- AverageExpression(obj, assays = "RNA", features = features, group.by = "Sample")
  avgExpr <- AverageExpression(obj, assays = "RNA", features = features, group.by = "Sample")
  avgExpr <- avgExpr$RNA %>% as.data.frame()
  chr <- std_chrm[ix]
  avgExpr$chr <- chr
  chrByExpr <- bind_rows(chrByExpr, avgExpr)
}

# Pivot longer for plotting purposes
chrByExprPlot <- chrByExpr %>% pivot_longer(!chr, names_to = "Sample", values_to = "normExpr")

pdf(paste0(plotDir, "/chrByExpression_RNAlognorm.pdf"), w = 16, h = 4)
ggplot(chrByExprPlot, 
       aes(x = factor(chr, level = std_chrm), y = normExpr, fill = Sample)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 0.9))+
  scale_y_log10()+
  viridis::scale_fill_viridis(discrete = T, name = "") + theme_bw()
dev.off()

saveRDS(chrByExpr, paste0(plotDir, "/chrByExprPlot_RNAlognorm.rds"))

##########################################
# Plot MBP and BCL2 over the samples
##########################################

pdf(paste0(plotDir, "/ViolinPlot_MBP_BCL2_Expression.pdf"))
VlnPlot(obj, features = c("MBP", "BCL2"), assay = "RNA", group.by = "Sample", log = TRUE)
dev.off()

##########################################
# Plot differential gene expression 
##########################################

# Distribution of log2fc across all genes
diff.genes <- FindAllMarkers(
  obj,
  ident.1 = "Trisomy18",
  ident.2 = "Control",
  group.by = "Sample",
  logfc.threshold = 0,
  min.pct = 0.3
)

pdf(paste0(plotDir, "/hist_diffGenes.pdf"))
hist(diff.genes$avg_log2FC)
dev.off()

# Distribution of differential genes for each chromosome
std_chrm <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
              "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
              "chr21", "chr22", "chrX", "chrY")

diffGenes_df <- data.frame()

for(ix in seq_along(std_chrm)) {
  features <- chrByGenes[ix] %>% unname() %>% unlist()
  diffGenes <- FindMarkers(obj, ident.1 = "Trisomy18", ident.2 = "Control",
                           group.by = "disease", features = features, logfc.threshold = 0, min.pct = 0.1)
  chr <- std_chrm[ix]
  diffGenes$chr <- chr
  diffGenes_df <- bind_rows(diffGenes_df, diffGenes)
}

pdf(paste0(plotDir, "/chrByDiffExpression.pdf"), w = 12, h = 4)
ggplot(diffGenes_df, 
       aes(x = factor(chr, level = std_chrm), y = avg_log2FC)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2)+
  #scale_y_log10()+
  theme_bw()
dev.off()

pdf(paste0(plotDir, "/chrByDiffExpression_padj.pdf"), w = 12, h = 4)
ggplot(diffGenes_df, 
       aes(x = factor(chr, level = std_chrm), y = -log(p_val_adj))) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2)+
  #scale_y_log10()+
  theme_bw()
dev.off()

saveRDS(diffGenes_df, paste0(plotDir, "/chrByDiffExpression.rds"))