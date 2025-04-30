suppressPackageStartupMessages({
  library(Seurat) 
  library(SeuratObject)
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(ggthemes)
  library(BPCells)
  library(celda)
  library(here)
})
scriptPath <- here("code/utils/")
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/seurat_helpers.R"))

# user inputs -------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
output.dir <- args[1]
organCode <- args[2]
blacklist.genes <- readRDS(args[3])
known.markers <- args[4]
iter <- args[5]
newiter <- args[6]
if (!is.na(args[7])){ # if arg given
  norm.flag <- as.logical(args[7])
} else{ # default no normalization if no arg given
  norm.flag <- F
}

setwd(output.dir)
obj <- readRDS(paste0("RNA_obj_clustered_",iter,".rds"))
clustered.plot.dir <- file.path(output.dir,"plots/clustered", newiter)
dir.create(clustered.plot.dir, recursive = T, showWarnings = F)

## decontX on SCT counts --------------------------------------------------------
counts <- GetAssayData(object = obj, slot = "counts") # using SCT corrected counts
clusters <- Idents(obj) %>% as.numeric()

# Run on only expressed genes
x <- counts[rowSums(counts)>0,]
message(sprintf("Running decontX on %s cells with %s non-zero genes...", dim(x)[2], dim(x)[1]))
decon <- decontX(x, z=clusters, verbose=TRUE, seed=1) # with cluster prior
saveRDS(decon, paste0("decontX_",newiter,".rds"))
newCounts <- decon$decontXcounts

# Add back unexpressed genes and sort according to original counts
newCounts <- rbind(newCounts, counts[rowSums(counts)==0,])[rownames(counts),]
obj$estConp <- decon$contamination # Estimated 'contamination proportion, 0 to 1'
obj[["decontX"]] <- CreateAssayObject(counts = as(round(newCounts), "sparseMatrix"))

# filter out any cells that have 0 UMI after background removal, they will mess up subsequent analysis
obj <- obj %>% subset(nCount_decontX > 0)

# plot estimated contamination on umap
FeaturePlot(obj, features="estConp", pt.size=1)
ggsave(file.path(clustered.plot.dir,"est_contamination_postdecontX.pdf"), width=7, height=7)

# violin plot of estimated contamination proportion 
cmap <- getColorMap(cmaps_BOR$stallion, n=length(unique(obj$Clusters)))
VlnPlot(obj, features="estConp", group.by="Clusters", pt.size=0, 
        cols=cmap)
ggsave(file.path(clustered.plot.dir, "vlnplot_est_contamination_cluster.pdf"), width=7, heigh=7)


# repeat steps 10-12 from the main script
## find consensus variable features --------------------------------------

# split into obj.list
sample.names <- unique(obj$Sample)
obj.list <- lapply(seq_along(sample.names),function(i){
  obj %>% subset(Sample==sample.names[i])
})
names(obj.list) <- sample.names

obj.list <- lapply(seq_along(obj.list), function(i){
  if (norm.flag){
    obj.list[[i]] %>% NormalizeData %>% FindVariableFeatures(nfeatures=3000)
  }else{
    print("WARNING: No normalization performed on decontX counts per sample!")
    FindVariableFeatures(obj.list[[i]], nfeatures=3000)
  }
  }) # don't do SCT again
names(obj.list) <- sample.names

consFeatures <- getConsensusVarFeatures(obj.list, nfeatures = 3000, blacklist.genes = blacklist.genes)
dir.create("consensus_features", recursive = T, showWarnings = F)
saveRDS(consFeatures, paste0("consensus_features/consensusFeatures_",newiter,".rds"))

# define variable features in merged object
VariableFeatures(obj) <- consFeatures

saveRDS(obj, paste0("RNA_obj_preSCT_",newiter,".rds")) # it's named preSCT to keep the naming consistent with iteration v1, but don't do SCT

## scaling, PCA, clustering ----------------------------------------------------------
DefaultAssay(obj) <- "decontX"
if (norm.flag){
  obj <- obj %>% NormalizeData %>% ScaleData %>% RunPCA(features = consFeatures) %>%
    FindNeighbors(dims=1:50) %>% FindClusters(resolution=0.3, random.seed=1) %>%
    RunUMAP(dims=1:50, min.dist=0.3, n.neighbors=50)
} else{
  print("WARNING: No normalization performed on decontX counts!")
  obj <- obj %>% ScaleData %>% RunPCA(features = consFeatures) %>%
    FindNeighbors(dims=1:50) %>% FindClusters(resolution=0.3, random.seed=1) %>%
    RunUMAP(dims=1:50, min.dist=0.3, n.neighbors=50)
}


obj$Clusters <- obj$seurat_clusters
saveRDS(obj, paste0("RNA_obj_clustered_",newiter,".rds"))

# find markers
cluster.markers <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=1.0, min.diff.pct = 0.2)
dir.create("cluster_markers", recursive = T, showWarnings = F)
write_tsv(cluster.markers, paste0("cluster_markers/cluster_markers_",newiter,".tsv"))


## plot cluster QC -------------------------------------------------------------------
plotClusterQC(obj, plotDir = clustered.plot.dir, subgroup=organCode, pointSize = 0.1)

# plot known marker genes
source(known.markers)
feature.list <- c(markers.endothelial, markers.epithelial, markers.immune, markers.stromal)
FeaturePlot(obj, slot = "data", reduction="umap",
            features = feature.list, raster = T, raster.dpi = c(1024,1024))
ggsave(paste0(clustered.plot.dir,"/compartment_features_UMAP_",organCode,".pdf"),
       width=20, height = 15)

DotPlot(obj, features=feature.list, group.by="Clusters")
ggsave(paste0(clustered.plot.dir,"/compartment_features_dotplot_",organCode,".pdf"), 
       width=2+length(feature.list)*0.8, height=8)

# plot sex genes per sample again
sex.genes <- c("XIST", "SRY", "UTY","TTTY10", "TTTY14")
DotPlot(object = obj, group.by="Sample", features = sex.genes)
ggsave(file.path(clustered.plot.dir,paste0("sample_sex_",newiter,".pdf")), width=10, height=5)

# plot organ specific feature sets
featureSets.valid <- featureSets[featureSets %in% rownames(obj)]
pdf(paste0(clustered.plot.dir,"/organ_features_UMAP_",organCode,".pdf"), height=7, width=7)
for (feat in featureSets.valid){
  print(FeaturePlot(obj, slot = "data", reduction="umap",
                    features = feat, raster = T, raster.dpi = c(1024,1024)))
}
dev.off()

DotPlot(obj, features=featureSets.valid, group.by="Clusters") + coord_flip()
ggsave(paste0(clustered.plot.dir,"/organ_features_dotplot_",organCode,".pdf"), 
       width=4+length(unique(obj$Clusters)) * 0.4, height=length(featureSets.valid)*0.3)


## cluster filtering plotting ------------------------------------------------------------------
# plot sample using decontX embedding, but color by per sample clusters
obj.list <- readRDS(paste0("RNA_objlist_preSCT_",iter,".rds")) # per sample SCT transform and clustering
for (sample in obj$Sample %>% unique){
  subobj <- obj[,colnames(obj.list[[sample]])]
  umapDF <- data.frame(Embeddings(object=subobj, reduction="umap"), 
                       persample_cluster=obj.list[[sample]]@meta.data[colnames(subobj),"Clusters"])
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]
  pdf(paste0(clustered.plot.dir, sprintf("/cluster_UMAP_colbypersampclust_%s.pdf", sample)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=cmaps_BOR$stallion, point_size=0.1))
  dev.off()
}
