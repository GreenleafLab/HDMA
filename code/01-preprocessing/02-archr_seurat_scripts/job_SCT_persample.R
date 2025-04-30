suppressPackageStartupMessages({
  library(Seurat) 
  library(SeuratObject)
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(ggthemes)
  library(here)
})

scriptPath <- here("code/utils")
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/seurat_helpers.R")) # added some SHAREseq functions in here

# user inputs -------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
iter <- args[1]
obj <- readRDS(args[2])
clustered.plot.dir <- args[3]
blacklist.genes.path <- args[4]

## SCT, cluster, umap per sample, find cluster markers --------------------------------------
# parallelization
plan("multisession", workers = 4)
options(future.globals.maxSize = 2*1000 * 1024^2)
options(future.seed=TRUE)

obj.list <- SplitObject(obj, split.by="Sample") # recreate obj.list from obj to get meta data

sample.names <- names(obj.list)
dir.create("cluster_markers", recursive=T, showWarnings=F)

obj.list <- lapply(seq_along(obj.list), function(i){
  obj <- seuratPreProcess(obj.list[[i]], vst.flavor="v2", dims=1:50, res=0.3, seed=1)
  obj$Clusters <- obj$seurat_clusters
  plotClusterQC(obj, plotDir = clustered.plot.dir, subgroup=sample.names[i], pointSize = 0.1)
  
  cluster.markers <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=1.0, min.diff.pct = 0.2)
  write_tsv(cluster.markers, paste0("cluster_markers/cluster_markers_",iter,"_",sample.names[i],".tsv"))
  return(obj)
})
names(obj.list) <- sample.names
saveRDS(obj.list, paste0("RNA_objlist_preSCT_",iter,".rds"))


## get consensus features -------------
blacklist.genes <- readRDS(blacklist.genes.path)
consFeatures <- getConsensusVarFeatures(obj.list, nfeatures = 3000, blacklist.genes = blacklist.genes)

# define variable features in merged object
VariableFeatures(obj) <- consFeatures

# save
dir.create("consensus_features", recursive=T, showWarnings=F)
saveRDS(consFeatures, paste0("consensus_features/consensusFeatures_",iter,".rds"))
saveRDS(obj, paste0("RNA_obj_preSCT_",iter,".rds"))