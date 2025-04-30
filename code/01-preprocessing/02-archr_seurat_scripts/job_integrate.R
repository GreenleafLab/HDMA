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

scriptPath <- here("code/utils")
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))
source(paste0(scriptPath, "/seurat_helpers.R")) # added some SHAREseq functions in here

# user inputs -------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
iter <- args[1] # name for new iteration
obj.list <- readRDS(args[2]) # rds file path to the individually SCT normalized obj list
output.dir <- args[3] # work directory
blacklist.genes.path <- args[4] # list of black list genes to exclude from consensus features

dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
setwd(output.dir)

## 11b. try integration (brain only) ----------------------
# too many sample specific cluster that are cross compartment

blacklist.genes <- readRDS(blacklist.genes.path)
consFeatures <- getConsensusVarFeatures(obj.list, nfeatures = 3000, blacklist.genes = blacklist.genes)

# run PCA on each sample using consensus features
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = consFeatures)
obj.list <- lapply(X = obj.list, FUN = RunPCA, features = consFeatures)

int.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                      anchor.features = consFeatures, dims = 1:30, reduction = "rpca", k.anchor = 20)
obj.combined.sct <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", dims = 1:30)
obj.combined.sct <- RunPCA(obj.combined.sct, verbose = FALSE)
obj.combined.sct <- RunUMAP(obj.combined.sct, reduction = "pca", dims = 1:30)
obj.combined.sct <- FindNeighbors(obj.combined.sct, reduction = "pca", dims = 1:30)
obj.combined.sct <- FindClusters(obj.combined.sct, resolution = 0.5)

saveRDS(obj.combined.sct,  paste0("RNA_obj_clustered_",iter,".rds"))
