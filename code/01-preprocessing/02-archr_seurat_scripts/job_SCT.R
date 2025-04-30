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
# iter <- "SCTdecontX_v1" # this is used to name outputs

args <- commandArgs(trailingOnly = T)
iter <- args[1]
obj <- readRDS(args[2])
output.dir <- args[3]
setwd(output.dir)

#optional args for plotting
if (length(args)==6){
  known.markers <- args[4]
  clustered.plot.dir <- args[5]
  organCode <- args[6]
}

## scaling, PCA, clustering ----------------------------------------------------------
# parallelization
# library(future)
# plan("multisession", workers = 4)
# options(future.globals.maxSize = 30*1000 * 1024^2)
# options(future.seed=TRUE)

# this takes a lot of memory
obj <- SCTransform(obj, vst.flavor="v2") %>% RunPCA(features = VariableFeatures(obj)) %>%
  FindNeighbors(dims=1:50) %>% FindClusters(resolution=0.3, random.seed=1) %>%
  RunUMAP(dims=1:50, min.dist=0.3, n.neighbors=50)

obj$Clusters <- obj$seurat_clusters
saveRDS(obj, paste0("RNA_obj_clustered_",iter,".rds"))

# find markers
cluster.markers <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=1.0, min.diff.pct = 0.2)
write_tsv(cluster.markers, paste0("cluster_markers/cluster_markers_",iter,".tsv"))

## (optional) plot cluster QC -------------------------------------------------------------------
if (length(args)==6){
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
  ggsave(file.path(clustered.plot.dir,paste0("sample_sex_",iter,".pdf")), width=10, height=5)
  
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
}