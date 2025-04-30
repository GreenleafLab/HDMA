# Preprocess scRNA using Seurat -----------------------------------------------------------------

# imports -----------------------------------------------------------------
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
input.dir <- here("output/01-preprocessing/01/##ORGAN##/")
output.dir <- here("output/01-preprocessing/02/##ORGAN##/rna_preprocess_output")
code.dir <- here("code/01-preprocessing/02-archr_seurat_scripts")
dir.create(output.dir, recursive = T, showWarnings = F)
setwd(output.dir)

organCode <- "##"
whitelist.dir <- here("output/01-preprocessing/02/shared/whitelist_r1")
blacklist.genes.path <- here("data/external/ENCODE_blacklist/blacklist.genes.rds")
meta.path <- here("output/01-preprocessing/02/shared/SampleMetadata.csv")
cutoffs.path <- here("output/01-preprocessing/02/shared/sample_filter_cutoffs_metadata.csv")
known.markers <- paste0(here("output/01-preprocessing/02/shared/markergenes/markergenes_"), organCode, ".R")

# post RNA preprocess final cell list (not sample meta data)
dir_meta <- here("output/01-preprocessing/02/shared/meta/")

# main --------------------------------------------------------------------
# each section is designed to be the minimal "automated" run unit with no pausing needed

## 1. read data, create seurat v4 object ----------------------------------
plot.dir <- file.path(output.dir,"plots/raw_qc")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

obj.list <- read_rna_data(input.dir, min.umi=300)
saveRDS(object = obj.list,
        file = "RNA_objlist_raw.rds")
#obj.list <- readRDS("RNA_objlist_raw.rds")

# plot QC on raw data
plot_rna_qc(obj.list, plot.dir, umi.cutoffs=1e3, gene.cutoffs=500, pct.mt.cutoffs=20)

## 2. determine cutoffs --------------------------------------------------------
#    manually inspect the output, fill in the cutoff meta file, then import cutoffs
print(cutoffs.path)
cutoffs <- read.csv(cutoffs.path, header=T, row.names=1)

# try cutoffs
sample.names <- names(obj.list)
tmp.obj.list <- lapply(seq_along(obj.list), function(i){
  obj <- obj.list[[i]]
  sample <- names(obj.list)[i]
  obj <- subset(obj, nCount_RNA>cutoffs[sample,"rna_nUMIs"] &
                  nFeature_RNA>cutoffs[sample,"rna_nGenes"] &
                  percent.mt<cutoffs[sample,"rna_pctMT"])
  return(obj)
})
names(tmp.obj.list) <- sample.names

# replot QC on filtered data
#    manually inspect plots, repeat step 2 as needed
filtered.plot.dir <- file.path(output.dir,"plots/filtered")
dir.create(filtered.plot.dir, recursive = TRUE, showWarnings = FALSE)
plot_rna_qc(tmp.obj.list, filtered.plot.dir, 
        umi.cutoffs=cutoffs[sample.names,"rna_nUMIs"], 
        gene.cutoff=cutoffs[sample.names,"rna_nGenes"], 
        pct.mt.cutoffs=cutoffs[sample.names,"rna_pctMT"]) 

## 3. save filtered object ------------------------------------------------------
obj.list <- tmp.obj.list
rm(tmp.obj.list)
saveRDS(object = obj.list, file = "RNA_objlist_filtered.rds")
#obj.list <- readRDS("RNA_objlist_filtered.rds")

rna.cells.passfilter <- lapply(obj.list, function(obj){colnames(obj)}) %>% unlist
umis <- lapply(obj.list, function(obj){obj$nCount_RNA}) %>% unlist
towrite <- data.frame(rna_cb=rna.cells.passfilter, rna_nUMIs=umis)
rownames(towrite) <- NULL
write.table(towrite, sep="\t",paste0(whitelist.dir,organCode,"_rna.txt"), 
            quote=F, row.names=F, col.names=T)

## 4. merge RNA and ATAC whitelist ------------------------------------------------
rna.cells.passfilter <- read.csv(paste0(whitelist.dir, organCode,"_rna.txt"), sep="\t")
atac.cells.passfilter <- read.csv(paste0(whitelist.dir, organCode,"_atac.txt"), sep="\t")
cells.passfilter <- merge(rna.cells.passfilter, atac.cells.passfilter, by.x="rna_cb", by.y="atac_cb") %>%
  dplyr::rename(cb=rna_cb) %>% mutate(nFrags_nUMIs_ratio=atac_nFrags/rna_nUMIs) 
rownames(cells.passfilter) <- cells.passfilter$cb
# plot nFrags vs nUMI plot per sample
whitelist.plot.dir <- file.path(output.dir,"plots/whitelist")
dir.create(whitelist.plot.dir, recursive = TRUE, showWarnings = FALSE)
plot_nfrag_numi(cells.passfilter, whitelist.plot.dir)

## 5. determine joint cutoffs -----------------------------------------------------
# manually inspect the output
# if needed, fill in the nFrags_nUMIs_ratio in cutoff meta file, then import cutoffs
print(cutoffs.path)
cutoffs <- read.csv(cutoffs.path, header=T, row.names=1)

# try cutoffs, manually inspect plots and repeat step 5 as needed
tmp <- cells.passfilter
tmp$sample <- lapply(cells.passfilter$cb, function(n){strsplit(n, split="#")[[1]][1]})
tmp <- lapply(unique(tmp$sample), function(n){
  sub <- tmp[tmp$sample==n,]
  cut <- cutoffs[n, "nFrags_nUMIs_ratio"]
  return(sub[sub$nFrags_nUMIs_ratio > cut,])
}) %>% do.call(rbind, .)

whitelist.filtered.plot.dir <- file.path(output.dir,"plots/whitelist_filtered")
dir.create(whitelist.filtered.plot.dir, recursive = TRUE, showWarnings = FALSE)
plot_nfrag_numi(tmp, whitelist.filtered.plot.dir)

## 6. save filtered whitelist -----------------------------------------------------
cells.passfilter <- cells.passfilter[tmp$cb,]
rm(tmp)
write.table(cells.passfilter, sep="\t", paste0(whitelist.dir, organCode, "_both.txt"), 
            quote=F, row.names=F)

## 7. read in cell names whitelisted from both ATAC and RNA ------------------------
cells.passfilter <- read.csv(paste0(whitelist.dir, organCode, "_both.txt"), sep="\t")
#obj.list <- readRDS("RNA_objlist_filtered.rds") 
sample.names <- names(obj.list)
obj.list <- lapply(seq_along(obj.list), function(i){
    obj <- obj.list[[i]]
    obj$passfilter <- colnames(obj) %in% cells.passfilter$cb
    return(subset(obj, subset=(passfilter==TRUE), return.null=T))
})
names(obj.list) <- sample.names

## 8. merge into a single object ----------------------------------------------------
obj <- merge(x=obj.list[[1]], y=obj.list[2:length(obj.list)], project=organCode)
obj$Sample <- lapply(rownames(obj@meta.data), function(n){strsplit(n, split="#")[[1]][1]}) %>% unlist

# plot sex genes per sample
sex.genes <- c("XIST", "SRY", "UTY","TTTY10", "TTTY14")
DotPlot(object = obj, group.by="Sample", features = sex.genes)
ggsave(file.path(whitelist.filtered.plot.dir,"sample_sex.pdf"), width=10, height=5)

## 9. add meta data ---------------------------------------------------------
# manually assign sex then import meta file
print(meta.path)
meta <- read.csv(meta.path)
rownames(meta) <- meta$SampleName

obj$sex <- meta[obj$Sample, "Sex"] 
obj$PCW <- paste0("PCW", meta[obj$Sample, "PCW"])
obj$PCD <- paste0("PCD", meta[obj$Sample, "PCD"])
obj$age <- obj$PCW

# add cell cycle information
obj <- CellCycleScoring(obj, s.features=cc.genes$s.genes, 
                        g2m.features=cc.genes$g2m.genes, set.ident=FALSE)


## 10. find consensus variable features --------------------------------------
iter <- "SCTdecontX_v1" # this is used to name outputs
clustered.plot.dir <- file.path(output.dir,"plots/clustered", iter)
dir.create(clustered.plot.dir, recursive = T, showWarnings = F)
dir.create("slurm_logs", recursive = T, showWarnings = F)

# submit as batch script, SCT -> cluster -> umap -> find marker genes per sample
# can submit the sbatch commands from blocks 10, 11 and 13 together, they will run sequentially based on job dependency 
saveRDS(obj, "tmp_RNA_obj_preSCTpersample.rds", compress = F)
command <- sprintf("sbatch -p sfgf,biochem,wjg --mem-per-cpu=200g -n 1 --time=11:00:00 --job-name=%s_SCTpersample --output=slurm_logs/slurm-%%j.out --wrap \"Rscript %s/job_SCT_persample.R %s %s %s %s\"",
                   organCode, code.dir, iter, "tmp_RNA_obj_preSCTpersample.rds", clustered.plot.dir, blacklist.genes.path)
out <- system(command, wait=F, intern=T)
presct_jobid <- str_split(out, "batch job ")[[1]][2] # get the job id

## 11. scaling, PCA, clustering ----------------------------------------------------------
# use slurm command, took ~4 hours and 270GB memory for 250k cells
command <- sprintf("sbatch --dependency=afterok:%s -p sfgf,biochem,wjg --mem-per-cpu=200g -n 1 --time=5:00:00 --job-name=%s_SCT --output=slurm_logs/slurm-%%j.out --wrap \"Rscript %s/job_SCT.R %s %s %s %s %s %s\"",
                   presct_jobid, organCode, code.dir, iter, paste0("RNA_obj_preSCT_",iter,".rds"), output.dir, known.markers, clustered.plot.dir, organCode)
out <- system(command, wait=F, intern=T)
sct_jobid <- str_split(out, "batch job ")[[1]][2] # get the job id

## 13. decontX on SCT counts --------------------------------------------------------
# slurm command, took ~10 hours and 80GB memory for 250k cells
# this script performs decontX, then repeats steps 10-12
newiter <- "SCTdecontXnorm_v2"
command <- sprintf("sbatch --dependency=afterok:%s -p sfgf,biochem,wjg --mem-per-cpu=100g -n 1 --time=11:00:00 --job-name=%s_decontX --output=slurm_logs/slurm-%%j.out --wrap \"Rscript %s/job_decontX.R %s %s %s %s %s %s T\"",
                   sct_jobid, organCode, code.dir, output.dir, organCode, blacklist.genes.path, known.markers, iter, newiter)
system(command, wait=F)

# manually check job status and proceed to next step once done

## 14. cluster filtering ------------------------------------------------------------------
# manually remove doublet-like clusters, then repeat steps 10-13
# for reproducibility, record notes on each cluster in a file cluster_annotations.csv 
# with columns cluster,celltype,celltype2,keep,note
# an example line in cluster_annotations.csv: 1,LU_end,vascular endothelial,T,PECAM1 marker gene
previter <- "SCTdecontX_v1" # SCT only no decontX iteration
newiter <- "SCTdecontXnorm_v2" # post decontX iteration
clustered.filtered.plot.dir <- paste0(output.dir,"/plots/clustered/",newiter, "_filtered")
dir.create(clustered.filtered.plot.dir, recursive=T, showWarnings=F)

obj <- readRDS(paste0("RNA_obj_clustered_",newiter,".rds"))

annot <- paste0("plots/clustered/SCTdecontX_v2/cluster_annotations.csv") %>% read.csv
rownames(annot) <- annot$cluster
obj$Clusters <- as.character(obj$seurat_clusters) # check this to make sure you are using the correct cluster IDs
obj$annotv1 <- annot[obj$Clusters,"celltype"]
obj$keepcell <- annot[obj$Clusters,"keep"] 
obj$Clusters <- obj$annotv1
plotClusterQC(obj, plotDir = paste0(clustered.filtered.plot.dir), subgroup=organCode, pointSize = 0.1)

obj$annotv2 <- annot[as.character(obj$seurat_clusters),"celltype2"] # this is more granular annotation
obj$Clusters <- obj$annotv2
plotClusterQC(obj, plotDir = paste0(clustered.filtered.plot.dir), subgroup=paste0(organCode,"_celltype2"), pointSize = 0.1)

## 14b. prepare object for next iteration, skip if no more iterations needed ------------
cellpf <- obj[,!is.na(obj$annotv1) & (obj$keepcell==TRUE)] %>% colnames
obj <- readRDS("tmp_RNA_obj_preSCTpersample.rds")
obj <- obj[,cellpf]

## 15. export round 2 cell white list -----------------------------------------------------
annot <- paste0("/path/to/your/cluster_annotations.csv") %>% read.csv

for (cluster in unique(annot$celltype)){
  sub <- annot[annot$celltype==cluster,]
  if (nrow(sub)>1){suffix <- 1:nrow(sub)}else{suffix <- ""}
  annot[rownames(sub),"L2_clusterID"] <- paste0(cluster,suffix)
}

for (cluster in unique(annot$celltype2)){
  sub <- annot[annot$celltype2==cluster,]
  if (nrow(sub)>1){suffix <- paste0(" ", 1:nrow(sub))}else{suffix <- ""}
  annot[rownames(sub),"L3_clusterName"] <- paste0(cluster,suffix)
}


obj$Clusters <- obj$seurat_clusters # make sure this is the cluster numbers, will be used to annotate ATAC
saveRDS(obj, "RNA_obj_clustered_final.rds")
towrite <- data.frame(cb=colnames(obj), 
                      L1_clusterID=obj$Clusters,
                      L1_clusterName=annot[obj$Clusters,"celltype"], # first level annotations
                      L2_clusterID=annot[obj$Clusters,"L2_clusterID"],
                      L2_clusterName=annot[obj$Clusters,"celltype2"], # second level annotations
                      L3_clusterName=annot[obj$Clusters,"L3_clusterName"]) # third level annotations

rownames(towrite) <- NULL

write.table(towrite, sep="\t", paste0(dir_meta, organCode, "_meta.txt"), quote=F, row.names=F, col.names=T)


# plot final marker set dotplot
source(known.markers) # contains a variable called markerSets
DotPlot(obj, features=rev(markerSets), group.by="Clusters") + coord_flip()
ggsave(paste0(clustered.filtered.plot.dir,"/marker_features_dotplot_",organCode,".pdf"), 
       width=2+length(unique(obj$Clusters)) * 0.4, height=length(markerSets)*0.3)

# reorder celltype3
#rownames(annot) <- annot$cluster
#obj$annotv3 <- annot[as.character(obj$seurat_clusters),"L3_clusterName"]
#obj$annotv3 <- factor(obj$annotv3, levels=c("")) 

DotPlot(obj, features=rev(markerSets), group.by="annotv3") + coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(clustered.filtered.plot.dir,"/marker_features_dotplot_",organCode,"_celltype3.pdf"), 
       width=2+length(unique(obj$Clusters)) * 0.4, height=2+length(markerSets)*0.3)

# a few more QC plots based on new cluster names
source(paste0(scriptPath, "/hdma_palettes.R"))
obj$age <- factor(obj$PCW, levels=names(cmap_pcw))
obj$annotv4 <- factor(paste0(organCode, as.character(obj$seurat_clusters)),
                      levels=paste0(organCode, levels(obj$seurat_clusters)))
obj$Clusters <- obj$annotv4
plotClusterQC(obj, plotDir = paste0(clustered.filtered.plot.dir), subgroup=paste0(organCode,"_celltype4"), pointSize = 0.1,
              gest_age_cmap = cmap_pcw)


obj$annotv5 <- paste0(obj$annotv4, ".", gsub(" ", "_", obj$annotv3))
dict <- unique(obj@meta.data[c("annotv3", "annotv4", "annotv5")])
rownames(dict) <- dict$annotv4
obj$annotv5 <- factor(obj$annotv5, levels=dict[levels(obj$annotv4), "annotv5"])
obj$Clusters <- obj$annotv5
plotClusterQC(obj, plotDir = paste0(clustered.filtered.plot.dir), subgroup=paste0(organCode,"_celltype5"), pointSize = 0.1,
              gest_age_cmap = cmap_pcw)


## rearrange the factor levels based on marker gene order
rownames(dict) <- dict$annotv3
obj$annotv4 <- factor(obj$annotv4, levels=dict[levels(obj$annotv3), "annotv4"])
DotPlot(obj, features=rev(markerSets), group.by="annotv4") + coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(clustered.filtered.plot.dir,"/marker_features_dotplot_",organCode,"_celltype4.pdf"), 
       width=2+length(unique(obj$Clusters)) * 0.4, height=2+length(markerSets)*0.3)

rownames(dict) <- dict$annotv3
obj$annotv5 <- factor(obj$annotv5, levels=dict[levels(obj$annotv3), "annotv5"])
DotPlot(obj, features=rev(markerSets), group.by="annotv5") + coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(clustered.filtered.plot.dir,"/marker_features_dotplot_",organCode,"_celltype5.pdf"), 
       width=2+length(unique(obj$Clusters)) * 0.4, height=2+length(markerSets)*0.3)

# one seurat object "*final*"

# "LU_meta.txt"

# cellbarcode, L1_clusterID, L1_clusterName, L2_clusterID, L2_clusterName, L3_clusterName
# T361_b11_Heart_PCW20#CL73_H01+A01+B09, 1, LU_end, LU_end1, vascular endothelial, vascular endothelial 1

## 15b. export cluster and sample level meta data -----------------------------------------------------
# proceed with this once done with ATAC peak calling, check the inputs
library(ArchR)
atac.proj.dir <- here("output/01-preprocessing/02/##ORGAN##/atac_preprocess_output/ATAC_obj_clustered_peaks_final/")
organ.full.name <- "##ORGAN##"

# cluster level meta
# Organ, organ code, cluster, celltype1-3, note, ncells per cluster, ngene, numi, pctmt, nfrags, TSS, FRiP
atac <- loadArchRProject(atac.proj.dir)
obj$frip <- getCellColData(atac)[colnames(obj),"FRIP"]
obj$nfrags <- getCellColData(atac)[colnames(obj),"nFrags"]
obj$tss <- getCellColData(atac)[colnames(obj),"TSSEnrichment"]

meta_df <- obj@meta.data %>% group_by(seurat_clusters) %>% 
  dplyr::summarise(ncell=n(),
                   ngene=median(nFeature_RNA),
                   numi=median(nCount_RNA),
                   pctmt=median(percent.mt),
                   nfrags=median(nfrags),
                   tss=median(tss),
                   frip=median(frip)
  ) %>% dplyr::rename(cluster=seurat_clusters)
annot <- merge(meta_df,annot,by="cluster")
cmeta <- data.frame(organ=organ.full.name,
                    organ_code=organCode,
                    L1_clusterID=annot$cluster, # e.g. 0
                    L0_clusterName=lapply(annot$celltype, function(n){strsplit(n, split="_")[[1]][2]}) %>% unlist, # e.g. epi
                    L1_clusterName=annot$celltype, # e.g. BR_epi
                    L2_clusterID=annot$L2_clusterID, # e.g. BR_epi1
                    L2_clusterName=annot$celltype2, # e.g. excitatory neurons
                    L3_clusterName=annot$L3_clusterName, # e.g. excitatory neurons1
                    ncell=annot$ncell,
                    median_numi=annot$numi,
                    median_ngene=annot$ngene,
                    median_pctmt=annot$pctmt,
                    median_nfrags=annot$nfrags,
                    median_tss=annot$tss,
                    median_frip=annot$frip,
                    note=annot$note
)

write.table(cmeta, sep="\t",paste0(dir_meta, organCode,"_meta_cluster.txt"), 
            quote=F, row.names=F, col.names=T)

# sample level meta
#Organ, organ code, sample id, pcw, pcd, sex, median ngene, numi, pctmt, nfrags, TSS, FRiP
smeta <- obj@meta.data %>% group_by(Sample) %>% 
  dplyr::summarise(organ=organ.full.name,
                   organ_code=organCode,
                   pcw=strsplit(unique(PCW), split="W")[[1]][2],
                   pcd=strsplit(unique(PCD), split="D")[[1]][2],
                   sex=unique(sex),
                   ncell=n(),
                   median_ngene=median(nFeature_RNA),
                   median_numi=median(nCount_RNA),
                   median_pctmt=median(percent.mt),
                   median_nfrags=median(nfrags),
                   median_tss=median(tss),
                   median_frip=median(frip)
  ) %>% dplyr::rename(sample=Sample)
write.table(smeta, sep="\t",paste0(dir_meta, organCode,"_meta_sample.txt"), 
            quote=F, row.names=F, col.names=T)

## (optional) 12. plot cluster QC -------------------------------------------------------------------
# plotting script in block 12 is now automatically run in job_SCT.R
# keeping these code here for quick reference and custom plotting if needed
obj <- readRDS(paste0("RNA_obj_clustered_",iter,".rds"))
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



## temp. joint embedding -----------------------------------------------------
# inputs
library(ArchR)
atac.proj.dir <- here("output/01-preprocessing/02/##ORGAN##/atac_preprocess_output/ATAC_obj_clustered_peaks_final/")
rna.proj.rds <- paste0(output.dir, "/RNA_obj_clustered_final.rds")

obj <- WNN(atac.proj.dir, rna.proj.rds)
saveRDS(obj, "RNA_obj_clustered_final_wnn.rds")

# plotting
plot.dir <- "plots/clustered/SCTdecontX_ITERATION_filtered/"
annot <- read_tsv(paste0(dir_meta, organCode,"_meta.txt")) %>% as.data.frame
rownames(annot) <- annot$cb
obj$celltype2 <- annot[colnames(obj),"L2_clusterName"]
obj$celltype <- annot[colnames(obj),"L1_clusterName"]

p1 <- ggboxplot(obj@meta.data, y="RNA.weight")
p2 <- ggboxplot(obj@meta.data, y="ATAC.weight")

p3 <- quickumap(obj, "wnn.umap", "Clusters", "WNN ATAC+RNA")
p4 <- quickumap(obj, "umap", "Clusters", "RNA only")
p5 <- quickumap(obj, "ATACUMAP", "Clusters", "ATAC only")

p6 <- quickumap(obj, "wnn.umap", "celltype", "WNN ATAC+RNA")
p7 <- quickumap(obj, "umap", "celltype", "RNA only")
p8 <- quickumap(obj, "ATACUMAP", "celltype", "ATAC only")

p9 <- quickumap(obj, "wnn.umap", "celltype2", "WNN ATAC+RNA")
p10 <- quickumap(obj, "umap", "celltype2", "RNA only")
p11 <- quickumap(obj, "ATACUMAP", "celltype2", "ATAC only")

p12 <- quickumap(obj, "wnn.umap", "Sample", "WNN ATAC+RNA")
p13 <- quickumap(obj, "umap", "Sample", "RNA only")
p14 <- quickumap(obj, "ATACUMAP", "Sample", "ATAC only")

pdf(paste0(plot.dir, "/cluster_UMAP_jointwnn_", organCode,".pdf"), width=7, height=7)
lapply(paste0("p", 1:14), function(n){get(n)})
dev.off()

p15 <- ggboxplot(obj@meta.data, y="RNA.weight", x="Sample", color="Sample") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p16 <- ggboxplot(obj@meta.data, y="RNA.weight", x="Clusters", color="Clusters") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p17 <- ggboxplot(obj@meta.data, y="ATAC.weight", x="Sample", color="Sample") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p18 <- ggboxplot(obj@meta.data, y="ATAC.weight", x="Clusters", color="Clusters") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf(paste0(plot.dir, "/cluster_UMAP_jointwnn_RNA_weights", organCode,".pdf"), width=7, height=7)
lapply(paste0("p", 15:18), function(n){get(n)})
dev.off()
