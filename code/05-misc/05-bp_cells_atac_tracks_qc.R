# Use BPCells to plot ATAC QC browser tracks efficiently -----------------------------------------------------------------

# imports -----------------------------------------------------------------
library(BPCells)
library(ArchR)
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(here)
})


# user inputs -------------------------------------------------------------
wd <- here::here("output/05-misc/05")
figout <- here::here("figures/05-misc/05")
hh
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
dir.create(figout, recursive = TRUE, showWarnings = FALSE)
setwd(wd)

scriptPath <- here::here("code/utils")
source(paste0(scriptPath, "/track_helpers_BL.R"))
source(paste0(scriptPath, "/track_helpers_SJ.R"))
source(paste0(scriptPath, "/trackplots.R")) # BPCells trackplot helpers in dev

# main ---------------------------------------------------------------------
global_bp_obj <- readRDS(here::here("output/05-misc/01/global_bp_obj.rds"))

global_bp_obj$cell_metadata <- global_bp_obj$cell_metadata %>% 
  dplyr::mutate(clust_name_new=paste0(organ_code, cluster_id))

for (gene in c("ACTB", "GAPDH", "PECAM1", "NEUROD1", "MYH7", "MYOD1", "PTPRC")){
  print(gene)
  pdf(file=paste0(figout, "/qc_plot_", gene, ".pdf"), width=8, height=6)
  for (org in unique(global_bp_obj$cell_metadata$organ)){
  # for (org in c("Adrenal", "Thymus")){ # test code
    submeta <- global_bp_obj$cell_metadata %>% dplyr::filter(organ==org)
    p1 <- trackplot_helper_v2c(
      gene=gene,
      clusters=submeta$clust_name_new, 
      fragments=global_bp_obj$frags %>% select_cells(submeta$cb), 
      cell_read_counts=submeta$nFrags,  
      transcripts=global_bp_obj$transcripts, 
      flank=2e4
    )
    print(p1)
  }
  dev.off()
}

pdf(file=paste0(figout, "/qc_plot_GAPDH_zoomedin.pdf"), width=6, height=9)
for (org in unique(global_bp_obj$cell_metadata$organ)){
  submeta <- global_bp_obj$cell_metadata %>% dplyr::filter(organ==org)
  ordered_name <- factor(submeta$clust_name_new, levels=unique(submeta$clust_name_new)[submeta$cluster_id %>% unique %>% order])
  p1 <- trackplot_helper_v2c(
    region="chr12:6530000-6560000",
    clusters=ordered_name, 
    fragments=global_bp_obj$frags %>% select_cells(submeta$cb), 
    cell_read_counts=submeta$nFrags,  
    transcripts=global_bp_obj$transcripts, 
    flank=2e4
  )
  print(p1)
}
dev.off()

pdf(file=paste0(figout, "/qc_plot_GAPDH_zoomedin_v2.pdf"), width=6, height=8)
for (org in unique(global_bp_obj$cell_metadata$organ)){
  submeta <- global_bp_obj$cell_metadata %>% dplyr::filter(organ==org)
  ordered_name <- factor(submeta$clust_name_new, levels=unique(submeta$clust_name_new)[submeta$cluster_id %>% unique %>% order])
  p1 <- trackplot_helper_v2c(
    region="chr12:6530000-6545000",
    clusters=ordered_name, 
    fragments=global_bp_obj$frags %>% select_cells(submeta$cb), 
    cell_read_counts=submeta$nFrags,  
    transcripts=global_bp_obj$transcripts, 
    flank=2e4
  )
  print(p1)
}
dev.off()
