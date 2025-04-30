#!/share/software/user/open/R/4.1.2/bin/Rscript

# Merge gex matrix feature set across organs --------------------------------------------------

# imports -----------------------------------------------------------------------------------------
#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(ggrepel)
  library(rhdf5)
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
wd <- here("output/01-preprocessing/03")
input.dir <- here("output/01-preprocessing/02/")
global.gexfeat.file <- paste0(wd, "/global_gexfeat_decontx_gr.rds")

# set working directory
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
pointSize <- 0.5
addArchRGenome("hg38")
addArchRThreads(8)

# rna gtf file used in RNA feature assignment
gtf.file <- here("data/reference/gtf/gencode.v41.annotation.BPfiltered.gtf")

# Manually define which organs to keep
toKeep <- read_tsv(here("code/02-global_analysis/01-organs_keep.tsv"))
organs <- toKeep$organs

# main ----------------------------------------------------------------------------------------------------

## 1. Create a common gex feature set from the global ArchR arrow files -----------------------------------------------
# these are genes considered expressed in at least 1 organ, excluding chrM and chrY
if (!file.exists(global.gexfeat.file)){
  all_arrow <- list()
  for (organ in organs){
    atac_proj <- loadArchRProject(file.path(input.dir, organ, "atac_preprocess_output", "ATAC_obj_clustered_peaks_final_decontx"))
    all_arrow[[organ]] <- getArrowFiles(atac_proj)
  }
  
  all_arrow <- all_arrow %>% unlist %>% unname
  
  # take a look at the number of features in each arrow file
  n_gexfeat <- lapply(1:length(all_arrow), function(i){
    ArchR:::.getFeatureDF(all_arrow[i], "GeneExpressionMatrix")$name %>% length
  }) %>% unlist
  # > n_gexfeat
  # [1] 25263 25263 25263 25263 32431 32431 32431 32431 32431 32431 32431 32431 32000 32000 32000 32000 32000 32000 32000 31980 31980 31980 31980 31980 31980 31980
  # [27] 31980 32673 32673 32673 32673 32673 32673 32673 33395 33395 33395 33395 33395 33395 33395 33395 33395 33395 33395 33395 33395 33395 33395 32793 32793 32793
  # [53] 32793 32793 32793 32793 32793 34122 34122 34122 32442 32442 32442 32965 32965 32965 32965 32965 32965 32965 33576 33576 33576 26115 26115 26115
  
  
  global_gexfeat <- lapply(1:length(all_arrow), function(i){
    ArchR:::.getFeatureDF(all_arrow[i], "GeneExpressionMatrix")$name
  }) %>% unlist %>% unique
  # > global_gexfeat %>% length
  # [1] 34788
  
  # get ranges for global gex genes
  # read RNA annotation
  message(sprintf("reading RNA annotations from %s", gtf.file))
  rna.annot <- rtracklayer::import(gtf.file)
  rna.annot <- rna.annot[rna.annot$type=="gene"]
  names(rna.annot) <- rna.annot$gene_name
  is.gene.related <- grep("gene_", colnames(mcols(rna.annot)))
  mcols(rna.annot) <- mcols(rna.annot)[,is.gene.related]
  
  # read genes from ArchR annotation
  atac.annot<- getArchRGenome(geneAnnotation = T)$gene
  names(atac.annot) <- atac.annot$symbol
  
  # for genes present in the ATAC annotation, use ATAC annotation genome ranges
  inatacmask <- global_gexfeat %in% names(atac.annot)
  message(sprintf("# genes in global gene set and in ArchR annotation: %s", sum(inatacmask)))
  # output: # genes in global gene set and in ArchR annotation: 20278
  ranges1 <- atac.annot[global_gexfeat[inatacmask]]
  
  # for genes not present in the ATAC annotation, use RNA annotation genome ranges
  # e.g. these are mostly lncRNAs, antisense genes
  message(sprintf("# genes in global gene set but not in ArchR annotation: %s", sum(!inatacmask)))
  message(paste0("Examples of genes not in ArchR annotation: "))
  message(paste0(global_gexfeat[!inatacmask] %>% head(20), ","))
  ranges2 <- rna.annot[global_gexfeat[!inatacmask]]
  
  # Rename/format/add some meta columns for consistency before merging
  genome(ranges2) <- "hg38"
  ranges2$symbol <- ranges2$gene_name
  ranges2$gene_name <- NULL
  ranges2$ensembl_gene_id <- ranges2$gene_id
  ranges2$gene_id <- NULL
  
  ranges1$gene_type <- rna.annot[names(ranges1)]$gene_type
  ranges1$ensembl_gene_id <- rna.annot[names(ranges1)]$gene_id
  
  global_gexfeat_ranges <- c(ranges1, ranges2) %>% sort
  
  saveRDS(global_gexfeat_ranges, global.gexfeat.file)
} else {
  message("global gex feature granges already exists")
}

## 2a. Add gene expression matrix to ArchRProject per organ -------------------------------------------------------------------

add_global_feat_wrapper <- function(organ, input.dir, global.gexfeat.file){
  message(organ)
  # Grab atac and rna paths
  atac_wd <- paste0(input.dir, "/", organ, "/atac_preprocess_output/")
  rna_wd <- paste0(input.dir, "/", organ, "/rna_preprocess_output/")
  atac_proj <- loadArchRProject(paste0(atac_wd, "/ATAC_obj_clustered_final_decontx/"))
  rna_proj <- paste0(rna_wd, "/RNA_obj_clustered_final.rds") %>% readRDS
  
  # save a copy of the atac proj to a new folder
  atac_proj_out <- paste0(atac_wd, "/ATAC_obj_clustered_final_common_gexfeat_decontx")
  atac_proj <- saveArchRProject(atac_proj, atac_proj_out, load=T)

  atac_proj <- addSeuratGeneExpressionMatrixFixedFeat(atac_proj, global.gexfeat.file, rna_proj, assay="decontX")
  
  saveArchRProject(atac_proj, dropCells = TRUE)
}

## exclude organs that already ran successfully
#organs <- organs[!organs %in% c("Adrenal", "Brain", "Eye")]


for (organ in organs){
  add_global_feat_wrapper(organ, input.dir, global.gexfeat.file)
}

## 2b. Replace GEX matrix in the global ArchR project Arrow Files ------------------------------------------------------------------------------------

replace_arrow_gex_mtx_wrapper <- function(input.dir, organ){
  atac_wd <- paste0(input.dir, "/", organ, "/atac_preprocess_output/")
  arrow_path <- paste0(atac_wd, "/ATAC_obj_clustered_final_common_gexfeat_decontx/ArrowFiles")
  arrow_files <- list.files(path = arrow_path, pattern = "\\.arrow$", full.names = TRUE)
  
  for (f in arrow_files){
    message(paste0("processing ", basename(f)))
    source_file <- f
    destination_file <- paste0(wd, "/allSamples_decontx/ArrowFiles/", basename(f))
    stopifnot("Destination file does not exist"= file.exists(destination_file))
    
    # Specify the paths to the source and destination groups
    source_group_path <- "GeneExpressionMatrix"
    destination_group_path <- "GeneExpressionMatrix"
    
    # Open the source HDF5 file in read-only mode
    source_h5file <- H5Fopen(source_file, "H5F_ACC_RDONLY")
    source_group_data <- h5read(source_h5file, name = source_group_path)
    H5Fclose(source_h5file)
    
    # Open the destination HDF5 file in read-write mode
    destination_h5file <- H5Fopen(destination_file, "H5F_ACC_RDWR")
    H5Ldelete(destination_h5file, destination_group_path) # delete existing group
    h5write(source_group_data, destination_h5file, name = destination_group_path)
    H5Fclose(destination_h5file)
    
    h5closeAll() # in case anything isn't closed properly
    message("success")
  }
}

for (organ in organs){
  replace_arrow_gex_mtx_wrapper(input.dir, organ)
}



