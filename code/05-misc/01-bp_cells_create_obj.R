# Use BPCells to plot browser tracks efficiently -----------------------------------------------------------------

# imports -----------------------------------------------------------------
library(BPCells)
library(ArchR)
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(here)
})


# user inputs -------------------------------------------------------------
wd <- here::here("output/05-misc/01")
plotdir <- here::here("figures/05-misc/01")
rawdata.dir <- here::here("output/01-preprocessing/01")
input.dir <- here::here("output/01-preprocessing/02")
global.peaks.dir <- here::here("output/01-preprocessing/03")
tissue_meta <- read_tsv(here::here("code/02-global_analysis/01-organs_keep.tsv"))

dir.create(wd, recursive = TRUE, showWarnings = FALSE)
setwd(wd)

scriptPath <- here::here("code/utils")
source(paste0(scriptPath, "/track_helpers_BL.R"))
source(paste0(scriptPath, "/track_helpers_SJ.R"))
source(paste0(scriptPath, "/trackplots.R")) # BPCells trackplot helpers in dev

# main ---------------------------------------------------------------------
## load ATAC fragment files -------------------------------------------------

# consider parallel::mclapply(list, function, mc.cores) for higher sample numbers

# get all fragment files first, exclude fail batches
frags <- Sys.glob(paste0(rawdata.dir, "/*/ATAC/*fragments.tsv.gz"))
frags <- frags[-grep("dummy|b14|b1/|b2/|_b12/|_b18/|Kidney", frags)]

# get all samples pf names
sample_names <- list.files(file.path(global.peaks.dir, "allSamples/ArrowFiles")) %>% 
  basename() %>% str_replace(".arrow","") %>% paste0(collapse="|")
frags <- frags[grep(sample_names,frags)]

# get all whitelisted cell names
global_proj <- loadArchRProject(file.path(global.peaks.dir, "allSamples_decontx"))
cell_whitelist <- getCellNames(global_proj)           

frags_all <- c()
for (f in frags){
  name <- basename(f) %>% str_replace(".fragments.tsv.gz", "")
  # Check if we already ran import
  if (!file.exists(name)) {
    message(paste0("creating BPcells folder ", name))
    frags_raw <- open_fragments_10x(f) %>%
      write_fragments_dir(paste0(name, "_tmp"))
    
    # add prefix to cell names and filter by whitelist
    frags_raw@cell_names <- paste0(name, "#", frags_raw@cell_names)
    frags_filt <- frags_raw %>% 
                  select_cells(frags_raw@cell_names[(frags_raw@cell_names %in% cell_whitelist)]) %>% 
                  write_fragments_dir(name)
  } else {
    message(paste0("reading BPcells folder ", name))
    frags_filt <- open_fragments_dir(name)
  }
  frags_all <- c(frags_all, frags_filt %>% select_chromosomes(paste0("chr", c(1:22, "X"))))
}
frags_all <- do.call(c, frags_all)

## load decontx RNA matrices -----------------------------------------------
# read consensus gene set
genes <-readRDS(paste0(global.peaks.dir, "/global_gexfeat_decontx_gr.rds"))

rna_all <- c()
for (organ in tissue_meta$organ){
  output <- paste0("RNA_", organ)
  if (!file.exists(output)){
    message(paste0("creating BPcells folder ", output))
    seurat <- readRDS(paste0(input.dir, "/", organ, "/rna_preprocess_output/RNA_obj_clustered_final.rds")) # load individual organ RNA objs

    message(paste0("joining with consensus genes ", output))
    output_mat <- Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(length(genes), ncol(seurat)), dimnames=list(names(genes), colnames(seurat)))
    output_mat <- as(output_mat, "IterableMatrix")
    keeper_genes <- intersect(rownames(seurat), names(genes))
    output_mat[keeper_genes,] <- seurat@assays$decontX@counts[keeper_genes,]
    rm(seurat) # save memory
    
    # # Check for all integers from decontX, and not tossing too many genes
    # stopifnot(all.equal(output_mat@x, floor(output_mat@x)))
    # stopifnot(length(keeper_genes) > .95*nrow(seurat))
    
    # Write output
    message(paste0("writing BPcells folder ", output))
    rna <- output_mat %>%
      convert_matrix_type("uint32_t") %>%
      transpose_storage_order() %>%
      write_matrix_dir(output)
  } else {
    message(paste0("reading BPcells folder ", output))
    rna <- open_matrix_dir(output)
  }
  rna_all <- c(rna_all, rna)
}
rna_all <- do.call(cbind, rna_all)

## load raw RNA matrices -----------------------------------------------
# read consensus gene set
genes <- readRDS(paste0(global.peaks.dir, "/global_gexfeat_gr.rds"))

rna_raw_all <- c()
for (organ in tissue_meta$organ){
  output <- paste0("RNA_raw_", organ)
  if (!file.exists(output)){
    message(paste0("creating BPcells folder ", output))
    seurat <- readRDS(paste0(input.dir, "/", organ, "/rna_preprocess_output/RNA_obj_clustered_final.rds")) # load individual organ RNA objs
    
    message(paste0("joining with consensus genes ", output))
    output_mat <- Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(length(genes), ncol(seurat)), dimnames=list(names(genes), colnames(seurat)))
    output_mat <- as(output_mat, "IterableMatrix")
    keeper_genes <- intersect(rownames(seurat), names(genes))
    output_mat[keeper_genes,] <- seurat@assays$RNA@counts[keeper_genes,]
    rm(seurat) # save memory
    
    # # Check for all integers from decontX, and not tossing too many genes
    # stopifnot(all.equal(output_mat@x, floor(output_mat@x)))
    # stopifnot(length(keeper_genes) > .95*nrow(seurat))
    
    # Write output
    message(paste0("writing BPcells folder ", output))
    rna <- output_mat %>%
      convert_matrix_type("uint32_t") %>%
      transpose_storage_order() %>%
      write_matrix_dir(output)
  } else {
    message(paste0("reading BPcells folder ", output))
    rna <- open_matrix_dir(output)
  }
  rna_raw_all <- c(rna_raw_all, rna)
}
rna_raw_all <- do.call(cbind, rna_raw_all)

## make metadata file -----------------------------------------------------

# cell level annot
organ.code.list <- tissue_meta$organcode
cell.annots <- lapply(1:length(organ.code.list), function(i){
  # read cluster id annotation
  annot <- read.csv(sprintf(here::here("output/01-preprocessing/02/shared/meta/%s_meta.txt"), organ.code.list[i]), sep="\t") %>% as.data.frame
  rownames(annot) <- annot$cb
  annot$L0_clusterID <- paste0(organ.code.list[i],"_",annot$L1_clusterID)
  annot$L3_clusterID <- paste0(organ.code.list[i],"_",annot$L3_clusterName)
  annot$organ <- organ.code.list[i]
  return(annot)
})
cell.annots <- dplyr::bind_rows(cell.annots)
cell.annots$nFrags <- getCellColData(global_proj)[cell.annots$cb, "nFrags"]
cell.annots <- cell.annots %>% dplyr::rename(organ_code=organ)

# update with latest L1/L2/L3 cell annotations
new_clust_annots <- read_csv(here::here("output/05-misc/03/TableS2_cluster_meta_qc.csv")) %>%
                  dplyr::select(Cluster, L1_annot, L2_annot, L3_annot, cluster_id, compartment, organ)
colnames(cell.annots) <- gsub("L", "archive_L", colnames(cell.annots))
df <- merge(new_clust_annots, cell.annots, by.x="Cluster", by.y="archive_L0_clusterID")
df <- df %>% mutate(Cluster_chrombpnet = ifelse(organ == "StomachEsophagus", paste0("Stomach_c", cluster_id), paste0(organ, "_c", cluster_id)))
cell.annots <- df %>% dplyr::select(cb, Cluster, Cluster_chrombpnet, organ, organ_code, nFrags, L1_annot, L2_annot, L3_annot, cluster_id, compartment, 
                           grep("archive", colnames(df), value=T))

## concatenate all frags + RNA --------------------------------------------
frags_all <- select_cells(frags_all, cell.annots$cb) 
rna_all <- rna_all[,cell.annots$cb]
rna_all <- log1p(multiply_cols(rna_all, 1/colSums(rna_all))*10000) # log norm RNA decontX'ed

rna_raw_all <- rna_raw_all[,cell.annots$cb]
rna_raw_all <- log1p(multiply_cols(rna_raw_all, 1/colSums(rna_raw_all))*10000) # log norm RNA raw

## get the transcript reference ----------------------------------------------
transcripts <- read_gencode_transcripts(
  "./references", 
  release="42", 
  transcript_choice="MANE_Select",
  annotation_set = "basic", 
  features="transcript" # Make sure to set this so we don't get exons as well
)

## save global bp object -------------------------------------------------------
rna_all <- write_matrix_dir(rna_all, "RNA_merged")
rna_raw_all <- write_matrix_dir(rna_raw_all, "RNA_raw_merged")
frags_all <- write_fragments_dir(frags_all, "ATAC_merged")

global_bp_obj <- list()
global_bp_obj$cell_metadata <- cell.annots
global_bp_obj$frags <- frags_all
global_bp_obj$rna <- rna_all
global_bp_obj$rna_raw <- rna_raw_all
global_bp_obj$transcripts <- transcripts
saveRDS(global_bp_obj, "global_bp_obj.rds")

## get p2g loops --------------------------------------------------
global_bp_obj <- readRDS("global_bp_obj.rds")
p2g_list <- list()

# add global p2g loops to BPCells object
global_proj <- loadArchRProject(file.path(global.peaks.dir, "/allSamples_decontx"))
loops <- get_p2g_loops(global_proj)
#global_bp_obj$p2g_loops_global <- GRanges(loops)
p2g_list[["global"]] <- GRanges(loops)
         
# add per organ p2g loops
for (organ in tissue_meta$organ){
  message(organ)
  atac_proj <- loadArchRProject(paste0(input.dir, "/", organ, "/atac_preprocess_output/ATAC_obj_clustered_peaks_final_decontx"))
  loops <- get_p2g_loops(atac_proj)
  p2g_list[[organ]] <- GRanges(loops)
}

for(organ in tissue_meta$organ){
  message(organ)
  loops <- global_bp_obj[[paste0("p2g_loops_", organ)]]
  p2g_list[[organ]] <- loops
}

global_bp_obj$loops_p2g <- p2g_list
saveRDS(global_bp_obj, "global_bp_obj.rds")

## get abc loops --------------------------------------------------
abc_list <- list()

for (organ in tissue_meta$organ){
  message(organ)
  abc_raw <- readRDS(paste0(here::here("output/04-enhancers/06/abc_raw_"), organ, ".rds"))
  loops <- abc_raw[,c("chr", "start", "end", "name", "class", "TargetGene","TargetGeneTSS", "ABC.Score", "CellType")] %>%
            dplyr::rename(peak_start=start, peak_end=end)
  loops$start <- pmin((loops$peak_start + loops$peak_end)%/%2, loops$TargetGeneTSS)
  loops$end <- pmax((loops$peak_start + loops$peak_end)%/%2, loops$TargetGeneTSS)
  abc_list[[organ]] <- GRanges(loops)
}

global_bp_obj$loops_abc <- abc_list
saveRDS(global_bp_obj, "global_bp_obj.rds")

## save ABC/P2G loops for WashU browser visualization ------------------------
color_fn <- colorRampPalette(cmap_cor)
color_column <- color_fn(100) # Generate a range of 100 colors

for (i in seq_along(tissue_meta$organ)){
  organ <- tissue_meta$organ[i]
  organcode <- tissue_meta$organcode[i]
  message(organ)
  org_abc <- global_bp_obj$loops_abc[[organ]]
  for (clust in unique(org_abc$CellType)){
    message(clust)
    clust_abc <- org_abc[org_abc$CellType == clust]
    hex_colors <- color_column[ceiling(clust_abc$ABC.Score * (length(color_column) - 1)) + 1]
    clust_abc$chr <- seqnames(clust_abc)
    df <- mcols(clust_abc) %>% as.data.frame %>% 
      dplyr::mutate(seqnames=chr, start=TargetGeneTSS, end=TargetGeneTSS+1,
                    peak=paste0(chr, ":", peak_start, "-", peak_end,",",ABC.Score),
                    name=paste0(TargetGene, "|", class), color=hex_colors) %>%
      dplyr::select(c(seqnames, start, end, peak, name, color))
    # clust_name <- paste0(organcode, "_", str_split(clust, "_", simplify=T)[2])
    clust_name <- clust
    outfile <- paste0(here::here("output/04-enhancers/06/"), clust_name, "_ABC_threshold0.013.longrange")
    write.table(df, file=outfile,
                quote=F, row.names=F, col.names=F, sep="\t")
  }
  
}

## add global peaks --------------------------------------------------
global_peaks <- readRDS(paste0(global.peaks.dir, "/peaks_all_df_annots.rds"))
global_bp_obj$peaks$global <- GRanges(global_peaks)
saveRDS(global_bp_obj, "global_bp_obj.rds")

## add motif hits from chrombpnet --------------------------------------------------
hits_dir <- here::here("output/03-chrombpnet/02-compendium/hits_unified_motifs/reconciled_per_celltype_peaks/") 
clusts <- list.files(hits_dir)

hits_ls <- list()
for (clust in clusts){
  hits_ls[[clust]] <- rtracklayer::import.bed(
    Sys.glob(glue(hits_dir, "/{clust}/counts_v0.23_a0.8_all/*reconciled.bed.gz")),
    extraCols = c("pattern_class" = "character"))
}
global_bp_obj$hits <- hits_ls
saveRDS(global_bp_obj, "global_bp_obj.rds")
