#!/usr/bin/env Rscript

# Purpose: load global ArchR project which contains the ATACdata for all cells
# passing QC, and extract fragments per cluster per sample.

# Set up -----------------------------------------------------------------------

# load ArchR (and associated libraries)
library(ArchR)
library(dplyr)
library(purrr)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(parallel)

# parse args
args <- commandArgs(trailingOnly = TRUE)

# set working directory
setwd(args[1])

# set out
outdir <- file.path(args[2], "fragments")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# misc options
addArchRGenome("hg38")

# set genome
genome <- BSgenome.Hsapiens.UCSC.hg38.masked

proj_dir <- trimws(readr::read_lines("../../../ROOT_DIR.txt"))

# Get fragments per cluster per sample -----------------------------------------

# load global ArchR Project (which is a directory)
atac_proj <- loadArchRProject(file.path(proj_dir, "output/01-preprocessing/03/allSamples"))

# get a list of cell names in each cluster
clusters <- data.frame("cellName" = getCellNames(atac_proj),
                       "cluster" = getCellColData(ArchRProj = atac_proj, "RNA_Clusters")[[1]])

# number of cells per cluster
table(clusters$cluster)

# get fragments for one arrow file
split_fragments_by_cluster <- function(arrow) {
  
  message("@ processing ", arrow)
  
  # specify a file indicating the arrow file has been processed; skip processing it
  # already exists
  done <- file.path(outdir, paste0(".", gsub(".arrow", ".done", basename(arrow))))
  
  if (file.exists(done)) {
    message("@@ already done ", basename(arrow), "; skipping.")
    return(NULL)
  }
  
  # returns granges
  frags <- getFragmentsFromArrow(arrow)
  
  # convert to char
  mcols(frags)$barcode <- as.character(mcols(frags)$RG)
  
  # get the cluster assignment for each barcode
  frags$cluster <- clusters$cluster[match(frags$barcode, clusters$cellName)]
  
  # split per cluster
  frags_per_cluster <- split(frags, frags$cluster)
  
  # write to bed files per cluster, cluster__sample.tsv
  out_suffix <- gsub(".arrow", ".tsv", basename(arrow))
  
  # export as a tsv, 1-based coordinates
  iwalk(frags_per_cluster,
       ~ as.data.frame(.x) %>%
         dplyr::select(chr = seqnames, start, end, barcode) %>% 
         data.table::fwrite(file.path(outdir, paste0(.y, "__", out_suffix)),
                            col.names = FALSE, sep = "\t"))
  
  # specify done.
  file.create(done)
   
}

# loop over arrow files, split each one by cluster
# and write to tsv files
mclapply(getArrowFiles(atac_proj), split_fragments_by_cluster, mc.cores = 4)

message("@ done.")
