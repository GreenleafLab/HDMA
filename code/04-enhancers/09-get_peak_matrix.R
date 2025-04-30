library(ArchR)
library(tidyverse)
library(here)

proj <- loadArchRProject(here::here("output/01-preprocessing/03/allSamples_decontx"))
vista_peak_matrix <- readRDS(here::here("output/04-enhancers/09/vista_peak_matrix_cell.rds"))
addArchRThreads(8)


peak_matrix <- list()
for (chr in seqlevels(proj@peakSet)){
  peak_matrix_part <- getMatrixFromProject(proj, useMatrix = "PeakMatrix", useSeqnames = chr)
  overlaps <- findOverlaps(peak_matrix_part %>% rowRanges, vista_peak_matrix@rowRanges, minoverlap=375)
  peak_matrix[[chr]] <- peak_matrix_part[overlaps@from %>% unique]
  rm(peak_matrix_part)
}

peak_matrix <- do.call(rbind, peak_matrix)
saveRDS(peak_matrix, here::here("output/04-enhancers/09/peak_matrix_cell_overlapvista.rds"))


