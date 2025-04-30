# Extract peaks and counts from ArchR projects
# For gchromvar inputs
# Also update immune cell metadata
library(tidyverse)
library(data.table)
library(ArchR)
library(GenomicRanges)
library(here)


out <- here::here("output/06-variants/00/")
figout <- here::here("figures/06-variants/00/")
dir.create(out, showWarnings = F, recursive=T)
dir.create(figout, showWarnings = F, recursive=T)
setwd(out)

dir.create("count_matrices", showWarnings = F)
dir.create("peaks", showWarnings = F)

tissues = c("Adrenal", "Heart", "Lung", "Skin", "Spleen", "StomachEsophagus",
            "Brain", "Eye", "Liver", "Muscle", "Thyroid", "Thymus")
tis2abbrev = data.frame(tis = tissues,
                        abbrev = c("AG", "HT", "LU", "SK", "SP", "ST", "BR", "EY", "LI", "MU", "TR", "TM"))

lapply(tissues, function(tis){
  cat(tis, "\n")
  meta = fread(paste0(here::here("output/01-preprocessing/02/shared/meta/"),
                      tis2abbrev[tis2abbrev$tis == tis,]$abbrev, "_meta.txt"), data.table = F) %>%
    select(L1_clusterID, L2_clusterName, L3_clusterName) %>% distinct()
  
  proj = loadArchRProject(paste0(here::here("output/01-preprocessing/02/"), tis, "/atac_preprocess_output/ATAC_obj_clustered_peaks_final_decontx"),
                          showLogo = F)
  coldat = as.data.frame(getCellColData(ArchRProj = proj, select = "L1_clusterID"))
  coldat$cb = row.names(coldat)
  coldat = left_join(coldat, meta, by = "L1_clusterID")
  
  proj = addCellColData(
    ArchRProj = proj,
    data = coldat$L2_clusterName,
    name = "L2_clusterName",
    cells = coldat$cb,
    force = T
  )
  proj = addCellColData(
    ArchRProj = proj,
    data = coldat$L3_clusterName,
    name = "L3_clusterName",
    cells = coldat$cb,
    force = T
  )
  
  grpSE = getGroupSE(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "L2_clusterName"
  )
  
  cts = grpSE@assays@data@listData$PeakMatrix
  fwrite(cts, paste0("count_matrices/",
                     tis, ".tsv"), row.names = F, sep = "\t")
  system(paste0("gzip -f count_matrices/",
                tis, ".tsv"))
  
  peaks = grpSE@elementMetadata %>%
    as.data.frame() %>%
    select(seqnames, start, end)
  # BEDs are 0-based
  peaks$start = as.numeric(peaks$start -1)
  fwrite(peaks, paste0("peaks/",
                       tis, ".bed"), row.names = F, sep = "\t", col.names = F)
  system(paste0("gzip -f peaks/",
                tis, ".bed"))
})

