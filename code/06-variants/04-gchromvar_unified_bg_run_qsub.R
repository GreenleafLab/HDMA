options(scipen = 999)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Matrix)
library(parallel)
library(here)

out <- here::here("output/06-variants/04/")
figout <- here::here("figures/06-variants/04/")
dir.create(out, showWarnings = F, recursive=T)
dir.create(figout, showWarnings = F, recursive=T)
setwd(out)

ch = import.chain(here::here("data/external/variants/hg38ToHg19.over.chain"))

args = commandArgs(trailingOnly = T)
tissue = args[1]
iteration = args[2]

dir.create(tissue, showWarnings = F)
load(paste0(out, "/hg19_inputs/", tissue, "_SE_scores.RData"))

bkgd_peaks = getBackgroundPeaks(SE, niterations = 500, w = 0.1, bs = 50)
wdev = computeWeightedDeviations(SE, scores, background_peaks = bkgd_peaks)
zscoreWeighted = t(assays(wdev)[["z"]])
mdf = reshape2::melt(zscoreWeighted)
mdf$iteration = iteration

fwrite(mdf, paste0(tissue, "/gchromvar_results_part", iteration, ".txt"), sep = "\t", row.names = F)


