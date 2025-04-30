#4 script wrapper
library(tidyverse)

setwd(here("code/logs"))

toKeep <- read_tsv(here("code/02-global_analysis/01-organs_keep.tsv"))
organs <- toKeep$organs

script <- here("code/01-preprocessing/02-archr_seurat_scripts/4_readd_gex_peak2gene_decontx.R")

lapply(organs, function(organ){
  com <- sprintf("sbatch -p wjg,sfgf,biochem --mem-per-cpu=200g --time=20:00:00 --job-name=p2g_%s --wrap \"Rscript %s %s\"", 
          organ, script, organ)
  system(com)
})

