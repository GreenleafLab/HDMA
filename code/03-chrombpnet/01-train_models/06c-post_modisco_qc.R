# Purpose: after QC-ing modisco reports, we exclude models where we
# seemed to be strongly affected by GC bias. These generally came from low-coverage
# cell types with <1M fragments. Therefore, we remove those from our models-to-keep file.

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)
library(glue)
library(purrr)
library(stringr)
library(rvest) # for parsing modisco HTML reports
library(universalmotif) # for working with motifs

hdma_path   <- here::here()
bias_params <- "bias_Heart_c0_thresh0.4"
out         <- here("output/03-chrombpnet/01-models/qc")

chrombpnet_models_keep <- read_tsv(file.path(out, "chrombpnet_models_keep.tsv"),
                                   col_names = c("Cluster", "Folds_keep", "Cluster_ID"))
length(unique(chrombpnet_models_keep$Cluster))

models_rm_post_modisco <- c("Adrenal_c4", "Adrenal_c5", "Adrenal_c6", "Eye_c19", "Eye_c20" ,"Thymus_c16", "Thyroid_c8", "Thyroid_c10", "Thyroid_c11", "Thyroid_c9")

n_frags <- read_tsv(here("output/01-preprocessing/03/fragmentsPerCluster.tsv")) %>% 
  dplyr::rename(total_frags = total_reads)

cluster_meta <- read_csv(here("output/05-misc/03/TableS2_cluster_meta_qc.csv")) %>% 
  left_join(n_frags, by = c("Cluster" = "RNA_Clusters"))

# check # frags
cluster_meta %>% filter(Cluster_ChromBPNet %in% models_rm_post_modisco) %>% dplyr::select(Cluster_ChromBPNet, ncell, total_frags) %>% arrange(desc(total_frags))
# # A tibble: 10 × 3
# Cluster_chrombpnet ncell total_frags
# <chr>              <dbl>       <dbl>
#   1 Eye_c19              519     3611353
# 2 Eye_c20              329     2798813
# 3 Thyroid_c8           200     1932991
# 4 Thymus_c16            98     1333100
# 5 Adrenal_c4           155      991281
# 6 Adrenal_c5           115      747083
# 7 Thyroid_c10          104      733861
# 8 Thyroid_c9           118      681542
# 9 Adrenal_c6            89      526401
# 10 Thyroid_c11           38      227150

# 189 models remaining
length(setdiff(chrombpnet_models_keep$Cluster, models_rm_post_modisco))

chrombpnet_models_keep %>% 
  filter(!(Cluster %in% models_rm_post_modisco)) %>% 
  write_tsv(file.path(out, "chrombpnet_models_keep2.tsv"))
