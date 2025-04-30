library(ArchR)
library(tidyverse)
library(here)

# check if ncells in ##_atac.txt matches ncells in atac_preprocess_output/filtered_r1
organ_list <- c("Adrenal","Brain", "Eye", "Heart", "Liver", "Lung", "Muscle",
                "Skin", "Spleen", "Thymus")
organ_codes <- c("AG", "BR", "EY", "HT", "LI", "LU", "MU", "SK", "SP", "TM")

for (i in 1:length(organ_list)){
  ncell_whitelist <- dim(read_table(sprintf(here::here("output/01-preprocessing/02/shared/whitelist_r1/%s_atac.txt"), organ_codes[i])))[1]
  ncell_filtered <- loadArchRProject(sprintf(here::here("output/01-preprocessing/02/%s/atac_preprocess_output/filtered_r1"), organ_list[i]))$cellNames %>% length  
  flag <- (ncell_whitelist == ncell_filtered)
  print(sprintf("%s: %s whitelist %s filtered %s", organ_codes[i],flag, ncell_whitelist, ncell_filtered))
}

for (i in 1:length(organ_list)){
  ncell_whitelist <- dim(read_table(sprintf(here::here("output/01-preprocessing/02/shared/shared/meta/%s_meta.txt"), organ_codes[i])))[1]
  ncell_filtered <- loadArchRProject(sprintf(here::here("output/01-preprocessing/02/%s/atac_preprocess_output/ATAC_obj_clustered_final"), organ_list[i]))$cellNames %>% length  
  flag <- (ncell_whitelist == ncell_filtered)
  print(sprintf("%s: %s rnafinal %s atacfinal %s", organ_codes[i],flag, ncell_whitelist, ncell_filtered))
}
