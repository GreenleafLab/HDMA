# Purpose: this script loops through all modisco reports and converts the HTMLs
# to TSV for easy parsing and re-use later.

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
out         <- file.path(hdma_path, "output/03-chrombpnet/01-models/modisco_tsv")
modisco_dir <- file.path(hdma_path, "output/03-chrombpnet/01-models/modisco")

dir.create(out, showWarnings = FALSE)

#' Convert the motifs.html output from modisco to a TSV, adding the name of the
#' model ("Cluster") as well as the model head used for interpretation ("Model_head")
#' 
#' @param key character, Name of model/cluster
#' @param key characger, Either "counts" or "profile"
modisco_report_to_tsv <- function(key, model_head = "counts") {
  
  message("@ ", key)
  
  modisco_report <- rvest::read_html(glue(
    "{modisco_dir}/{bias_params}/{key}/{model_head}_modisco_report/motifs.html"))
  
  # read in MEME format (only for pos patterns though)
  meme <- universalmotif::read_meme(glue(
    "{modisco_dir}/{bias_params}/{key}/{key}_memedb.{model_head}.txt"))
  
  # get consensus sequences using universal motif
  meme_df <- data.frame(pattern = map_chr(meme, ~ glue('{key}__pos_patterns.{.x["name"]}')),
                        consensus_fwd = map_chr(meme, ~ .x["consensus"]),
                        consensus_rev = map_chr(meme, ~ universalmotif::motif_rc(.x)["consensus"]))
  
  
  col_names <- modisco_report %>% html_elements("thead") %>% html_table() %>% getElement(1) %>% colnames()
  modisco_df <- modisco_report %>% html_elements("tbody") %>% html_table() %>% 
    getElement(1) %>% 
    set_colnames(col_names) %>% 
    mutate(pattern = paste0(key, "__", pattern)) %>% 
    tibble::add_column("Cluster" = key, "Model_head" = model_head, .before = 1) %>% 
    left_join(meme_df, by = "pattern") %>% 
    dplyr::relocate(consensus_fwd, .after = 3) %>% 
    dplyr::relocate(consensus_rev, .after = consensus_fwd)
  
  data.table::fwrite(modisco_df, file = glue("{out}/{key}_modisco_report.counts.tsv"), sep = "\t")
  
  return(modisco_df)
  
}

chrombpnet_models_keep <- read_tsv(file.path(hdma_path, "output/03-chrombpnet/01-models/qc/chrombpnet_models_keep.tsv"),
                                   col_names = c("Cluster", "Folds_keep", "Cluster_ID"))
length(unique(chrombpnet_models_keep$Cluster))

counts_modisco_reports <- map_dfr(chrombpnet_models_keep$Cluster, ~ modisco_report_to_tsv(.x))

write_tsv(counts_modisco_reports, glue("{out}/all_modisco_report.motifs.counts.tsv"))

message("@ done.")