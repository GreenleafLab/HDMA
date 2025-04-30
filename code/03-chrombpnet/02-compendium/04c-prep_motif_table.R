# Purpose: Given the clustered motifs, we prepare a table where each row is
# a clustered output pattern. We aggregate info about the cell types from which
# the pattern was derived, the number of seqlets supporting it, the most
# similar motifs from the TOMTOM matching, the family of the best-matching TFs
# from HOCOMOCO, and the CWM motif logos.

proj_dir <- here::here()
kundaje_dir <- trimws(readr::read_lines("../../AK_PROJ_DIR.txt"))
doc_id <- "04c"
figout <- here("figures/03-chrombpnet/02-compendium/", doc_id, "/"); dir.create(figout, recursive = TRUE, showWarnings = FALSE)
out <- here("output/03-chrombpnet/02-compendium/modisco_compiled_anno")

# setup ------------------------------------------------------------------------
set.seed(100)

# libraries --------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)
library(glue)
library(purrr)
library(stringr)
library(cowplot)

script_path <- here("code/utils/")
source(file.path(script_path, "plotting_config.R")) #Assumes these files are on same level if in a subdirectory in hdma/utils
source(file.path(script_path, "hdma_palettes.R"))
source(file.path(script_path, "sj_scRNAseq_helpers.R"))

ggplot2::theme_set(theme_BOR())


# load data --------------------------------------------------------------------

# Cluster metadata
cluster_meta <- read_csv(here("output/05-misc/03/TableS2_cluster_meta_qc.csv"))

chrombpnet_models_keep <- read_tsv("../output/01-models/qc/chrombpnet_models_keep2.tsv",
                                   col_names = c("Cluster", "Folds_keep", "Cluster_ID"))

cluster_order <- cluster_meta %>% arrange(organ, L1_annot) %>% pull(Cluster_ChromBPNet)
organs <- unique(cluster_meta$organ) %>% sort()
organs[organs == "StomachEsophagus"] <- "Stomach"

head(chrombpnet_models_keep)


# Load merged modisco patterns

# indiv reports
counts_modisco_reports <- read_tsv(here("output/03-chrombpnet/01-models/modisco_tsv/all_modisco_report.motifs.counts.tsv")) %>%
  filter(Cluster %in% chrombpnet_models_keep$Cluster)
dim(counts_modisco_reports)

# merged reports
modisco_merged <- map_dfr(list.files(here("output/03-chrombpnet/02-compendium/modisco_merged/"), recursive = TRUE, pattern = "merged_modisco.tsv", full.names = TRUE),
                          ~ read_tsv(.x))

modisco_merged <- modisco_merged %>%
  mutate(component_pattern = paste0(component_celltype, "__", pattern_class, ".", pattern))
dim(modisco_merged)

# save merged modisco
write_tsv(modisco_merged, file = glue("{out}/modisco_merged_reports.tsv"))


# compiled reports with tomtom matches againt Vierstra database
modisco_compiled <- read_tsv(here("output/03-chrombpnet/02-compendium/modisco_compiled/modisco_compiled.tsv")) %>%
  # fix logo names
  mutate(modisco_cwm_fwd = ifelse(grepl("^pos", pattern), gsub("pos.", "pos_patterns.pos.", modisco_cwm_fwd), gsub("neg.", "neg_patterns.neg.", modisco_cwm_fwd)),
         modisco_cwm_rev = ifelse(grepl("^pos", pattern), gsub("pos.", "pos_patterns.pos.", modisco_cwm_rev), gsub("neg.", "neg_patterns.neg.", modisco_cwm_rev)))


# Load external motif and TF info

# Vierstra annotations:
vierstra_metadata <- read_tsv(file.path(kundaje_dir, "refs/Vierstra_motifs/metadata.tsv"))

# NOTE: there are a few motifs in the Vierstra meme db that are not in the metadata.
motifs_to_add <- map(0:9, ~ setdiff(modisco_compiled[[glue("match{.x}")]], vierstra_metadata$motif_id)) %>%
  unlist() %>%
  unique()
motifs_to_add

# they seem to be contained in the source id column
map(0:9, ~ setdiff(modisco_compiled[[glue("match{.x}")]], c(vierstra_metadata$motif_id, vierstra_metadata$source_id))) %>%
  unlist() %>%
  unique()

# let's create new rows for those
motifs_to_add_rows <- vierstra_metadata %>% filter(source_id %in% setdiff(motifs_to_add, "NA")) %>% mutate(motif_id = source_id)
vierstra_metadata <- bind_rows(vierstra_metadata, motifs_to_add_rows)

map(0:9, ~ setdiff(modisco_compiled[[glue("match{.x}")]], vierstra_metadata$motif_id)) %>%
  unlist() %>%
  unique()

write_tsv(vierstra_metadata, file = glue("{out}/vierstra_metadata_updated.tsv"))

# HOCOMOCO annotations:
hocomoco <- read_tsv(here("output/02-global_analysis/04/H12CORE_annotation.tsv"))
hocomoco_unique <- hocomoco %>%
  distinct(masterlist_info.species.HUMAN.gene_symbol, .keep_all = TRUE) %>%
  mutate(
    tf = masterlist_info.species.HUMAN.gene_symbol,
    tf_family = case_when(
      !is.na(masterlist_info.tfclass_subfamily) ~ paste0(masterlist_info.tfclass_family, "_", masterlist_info.tfclass_subfamily),
      is.na(masterlist_info.tfclass_subfamily) ~ paste0(masterlist_info.tfclass_family)
    )) %>%
  dplyr::select(tf, tf_family)


# data wrangling ---------------------------------------------------------------

# First, for each merged pattern, let's get a summary of which input patterns produced it,
# and which cell types and organs they come from.
modisco_merged_summarized <- modisco_merged %>%
  left_join(cluster_meta %>% dplyr::select(organ, Cluster_chrombpnet, L2_annot), by = c("component_celltype" = "Cluster_chrombpnet")) %>%
  arrange(merged_pattern, organ, component_celltype) %>%
  group_by(merged_pattern) %>%
  summarize(total_n_seqlets = sum(n_seqlets),
            n_component_patterns = n(),
            n_component_celltypes = length(unique(component_celltype)),
            n_component_organs = length(unique(organ)),
            component_organs = glue_collapse(unique(organ), sep = ","),
            component_celltypes = glue_collapse(unique(L2_clusterName), sep = ",")) %>%
  dplyr::rename(pattern = merged_pattern)

head(modisco_merged_summarized)
dim(modisco_merged_summarized)

# Now let's add this to the compiled report:
all(modisco_compiled$pattern %in% modisco_merged_summarized$pattern)

modisco_compiled2 <- modisco_compiled %>%
  left_join(modisco_merged_summarized, by = "pattern") %>%
  dplyr::select(pattern, matches("n_"), component_organs, component_celltypes, modisco_cwm_fwd, modisco_cwm_rev, query_consensus, everything())

# And add URLs to logos:
s3_url <- "https://human-dev-multiome-atlas.s3.amazonaws.com/"

get_url <- function(logo_path) {
  if (!is.na(logo_path)) {
    paste0('=IMAGE("', gsub('\\./', s3_url, logo_path), '",4,100,250)')
  } else {
    NA_character_ # Return NA as character
  }
}

for (i in 0:2) {
  col_name <- glue("match{i}_logo")
  modisco_compiled2[, col_name] <- map_chr(modisco_compiled2[[col_name]], ~ get_url(.x))
}

modisco_compiled2$modisco_cwm_fwd <- map_chr(modisco_compiled2$modisco_cwm_fwd, ~ get_url(.x))
modisco_compiled2$modisco_cwm_rev <- map_chr(modisco_compiled2$modisco_cwm_rev, ~ get_url(.x))

cat(modisco_compiled2[1, ]$modisco_cwm_fwd)
cat(modisco_compiled2[1, ]$modisco_cwm_rev)
cat(modisco_compiled2[1, ]$match0_logo)


# Add external annotations to compiled report, for the first three matches
annotate_motif_vierstra <- function(motif_id) {
  if (motif_id %in% vierstra_metadata$motif_id) {
    paste0(vierstra_metadata[vierstra_metadata$motif_id == motif_id, ]$tf_name, ": ", motif_id)
  } else {
    motif_id
  }
}

modisco_compiled3 <- modisco_compiled2

for (i in 0:9) {
  col_name <- glue("match{i}")
  modisco_compiled3[, col_name] <- map_chr(modisco_compiled3[[col_name]], ~ annotate_motif_vierstra(.x))
}

# for a given TF:motif_id pair, find the other TFs in the same family
# and the other TFs in the same archetype
expand_tf_set <- function(i) {
  if (!is.na(i)) {
    inputs <- str_split(i, ": ") %>% unlist()
    
    if (inputs[1] %in% hocomoco_unique$tf) {
      tf_family_i <- hocomoco_unique %>% filter(tf == inputs[1]) %>% pull(tf_family)
      tfs_in_family <- hocomoco_unique %>% filter(tf_family == tf_family_i) %>% pull(tf) %>%
        unique() %>% sort() %>% glue_collapse(sep = ",")
    } else {
      tfs_in_family <- "not found in HOCOMOCO"
    }
    
    archetype_i <- vierstra_metadata %>% filter(motif_id == inputs[2]) %>% pull(cluster)
    tfs_in_archetype <- vierstra_metadata %>% filter(cluster == archetype_i) %>%
      pull(tf_name) %>% unique() %>% sort() %>% glue_collapse(sep = ",")
    tfs_in_archetype <- paste0(archetype_i, ": ", tfs_in_archetype)
    
  } else {
    tfs_in_family <- "no motif match"
    tfs_in_archetype <- "no motif match"
  }
  
  return(list("TFs_in_family" = tfs_in_family,
              "TFs_in_archetype" = tfs_in_archetype))
}

modisco_compiled4 <- modisco_compiled3

expanded_tfs <- map(modisco_compiled4$match0, ~ expand_tf_set(.x)) %>% purrr::transpose()
modisco_compiled4 <- modisco_compiled4 %>%
  tibble::add_column(match0_TFs_in_family = unlist(expanded_tfs$TFs_in_family),  .after = "match0_logo") %>%
  tibble::add_column(match0_TFs_in_archetype = unlist(expanded_tfs$TFs_in_archetype), .after = "match0_TFs_in_family")


# save the table. --------------------------------------------------------------
modisco_compiled4 %>%
  mutate(pattern_class = ifelse(grepl("^pos", pattern), "pos", "neg")) %>%
  arrange(desc(pattern_class), pattern) %>%
  dplyr::relocate(pattern_class, .before = pattern) %>%
  tibble::rowid_to_column(var = "idx") %>%
  dplyr::select(-matches("e_val")) %>%
  write_tsv(file = glue("{out}/modisco_compiled_annotated.tsv"), escape = "none")

# done.