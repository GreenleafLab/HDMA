# Purpose: Define the set of composite motifs to test for syntax/synergy.
# When there are several composites for the same pair of motifs (e.g. GATA GATA)
# which only differ by spacing/orientation, we only need to test one, since
# we'll exhaustively evaluate all arrangements. We also define which cell type
# to test the sequences in, based on the cell type with the most hits for each composite.


# set up -----------------------------------------------------------------------
library(here)
kundaje_dir <- trimws(readr::read_lines("../../AK_PROJ_DIR.txt"))
doc_id      <- "04"
out         <- here( "output/03-chrombpnet/03-syntax/", doc_id); dir.create(out, recursive = TRUE)
figout      <- here( "figures/03-chrombpnet/03-syntax", doc_id, "/"); dir.create(figout, recursive = TRUE)


# libraries --------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)
library(glue)
library(purrr)
library(stringr)


# Here, we look at the unique composite motifs (same granular annotations, but choosing)
# one variant per annotation e.g. if there are several GATA_MEF2 motifs:
  
composites_to_test <- motifs_compiled_unique %>%
  filter(category %in% c("homocomposite", "heterocomposite")) %>%
  filter(!grepl("half", motif_name) & !grepl("unresolved", motif_name)) %>%
  mutate(motif_name2 = motif_name) %>% 
  separate(motif_name2, into = c("idx_uniq", "annotation_varidx"), sep = "\\|") %>%
  separate(annotation_varidx, into = c("annotation", "variant_idx"), sep = "#") %>% 
  group_by(annotation) %>% 
  # there are many variants, choose the one with the most hits
  slice_max(order_by = total_hits, n = 1) %>% 
  dplyr::select(motif_name, category, idx_uniq, annotation_broad, annotation, variant_idx, total_hits) %>% 
  separate(annotation, into = c("component_motifA", "component_motifB"), sep = "_")

# for each composite, figure out which cell type has the most hits
composites_to_test <- hits_all_anno %>%
  dplyr::select(Cluster, motif_name, n_hits_per_cluster = n) %>% 
  right_join(composites_to_test, by = "motif_name") %>% 
  group_by(motif_name) %>%
  slice_max(order_by = n_hits_per_cluster, n = 1, with_ties = FALSE) %>%
  ungroup()

# here, we look at the granular annotations for base motifs, and select the most frequent between variants
# i.e. if there's an GATA#1, GATA#2, GATA#3; we'll choose the most frequent (by hits)
granular_base <- motifs_compiled_unique %>%
  filter(category %in% c("base", "base_with_flank")) %>%
  group_by(annotation_broad) %>%
  mutate(motif_name2 = motif_name) %>%
  separate(motif_name2, into = c("idx_uniq", "annotation_varidx"), sep = "\\|") %>%
  separate(annotation_varidx, into = c("annotation", "variant_idx"), sep = "#") %>%
  dplyr::select(motif_name, idx_uniq, annotation, variant_idx, total_hits, category, annotation_broad) %>%
  group_by(annotation) %>%
  slice_max(order_by = total_hits, n = 1, with_ties = FALSE)

# now, for each component motif, add the most frequent version among the granular_base df
composites_to_test <- composites_to_test %>%
  left_join(granular_base %>% dplyr::select(component_motifA_motif_name = motif_name, annotation), by = c("component_motifA" = "annotation")) %>%
  left_join(granular_base %>% dplyr::select(component_motifB_motif_name = motif_name, annotation), by = c("component_motifB" = "annotation")) %>%
  mutate(test_spacing = ifelse(!is.na(component_motifA_motif_name) & !is.na(component_motifB_motif_name), "Y", "N")) 

# now, for each component, get the high affinity sequence corresponding to the
# trimmed CWM of each component.
composites_to_test$high_affinity_seqA <- map_chr(
  composites_to_test$component_motifA_motif_name, ~ ifelse(is.na(.x), NA, get_high_affinity_seq(trim_cwm(cwm_list[[.x]], flank = 0))))

composites_to_test$high_affinity_seqB <- map_chr(
  composites_to_test$component_motifB_motif_name, ~ ifelse(is.na(.x), NA, get_high_affinity_seq(trim_cwm(cwm_list[[.x]], flank = 0))))

# add the motif CWM images
get_url <- function(i) {
  
  s3_url <- "https://human-dev-multiome-atlas.s3.amazonaws.com/trimmed_logos/"
  
  if (!is.na(i)) {
    pattern_class <- motifs_compiled_unique %>% filter(motif_name == i) %>% pull(pattern_class)
    pattern <- motifs_compiled_unique %>% filter(motif_name == i) %>% pull(pattern)
    key <- paste0(pattern_class, ".", pattern, ".cwm.fwd.png")
    
    paste0('=IMAGE("', s3_url, key, '",4,100,250)')
    
  } else return(NA)
  
}

composites_to_test$composite_cwm_fwd <- map_chr(composites_to_test$motif_name, ~ get_url(.x))
composites_to_test$componentA_cwm_fwd <- map_chr(composites_to_test$component_motifA_motif_name, ~ get_url(.x))
composites_to_test$componentB_cwm_fwd <- map_chr(composites_to_test$component_motifB_motif_name, ~ get_url(.x))

composites_to_test <- composites_to_test %>% 
  dplyr::relocate(composite_cwm_fwd, .after = motif_name) %>% 
  dplyr::relocate(Cluster, .before = n_hits_per_cluster) %>% 
  dplyr::relocate(componentA_cwm_fwd, high_affinity_seqA, .after = component_motifA_motif_name) %>% 
  dplyr::relocate(componentB_cwm_fwd, high_affinity_seqB, .after = component_motifB_motif_name)

composites_to_test %>%
  write_tsv(glue("{out}/composite_motif_spacing_tests.tsv"), escape = "none")