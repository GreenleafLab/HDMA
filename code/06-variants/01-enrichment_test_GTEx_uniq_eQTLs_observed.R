# Test if motifs are enriched in eQTL variants from
# GTEx v8 eQTL gene-variant pairs

library(data.table)
library(tidyverse)
library(arrow)
library(GenomicRanges)
library(parallel)
library(here)

out <- here::here("output/06-variants/01/")
figout <- here::here("figures/06-variants/01/")
dir.create(out, showWarnings = F, recursive=T)
dir.create(figout, showWarnings = F, recursive=T)
setwd(out)

#### File paths and functions --------------------------------------------------
# eQTL variants
# Each file is parquet
ct.eqtl.files = list.files(here::here("data/external/variants/GTEx_eQTL_SuSiE"),
                           full.names = T)
ct.eqtl.files = ct.eqtl.files[grep("README.txt", ct.eqtl.files, invert = T)]

# Motifs
motif_hits_list = list.files(here::here("output/03-chrombpnet/02-compendium/hits_unified_motifs/reconciled_per_celltype_peaks/"),
                             pattern = "reconciled.bed.gz$", full.names = T, recursive=T)

# Function to count variants in motifs
count_hits = function(ct.eqtls, ct.hits){
  ct.eqtls.gr = makeGRangesFromDataFrame(ct.eqtls,
                                         seqnames.field = "chr",
                                         start.field = "pos",
                                         end.field = "pos",
                                         keep.extra.columns = T)
  ct.hits.gr = makeGRangesFromDataFrame(ct.hits,
                                        seqnames.field = "chr",
                                        start.field = "start",
                                        end.field = "end",
                                        keep.extra.columns = T)
  overlaps = findOverlaps(ct.eqtls.gr, ct.hits.gr)
  eqtls.in.motifs.df = ct.eqtls.gr[queryHits(overlaps)]
  eqtls.in.motifs.df = cbind.data.frame(
    mcols(eqtls.in.motifs.df),
    mcols(ct.hits.gr[subjectHits(overlaps)])
  )
  cts.df = eqtls.in.motifs.df %>%
    group_by(motif) %>%
    summarize(cts = n()) %>%
    as.data.frame()
  return(cts.df)
}

add_missing = function(results.df, ct.hits.df){
  all_motifs = unique(ct.hits.df$motif)
  not_in_results = all_motifs[!all_motifs %in% results.df$motif]
  missing = data.frame(motif = not_in_results,
                       cts = 0)
  results.df = rbind(results.df, missing)
  return(results.df)
}


#### Assemble eQTL variants per organ ------------------------------------------
ct.eqtls = lapply(ct.eqtl.files, function(eqtlfile){
  organ = word(eqtlfile, 2, 2, sep = "\\.")
  df = read_parquet(eqtlfile)
  # Assign effect on gene expression based on allelic fold change (aFC)
  # Positive is increase, negative is decrease
  df$direction = ifelse(df$afc < 0, "down", "up")
  df$chr = word(df$variant_id, 1, 1, sep = "_")
  df$pos = as.numeric(word(df$variant_id, 2, 2, sep = "_"))
  df$id = paste0(df$variant_id, "_", df$phenotype_id)
  df = select(df, chr, pos, id, direction, pip, afc)
  df$organName = organ
  df
}) %>% bind_rows()

ct.eqtls$organ = word(ct.eqtls$organName, 1, 1, sep = "_")
# Group Esophagus together with Stomach
ct.eqtls$organ = gsub("Esophagus", "Stomach", ct.eqtls$organ)
# For repeated ids in each organ, pick the one with the higher SuSiE PIP
ct.eqtls = ct.eqtls %>%
  group_by(organ, chr, pos, id) %>%
  arrange(desc(pip)) %>%
  slice_head(n = 1)

##### Plot some stats ----------------------------------------------------------
# Number of gene-variant pairs tested per organ
tbl.eqtls = janitor::tabyl(ct.eqtls, var1 = organ, var2 = direction)
tbl.eqtls = tbl.eqtls %>%
  pivot_longer(cols = c("down", "up"), names_to = "eQTL_effect_dir", values_to = "n") %>%
  as.data.frame()
tbl.eqtls$eQTL_effect_dir = gsub("down", "Decrease expression", tbl.eqtls$eQTL_effect_dir)
tbl.eqtls$eQTL_effect_dir = gsub("up", "Increase expression", tbl.eqtls$eQTL_effect_dir)
tbl.eqtls$eQTL_effect_dir = factor(tbl.eqtls$eQTL_effect_dir, levels = c("Increase expression", "Decrease expression"))

p.eqtl.pairs = ggplot(tbl.eqtls, aes(y = forcats::fct_rev(organ), x = n, fill = forcats::fct_rev(eQTL_effect_dir))) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n), hjust = -0.2, position = position_dodge(width = 0.85), size = 2) +
  theme_classic(base_size = 8) +
  labs(y = "", x = "Number of unique gene-variant pairs", fill = "Variant effect on gene") +
  scale_fill_manual(guide = guide_legend(reverse = T), values = c("royalblue", "indianred1")) +
  scale_x_continuous(expand = c(0,0), limits = c(0,30000)) +
  theme(panel.grid.major.x = element_line(color = "lightgrey", linewidth = 0.25),
        panel.grid.minor.x = element_line(color = "lightgrey", linewidth = 0.25),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(figout, "/eQTL_uniq_gene_variant_pairs_organ_distribution.pdf"), p.eqtl.pairs, height = 3, width = 5, units = "in")


#### Observed counts -----------------------------------------------------------
# Work on each organ separately
enrich.results = mclapply(
  c("Adrenal", "Brain", "StomachEsophagus", "Heart", "Liver", "Lung",
    "Muscle", "Skin", "Spleen", "Thyroid"), function(tissue){
      cat(tissue, "\n")
      
      if(tissue == "StomachEsophagus"){
        tissue.grep = "Stomach|Esophagus"
      } else {
        tissue.grep = tissue
      }
      ct.hits.files = motif_hits_list[grep(tissue.grep, motif_hits_list)]
      ct.hits = lapply(ct.hits.files, function(ctfile){
        df = fread(ctfile, data.table = F, col.names = c("chr", "start", "end", "motif", "score", "strand", "type"))
        df = select(df, chr, start, end, motif, type)
        df
      }) %>% bind_rows() %>% distinct()
      ct.hits.pos = ct.hits[ct.hits$type == "pos_patterns",]
      ct.hits.neg = ct.hits[ct.hits$type == "neg_patterns",]
      
      ct.eqtl.files = ct.eqtl.files[grep(tissue.grep, ct.eqtl.files)]
      ct.eqtls = lapply(ct.eqtl.files, function(eqtlfile){
        df = read_parquet(eqtlfile)
        df$direction = ifelse(df$afc < 0, "down", "up")
        df$chr = word(df$variant_id, 1, 1, sep = "_")
        df$pos = as.numeric(word(df$variant_id, 2, 2, sep = "_"))
        df$id = paste0(df$variant_id, "_", df$phenotype_id)
        df = select(df, chr, pos, id, direction, pip, afc)
        df
      }) %>% bind_rows() %>% as.data.frame()
      # For repeated ids, pick the one with the higher SuSiE PIP
      ct.eqtls = ct.eqtls %>%
        group_by(chr, pos, id) %>%
        arrange(desc(pip)) %>%
        slice_head(n = 1)
      
      ct.eqtls.up = ct.eqtls[ct.eqtls$direction == "up",]
      ct.eqtls.down = ct.eqtls[ct.eqtls$direction == "down",]
      
      # Each motif is tested separately
      results.pos.up = count_hits(ct.eqtls.up, ct.hits.pos)
      results.pos.up = add_missing(results.pos.up, ct.hits.pos)
      results.pos.up$test_type = "pos_up"
      results.pos.down = count_hits(ct.eqtls.down, ct.hits.pos)
      results.pos.down = add_missing(results.pos.down, ct.hits.pos)
      results.pos.down$test_type = "pos_down"
      results.neg.up = count_hits(ct.eqtls.up, ct.hits.neg)
      results.neg.up = add_missing(results.neg.up, ct.hits.neg)
      results.neg.up$test_type = "neg_up"
      results.neg.down = count_hits(ct.eqtls.down, ct.hits.neg)
      results.neg.down = add_missing(results.neg.down, ct.hits.neg)
      results.neg.down$test_type = "neg_down"
      
      results = rbind(results.pos.up, results.pos.down, results.neg.up, results.neg.down)
      results$tissue = tissue
      results
    }, mc.cores = 2) %>% bind_rows()
fwrite(enrich.results, "eQTLs_motifsv3_uniq_pairs_with_0_counts.txt", row.names = F, sep= "\t")

# Shuffled counts are in a separate script for qsub
