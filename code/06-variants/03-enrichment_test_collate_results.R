# Collate results for eQTL enrichment test
# Downstream analysis
# Plotting

library(data.table)
library(tidyverse)
library(parallel)
library(here)

out <- here::here("output/06-variants/03/")
figout <- here::here("figures/06-variants/03/")
dir.create(out, showWarnings = F, recursive=T)
dir.create(figout, showWarnings = F, recursive=T)
setwd(out)

#### Collate results -----------------------------------------------------------
enrich.results = fread(here::here("output/06-variants/01/eQTLs_motifsv3_uniq_pairs_with_0_counts.txt"), data.table = F)
enrich.results$shuf.idx = 0 # index is shuffle round. 0 for observed.

shuffle.res.list = list.files(here::here("output/06-variants/02/"), pattern = "results_100k_", full.names = T)
shuffle.results = lapply(shuffle.res.list, function(srl){
  cat(srl, "\n")
  df = fread(srl, data.table = F)
  df
}) %>% bind_rows()
gc()
sr = bind_rows(enrich.results, shuffle.results)

# Save different test type results separately
fwrite(sr[sr$test_type == "pos_up",], "v3_uniq_pairs_shuf_afc_results_100k_pos_up.txt", sep = "\t", row.names = F)
fwrite(sr[sr$test_type == "pos_down",], "v3_uniq_pairs_shuf_afc_results_100k_pos_down.txt", sep = "\t", row.names = F)
fwrite(sr[sr$test_type == "neg_up",], "v3_uniq_pairs_shuf_afc_results_100k_neg_up.txt", sep = "\t", row.names = F)
fwrite(sr[sr$test_type == "neg_down",], "v3_uniq_pairs_shuf_afc_results_100k_neg_down.txt", sep = "\t", row.names = F)

system("gzip -f v3_uniq_pairs_shuf_afc_results_100k_pos_up.txt")
system("gzip -f v3_uniq_pairs_shuf_afc_results_100k_pos_down.txt")
system("gzip -f v3_uniq_pairs_shuf_afc_results_100k_neg_up.txt")
system("gzip -f v3_uniq_pairs_shuf_afc_results_100k_neg_down.txt")

rm(sr, shuffle.results) ; gc()

#### Calculate metrics ---------------------------------------------------------
# Get p values and fold changes
for(testType in c("pos_up", "pos_down", "neg_up", "neg_down")){
  sr = fread(paste0("v3_uniq_pairs_shuf_afc_results_100k_", testType, ".txt.gz"), data.table = F, sep = "\t")
  df.split = split(sr, sr$tissue)
  results = lapply(df.split, function(dfs){
    cat(dfs$tissue[1], "\n")
    dfs = as.data.table(dfs)
    dfs$shuf.idx = paste0("shuf", dfs$shuf.idx)
    dfs.wide = dcast(dfs,
                     tissue + motif + test_type ~ shuf.idx,
                     value.var = "cts",
                     fill = 0) %>%
      dplyr::rename(cts = shuf0) %>%
      as.data.frame()
    
    # p values is number of times shuffled counts is greater than observed counts
    # divided by number of shuffles (100,000)
    dfs.wide = dfs.wide %>%
      rowwise() %>%
      mutate(gt_obs = sum(c_across(starts_with("shuf")) > cts))
    dfs.wide$p = dfs.wide$gt_obs/100000
    dfs.wide$fdr = p.adjust(dfs.wide$p, method = "BH")
    dfs.wide$neglogfdr = -log(dfs.wide$fdr, 10)
    # neglogfdr capped at 5 because we did 100,000 shuffles
    dfs.wide[dfs.wide$neglogfdr == Inf,]$neglogfdr = 5
    dfs = dfs.wide %>%
      dplyr::rename(shuf0 = cts) %>%
      pivot_longer(cols = starts_with("shuf"), names_to = "shuf.idx", values_to = "cts")
    
    dfs.means = dfs %>%
      filter(shuf.idx != "shuf0") %>%
      group_by(motif) %>%
      mutate(mean_shuf = mean(cts)) %>%
      select(tissue, motif, mean_shuf) %>%
      unique()
    dfs = left_join(dfs, dfs.means, by = c("tissue", "motif"))
    
    # Enrichment score is fold change of observed counts over mean of shuffled counts
    # If 0 counts, set to 0 to avoid division by 0
    dfs = dfs %>%
      mutate(fc = ifelse(cts == 0 & mean_shuf == 0, 0, cts/mean_shuf))
    dfs.ci = dfs %>%
      filter(shuf.idx != "shuf0") %>%
      group_by(motif) %>%
      mutate(ci_low_fc = quantile(fc, 0.025),
             ci_high_fc = quantile(fc, 0.975),
             ci_low_cts = quantile(cts, 0.025),
             ci_high_cts = quantile(cts, 0.975)) %>%
      select(tissue, motif, ci_low_fc, ci_high_fc, ci_low_cts,ci_high_cts) %>%
      unique()
    
    dfs = left_join(dfs, dfs.ci, by = c("tissue", "motif"))
    
    # Consider enriched if:
    # FDR < 0.05
    # Observed count is above 95% quantile of shuffled counts
    dfs.enrich = dfs %>%
      filter(shuf.idx == "shuf0") %>%
      group_by(motif) %>%
      mutate(is.enrich.fc = ifelse((fdr < 0.05) & (fc > ci_high_fc) & (ci_high_fc != 0), 1, 0),
             is.enrich.cts = ifelse((fdr < 0.05) & (cts > ci_high_cts) & (ci_high_cts != 0), 1, 0)) %>%
      select(tissue, motif, is.enrich.fc, is.enrich.cts) %>%
      unique()
    dfs = left_join(dfs, dfs.enrich, by = c("tissue", "motif"))
    dfs = dfs[dfs$shuf.idx == "shuf0",]
    dfs$test_type = testType
    dfs
  }) %>% bind_rows()
  results$test_type = testType
  
  fwrite(results, paste0("v3_shuf_afc_100k_with_metrics_", testType, ".txt"),
         sep = "\t", row.names = F)
  system(paste0("gzip -f v3_shuf_afc_100k_with_metrics_", testType, ".txt"))
  
}

#### Fisher test ---------------------------------------------------------------
# By tissue and test type
pos.df = lapply(c("pos_up", "pos_down"), function(testType){
  df = fread(paste0("v3_uniq_pairs_shuf_afc_100k_with_metrics_", testType, ".txt.gz"),
             sep = "\t", data.table = F, nThread = 4)
  df$tissue = gsub("StomachEsophagus", "Stomach", df$tissue)
  df
}) %>% bind_rows()

neg.df = lapply(c("neg_up", "neg_down"), function(testType){
  df = fread(paste0("v3_uniq_pairs_shuf_afc_100k_with_metrics_", testType, ".txt.gz"),
             sep = "\t", data.table = F, nThread = 4)
  df$tissue = gsub("StomachEsophagus", "Stomach", df$tissue)
  df
}) %>% bind_rows()

##### Activating motifs --------------------------------------------------------
pos_up_data = pos.df %>%
  filter(test_type == "pos_up") %>%
  mutate(is.sig = ifelse((fdr < 0.05) & (is.enrich.cts == 1), 1, 0)) %>%
  summarise(sig.enriched = sum(is.sig), total = n())
pos_up_enriched = pos_up_data$sig.enriched
pos_up_total = pos_up_data$total

pos_down_data = pos.df %>%
  filter(test_type == "pos_down") %>%
  mutate(is.sig = ifelse((fdr < 0.05) & (is.enrich.cts == 1), 1, 0)) %>%
  summarise(sig.enriched = sum(is.sig), total = n())
pos_down_enriched = pos_down_data$sig.enriched
pos_down_total = pos_down_data$total

# Contingency table for activating motifs
contingency_table_pos = matrix(c(pos_up_enriched, pos_down_enriched,
                                 pos_up_total - pos_up_enriched, pos_down_total - pos_down_enriched),
                               nrow = 2,
                               dimnames = list(c("Pos_up", "Pos_down"),
                                               c("Enriched", "Not_Enriched")))

#          Enriched Not_Enriched
# Pos_up          2         3057
# Pos_down       42         3017

fisher.test(contingency_table_pos)
# p-value = 9.906e-11
# odds ratio 0.04700608

##### Reducing motifs ----------------------------------------------------------
neg_up_data = neg.df %>%
  filter(test_type == "neg_up") %>%
  mutate(is.sig = ifelse((fdr < 0.05) & (is.enrich.cts == 1), 1, 0)) %>%
  summarise(sig.enriched = sum(is.sig), total = n())
neg_up_enriched = neg_up_data$sig.enriched
neg_up_total = neg_up_data$total

neg_down_data = neg.df %>%
  filter(test_type == "neg_down") %>%
  mutate(is.sig = ifelse((fdr < 0.05) & (is.enrich.cts == 1), 1, 0)) %>%
  summarise(sig.enriched = sum(is.sig), total = n())
neg_down_enriched = neg_down_data$sig.enriched
neg_down_total = neg_down_data$total

# Contingency table for reducing motifs
contingency_table_neg = matrix(c(neg_up_enriched, neg_down_enriched,
                                 neg_up_total - neg_up_enriched, neg_down_total - neg_down_enriched),
                               nrow = 2,
                               dimnames = list(c("neg_up", "neg_down"),
                                               c("Enriched", "Not_Enriched")))

#          Enriched Not_Enriched
# neg_up         11          133
# neg_down        0          144

fisher.test(contingency_table_neg)
# p-value = 0.0008009
# odds ratio = Inf

#### Plot ----------------------------------------------------------------------
# Colour schemes
load(here::here("code/utils/color_scheme.RData"))

cmap_category = c("base"            = "firebrick1",
                  "base_with_flank" = "firebrick4",
                  "homocomposite"   = "darkorchid4",
                  "heterocomposite" = "royalblue3",
                  "unresolved"      = "black",
                  "repeat"          = "gray30",
                  "partial"         = "gray50",
                  "exclude"         = "gray90")

all_motifs = fread(here::here("output/03-chrombpnet/03-syntax/01/motifs_compiled_unique_short.tsv"), header = T,
                   data.table = F, sep = "\t")[,c("motif_name", "idx_uniq", "pattern_class",
                                                  "annotation", "annotation_broad", "category")]

##### Prepare for plotting -----------------------------------------------------
plt.df = lapply(c("pos_up", "pos_down", "neg_up", "neg_down"), function(testType){
  df = fread(paste0("v3_uniq_pairs_shuf_afc_100k_with_metrics_", testType, ".txt.gz"),
             sep = "\t", data.table = F, nThread = 4)
  df = df[(df$is.enrich.cts == 1) & (df$fdr < 0.05),]
  df$tissue = gsub("StomachEsophagus", "Stomach", df$tissue)
  df
}) %>% bind_rows()
plt.df = left_join(plt.df, all_motifs, by = c("motif"="motif_name"))

# Group by broad motif annotations for plotting
plt.df2 = plt.df %>%
  group_by(tissue, test_type, annotation_broad) %>%
  mutate(total_fc = sum(fc)) %>%
  ungroup() %>%
  group_by(test_type, annotation_broad) %>%
  mutate(relative_score = total_fc/sum(total_fc)) %>%
  as.data.frame()
# For unresolved motifs, the annotation_broad is replaced with annotation to avoid
# having all of them being called "unresolved" since many of them are 
# very different from each other
plt.df2[plt.df2$annotation_broad == "unresolved",]$annotation_broad = plt.df2[plt.df2$annotation_broad == "unresolved",]$annotation

# Arrange motifs by decreasing entropy of relative enrichment scores
calculate_entropy = function(scores) {
  -sum(scores * log(scores), na.rm = T)
}

entropy_df = plt.df2 %>%
  group_by(test_type, annotation_broad) %>%
  summarise(entropy = calculate_entropy(relative_score)) %>%
  ungroup()

plt.df2 = left_join(plt.df2, entropy_df, by = c("test_type", "annotation_broad"))

plt.df3 = plt.df2 %>%
  mutate(split_fc = ifelse(test_type %in% c("pos_up", "neg_up"), fc, -fc)) %>%
  select(tissue, motif, test_type, neglogfdr, fc, pattern_class, annotation_broad,
         category, entropy, split_fc, relative_score)
plt.df3$eqtl_dir = ifelse(plt.df3$split_fc < 0,
                          "Enrichment in\nnegative eQTL variants",
                          "Enrichment in\npositive eQTL variants")

##### Save results for sup table -----------------------------------------------
sup.table = plt.df2 %>%
  mutate(split_fc = ifelse(test_type %in% c("pos_up", "neg_up"), fc, -fc))
sup.table$eqtl_dir = ifelse(sup.table$split_fc < 0,
                            "Enrichment in\nnegative eQTL variants",
                            "Enrichment in\npositive eQTL variants")
sup.table = sup.table %>%
  select(organ = tissue, motif, motif_group = annotation_broad, pattern_class, category, 
         observed_counts = cts, mean_shuf_counts = mean_shuf, enrichment_score = fc,
         n_shuf_gt_obs = gt_obs, p, neglogfdr, enriched_in = eqtl_dir)
sup.table$enriched_in = gsub("Enrichment in\n", "", sup.table$enriched_in)
sup.table$enriched_in = gsub(" ", "_", sup.table$enriched_in)
sup.table = sup.table %>%
  arrange(enriched_in, desc(neglogfdr), desc(enrichment_score), organ, motif_group, motif)
sup.table$pattern_class = gsub("pos_patterns", "activating", sup.table$pattern_class)
sup.table$pattern_class = gsub("neg_patterns", "repressive", sup.table$pattern_class)
fwrite(sup.table, "eQTL_uniq_gene_variant_pairs_enriched_fdr0.05.tsv", sep = "\t", row.names = F)

##### Enrichment scores --------------------------------------------------------
enrich.up = filter(plt.df3, eqtl_dir == "Enrichment in\nnegative eQTL variants")
missing.up = plt.df3[!plt.df3$motif %in% enrich.up$motif,][,c("motif", "pattern_class", "annotation_broad", "entropy", "category", "relative_score")]
if(length(missing.up) != 0){
  enrich.up = bind_rows(enrich.up, missing.up)
}
# Sort
enrich.up = enrich.up %>%
  arrange(desc(pattern_class), relative_score)
enrich.up$annotation_broad = factor(enrich.up$annotation_broad, levels = rev(unique(enrich.up$annotation_broad)))

enrich.down = filter(plt.df3, eqtl_dir == "Enrichment in\npositive eQTL variants")
missing.down = plt.df3[!plt.df3$motif %in% enrich.down$motif,][,c("motif", "pattern_class", "annotation_broad", "entropy", "category", "relative_score")]
if(length(missing.down) != 0){
  enrich.down = bind_rows(enrich.down, missing.down)
}
# Sort for consistency
enrich.down$annotation_broad = factor(enrich.down$annotation_broad, levels = rev(unique(enrich.up$annotation_broad)))

p1 = ggplot(enrich.up,
            aes(x = split_fc, y = annotation_broad, color = neglogfdr)) +
  geom_jitter(width = 0, height = 0.3, size = 1.5, show.legend = F) +
  geom_hline(yintercept = seq(1.5, length(unique(enrich.up$annotation_broad))-0.5, 1),
             lwd = 0.2, color = "black") +
  geom_hline(yintercept = 5.5, lwd = 1, color = "black") + # To visually separate activating and repressive motifs
  scale_x_continuous(breaks = c(-2,-1.5,-1), limits = c(-2.15, -0.95), labels = c(2,1.5,1), expand = c(0,0)) +
  scale_color_gradient(low = "navy", high = "gold",
                       na.value = "white", limits = c(0,5)) +
  scale_y_discrete(position = "right") +
  labs(title = "Enriched in\nnegative eQTL variants", y = "", x = "Enrichment\nscore") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    panel.border = element_rect(linewidth = 1, color = "black", fill = "transparent"),
    panel.grid.major.x = element_line(linewidth = 0.3, color = "lightgrey"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 0, b = 0, r = -10, l = -10)
  )

p5 = ggplot(enrich.down,
            aes(x = split_fc, y = annotation_broad, color = neglogfdr)) +
  geom_jitter(width = 0, height = 0.3, size = 1.5, show.legend = F) +
  scale_x_continuous(breaks = c(1,1.5,2), limits = c(0.95,2.15), labels = c(1,1.5,2), expand = c(0,0)) +
  geom_hline(yintercept = seq(1.5, length(unique(enrich.up$annotation_broad))-0.5, 1),
             lwd = 0.2, color = "black") +
  geom_hline(yintercept = 5.5, lwd = 1, color = "black") +
  scale_color_gradient(low = "navy", high = "gold",
                       na.value = "white", limits = c(0,5)) +
  labs(title = "Enriched in\npositive eQTL variants", y = "", x = "Enrichment\nscore") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    panel.border = element_rect(linewidth = 1, color = "black", fill = "transparent"),
    panel.grid.major.x = element_line(linewidth = 0.3, color = "lightgrey"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 0, b = 0, r = -10, l = -10)
  )

##### Organ proportions --------------------------------------------------------
# Proportion of the significant hits by organ
annotation_order = rev(unique(enrich.up$annotation_broad))

score_df = plt.df3 %>%
  dplyr::select(tissue, annotation_broad, test_type, relative_score) %>%
  distinct()

eQTL_up = score_df[score_df$test_type %in% c("pos_up", "neg_up"),]
motifs.missing.in.up = score_df[!score_df$annotation_broad %in% eQTL_up$annotation_broad,]
if(nrow(motifs.missing.in.up) != 0){
  motifs.missing.in.up$relative_score = 0
  eQTL_up = rbind(eQTL_up, motifs.missing.in.up)
}
eQTL_up$annotation_broad = factor(eQTL_up$annotation_broad, levels = annotation_order)

eQTL_down = score_df[score_df$test_type %in% c("pos_down", "neg_down"),]
motifs.missing.in.down = score_df[!score_df$annotation_broad %in% eQTL_down$annotation_broad,]
if(nrow(motifs.missing.in.down) != 0){
  motifs.missing.in.down$relative_score = 0
  eQTL_down = rbind(eQTL_down, motifs.missing.in.down)
}
eQTL_down$annotation_broad = factor(eQTL_down$annotation_broad, levels = annotation_order)

p2 = ggplot(eQTL_down, aes(x = relative_score, y = annotation_broad, fill = tissue)) +
  scale_fill_manual(values = cmap_organ) +
  geom_col(width = 0.9, position = "fill", show.legend = F) +
  labs(x = "% Organ\ndistribution", y = "", fill = "Organ") +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,1)) +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "mm"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 0, b = 0, r = -50, l = -50))

p4 = ggplot(eQTL_up, aes(x = relative_score, y = annotation_broad, fill = tissue)) +
  scale_fill_manual(values = cmap_organ) +
  geom_col(width = 0.9, position = "fill", show.legend = F) +
  labs(x = "% Organ\ndistribution", y = "", fill = "Organ") +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,1)) +
  scale_y_discrete(position = "right") +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "mm"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 0, b = 0, r = -50, l = -50))

# Broad annotation names
names.df = plt.df3 %>%
  select(annotation_broad, pattern_class) %>%
  distinct()
names.df$value = 1
names.df$annotation_broad = factor(names.df$annotation_broad, levels = annotation_order)

p3 = ggplot(names.df, aes(x = value, y = annotation_broad, fill = pattern_class)) +
  geom_col(width = 0.9, show.legend = F) +
  geom_text(aes(x = value / 2, label = annotation_broad), hjust = 0.5, size = 3.5) +
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_y_discrete(position = "right") +
  theme_bw() +
  scale_fill_manual(values = c("#FF999A", "#99C199")) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, b = 0, r = 0, l = 2)
  )

# Legends
p.legend = ggplot(enrich.up, aes(x = split_fc, y = annotation_broad, color = neglogfdr)) +
  geom_jitter() +
  labs(color = "-log(FDR)") +
  scale_color_gradient(low = "navy", high = "gold",
                       na.value = "white", limits = c(0,5))

legend.fdr = cowplot::get_legend(p.legend)

p.organ = ggplot(rbind(eQTL_down, eQTL_up), aes(x = relative_score, y = annotation_broad, fill = tissue)) +
  scale_fill_manual(values = cmap_organ) +
  geom_col(position = "fill") +
  labs(fill = "Organ")
legend.organ = cowplot::get_legend(p.organ)

p6 = cowplot::plot_grid(legend.fdr, legend.organ, ncol = 1, align = "hv", axis = "trbl")

layout = "AAABCCCDEEE"
svg(paste0(figout, "/eQTL_enrichment_uniq_pairs.svg"), height = 7, width = 8.2)
p = p1+p2+p3+p4+p5+plot_layout(design = layout)
dev.off()

ggsave(filename = paste0(figout, "/eQTL_enrichment_legend.pdf"),
       p6, height = 6, width = 2, units = "in")


#### Variants in each motif ----------------------------------------------------
# What variants hit each motif?
search_df = plt.df2 %>%
  select(tissue, motif, test_type) %>%
  distinct()

lapply(1:nrow(search_df), function(i){
  motifname =  gsub("|", "_", search_df$motif[i], fixed = T)
  motifname =  gsub(":", "_", motifname, fixed = T)
  motifname =  gsub(",", "_", motifname, fixed = T)
  motifname =  gsub("#", "_", motifname, fixed = T)
  motifname =  gsub("/", "_", motifname, fixed = T)
  cat(i, motifname, "\n")
  outfilename = paste0("eQTL_genes_hit/", motifname, "_", search_df$tissue[i], "_", search_df$test_type[i], ".tsv")
  pattern_type = paste0(word(search_df$test_type[i], 1, 1, sep = "_"), "_patterns")
  eqtl_dir = word(search_df$test_type[i], 2, 2, sep = "_")
  find_hit_genes(search_df$tissue[i], search_df$motif[i], pattern_type, eqtl_dir) %>%
    fwrite(outfilename, sep = "\t", row.names = F)
})
rm(search_df)

# Collate
filelist = list.files("eQTL_genes_hit", full.names = T)
n_gene_pairs = lapply(filelist, function(f){
  cat(f, "\n")
  test_type = gsub(".tsv", "", word(f, -2, -1, sep = "_"))
  df = fread(f, sep = "\t", data.table = F)
  df$test_type = test_type
  df
}) %>% bind_rows()

# Save as sup table
n_gene_pairs = left_join(n_gene_pairs, distinct(plt.df2[,c("motif", "pattern_class", "annotation_broad")]))
n_gene_pairs = n_gene_pairs %>%
  select(organ, eQTL_direction = test_type, motif_pattern_class = pattern_class, motif, motif_group = annotation_broad,
         variant, gene_ID)
n_gene_pairs$motif_pattern_class = gsub("pos_patterns", "activating", n_gene_pairs$motif_pattern_class)
n_gene_pairs$motif_pattern_class = gsub("neg_patterns", "repressive", n_gene_pairs$motif_pattern_class)
n_gene_pairs$eQTL_direction = gsub("pos_up", "positive", n_gene_pairs$eQTL_direction)
n_gene_pairs$eQTL_direction = gsub("pos_down", "negative", n_gene_pairs$eQTL_direction)
n_gene_pairs$eQTL_direction = gsub("neg_up", "positive", n_gene_pairs$eQTL_direction)
n_gene_pairs$eQTL_direction = gsub("neg_down", "negative", n_gene_pairs$eQTL_direction)
n_gene_pairs = n_gene_pairs %>%
  arrange(organ, motif_pattern_class, eQTL_direction, motif)
fwrite(n_gene_pairs, "eQTL_enriched_variant_gene_pairs.tsv", sep = "\t", row.names = F)

##### Plot for sup figure ------------------------------------------------------
plt_n_gene_pairs = n_gene_pairs %>%
  mutate(var_gene_pair = paste0(variant, "_", gene_ID)) %>%
  select(-organ, -variant, -gene_ID, -eQTL_direction, -motif_pattern_class) %>%
  distinct()
tbl.gene.pairs = janitor::tabyl(plt_n_gene_pairs, motif)
colnames(tbl.gene.pairs) = c("motif", "n_uniq_gene_vars", "percent")
tbl.gene.pairs$percent = NULL

plt.df4 = plt.df2
plt.df4 = plt.df4 %>%
  group_by(motif) %>%
  mutate(totalCts = sum(cts)) %>%
  ungroup() %>%
  as.data.frame()
plt.df4 = left_join(plt.df4, tbl.gene.pairs, by = "motif")
plt.df4$label = paste0(plt.df4$totalCts, " (", plt.df4$n_uniq_gene_vars, ")")
# Sort for consistency
plt.df4$annotation_broad = factor(plt.df4$annotation_broad, levels = rev(unique(enrich.up$annotation_broad)))
plt.df4 = arrange(plt.df4, annotation_broad, entropy)
plt.df4$motif = factor(plt.df4$motif, levels = unique(plt.df4$motif))
label.df = plt.df4 %>%
  select(motif, label, totalCts) %>%
  distinct() %>%
  dplyr::rename(cts = totalCts)
label.df$motif = factor(label.df$motif, levels = unique(label.df$motif))

p.var.motif = ggplot() +
  geom_bar(data = plt.df4, aes(x = cts, y = motif, fill = tissue),
           position = "stack", stat = "identity") +
  geom_text(data = label.df, aes(x = cts, y = motif, label = label),
            hjust = -0.1, size = 2) +
  scale_fill_manual(values = cmap_organ) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,1000,100),
                     limits = c(0, 1000)) +
  labs(title = "Motif instances significantly enriched with eQTL variants",
       fill = "Organ", y = "", x = "Number of eQTL variants in motif instances") +
  theme_bw(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_blank()
  )
ggsave(paste0(figout, "/num_eQTL_vars_in_motifs.pdf"), p.var.motif, height = 4, width = 6, units = "in")
