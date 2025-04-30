# Plot variant scoring results for all fetal-only hits
library(data.table)
library(BPCells)
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
library(rtracklayer)
library(cowplot)
library(trqwe)
library(glue)
library(here)


scriptPath <- here::here("code/utils")
source(paste0(scriptPath, "/track_helpers_BL.R"))
source(paste0(scriptPath, "/track_helpers_SJ.R"))

out <- here::here("output/06-variants/07a/")
figout <- here::here("figures/06-variants/07a/")
dir.create(out, showWarnings = F, recursive=T)
dir.create(figout, showWarnings = F, recursive=T)

global_bp_obj = readRDS(here::here("output/05-misc/01/global_bp_obj.rds"))
global_bp_obj$frags@dir = here::here("output/05-misc/01/ATAC_merged")
global_bp_obj$rna@dir = here::here("output/05-misc/01/RNA_merged")

hits = fread(here::here("output/06-variants/05/v3motifs_w_causal_snps_all_tissues_disease_only_pipsum_5_fdr_0.05_fetalAdult.txt"),
             data.table = F)
hits = hits[hits$foundInMin2Adults == "fetalOnly",]
# Add trait names
meta_trait = fread(here::here("data/external/variants/standardized_trait_list.txt"),
                   data.table = F)
hits = merge(hits, meta_trait[,c("meta_id", "newname")])
hits$traitname = word(hits$newname, 3, -1, sep = "_")
hits$traitname = word(hits$traitname, 1, 1, sep = "\\.")
hits$newname = NULL

organ2abbr = data.frame(organ = c("Adrenal", "Brain", "Eye", "Heart", "Liver", "Lung", "Muscle",
                                  "Skin", "Spleen", "Stomach", "Thyroid", "Thymus"),
                        abbrev = c("AG", "BR", "EY", "HT", "LI", "LU", "MU", "SK", "SP", "ST", "TR", "TM"))


# Representative set of hits
# For each MeSH ID and rs ID in a cell type, pick one representative entry to avoid redundant results
# Sorted by decreasing FDR then by decreasing gchromvar Z score
# GD03523	D009203	Myocardial infarction and
# AT000	D012140	J40 J47 Chronic lower respiratory diseases
# are also dropped for being similar to GCST90038616 D001249	Asthma and GCST90132314	D003324	Coronary artery disease
variant.list = fread(here::here("output/06-variants/05/fetal_only_variant_scoring_representative.txt"),
                     data.table = F)
variant.list$cluster = paste0(variant.list$organ, "_c", variant.list$selected_L1_cluster)
variant.list$traitname = paste0(variant.list$meta_id, "_", variant.list$mesh_id, "_", variant.list$trait,
                                "_", variant.list$rsid, "_", variant.list$cluster)
variant.list$traitname = gsub(" ", "_", variant.list$traitname)

file.list = list.files(here::here("output/06-variants/06/"),
                       pattern = "peakcentered.tsv", full.names = T)
file.list = data.frame(path = file.list)
file.list$traitname = word(basename(file.list$path), 1, -2, sep = "_")
file.list = file.list[file.list$traitname %in% variant.list$traitname,]

indiv_theme = theme(plot.margin = unit(c(0, 0, 0, 0), "pt"),
                    legend.position = "none")

lapply(1:nrow(file.list), function(i){
  cat(i, "\n")
  fl = file.list$path[i]
  filename = gsub(".tsv", "", basename(fl))
  trait = word(basename(fl), 1, -5, sep = "_")
  metaid = word(trait, 1, 1, sep = "_")
  meshid = word(trait, 2, 2, sep = "_")
  RSID = word(basename(fl), -4, -4, sep = "_")
  organ = word(basename(fl), -3, -3, sep = "_")
  accpred = fread(fl, data.table = F)
  
  plotTitle = paste0(
    gsub("_", " ", accpred$trait), "\n",
    accpred$variant, "\n",
    organ, " ", gsub(".", " ", accpred$celltype, fixed = T), " (cluster ", gsub("c", "", word(accpred$cluster, -1, -1, sep = "_")), ")\n\n",
    "Mean ALT/REF log2FC: ", accpred$l2fc_mean, "\n", "Mean ALT/REF diff in local counts: ", accpred$diff_local_mean
  )
  
  varscore_bed = rtracklayer::import.bed(paste0(here::here("output/06-variants/06/"),
                                                "/", filename, ".bed"),
                                         extraCols = c(
                                           "pred_ref" = "numeric",
                                           "pred_alt" = "numeric",
                                           "contrib_sum_ref" = "numeric",
                                           "contrib_sum_alt" = "numeric",
                                           "contrib_ref_A" = "numeric",
                                           "contrib_ref_C" = "numeric",
                                           "contrib_ref_G" = "numeric",
                                           "contrib_ref_T" = "numeric",
                                           "contrib_alt_A" = "numeric",
                                           "contrib_alt_C" = "numeric",
                                           "contrib_alt_G" = "numeric",
                                           "contrib_alt_T" = "numeric"
                                         ))
  
  variant.df = variant.list %>%
    filter(meta_id == metaid & mesh_id == meshid & organ == organ & rsid == RSID & cluster == accpred$cluster)
  variant.chr = variant.df$snpChr.hg38
  variant.pos = as.numeric(variant.df$snpPos.hg38)
  
  rangeStart = max(accpred$peak_start, variant.pos-200)
  rangeEnd = min(accpred$peak_end, variant.pos+200)
  
  if((variant.pos+30) > accpred$peak_end){
    highlightStartSmall = accpred$peak_end-61
    highlightEndSmall = accpred$peak_end
  } else {
    highlightStartSmall = variant.pos-30
    highlightEndSmall = variant.pos+30
  }
  
  region = paste0(variant.chr, ":", rangeStart, "-", rangeEnd)
  region_zoom = paste0(variant.chr, ":", highlightStartSmall, "-", highlightEndSmall)
  
  preds_gr = varscore_bed %>% 
    gr_select(c("pred_ref", "pred_alt"), c("score1", "score2"))
  
  # Accessibility predictions ----
  p.accpred = trackplot_bw_paired(bw = preds_gr, region = region, alpha = 0.8,
                                  facet_label = "", track_label = "",
                                  score_cmap = c("ref" = "darkgrey", "alt" = "red3")) +
    highlight_region(region_zoom, color = "skyblue", alpha = 0.2) +
    ggtitle(plotTitle) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, hjust = 0.5))
  
  # REF and ALT contributions ----
  contr_gr = varscore_bed %>% 
    gr_select(c("contrib_sum_ref", "contrib_sum_alt"), c("score1", "score2"))
  
  regionzoom_gr = str_to_gr(region_zoom)
  contr_gr = contr_gr[contr_gr %over% regionzoom_gr]
  
  ymax = max(c(contr_gr$score1, contr_gr$score2)) * 1.05
  ymin = min(c(contr_gr$score1, contr_gr$score2)) * 1.05
  ymax_accuracy = 10^as.integer(log10(0.01 * abs(ymax)))
  ymin_accuracy = 10^as.integer(log10(0.01 * abs(ymin)))
  
  range_label = sprintf("[%s-%s]", 
                        scales::label_comma(accuracy = ymin_accuracy, big.mark=" ")(ymin),
                        scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))
  
  contrib_ref = varscore_bed %>%
    gr_select(c("contrib_ref_A", "contrib_ref_C", "contrib_ref_G", "contrib_ref_T"),
              c("A", "C", "G", "T"))
  
  contrib_alt = varscore_bed %>%
    gr_select(c("contrib_alt_A", "contrib_alt_C", "contrib_alt_G", "contrib_alt_T"),
              c("A", "C", "G", "T"))
  
  if((variant.pos+30) > accpred$peak_end){
    # In case peak ends early
    p.ref = trackplot_contribs2(contrib_ref, region_zoom, facet_label = "", track_label = NULL) +
      scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
      annotate("text", x = 1, y = ymax, label = range_label, vjust = 1.5, hjust = -0.1, size = 11*.8/ggplot2::.pt) +
      highlight_relative_region(44.5, 45.5, color = "gold", alpha = 0.25)
    p.alt = trackplot_contribs2(contrib_alt, region_zoom, facet_label = "", track_label = NULL) +
      scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
      annotate("text", x = 1, y = ymax, label = range_label, vjust = 1.5, hjust = -0.1, size = 11*.8/ggplot2::.pt) +
      highlight_relative_region(44.5, 45.5, color = "gold", alpha = 0.25)
  } else {
    p.ref = trackplot_contribs2(contrib_ref, region_zoom, facet_label = "", track_label = NULL) +
      scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
      annotate("text", x = 1, y = ymax, label = range_label, vjust = 1.5, hjust = -0.1, size = 11*.8/ggplot2::.pt) +
      highlight_relative_region(30.5, 31.5, color = "gold", alpha = 0.25)
    p.alt = trackplot_contribs2(contrib_alt, region_zoom, facet_label = "", track_label = NULL) +
      scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
      annotate("text", x = 1, y = ymax, label = range_label, vjust = 1.5, hjust = -0.1, size = 11*.8/ggplot2::.pt) +
      highlight_relative_region(30.5, 31.5, color = "gold", alpha = 0.25)
  }
  
  
  # Motifs ----
  sub.hits = hits %>%
    filter(rsid == RSID, organ == organ, meta_id == metaid, mesh_id == meshid)
  celltype_name = gsub(".", " ", sub.hits$celltype, fixed = T)
  
  abbr = organ2abbr[organ2abbr$organ == organ,]$abbrev
  submeta = global_bp_obj$cell_metadata %>% 
    dplyr::filter(L2_clusterName %in% c(celltype_name)) %>%
    dplyr::filter(organ == abbr)
  
  hits_input = unique(global_bp_obj$hits[[accpred$cluster]])
  
  p.motifs = trackplot_helper_v2c(
    region = region_zoom,
    clusters = submeta$L2_clusterID, 
    transcripts = global_bp_obj$transcripts, 
    fragments = global_bp_obj$frags %>% select_cells(submeta$cb), 
    cell_read_counts = submeta$nFrags,
    annot = hits_input,
    annot_color = "pattern_class",
    annot_labelby = "name",
    annot_strand = T,
    annot_label = ""
  ) + theme(legend.position = "None",
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  p.mot = p.motifs[[4]] +
    scale_color_manual(values = rep("#A6A8AB",2))
  
  # Zoom gap
  gap1 = data.frame(pos = c(rangeStart, seq(highlightStartSmall, highlightEndSmall, by = 5), rangeEnd),
                    score = 0)
  gap1[gap1$pos >= highlightStartSmall & gap1$pos <= highlightEndSmall,]$score = 1
  
  p.gap1 = ggplot(gap1, aes(x = pos, y = score)) +
    geom_area(fill = "skyblue", alpha = 0.25, show.legend = F) +
    scale_x_continuous(expand = c(0, 0), limits = c(rangeStart, rangeEnd)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
  
  p.out = trackplot_combine(list(
    BPCells:::wrap_trackplot(p.accpred, unit(0.33, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(p.gap1, unit(0.09, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(p.ref, unit(0.17, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(p.alt, unit(0.17, "null")) + indiv_theme,
    BPCells:::wrap_trackplot(p.mot, unit(0.2, "null")) + indiv_theme))
  
  cat(plotTitle, "\n")
  
  if(i < 10){
    num = paste0(0, i)
  } else {
    num = i
  }
  
  pdf(paste0(figout, "/", num, "_variant_scoring_", filename, ".pdf"),
      height = 3.5, width = 5.2)
  print(p.out)
  dev.off()
})

