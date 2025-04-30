# Plot multiple tracks with BPCells
library(ArchR)
library(BPCells)
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
library(BSgenome.Hsapiens.UCSC.hg38)
library(cowplot)
library(ggpubr)
library(rtracklayer)
library(ggrepel)
library(ggcoverage)
options(scipen = 999)
library(here)

scriptPath <- here::here("code/utils")
source(paste0(scriptPath, "/track_helpers_BL.R"))
source(paste0(scriptPath, "/track_helpers_SJ.R"))
load(here::here("code/utils/color_scheme.RData"))

out <- here::here("output/06-variants/07b/")
figout <- here::here("figures/06-variants/07b/")
dir.create(out, showWarnings = F, recursive=T)
dir.create(figout, showWarnings = F, recursive=T)
setwd(out)

rsid = "rs12740374"
gwas_id = "GCST90132314"

## Organ name match
tissue2abbrev = data.frame(tissue = c("Adrenal", "Brain", "Eye", "Heart", "Liver",
                                      "Lung", "Muscle", "Skin", "Spleen",
                                      "StomachEsophagus", "Thymus", "Thyroid"),
                           abbrev = c("AG", "BR", "EY", "HT", "LI", "LU", "MU",
                                      "SK", "SP", "ST", "TM", "TR"))

## Hits
hits = fread(here::here("output/06-variants/05/v3motifs_w_causal_snps_all_tissues_disease_only_pipsum_5_fdr_0.05_fetalAdult.txt"),
             data.table = F)
# Fetal only entries
hits = hits[hits$foundInMin2Adults == "fetalOnly",]
# Add trait names
meta_trait = fread(here::here("data/external/variants/standardized_trait_list.txt"),
                   data.table = F)
hits = merge(hits, meta_trait[,c("meta_id", "newname")])
hits$traitname = word(hits$newname, 3, -1, sep = "_")
hits$traitname = word(hits$traitname, 1, 1, sep = "\\.")
hits$newname = NULL

global_bp_obj = readRDS(here::here("output/05-misc/01/global_bp_obj.rds"))
global_bp_obj$frags@dir = here::here("output/05-misc/01/ATAC_merged")
global_bp_obj$rna@dir = here::here("output/05-misc/01/RNA_merged")

# Subset for the hit we want to plot
sub.hits = hits[hits$rsid == rsid,]
tis = sub.hits$organ[1]
organ_color = cmap_organ[names(cmap_organ)==tis][[1]]
celltype_name = sub.hits$celltype[1]
celltype_name = gsub(".", " ", celltype_name, fixed = T)
submeta = global_bp_obj$cell_metadata %>% 
  dplyr::filter(L2_clusterName %in% c(celltype_name)) %>%
  dplyr::filter(organ == "MU")
# Only plot peaks from MU_endothelial 2 (Muscle c6), which was used for variant scoring
submeta = submeta[submeta$L3_clusterID == "MU_endothelial 2",]

# First show genes and GWAS SNPs
minStart = sub.hits$snpPos.hg38[1]-50000
maxEnd = sub.hits$snpPos.hg38[1]+25000

highlightStart = sub.hits$snpPos.hg38[1]-10000
highlightEnd = sub.hits$snpPos.hg38[1]+10000


#### Genes ####
p.acc = trackplot_helper_v2c(
  region = paste0(sub.hits$snpChr.hg38[1], ":", minStart, "-", maxEnd),
  highlights = gr_to_highlight(GRanges(paste0(sub.hits$snpChr.hg38[1], ":", highlightStart, "-", highlightEnd))),
  clusters = submeta$L3_clusterID, 
  fragments = global_bp_obj$frags %>% select_cells(submeta$cb), 
  cell_read_counts = submeta$nFrags, 
  transcripts = global_bp_obj$transcripts
) + theme(legend.position = "None")

p.sbar = p.acc[[1]]

large_highlight_rect = as.data.frame(gr_to_highlight(GRanges(paste0(sub.hits$snpChr.hg38[1], ":", highlightStart, "-", highlightEnd))))

p.annot = p.acc[[3]] +
  theme(strip.text.y.left = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

#### GWAS SNPs + highlight causal SNP ####
## For lifting GWAS SNPS
ch = import.chain(here::here("data/external/variants/hg19ToHg38.over.chain"))

getGWAS = function(mi){
  gwasfiles = list.files(here::here("data/external/variants/20240710_sumstats_standardized/sumStats_input"), full.names = T)
  gwasfiles = gwasfiles[grep(mi, gwasfiles)]
  gwas = fread(gwasfiles, data.table = F, sep = "\t", nThread = 8)
  gwas = gwas %>%
    dplyr::select(chrom = CHR, start = BP, P)
  gwas$score = -log(gwas$P, 10)
  gwas$end = gwas$start
  gwas = gwas %>% 
    select(chrom, start, end, score)
  gwas$chrom = paste0("chr", gwas$chrom)
  gwas.gr = makeGRangesFromDataFrame(gwas, keep.extra.columns = T)
  gwas.gr.hg38 = unlist(liftOver(gwas.gr, ch))
  gwas = as.data.frame(gwas.gr.hg38)
  gwas = gwas %>%
    filter(seqnames == sub.hits$snpChr.hg38[1] & start >= minStart & end <= maxEnd) %>%
    select(-width, -strand)
  gwas$tohighlight = "N"
  causal_snp  = filter(gwas, start >= highlightStart & end <= highlightEnd) %>%
    filter(score == max(score))
  gwas[gwas$start == causal_snp$start & causal_snp$end == causal_snp$end & gwas$score == causal_snp$score,]$tohighlight = "Y"
  gwas$meta_id = mi
  return(gwas)
}

gwas.snps = getGWAS(gwas_id)
maxlimit = (ceiling(max(gwas.snps$score) / 5) * 5) + 0.5
gwas.snps$label = rsid

gc()

## Get causaldb info table
cdb_stats = fread(here::here("data/external/variants/20240508_causaldb_wMeta.txt.gz"),
                  data.table = F, nThread = 16)
snps = cdb_stats[cdb_stats$meta_id == gwas_id,]
snps = snps %>%
  select(chr, start = bp, score = susie, rsid) %>%
  mutate(end = start)
snps$chr = paste0("chr", snps$chr)
snps.gr = makeGRangesFromDataFrame(snps, keep.extra.columns = T)
snps.gr.hg38 = unlist(liftOver(snps.gr, ch))
rm(snps.gr) ; gc()
snps = as.data.frame(snps.gr.hg38)
snps = snps %>%
  select(seqnames, start, end, rsid) %>%
  mutate(is.causal = "Y") %>%
  distinct()

gwas.snps = left_join(gwas.snps, snps, by = c("seqnames", "start", "end"))
gwas.snps[is.na(gwas.snps$is.causal),]$is.causal = "N"

p.gwas = ggplot() +
  geom_rect(data = large_highlight_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "skyblue", alpha = 0.2) +
  geom_point(data = gwas.snps, aes(x = start, y = score), color = "grey50",
             shape = 19, size = 0.8, show.legend = F, alpha = 0.5) +
  geom_point(data = gwas.snps[gwas.snps$tohighlight == "Y",], aes(x = start, y = score), color = "orangered1",
             shape = 19, size = 0.8, show.legend = F) +
  scale_x_continuous(expand = c(0, 0), limits = c(minStart, maxEnd)) + 
  ylim(0, maxlimit) +
  labs(y = bquote(-log[10](p)), x = "") +
  BPCells:::trackplot_theme() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 11, angle = 90, vjust = 0.5))


#### First zoom gap ####
# Between GWAS and peaks
gap1 = data.frame(pos = c(minStart, seq(highlightStart, highlightEnd, by = 5), maxEnd),
                  score = 0)
gap1[gap1$pos >= highlightStart & gap1$pos <= highlightEnd,]$score = 1

p.gap1 = ggplot(gap1, aes(x = pos, y = score)) +
  geom_area(fill = "skyblue", alpha = 0.2, show.legend = F) +
  scale_x_continuous(expand = c(0, 0), limits = c(minStart, maxEnd)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

#### Predicted accessibility ####
# Upon installation of the variant, what are the effects?

# 0-based coordinates
accpred = fread(paste0("/illumina/scratch/deep_learning/yng7/gChromvar/20250407-fetal_only_variant_scoring/",
                       "GCST90132314_D003324_Coronary_artery_disease_rs12740374_Muscle_c6_peakcentered.tsv"), 
                data.table = F, sep = "\t")
print(accpred$l2fc_mean) # 0.244
print(accpred$diff_local_mean) # 13.69

# Center at the variant, then expand 250 bases on both sides unless peak ends early
highlightStartMid = max(accpred$peak_start, sub.hits$snpPos.hg38[1]-250)
highlightEndMid = min(accpred$peak_end, sub.hits$snpPos.hg38[1]+250)

flank = 50
highlightStartSmall = sub.hits$snpPos.hg38[1]-flank
highlightEndSmall = sub.hits$snpPos.hg38[1]+flank
small_highlight_rect = as.data.frame(gr_to_highlight(GRanges(paste0(sub.hits$snpChr.hg38[1], ":", highlightStartSmall, "-", highlightEndSmall))))

# 0-based coordinates
varscore_bed = rtracklayer::import.bed(paste0("/illumina/scratch/deep_learning/yng7/gChromvar/20250407-fetal_only_variant_scoring/",
                                              "GCST90132314_D003324_Coronary_artery_disease_rs12740374_Muscle_c6_peakcentered.bed"), 
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

region = paste0(accpred$chr, ":", highlightStartMid, "-", highlightEndMid)
region_zoom = paste0(sub.hits$snpChr.hg38[1], ":", highlightStartSmall, "-", highlightEndSmall)

preds_gr = varscore_bed %>% 
  gr_select(c("pred_ref", "pred_alt"), c("score1", "score2"))
p.accpred = trackplot_bw_paired(bw = preds_gr, region = region, alpha = 0.8,
                                facet_label = "Pred. accessibility", track_label = "",
                                score_cmap = c("ref" = "darkgrey", "alt" = "red3")) +
  highlight_region(region_zoom, color = "skyblue", alpha = 0.2) +
  theme(legend.position = "none")


#### Contributions for REF and ALT ####
library(glue)
contr_gr = varscore_bed %>% 
  gr_select(c("contrib_sum_ref", "contrib_sum_alt"), c("score1", "score2"))

# Set y axis limits to be the same for ref and alt contrib tracks
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

p.ref = trackplot_contribs2(contrib_ref, region_zoom, facet_label = "Contribs. REF", track_label = NULL) +
  highlight_relative_region(flank + 0.5, flank + 1.5, alpha = 0.2, color = "gold") +
  scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
  annotate("text", x = 1, y = ymax, label = range_label, vjust = 1.5, hjust = -0.1, size = 11*.8/ggplot2::.pt) +
  theme(legend.position = "none")

p.alt = trackplot_contribs2(contrib_alt, region_zoom, facet_label = "Contrib. ALT", track_label = NULL) +
  highlight_relative_region(flank + 0.5, flank + 1.5, alpha = 0.2, color = "gold") +
  scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
  annotate("text", x = 1, y = ymax, label = range_label, vjust = 1.5, hjust = -0.1, size = 11*.8/ggplot2::.pt) +
  theme(legend.position = "none")


#### Fetal peaks ####
minStart = sub.hits$snpPos.hg38[1]-10000
maxEnd = sub.hits$snpPos.hg38[1]+10000
mid_highlight_rect = as.data.frame(gr_to_highlight(GRanges(paste0(sub.hits$snpChr.hg38[1], ":", highlightStartMid, "-", highlightEndMid))))

p.acc2 = trackplot_helper_v2c(
  region = paste0(sub.hits$snpChr.hg38[1], ":", minStart, "-", maxEnd),
  cmap = organ_color, 5,
  clusters = submeta$L3_clusterID, 
  fragments = global_bp_obj$frags %>% select_cells(submeta$cb), 
  cell_read_counts = submeta$nFrags, 
  transcripts=global_bp_obj$transcripts
) + theme(legend.position = "None")

p.peaks = p.acc2[[2]] +
  geom_rect(data = mid_highlight_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "skyblue", alpha = 0.2) +
  theme(strip.text.y.left = element_blank(),
        axis.title.y.left = element_text(size = 11),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey80", linewidth = 0.2)) +
  labs(y = "Normalized\ninsertions\n(RPKM)")


#### Adult peaks ####
# ENCODE snATAC-seq pseudoreplicated peaks
filelist = fread(here::here("data/external/peaksets_adult/bedfile_annotations.tsv"), data.table = F)
filelist = filelist %>% 
  filter(tissue == "muscle organ" & peak_type == "pseudoreplicated peaks" & grepl("endothe", celltype))
filelist$encode_ids = gsub(".bigBed.bed", "", filelist$file)

adult.peaks = lapply(filelist$encode_ids, function(id){
  bed.df = fread(paste0(here::here("data/external/peaksets_adult/original/"), id, ".bigBed.bed"),
                 data.table = F,
                 col.names = c("chr", "start", "end", "peak_name", "score", "drop", "signal", "neglogp", "neglogq", "peak_num")) %>%
    mutate(encode_id = id) %>%
    select(encode_id, chr, start, end, signal, neglogp) %>%
    filter(start >= minStart & end <= maxEnd & chr == sub.hits$snpChr.hg38[1]) %>%
    filter(neglogp > 2)
  if(nrow(bed.df) == 0){
    data.frame(encode_id = id,
               chr = sub.hits$snpChr.hg38[1],
               start = 0,
               end = 0,
               signal = 0,
               neglogp = 0)
  } else { bed.df }
}) %>% bind_rows()

plotAdultPeaks = function(dat, plotColor){
  p = ggplot() +
    geom_rect(data = dat, aes(xmin = start, xmax = end, ymin = 0, ymax = signal), fill = plotColor) +
    geom_rect(data = mid_highlight_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "skyblue", alpha = 0.2) + 
    scale_x_continuous(expand = c(0, 0), limits = c(minStart, maxEnd)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(dat$signal))+1),
                       breaks = c(0, ceiling(max(dat$signal)))) +
    BPCells:::trackplot_theme() +
    labs(x = "", y = "Score") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 11, angle = 90, vjust = 0.5),
          strip.text = element_blank(),
          panel.border = element_rect(color = "grey80", linewidth = 0.2)) +
    facet_wrap(~encode_id, ncol = 1)
  return(p)
}

p.adult = plotAdultPeaks(adult.peaks, organ_color)


#### Second zoom gap ####
# Between peaks and accessibility scores
gap2 = data.frame(pos = c(minStart, seq(highlightStartMid, highlightEndMid, by = 5), maxEnd),
                  score = 0)
gap2[gap2$pos >= highlightStartMid & gap2$pos <= highlightEndMid,]$score = 1

p.gap2 = ggplot(gap2, aes(x = pos, y = score)) +
  geom_area(fill = "skyblue", alpha = 0.2, show.legend = F) +
  scale_x_continuous(expand = c(0, 0), limits = c(minStart, maxEnd)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


#### Third zoom gap ####
# Between accessbility/contribution scores and the rest
gap3 = data.frame(pos = c(highlightStartMid, seq(highlightStartSmall, highlightEndSmall, by = 5), highlightEndMid),
                  score = 0)
gap3[gap3$pos >= highlightStartSmall & gap3$pos <= highlightEndSmall,]$score = 1

p.gap3 = ggplot(gap3, aes(x = pos, y = score)) +
  geom_area(fill = "skyblue", alpha = 0.2, show.legend = F) +
  scale_x_continuous(expand = c(0, 0), limits = c(highlightStartMid, highlightEndMid)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

#### Motif hits ####
# Remove unresolved motifs and collapse hits from subclusters
hits_input = unique(global_bp_obj$hits$Muscle_c6)
hits_input = hits_input[grep("unresolved", hits_input$name, invert = T)]
p.motifs = trackplot_helper_v2c(
  region = paste0(sub.hits$snpChr.hg38[1], ":", highlightStartSmall, "-", highlightEndSmall),
  clusters = submeta$L2_clusterID, 
  fragments = global_bp_obj$frags %>% select_cells(submeta$cb), 
  cell_read_counts = submeta$nFrags, 
  transcripts = global_bp_obj$transcripts, 
  annot = hits_input,
  annot_color = "pattern_class",
  annot_labelby = "name",
  annot_strand = T,
  annot_label = "Motif calls"
) + theme(legend.position = "None")

p.mot = p.motifs[[4]] +
  scale_color_manual(values = rep("#A6A8AB",2)) +
  theme(strip.text.y.left = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


#### SNP label ####
pos_highlight_rect = data.frame(xmin = sub.hits$snpPos.hg38[1]-0.5, xmax = sub.hits$snpPos.hg38[1]+0.5,
                                ymin = -Inf, ymax = Inf)
snp.loc = data.frame(start = sub.hits$snpPos.hg38[1],
                     value = 1)
p.snp = ggplot() +
  geom_point(data = snp.loc, aes(x = start, y = value),
             shape = 19, size = 0.8, color = "orangered1") +
  scale_x_continuous(expand = c(0, 0), limits = c(highlightStartSmall, highlightEndSmall)) + 
  geom_rect(data = pos_highlight_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gold", alpha = 0.2) +
  labs(x = "", y = "") + 
  ylim(0.9, 1.1) +
  BPCells:::trackplot_theme() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

#### phyloP ####
system(paste0('bigWigToBedGraph -chrom=', sub.hits$snpChr.hg38[1],
              ' -start=', highlightStartSmall-1, ' -end=', highlightEndSmall, 
              ' ', here::here("data/external/variants/hg38.phyloP447wayPrimates.bw"), ' plotting_files_',
              rsid, '/phylopout.bedGraph'))

phylop = read.delim(paste0("plotting_files_", rsid, "/phylopout.bedGraph"), header = F)
colnames(phylop) = c("chr", "start", "end", "score")

# Need a value for each position
# Function to expand ranges into individual rows
expand_ranges = function(df) {
  expanded_df = do.call(rbind, lapply(1:nrow(df), function(i) {
    chr = df$chr[i]
    score = df$score[i]
    start = df$start[i]
    end = df$end[i]
    
    # Create a sequence for each base pair in the range
    data.frame(
      chr = chr,
      start = start:(end - 1),
      end = (start + 1):end,
      score = score
    )
  }))
  return(expanded_df)
}

expanded_df = expand_ranges(phylop)

phylop_pos = expanded_df
phylop_pos$score = ifelse(phylop_pos$score < 0, 0, phylop_pos$score)
phylop_neg = expanded_df
phylop_neg$score = ifelse(phylop_neg$score >= 0, 0, phylop_neg$score)

p.phylop = ggplot() +
  geom_col(data = phylop_pos, aes(x = end, y = score), fill = "royalblue", color = "white", linewidth = 0.3) +
  geom_col(data = phylop_neg, aes(x = end, y = score), fill = "indianred", color = "white", linewidth = 0.3) +
  scale_x_continuous(expand = c(0, 0), limits = c(highlightStartSmall+0.5, highlightEndSmall-0.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(min(phylop$score)-0.1, max(phylop$score)+0.1),
                     breaks = c(-1,0,1)) +
  geom_rect(data = pos_highlight_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gold", alpha = 0.2) +
  BPCells:::trackplot_theme() +
  labs(x = "", y = "") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())


#### Assemble facets ####
indiv_theme = theme(plot.margin = unit(c(0, 0, 0, 0), "pt"),
                    legend.position = "none")
# For second zoom gap height
# Each adult peak track gets 0.0125; there are 5
second_zoom_height = 1 - sum(0.03, 0.05, 0.085, 0.04, 0.08, 0.25, 0.07, 0.03, 0.03, 0.08, 0.045, 0.045, 0.06)

p.plots = trackplot_combine(list(BPCells:::wrap_trackplot(p.sbar, unit(0.03, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.annot, unit(0.05, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.gwas, unit(0.085, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.gap1, unit(0.04, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.peaks, unit(0.08, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.adult, unit(0.25, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.gap2, unit(second_zoom_height, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.accpred, unit(0.07, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.gap3, unit(0.03, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.snp, unit(0.03, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.phylop, unit(0.08, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.ref, unit(0.045, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.alt, unit(0.045, "null")) + indiv_theme,
                                 BPCells:::wrap_trackplot(p.mot, unit(0.06, "null")) + indiv_theme))

ggsave(paste0(figout, "/tracks_", rsid, "_muscle_endothelial_CAD_v4.pdf"), 
       p.plots, height = 7, width = 6.5)
