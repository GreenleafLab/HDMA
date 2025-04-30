# Purpose: after producing filtered and reconciled set of hit calls, here we 
# compute some annotations on a per-hit level, namely:
# - genomic annotations done using ArchR:::.fastAnnoPeaks, which annotates peaks
#   according to the genomic features they overlap (distal, intronic, exonic, promoter),
#   distance to TSS, GC content, etc
#   see https://github.com/GreenleafLab/ArchR/blob/d9e741c980c7c64e5348c97a74d146cc95f8ba76/R/ReproduciblePeakSet.R#L590
# - distance to peak summit

# libraries
library(here)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(ggplot2)
library(cowplot)

script_path <- here("code/utils/")
source(file.path(script_path, "plotting_config.R"))
source(file.path(script_path, "hdma_palettes.R"))

ggplot2::theme_set(theme_BOR())


# 0. GET ARGS ------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
print(args)

peaks_bed    = args[1]
hits_tsv_in  = args[2]
hits_tsv_out = args[3]
out_dir      = args[4]



# 1. LOAD DATA -----------------------------------------------------------------

message("@ loading data...")

hits <- data.table::fread(hits_tsv_in, data.table = FALSE, header = TRUE)

# if these exist already, remove so we can overwrite
cols_drop <- c("distToGeneStart", "nearestGene", "peakType", "distToTSS", "nearestTSS", "GC", "N", "distToPeakSummit")

if (any(cols_drop %in% hits)) {
  
  hits <- hits[, !(names(hits) %in% cols_drop)]
  
}

# convert 0-based to 1-based (typically done automatically by rtracklayer)
hits$start <- hits$start + 1
hits <- GRanges(hits)
print(head(hits))
print(length(hits))

peaks <- rtracklayer::import.bed(peaks_bed, extraCols = c("a" = "character",
                                                          "b" = "character",
                                                          "c" = "character",
                                                          "d" = "character",
                                                          "e" = "character",
                                                          "f" = "character",
                                                          "summit" = "numeric"))
print(head(peaks))
print(length(peaks))



# 2. COMPUTE GENOMIC ANNO ------------------------------------------------------

message("@ computing genomic annotation...")
# First we want to get an annotation per hit that says what type of genomic
# location it occurs in (promoter, intronic, exonic, distal) as per the ArchR definitions.

# set default genome
addArchRGenome("hg38")

# get gene annotation

# Avoid - this uses a very old v24 Gencode annotation packaged with ArchR
# gene_annotation <- getArchRGenome(geneAnnotation = TRUE, genomeAnnotation = FALSE) %>% as.list

# Get the latest anno
gene_annotation <- ArchR::createGeneAnnotation(OrgDb = org.Hs.eg.db, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)

# > gene_annotation$TSS@metadata %>% unlist() %>% tibble::enframe()
# # A tibble: 18 × 2
# name                                                value                                       
# <chr>                                               <chr>                                       
#   1 genomeInfo.Db type                                  TxDb                                        
# 2 genomeInfo.Supporting package                       GenomicFeatures                             
# 3 genomeInfo.Data source                              UCSC                                        
# 4 genomeInfo.Genome                                   hg38                                        
# 5 genomeInfo.Organism                                 Homo sapiens                                
# 6 genomeInfo.Taxonomy ID                              9606                                        
# 7 genomeInfo.UCSC Table                               knownGene                                   
# 8 genomeInfo.UCSC Track                               GENCODE V38                                 
# 9 genomeInfo.Resource URL                             http://genome.ucsc.edu/                     
#   10 genomeInfo.Type of Gene ID                          Entrez Gene ID                              
# 11 genomeInfo.Full dataset                             yes                                         
# 12 genomeInfo.miRBase build ID                         NA                                          
# 13 genomeInfo.Nb of transcripts                        258145                                      
# 14 genomeInfo.Db created by                            GenomicFeatures package from Bioconductor   
# 15 genomeInfo.Creation time                            2021-10-19 10:58:00 -0700 (Tue, 19 Oct 2021)
# 16 genomeInfo.GenomicFeatures version at creation time 1.45.2                                      
# 17 genomeInfo.RSQLite version at creation time         2.2.8                                       
# 18 genomeInfo.DBSCHEMAVERSION                          1.2    


gene_annotation <- gene_annotation %>% as.list()

# annotate peaks
# https://github.com/GreenleafLab/ArchR/blob/d9e741c980c7c64e5348c97a74d146cc95f8ba76/R/ReproduciblePeakSet.R#L590

# about 10 min
Sys.time()
hits_anno <- ArchR:::.fastAnnoPeaks(hits, BSgenome = BSgenome.Hsapiens.UCSC.hg38, geneAnnotation = gene_annotation, promoterRegion = c(2000, 2000))
Sys.time()
print(length(hits_anno))



# 3. COMPUTE DIST TO PEAK SUMMIT -----------------------------------------------

message("@ computing peak summit distances...")

# Secondly, we also want to get the distance to peak summit. 
# Since ChromBPNet-formatted peaks are defined as 1000bp around peak summit, we can 
# take the center of ChromBPNet peaks or the start+500 to get the peak summit.

# we should be able to reuse the functionality from ArchR:::.fastAnnoPeaks
peak_summits <- GenomicRanges::resize(peaks, 1, "center")
hit_centers <- GenomicRanges::resize(hits, 1, "center")

# get distances
dist_summits <- distanceToNearest(hit_centers, peak_summits, ignore.strand = TRUE)
print(length(dist_summits))

# add to our GRanges object
mcols(hits_anno)$distToPeakSummit <- mcols(dist_summits)$distance
print(head(hits_anno))



# 4. SAVE TSV OUT --------------------------------------------------------------

message("@ writing out...")

# convert from 1-based coords (GRanges) to 0-based coordinates (BED)
hits_anno_df <- data.frame(hits_anno)
hits_anno_df$start <- hits_anno_df$start - 1
print(head(hits_anno_df))

# write out to TSV
write_tsv(hits_anno_df, file = hits_tsv_out)



# 5. SAVE SUMMARY PLOTS & METRICS ----------------------------------------------
# precalculate a few metrics so that hits don't need to be loaded in memory
# for this later

message("@ calculating summaries...")

# read in hits per motif
hits_per_motif <- read_tsv(paste0(out_dir, "/hits_per_motif.tsv")) %>% 
  # in case we're re-running, then only select the columns in the original
  # file output by the reconciling step
  dplyr::select(motif_name, pattern_class, count)

# number of hits to each element type per motif
hits_anno_summary1 <- hits_anno_df %>%
  group_by(motif_name, peakType) %>%
  count() %>% 
  pivot_wider(names_from = peakType, values_from = n)

# average distances
hits_anno_summary2 <- hits_anno_df %>%
  group_by(motif_name) %>%
  summarise(mean_distToTSS = round(mean(distToTSS), 2),
            median_distToTSS = median(distToTSS),
            mean_distToPeakSummit = round(mean(distToPeakSummit), 2),
            median_distToPeakSummit = median(distToPeakSummit))

# overwrite hits_per_motif with anno info
hits_per_motif <- hits_per_motif %>% 
  left_join(hits_anno_summary1, by = "motif_name") %>% 
  left_join(hits_anno_summary2, by = "motif_name") %>% 
  arrange(desc(count))

print(head(hits_per_motif))

hits_per_motif %>% 
  write_tsv(paste0(out_dir, "/hits_per_motif.tsv"))

# binned distances to peak summits
summit_dist_df <- hits_anno_df %>% 
  mutate(bin_dist = cut(distToPeakSummit,
                        breaks = seq(0, 500, by = 10), 
                        labels = seq(0, 499, by = 10), include.lowest = TRUE)) %>%
  group_by(motif_name, bin_dist) %>% 
  count()

# save binned counts
write_tsv(summit_dist_df, file = paste0(out_dir, "/motif_dist_to_summit.tsv"))


message("@ making plots...")

# rank motifs by distance to TSS
tss_dist_order <- hits_anno_df %>%
  dplyr::group_by(motif_name) %>%
  dplyr::summarize(median_dist = median(distToTSS)) %>%
  arrange(dplyr::desc(median_dist)) %>%
  pull(motif_name)

# most frequent motifs
top_motifs <- hits_anno_df %>% 
  group_by(motif_name) %>% 
  dplyr::count() %>% 
  dplyr::arrange(desc(n)) %>% 
  ungroup() %>% 
  slice_max(order_by = n, n = 50) %>% 
  pull(motif_name)

# genomic annotation barplot
p1 <- hits_anno_df %>% 
  filter(motif_name %in% top_motifs) %>% 
  mutate(motif_name = factor(motif_name, levels = tss_dist_order)) %>% 
  ggplot(aes(x = motif_name)) +
  geom_bar(aes(fill = peakType), position = "fill") +
  scale_fill_manual(values = cmap_peaktype) +
  coord_flip() +
  ylab("# hits") + xlab("motif") +
  ggtitle("Breakdown of hits by \ngenomic annotation") +
  theme(legend.position = "bottom") +
  theme(axis.text.y = element_text(size = 7))

# hit counts
p2 <- hits_anno_df %>% 
  filter(motif_name %in% top_motifs) %>% 
  mutate(motif_name = factor(motif_name, levels = tss_dist_order)) %>% 
  ggplot(aes(x = motif_name)) +
  geom_bar(position = "stack", fill = "black") +
  coord_flip() +
  ylab("# hits") + xlab(NULL) +
  scale_y_continuous(labels = scales::comma) +
  ggtitle("# hits for the top 50 \nmost frequent motifs") +
  theme(axis.text.y = element_blank())

# distance to TSS
p3 <- hits_anno_df %>% 
  mutate(motif_name = factor(motif_name, levels = tss_dist_order)) %>% 
  filter(motif_name %in% top_motifs) %>%
  ggplot(aes(x = motif_name, y = distToTSS)) +
  geom_boxplot(aes(fill = pattern_class), outliers = FALSE, alpha = 0.4) +
  scale_fill_manual(values = c("pos_patterns" = "darkgreen", "neg_patterns" = "red")) +
  coord_flip() +
  xlab(NULL) + ylab("distance to nearest TSS (bp)") +
  scale_y_continuous(labels = scales::comma) +
  ggtitle("Distribution of distances to nearest \nTSS (top 50 most frequent motifs)") +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "bottom")

cowplot::plot_grid(p1, p2, p3, align = "h", rel_widths = c(0.4, 0.3, 0.3), nrow = 1)
ggsave(paste0(out_dir, "/genomic_localization.pdf"), width = 13, height = 8)

# distance to summit
summit_dist_order <- hits_anno_df %>%
  dplyr::group_by(motif_name) %>%
  dplyr::summarize(median_dist = median(distToPeakSummit)) %>%
  arrange(dplyr::desc(median_dist)) %>%
  pull(motif_name)

p4 <- hits_anno_df %>% 
  mutate(motif_name = factor(motif_name, levels = summit_dist_order)) %>% 
  filter(motif_name %in% top_motifs) %>%
  ggplot(aes(x = motif_name, y = distToPeakSummit)) +
  geom_boxplot(aes(fill = pattern_class), outliers = FALSE, alpha = 0.4) +
  scale_fill_manual(values = c("pos_patterns" = "darkgreen", "neg_patterns" = "red")) +
  coord_flip(ylim = c(0, 500)) +
  xlab(NULL) + ylab("distance to peak summit (bp)") +
  ggtitle("Distribution of distances to peak summit \n(top 50 most frequent motifs)") +
  theme(legend.position = "bottom") +
  theme(axis.text.y = element_text(size = 7))

p4
ggsave(paste0(out_dir, "/distance_to_summits.pdf"), width = 6, height = 8)


message("@ done.")

