# Purpose: after running NucleoATAC, load the nucleosome positions
# and the annotated hit calls and plot the binned distribution of hits
# for each motif relative to nucleosome dyads, as well as the binned distribution
# of pairwise distances between dyads.

library(here)
library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(ggplot2)
library(glue)
# library(cowplot)
# library(pheatmap)

script_path <- here("code/utils/")
source(file.path(script_path, "plotting_config.R"))
source(file.path(script_path, "hdma_palettes.R"))
source(file.path(script_path, "sj_scRNAseq_helpers.R"))

ggplot2::theme_set(theme_BOR())


# 0. GET ARGS ------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
print(args)

occ_bed              = args[1]
hits_tsv             = args[2]
out_prefix           = args[3]
dyad_dist_counts_tsv = args[4]



# 1. LOAD DATA -----------------------------------------------------------------

message("@ loading data...")

hits <- data.table::fread(hits_tsv, data.table = FALSE, header = TRUE)

# convert 0-based to 1-based (done automatically by rtracklayer)
hits$start <- hits$start + 1
hits <- GRanges(hits)
print(head(hits))
print(length(hits))


# Columns for occpeaks.bed
# (1) chrom
# (2) dyad position (0-based)
# (3) dyad position (1-based)
# (4) occupancy
# (5) occupancy lower bound
# (6) occupancy upper bound
# (7) of reads in 121 bp window
occ <- data.table::fread(occ_bed, data.table = FALSE, header = TRUE,
                         col.names = c("chr", "start", "end", "occupancy",
                                       "occupancy_lower", "occupancy_upper",
                                       "n_reads"))

# convert 0-based to 1-based
occ$start <- occ$start + 1
occ <- GRanges(occ)

head(occ)



# load dyad distances
dyad_dist_counts <- read_tsv(dyad_dist_counts_tsv)


# 2. CALC DISTANCES ------------------------------------------------------------

message("@ calculating distances...")

# Reuse the functionality from ArchR:::.fastAnnoPeaks
hit_centers <- GenomicRanges::resize(hits, 1, "center")
dist_occ <- distanceToNearest(hit_centers, occ, ignore.strand = TRUE)
print(length(dist_occ))

# If some hits do not have nearest matches (occurs when nucleosomes called on
# chrY), then drop those hits
if (length(unique(queryHits(dist_occ))) < length(hits) ) {
  
  hits <- hits[base::intersect(seq_along(hits), queryHits(dist_occ))]
  
}

mcols(hits)$distToDyad <- mcols(dist_occ)$distance

# convert from 1-based to 0-based
hits_df <- data.frame(hits)
hits_df$start <- hits_df$start - 1


# # Calculate pairwise distances between one nucleosome call and all others
# #
# # @param query_idx Numeric, index of region to use as query in \code{regions}
# # @param region GRanges containing all other regions
# pairwise_distances.per_region <- function(query_idx, all_regions, max_dist = 250, k_nearest = 20) {
#   
#   # get distance between one dyad and all others
#   subject_idx <- unlist(nearestKNeighbors(all_regions[query_idx], all_regions, k = k_nearest))
#   dist_per_region <- distance(all_regions[query_idx], all_regions[subject_idx], ignore.strand = TRUE)
#   
#   # keep the counts for nucleosomes within max-dist of another,
#   # and drop the distance to self (=0)
#   dist_per_region <- dist_per_region[which(dist_per_region < max_dist & dist_per_region > 0)]
#   
#   if (length(dist_per_region) == 0) return(NULL)
#   else return(dist_per_region)
#   
# }
# 
# # our calls are sparse. So we can safely calculate distances for each dyad call with only the nearest K dyads.
# start.time <- Sys.time()
# dist_occ_pairwise <- map(seq_along(occ), ~ pairwise_distances.per_region(.x, occ)) %>% 
#   unlist()
# end.time <- Sys.time()
# end.time - start.time




# 3. PLOT ----------------------------------------------------------------------
message("@ plotting...")

nuc_dist_df <- hits_df %>%
  # focus on hits within 250bp of a nucleosome dyad call
  filter(distToDyad < 250) %>% 
  # group distances into 25 bins representing 10 bp increments, and count total hits per bin
  mutate(bin_dist = cut(distToDyad,
                        breaks = seq(0, 250, by = 10),
                        labels = seq(10, 251, by = 10),
                        include.lowest = TRUE)) %>%
  group_by(motif_name, bin_dist) %>% 
  count()

# save binned counts
write_tsv(nuc_dist_df, file = glue("{out_prefix}.motif_dist_to_dyad.tsv"))

# calculate an order by median distance to dyad
median_dist_to_dyad <- hits_df %>%
  filter(distToDyad < 250) %>% 
  group_by(motif_name) %>%
  summarize(n_hits = n(),
            min_dist = min(distToDyad),
            median_dist = median(distToDyad)) %>%
  slice_max(order_by = n_hits, n = 50) %>% 
  arrange(median_dist)

# make the matrix
# nuc_dist_mat <- nuc_dist_df %>% 
#   pivot_wider(names_from = bin_dist, values_from = n) %>% 
#   filter(motif_name %in% median_dist_to_dyad$motif_name) %>% 
#   tibble::column_to_rownames(var = "motif_name") %>% 
#   as.matrix()

# pheatmap(nuc_dist_mat[rev(median_dist_to_dyad$motif_name), ],
#          scale = "row",
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          border_color = NA,
#          color = rdbu2,
#          fontsize_col = 9,
#          fontsize_row = 9,
#          filename = "nucleoatac_dist.pdf",
#          cellwidth = 10,
#          cellheight = 10)

# plot
nuc_dist_df %>% 
  filter(motif_name %in% median_dist_to_dyad$motif_name) %>% 
  group_by(motif_name) %>% 
  # zscore counts per motif (row-wise)
  mutate(zscore = scale(n)) %>% 
  mutate(motif_name = factor(motif_name, levels = median_dist_to_dyad$motif_name)) %>% 
  ggplot(aes(x = bin_dist, y = motif_name)) +
  geom_tile(aes(fill = zscore)) +
  scale_fill_gradientn(colours = rdbu2,
                       # set midpoint to 0
                       # https://github.com/tidyverse/ggplot2/issues/3738#issuecomment-1336166757
                       rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
  # scale_x_continuous(breaks = seq(10, 250, by = 10)) +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 9)) +
  xlab("binned distance from nucleosome dyad (bp)") +
  ylab("motif") +
  ggtitle("Distribution of distances from hits to nucleosome dyads\n(top 50 most frequent motifs)")

ggsave(glue("{out_prefix}.nuc_pos.pdf"), width = 8, height = 7.5)


# plot 
p1 <- dyad_dist_counts %>% 
  dplyr::select(bin_dist = Bin, n = Count) %>% 
  mutate(bin_dist = as.numeric(bin_dist)) %>% 
  ggplot(aes(x = bin_dist, y = 1)) +
  geom_tile(aes(fill = n)) +
  scale_fill_gradientn(colours = ylrd) +
  scale_x_continuous(breaks = seq(10, 250, by = 10)) +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank()) +
  xlab("binned distance between nucleosome dyads (bp)") +
  ylab(NULL) +
  ggtitle("Distribution of pairwise distances between \nnuc dyad positions")

p2 <- dyad_dist_counts %>% 
  dplyr::select(bin_dist = Bin, n = Count) %>% 
  mutate(zscore = scale(n)) %>%
  mutate(bin_dist = as.numeric(bin_dist)) %>% 
  ggplot(aes(x = bin_dist, y = 1)) +
  geom_tile(aes(fill = zscore)) +
  scale_fill_gradientn(colours = rdbu2) +
  scale_x_continuous(breaks = seq(10, 250, by = 10)) +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank()) +
  xlab("binned distance between nucleosome dyads (bp)") +
  ylab(NULL) +
  ggtitle("Distribution of pairwise distances between \nnuc dyad positions")

ggsave(plot = cowplot::plot_grid(p1, p2, nrow = 2),
       filename = glue("{out_prefix}.dyad_binned_dist.pdf"), width = 8, height = 4)

message("@ done.")


