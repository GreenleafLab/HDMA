# Script for qsub for shuffled counts
library(data.table)
library(tidyverse)
library(arrow)
library(GenomicRanges)
library(parallel)
library(here)

out <- here::here("output/06-variants/02/")
figout <- here::here("figures/06-variants/02/")
dir.create(out, showWarnings = F, recursive=T)
dir.create(figout, showWarnings = F, recursive=T)
setwd(out)

# eQTLs
eQTL_list = list.files(here::here("data/external/variants/GTEx_eQTL_SuSiE"),
                       full.names = T)

# Motif calls
motif_hits_list = list.files(here::here("output/03-chrombpnet/02-compendium/hits_unified_motifs/reconciled_per_celltype_peaks/"),
                             pattern = "reconciled.bed.gz$", full.names = T, recursive=T)

# Options
args = commandArgs(trailingOnly = T)

# 1: initial seed
initial_seed = as.numeric(args[1])
# 2: tissue
tissue = args[2]
# 3: set number
setNum = as.numeric(args[3])

set_df = read.delim("sets.txt", header = F)
colnames(set_df) = c("set", "idx_start", "idx_end", "initial_seed")
set_df = set_df[set_df$set == setNum,]

# sets.txt looks like this:
# set idx_start idx_end initial_seed
#   1         1    1000         1913
#   2      1001    2000         2914
#   3      2001    3000         3915
#   4      3001    4000         4916
#   5      4001    5000         5917
#   6      5001    6000         6918
# ...

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


# Shuffle ----------------------------------------------------------------------
# Prepare
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

ct.eqtl.files = eQTL_list[grep(tissue.grep, eQTL_list)]
ct.eqtls = lapply(ct.eqtl.files, function(eqtlfile){
  df = read_parquet(eqtlfile)
  df$chr = word(df$variant_id, 1, 1, sep = "_")
  df$pos = as.numeric(word(df$variant_id, 2, 2, sep = "_"))
  df$id = paste0(df$variant_id, "_", df$phenotype_id)
  df = select(df, chr, pos, id, pip, afc)
  df
}) %>% bind_rows()

ct.eqtls = ct.eqtls %>%
  group_by(chr, pos, id) %>%
  arrange(desc(pip)) %>%
  slice_head(n = 1)

# Shuffle
shuf.results = lapply(set_df$idx_start[1]:set_df$idx_end[1], function(i){
  cat(i, "\n")
  seed = initial_seed + i
  set.seed(seed)
  df.shuffle = ct.eqtls
  df.shuffle$afc = sample(df.shuffle$afc)
  df.shuffle$direction = ifelse(df.shuffle$afc < 0, "down", "up")
  
  df.shuffle.up = df.shuffle[df.shuffle$direction == "up",]
  df.shuffle.down = df.shuffle[df.shuffle$direction == "down",]
  
  # Each motif is tested separately
  results.pos.up = count_hits(df.shuffle.up, ct.hits.pos)
  results.pos.up = add_missing(results.pos.up, ct.hits.pos)
  if(nrow(results.pos.up) != 0){
    results.pos.up$test_type = "pos_up"
  }
  
  results.pos.down = count_hits(df.shuffle.down, ct.hits.pos)
  results.pos.down = add_missing(results.pos.down, ct.hits.pos)
  if(nrow(results.pos.down) != 0){
    results.pos.down$test_type = "pos_down"
  }
  results.neg.up = count_hits(df.shuffle.up, ct.hits.neg)
  results.neg.up = add_missing(results.neg.up, ct.hits.neg)
  if(nrow(results.neg.up) != 0){
    results.neg.up$test_type = "neg_up"
  }
  results.neg.down = count_hits(df.shuffle.down, ct.hits.neg)
  results.neg.down = add_missing(results.neg.down, ct.hits.neg)
  if(nrow(results.neg.down) != 0){
    results.neg.down$test_type = "neg_down"
  }
  
  # Bind even if dataframe is empty
  results = bind_rows(list(results.pos.up, results.pos.down, results.neg.up, results.neg.down))
  results$shuf.idx = i
  results
}) %>% bind_rows()

shuf.results$tissue = tissue

fwrite(shuf.results, paste0(out, "/results_100k_", tissue, "_set", setNum, ".txt"),
       sep = "\t", row.names = F)
system(paste0("gzip -f ", out, "/results_100k_", tissue, "_set", setNum, ".txt"))
