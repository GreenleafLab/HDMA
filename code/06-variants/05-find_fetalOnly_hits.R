options(scipen = 999)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(here)


out <- here::here("output/06-variants/05/")
figout <- here::here("figures/06-variants/05/")
dir.create(out, showWarnings = F, recursive=T)
dir.create(figout, showWarnings = F, recursive=T)
setwd(out)

dir.create(paste0(out, "/motifs_w_causal_snps"), showWarnings = F)

#### Find motifs hit by causal SNPs --------------------------------------------
motifdir = here::here("output/03-chrombpnet/02-compendium/hits_unified_motifs/reconciled_per_celltype_peaks/")
metadir = here::here("output/01-preprocessing/02/shared/meta/")
tissue2abbr = data.frame(tissue = c("Adrenal", "Brain", "Eye", "Heart", "Liver", "Lung", "Muscle",
                                    "Skin", "Spleen", "Stomach", "Thyroid", "Thymus"),
                         abbrev = c("AG", "BR", "EY", "HT", "LI", "LU", "MU", "SK", "SP", "ST", "TR", "TM"))
# CAUSALdb SNPs
cdb.full = fread(here::here("data/external/variants/20240508_causaldb_wMeta.txt.gz"),
                 nThread = 12, data.table = F)
# For liftover
ch38to19 = import.chain(here::here("data/external/variants/hg38ToHg19.over.chain"))
ch19to38 = import.chain(here::here("data/external/variants/hg19ToHg38.over.chain"))

# Get gchromvar results
is.sig = fread(here::here("output/06-variants/04/sig_results_all_tissues_disease_only_pipsum_5_fdr_0.05.txt.gz"),
               sep = "\t", data.table = F)
traits = unique(is.sig$trait) # 68

results.df = mclapply(traits, function(traitName){
  cat(traitName, "\n")
  meta_id = word(traitName, 1, 1, sep = "_")
  
  # Get causal SNPs for trait
  causalSNPs = cdb.full[cdb.full$meta_id == meta_id,]
  
  causalSNPs = causalSNPs[causalSNPs$susie >= 0.8,] %>%
    select(chr, bp, susie, trait, primary) %>%
    arrange(chr, bp)
  # For variants with more than 1 susie score,
  # prioritize the one where primary = 1
  causalSNPs = causalSNPs %>%
    group_by(chr, bp, trait) %>%
    slice_max(primary) %>%  # Pick the row with the maximum primary value
    ungroup() %>%
    select(-primary) %>%
    as.data.frame()
  causalSNPs$chr = paste0("chr", causalSNPs$chr)
  causalSNPs$snp.chr.hg19 = causalSNPs$chr
  causalSNPs$snp.pos.hg19 = causalSNPs$bp
  
  snps.gr = makeGRangesFromDataFrame(causalSNPs,
                                     seqnames.field = "chr",
                                     start.field = "bp",
                                     end.field = "bp",
                                     keep.extra.columns = T)
  
  trait.results = filter(is.sig, trait == traitName)
  
  motifs.w.causal.snps = lapply(1:nrow(trait.results), function(i){
    cat(i, "\n")
    meta = fread(paste0(metadir, "/",
                        tissue2abbr[tissue2abbr$tissue == trait.results$tissue[i],]$abbrev,
                        "_meta_cluster.txt"), data.table = F) %>%
      select(L1_clusterID, L2_clusterName) %>% distinct()
    meta$L2_clusterName = gsub(" ", ".", meta$L2_clusterName, fixed = T)
    meta$L2_clusterName = gsub("-", ".", meta$L2_clusterName, fixed = T)
    meta$L2_clusterName = gsub("+", ".", meta$L2_clusterName, fixed = T)
    meta$L2_clusterName = gsub("/", ".", meta$L2_clusterName, fixed = T)
    meta$L2_clusterName = gsub("(", ".", meta$L2_clusterName, fixed = T)
    meta$L2_clusterName = gsub(")", ".", meta$L2_clusterName, fixed = T)
    
    cat(trait.results$celltype[i], "\n")
    cat(meta[meta$L2_clusterName == trait.results$celltype[i],]$L2_clusterName, "\n")
    clusterIDs = meta[meta$L2_clusterName == trait.results$celltype[i],]$L1_clusterID
    motifs = lapply(clusterIDs, function(cid){
      if(file.exists(paste0(motifdir, "/",
                            trait.results$tissue[i], "_c", cid, "/counts_v0.23_a0.8_all/",
                            trait.results$tissue[i], "_c", cid,
                            "__hits_unified.counts_v0.23_a0.8_all.reconciled.bed.gz"))){
        mot = fread(paste0(motifdir, "/",
                           trait.results$tissue[i], "_c", cid, "/counts_v0.23_a0.8_all/",
                           trait.results$tissue[i], "_c", cid,
                           "__hits_unified.counts_v0.23_a0.8_all.reconciled.bed.gz"),
                    data.table = F)[,c(1:4)]
        colnames(mot) = c("chr", "start", "end", "motifName")
        mot$motifChr.hg38 = mot$chr
        mot$motifStart.hg38 = as.numeric(mot$start) + 1
        mot$motifEnd.hg38 = mot$end
        mot
      } else { cat(paste0("File for cluster ", cid, " not found.\n"))}
    }) %>% bind_rows() %>% distinct()
    
    if(nrow(motifs) != 0){
      motifs$tissue = trait.results$tissue[i]
      motifs$celltype = trait.results$celltype[i]
      motifs.gr = makeGRangesFromDataFrame(motifs,
                                           seqnames.field = "chr",
                                           start.field = "start",
                                           end.field = "end",
                                           starts.in.df.are.0based = T,
                                           keep.extra.columns = T)
      motifs.gr.hg19 = unlist(liftOver(motifs.gr, ch38to19))
      motifs.hg19 = as.data.frame(motifs.gr.hg19)
      motifs.hg19[,c("motifChr.hg19", "motifStart.hg19", "motifEnd.hg19")] = motifs.hg19[,c("seqnames", "start", "end")]
      motifs.hg19 = motifs.hg19 %>% 
        select("motifChr.hg38", "motifStart.hg38", "motifEnd.hg38",
               "motifChr.hg19", "motifStart.hg19", "motifEnd.hg19",
               "seqnames", "start", "end",
               "motifName", "tissue", "celltype")
      motifs.gr.hg19 = makeGRangesFromDataFrame(motifs.hg19,
                                                seqnames.field = "seqnames",
                                                start.field = "start",
                                                end.field = "end",
                                                keep.extra.columns = T)
      
      # Find motifs hit by causal SNPs
      hits = findOverlaps(snps.gr, motifs.gr.hg19)
      if(length(hits) != 0){
        motif.snp.hits = mcols(motifs.gr.hg19[subjectHits(hits)])
        motif.snp.hits = cbind.data.frame(
          motif.snp.hits,
          mcols(snps.gr[queryHits(hits)])
        )
        motif.snp.hits
      }
    }
  }) %>% bind_rows()
  
  if(nrow(motifs.w.causal.snps) != 0){
    motifs.w.causal.snps = motifs.w.causal.snps %>% 
      distinct() %>%
      arrange(motifName)
    rownames(motifs.w.causal.snps) = 1:nrow(motifs.w.causal.snps)
    
    # Lift snp coordinates to hg38
    df.hg19 = motifs.w.causal.snps %>%
      mutate(seqnames = snp.chr.hg19,
             start = as.numeric(snp.pos.hg19) -1,
             end = snp.pos.hg19)
    df.hg19.gr = makeGRangesFromDataFrame(df.hg19,
                                          seqnames.field = "seqnames",
                                          start.field = "start",
                                          end.field = "end",
                                          keep.extra.columns = T)
    df.hg38.gr = unlist(liftOver(df.hg19.gr, ch19to38))
    df.hg38 = as.data.frame(df.hg38.gr) %>%
      dplyr::rename(snpChr.hg38 = seqnames, snpPos.hg38 = end,
                    snpChr.hg19 = snp.chr.hg19, snpPos.hg19 = snp.pos.hg19,
                    organ = tissue) %>%
      select(-width, -strand, -start) %>%
      relocate(snpChr.hg38, snpPos.hg38, .after = trait) %>%
      mutate(chrnum = gsub("chr", "", snpChr.hg38)) %>%
      arrange(trait, organ, celltype, chrnum, snpPos.hg38) %>%
      select(-chrnum)
    
    # Add additional info like Z score and p values
    df.hg38 = left_join(df.hg38, trait.results[, !names(trait.results) %in% "trait"],
                        by = c("organ"="tissue", "celltype"))
    df.hg38$fullTraitName = trait.results$trait[1]
    df.hg38 = df.hg38 %>%
      dplyr::rename(meta_id = studyid) %>%
      relocate(meta_id, mesh_id, trait, motifName, organ, celltype,
               susie, zscore, pvalue, bonf, fdr)
    fwrite(df.hg38, paste0("motifs_w_causal_snps/", meta_id, "_pipsum_5_fdr_0.05.txt"), row.names = F, sep = "\t")
    df.hg38
  }
}, mc.cores = 2) %>% bind_rows()

fwrite(results.df, "v3motifs_w_causal_snps_all_tissues_disease_only_pipsum_5_fdr_0.05.txt",
       row.names = F, sep = "\t")

## Assign if fetal only (peaksets from unpublished ENCODE data) -----------------------------------
results.df = fread("v3motifs_w_causal_snps_all_tissues_disease_only_pipsum_5_fdr_0.05.txt", data.table = F)
results.df$idx = 1:nrow(results.df)

# Wrote out file for manual matching of cell types
# fetal2adult.matching = select(results.df, organ, celltype) %>%
#   distinct() %>%
#   arrange(organ, celltype)
# fwrite(fetal2adult.matching, "fetal2adult_matching.txt", sep = "\t")

# Saved as same file name
fetal2adult.matching = fread("fetal2adult_matching.txt", sep = "\t", data.table = F)

results.df = left_join(results.df, fetal2adult.matching, by = c("organ", "celltype"))

# Only those which have matching adult bed files can be checked
results.df$foundIn = "fetalOnly"
results.df$foundInXpeaksets = 0
results.df$totalAdultPeaksets = 0

adult.peaks.dir = here::here("data/external/peaksets_adult/original")

for(i in 1:nrow(results.df)){
  if(results.df$adult_celltype[i] != "unassigned"){
    cat(i, "\n")
    adultBED_meta = read.delim(paste0(here::here("data/external/peaksets_adult/meta_adultBEDs/meta_adultBEDs_"),
                                      results.df$organ[i], ".tsv"),
                               header = T)
    adult_bedfiles = adultBED_meta[adultBED_meta$adult_celltype == results.df$adult_celltype[i],]
    
    beds.grlist = lapply(unique(adult_bedfiles$file), function(y){
      bed.df = fread(paste0(adult.peaks.dir, "/", y), data.table = F,
                     col.names = c("seqnames", "start", "end", "name", "score",
                                   "drop", "signalValue", "pValue", "qValue", "peak"))
      # Drop those with insignificant qvalues (-log(0.01,10))
      bed.df = bed.df[bed.df$qValue > 2,]
      bed.gr = makeGRangesFromDataFrame(bed.df,
                                        seqnames.field = "seqnames",
                                        start.field = "start",
                                        end.field = "end",
                                        starts.in.df.are.0based = T)
      bed.gr
    })
    # Find causal snps within adult peaks
    for(j in 1:length(beds.grlist)){
      adult.gr = unique(beds.grlist[[j]])
      snp.gr = makeGRangesFromDataFrame(results.df[i,],
                                        seqnames.field = "snpChr.hg38",
                                        start.field = "snpPos.hg38",
                                        end.field = "snpPos.hg38",
                                        keep.extra.columns = T)
      hits = findOverlaps(snp.gr, adult.gr)
      if(length(hits) != 0){
        results.df[i,]$foundIn = "both"
        results.df[i,]$foundInXpeaksets = results.df[i,]$foundInXpeaksets + 1
      }
    }
    results.df[i,]$totalAdultPeaksets = length(beds.grlist)
  }
}

results.df[results.df$adult_celltype == "unassigned",]$foundIn = "notChecked"
table(results.df$foundIn)
# both  fetalOnly notChecked 
#   91         17         74
nrow(results.df) # 182

# Only if the snp in motif is found in at least 2 adult peaksets, then it's considered to be found in adult
results.df = results.df %>% 
  mutate(foundInMin2Adults = ifelse(foundInXpeaksets > 1, "both", "fetalOnly"))
results.df[results.df$adult_celltype == "unassigned",]$foundInMin2Adults = "notChecked"
table(results.df$foundInMin2Adults)
# both  fetalOnly notChecked 
#   80         28         74

fwrite(results.df, "v3motifs_w_causal_snps_all_tissues_disease_only_pipsum_5_fdr_0.05_fetalAdult.txt",
       row.names = F, sep = "\t")

# Append RS IDs
cdb.full = fread(here::here("data/external/variants/20240508_causaldb_wMeta.txt.gz"),
                 nThread = 12, data.table = F)
results.df = fread("v3motifs_w_causal_snps_all_tissues_disease_only_pipsum_5_fdr_0.05_fetalAdult.txt",
                   data.table = F, sep = "\t")
results.df$chr = as.numeric(gsub("chr", "", results.df$snpChr.hg19))
results.df = left_join(results.df, cdb.full[,c("meta_id", "mesh_id", "chr", "bp", "rsid", "susie")],
                       by = c("meta_id", "mesh_id", "chr", "susie", "snpPos.hg19" = "bp"))
results.df$chr = NULL
results.df = distinct(results.df)
results.df = results.df %>%
  relocate(rsid, .after = snpPos.hg19)

length(unique(results.df$motifName)) # 37
length(unique(results.df$meta_id)) # 37
length(unique(results.df$rsid)) # 43
length(unique(results.df$mesh_id)) # 25

fwrite(results.df,
       "v3motifs_w_causal_snps_all_tissues_disease_only_pipsum_5_fdr_0.05_fetalAdult.txt",
       row.names = F, sep = "\t")
