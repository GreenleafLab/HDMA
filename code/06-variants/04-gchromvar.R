options(scipen = 999)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Matrix)
library(parallel)
library(ComplexHeatmap)
library(here)

out <- here::here("output/06-variants/04/")
figout <- here::here("figures/06-variants/04/")
dir.create(out, showWarnings = F, recursive=T)
dir.create(figout, showWarnings = F, recursive=T)
setwd(out)


ch = import.chain(here::here("data/external/variants/hg38ToHg19.over.chain"))

#### Create bed files of scores ####
dir.create(paste0(out, "/score_bed_files"), showWarnings = F)
dir.create(paste0(out, "/hg19_inputs"), showWarnings = F)

cdb.full = fread(here::here("data/external/variants/20240508_causaldb_wMeta.txt.gz"),
                 nThread = 8, data.table = F)
# Split by trait
cdb = cdb.full %>%
  mutate(start = as.numeric(bp)-1,
         trait_name = paste0(meta_id, "_", mesh_id, "_", trait)) %>%
  dplyr::select(chr, start, end = bp, block_id, susie, meta_id, mesh_id, trait_name)
cdb$chr = paste0("chr", cdb$chr)
cdb$block_id = paste0("region", cdb$block_id)
rm(cdb.full) ; gc()

convert_symbols = function(input_string) {
  result_string = gsub("'s", "s", input_string)
  result_string = gsub("[^a-zA-Z0-9]", " ", result_string)
  result_string = gsub("\\s+", "_", result_string)
  result_string = gsub("_$", "", result_string)
  return(result_string)
}

cdb.split = mclapply(split(cdb, cdb$meta_id), function(trait_subset){
  trait_subset$chrnum = as.numeric(gsub("chr", "", trait_subset$chr))
  trait_subset$start = as.numeric(trait_subset$start)
  trait_subset$end = as.numeric(trait_subset$end)
  trait_subset = arrange(trait_subset, chrnum, start, end)
  trait_subset$trait_name = convert_symbols(trait_subset$trait_name)
  filename = trait_subset$trait_name[1]
  trait_subset = trait_subset[,1:5]
  fwrite(trait_subset, paste0(out, "/score_bed_files/",
                              filename, ".bed"), sep = "\t", row.names = F, col.names = F)
}, mc.cores = 24)
rm(cdb.split) ; gc()

# Peaks at L2 level
prepare_SE_scores = function(tissue){
  cat("Formatting peaks...\n")
  tissue_peaks = fread(paste0(here::here("output/06-variants/00/peaks/"), tissue, ".bed.gz"),
                       data.table = F, sep = "\t", col.names = c("chr", "start", "end"))
  # Append counts
  cts_matrix = fread(paste0(here::here("output/06-variants/00/count_matrices/"), tissue, ".tsv.gz"),
                     data.table = F, sep = "\t")
  pks_with_cts = cbind(tissue_peaks, cts_matrix)
  pks_with_cts$site = paste0(pks_with_cts$chr, "_", pks_with_cts$start, "_", pks_with_cts$end)
  pks_with_cts$idx = NULL
  pks_with_cts.gr = makeGRangesFromDataFrame(pks_with_cts,
                                             seqnames.field = "chr",
                                             start.field =  "start",
                                             end.field = "end",
                                             starts.in.df.are.0based = T,
                                             keep.extra.columns = T)
  
  cat("Lifting peaks to hg19...\n")
  # liftOver from hg38 to hg19 since causalDB is in hg19
  pks_with_cts.hg19 = unlist(liftOver(pks_with_cts.gr, ch))
  pks_with_cts.hg19 = pks_with_cts.hg19[!(duplicated(pks_with_cts.hg19$site) | duplicated(pks_with_cts.hg19$site, fromLast = T)),]
  
  names(pks_with_cts.hg19) = pks_with_cts.hg19$site
  pks_with_cts.hg19.df = as.data.frame(pks_with_cts.hg19)
  
  counts = pks_with_cts.hg19.df[,-c(1:5)]
  counts$site = NULL
  counts = as.matrix(counts)
  
  cat("Creating SummarizedExperiment object...\n")
  SE = SummarizedExperiment(assays = list(counts = counts),
                            rowData = pks_with_cts.hg19,
                            colData = DataFrame(names = colnames(counts)))
  SE = addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
  
  cat("Importing bed scores for traits...\n")
  files = list.files(path = paste0(out, "/score_bed_files"), full.names = T)
  scores = importBedScore(rowRanges(SE), files, colidx = 5) # Score is in the 5th column
  
  cat("Saving...\n")
  save(SE, scores, file = paste0(out, "/hg19_inputs/", tissue, "_SE_scores.RData"))
}

tissues = c("Adrenal", "Brain", "Eye", "Heart", "Liver", "Lung", "Muscle",
            "StomachEsophagus", "Skin", "Spleen", "Thyroid", "Thymus")

cat("Preparing inputs...\n")
for(tis in tissues){
  cat(tis, "\n")
  prepare_SE_scores(tis)
}

# Since running gchromvar requires 50 iterations each, do qsub for each iteration
# Then combine results afterwards
# See 04-gchromvar_unified_bg_run_qsub.R


## Collate results --------------------------------------------------------------------
# Piece togther results of individual iterations
# Use sum of posterior probability scores (PIP) cut-off of 5
# This removes underpowered traits
# and those with a lot of -1 PIP scores, i.e. traits with a lot of variants
# that are not given scores in CAUSALdb due to lack of LD info
collate_results = function(tissue, nRuns, pipsum_cutoff){
  cat(tissue, "\n")
  load(paste0(output, "/hg19_inputs/", tissue, "_SE_scores.RData"))
  
  gchromVAR.df = lapply(1:nRuns, function(i) {
    cat(i, "\n")
    mdf = fread(paste0(tissue, "/gchromvar_results_part", i, ".txt"), data.table = F)
    mdf
  }) %>%
    rbindlist() %>%
    as.data.frame() %>%
    group_by(Var1, Var2) %>%
    summarize(zscore = mean(value))
  
  # Filter traits
  trait_pipsums = colSums(scores@assays@data@listData$weights) %>%
    sort()
  gchromVAR.df = gchromVAR.df %>%
    dplyr::filter(Var2 %in% names(trait_pipsums[trait_pipsums > pipsum_cutoff]))
  nCelltypes = length(unique(gchromVAR.df$Var1))
  
  # Calculate P values and do multiple testing correction
  gchromVAR.df = gchromVAR.df %>%
    dplyr::mutate(pvalue = pnorm(zscore, lower.tail = F),
                  bonf = pvalue * dim(scores)[2] * nCelltypes)
  gchromVAR.df$fdr = p.adjust(gchromVAR.df$pvalue, method = "BH")
  gchromVAR.df = gchromVAR.df %>%
    dplyr::select(Var2, Var1, zscore, pvalue, bonf, fdr) %>%
    dplyr::rename(trait = Var2, celltype = Var1) %>%
    arrange(-zscore)
  
  fwrite(gchromVAR.df, paste0("gchromvar_results_",
                              tissue, "_pipsum", pipsum_cutoff, ".txt"),
         sep = "\t", row.names = F)
  system(paste0("gzip -f gchromvar_results_",
                tissue, "_pipsum", pipsum_cutoff, ".txt"))
}

# Do 50 iterations
for(tis in tissues){
  collate_results(tis, 50, 5)
}

## Collate all sig results -----------------------------------------------------
# Load list for subset of disease traits relevant to fetal atlas organs
disease_list = fread(here::here("data/external/variants/disease_traits_subset.tsv"), sep = "\t", data.table = F)
disease_ids = unique(disease_list$study_id)

pipsum = 5
fdr_cutoff = 0.05

res.df = lapply(tissues, function(tis){
  df = fread(paste0("gchromvar_results_", tis, "_pipsum", pipsum, ".txt.gz"), sep = "\t", data.table = F)
  df$tissue = tis
  df
}) %>% bind_rows()
res.df = res.df[!is.na(res.df$pvalue),]

# Renaming for standardization
res.df$tissue = gsub("StomachEsophagus", "Stomach", res.df$tissue)

is.sig.fdr = res.df[res.df$fdr < fdr_cutoff,]
cat("FDR", fdr_cutoff, ":", nrow(is.sig.fdr), "\n") # FDR 0.05 : 7088
fwrite(is.sig.fdr, paste0("sig_results_all_tissues_pipsum_", pipsum, "_fdr_", fdr_cutoff, ".txt"),
       row.names = F, sep = "\t")
system(paste0("gzip -f sig_results_all_tissues_pipsum_", pipsum, "_fdr_", fdr_cutoff, ".txt"))

# Filter for disease only
is.sig.disease = is.sig.fdr %>%
  mutate(mesh_id = word(trait, 2, 2, sep = "_"),
         studyid = word(trait, 1, 1, sep = "_")) %>%
  filter(studyid %in% disease_ids)
cat("Disease only FDR", fdr_cutoff, ":", nrow(is.sig.disease), "\n") # Disease only FDR 0.05 : 303
fwrite(is.sig.disease, paste0("sig_results_all_tissues_disease_only_pipsum_", pipsum, "_fdr_", fdr_cutoff, ".txt"),
       row.names = F, sep = "\t")
system(paste0("gzip -f sig_results_all_tissues_disease_only_pipsum_", pipsum, "_fdr_", fdr_cutoff, ".txt"))


#### Some tidying up for sup table ---------------------------------------------
sig.results = fread("sig_results_all_tissues_pipsum_5_fdr_0.05.txt.gz",
                    data.table = F)
sig.results$uniq_celltype = paste0(sig.results$tissue, "_", sig.results$celltype)
sig.results$mesh_id = word(sig.results$trait, 2, 2, sep = "_")
length(unique(sig.results$trait)) # 1483
length(unique(sig.results$uniq_celltype)) # 131
length(unique(sig.results$mesh_id)) # 349
# Save for supplementary table
sig.results.tab = sig.results %>%
  mutate(meta_ID = word(trait, 1, 1, sep = "_"),
         trait_name = word(trait, 3, -1, sep = "_")) %>%
  select(meta_ID, MeSH_ID = mesh_id, trait_name, celltype, zscore, pvalue, fdr, tissue)
fwrite(sig.results.tab, "gchromvar_all_tissues_pipsum_5_fdr_0.05.txt",
       sep = "\t", row.names = F)

sig.results = fread("sig_results_all_tissues_disease_only_pipsum_5_fdr_0.05.txt.gz",
                    data.table = F)
sig.results$uniq_celltype = paste0(sig.results$tissue, "_", sig.results$celltype)
sig.results$mesh_id = word(sig.results$trait, 2, 2, sep = "_")
length(unique(sig.results$trait)) # 68
length(unique(sig.results$uniq_celltype)) # 79
length(unique(sig.results$mesh_id)) # 36

#### Plotting ------------------------------------------------------------------
load(here::here("code/utils/color_scheme.RData"))
plotTheme = theme(plot.title = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5),
                  panel.grid = element_blank())
tissues = c("Adrenal", "Brain", "Eye", "Heart", "Liver", "Lung", "Muscle",
            "Skin", "Spleen", "Stomach", "Thyroid", "Thymus")

##### Condense results (for main figure plot) ----------------------------------
sf = "sig_results_all_tissues_disease_only_pipsum_5_fdr_0.05.txt.gz"
is.sig.disease = fread(sf, sep = "\t", data.table = F)
is.sig.disease$meta_id = word(is.sig.disease$trait, 1, 1, sep = "_")

# Append disease-relevant organ information (was manually assigned)
trait2organ = fread(here::here("data/external/variants/standardized_trait_list.txt"),
                    sep = "\t", data.table = F) %>%
  select(meta_id, mesh_id = MeSH_id, relevantOrgan = relevantTissue) %>%
  distinct()

is.sig.disease = left_join(is.sig.disease, trait2organ, by = c("meta_id", "mesh_id"))
# Condense to 1 study per trait name by picking the one with the highest mean Z score
oneStudyPerTrait = is.sig.disease
oneStudyPerTrait = oneStudyPerTrait[grep("Age_hay_fever_rhinitis_or_eczema_diagnosed",
                                         oneStudyPerTrait$trait, invert = T),]

# Some editing of trait name labels
{
  oneStudyPerTrait$trait_label = word(oneStudyPerTrait$trait, 3, -1, sep = "_")
  oneStudyPerTrait$trait_label = gsub("_", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("I10 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("I15 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("I20 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("I25 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("I30 I52 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("C43 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("C44 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("I48 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("E78 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("J40 J47 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("J45 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("E03 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("M72 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("H25 H28 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("E70 E90 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("E00 E07 ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("  ", " ", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("  ", " ", oneStudyPerTrait$trait_label)
  # Rename or truncate some long labels
  oneStudyPerTrait$trait_label = gsub("^hypertension", "Hypertension", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("thyroid problem not cancer", "Noncancerous thyroid problem", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Non cancer illness code self reported eczema dermatitis", "Eczema or dermatitis", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub(" and other lipidemias", "", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub(" and other lipidaemias", "", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Doctor diagnosed ", "", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub(" MTAG", "", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub(" stage 1 and 2", "", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub(" multi trait analysis", "", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Glaucoma primary open angle", "Primary open angle glaucoma", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("hayfever allergic rhinitis", "Allergic rhinitis", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("hayfever or allergic rhinitis", "Allergic rhinitis", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Hayfever or allergic rhinitis", "Allergic rhinitis", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Disorders of lens", "Lens disorders", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("D006973 hypertension", "D006973 Hypertension", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Other malignant neoplasms of skin", "Other malignant skin neoplasms", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Melanoma and other malignant neoplasms of skin", "Melanoma", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub(" or family history of alzheimers disease", "", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Alzheimers disease", "Alzheimer's disease", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Other forms of heart disease", "Other heart diseases", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Contracture of palmar fascia Dupuytrens disease", "Dupuytren's disease", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Other hypothyroidism", "Hypothyroidism", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("asthma", "Asthma", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("hypothyroidism", "Hypothyroidism", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Hypothyroidism myxoedema", "Hypothyroidism", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("Other non epithelial cancer of skin", "Non-epithelial skin cancer", oneStudyPerTrait$trait_label)
  oneStudyPerTrait$trait_label = gsub("^ ", "", oneStudyPerTrait$trait_label)
}

meta.to.keep = oneStudyPerTrait %>%
  group_by(trait_label) %>%
  mutate(avg_z = mean(zscore)) %>%
  ungroup() %>%
  dplyr::arrange(desc(avg_z)) %>%
  group_by(trait_label) %>%
  slice_head(n = 1)
oneStudyPerTrait = oneStudyPerTrait[oneStudyPerTrait$studyid %in% meta.to.keep$studyid,]
oneStudyPerTrait = dplyr::rename(oneStudyPerTrait, Organ = tissue)

# Add broader cell type annotations
broad_celltype = fread(here::here("output/05-misc/03/TableS2_cluster_meta_qc.csv"),
                       data.table = F, sep = ",") %>%
  select(Organ = organ, celltype = L2_annot, broad_cell_type = L3_annot, organ_color, compartment_color) %>%
  distinct()
broad_celltype$Organ = gsub("StomachEsophagus", "Stomach", broad_celltype$Organ)
broad_celltype$celltype = word(broad_celltype$celltype, 2, -1, sep = " ")
broad_celltype$celltype = gsub(" ", ".", broad_celltype$celltype, fixed = T)
broad_celltype$celltype = gsub("(", ".", broad_celltype$celltype, fixed = T)
broad_celltype$celltype = gsub(")", ".", broad_celltype$celltype, fixed = T)
broad_celltype$celltype = gsub("/", ".", broad_celltype$celltype, fixed = T)
broad_celltype$celltype = gsub("+", ".", broad_celltype$celltype, fixed = T)
oneStudyPerTrait = left_join(oneStudyPerTrait, broad_celltype, by = c("Organ", "celltype"))

# If there is a cell type that only belongs to one organ, group it as "specialized"
check_broad_if_specialized = broad_celltype %>%
  select(Organ, broad_cell_type) %>%
  distinct()
tbl.broad_ct_n = janitor::tabyl(check_broad_if_specialized, var1 = broad_cell_type)
specialized = tbl.broad_ct_n[tbl.broad_ct_n$n == 1,]
oneStudyPerTrait[oneStudyPerTrait$broad_cell_type %in% specialized$broad_cell_type,]$broad_cell_type = "specialized"

# Some editing of broad cell type labels
capFirstLetter = function(x) {
  splitphrase = strsplit(x, " ")[[1]]
  if(length(splitphrase) > 1){
    firstword = str_to_title(splitphrase[1])
    remainder = paste(splitphrase[2:length(splitphrase)], collapse = " ")
    cappedString = paste0(firstword, " ", remainder)
    return(cappedString)
  } else {
    cappedString = str_to_title(x)
    return(cappedString)
  }
}
oneStudyPerTrait$broad_cell_type = unlist(lapply(oneStudyPerTrait$broad_cell_type, function(x) capFirstLetter(x)))

# Edit cell type names
{
  oneStudyPerTrait$Cell_type = gsub(".", " ", oneStudyPerTrait$celltype, fixed = T)
  oneStudyPerTrait$Cell_type = gsub("  ", " ", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub(" $", "", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("macrophages", "macrophage", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("Macrophage", "macrophage", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("fibroblasts", "fibroblast", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("Fibroblast", "fibroblast", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("Pericyte", "pericyte", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("keratinocytes", "keratinocyte", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("smooth muscle$", "smooth muscle cell", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("cells$", "cell", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("NK T", "NK/T", oneStudyPerTrait$Cell_type, fixed = T)
  oneStudyPerTrait$Cell_type = gsub("Kupffer", "Kupffer cell", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("^endothelial$", "endothelial cell", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("Endothelial", "endothelial", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("vascular endothelial", "vascular endothelial cell", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("stellate", "stellate cell", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("endocardial", "endocardial cell", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("epicardial", "epicardial cell", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("MYH3 smooth muscle", "MYH3+ smooth muscle", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("aCM", "atrial cardiomyocyte", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("vCM", "ventricle cardiomyocyte", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("cycling vCM", "Cycling ventricle cardiomyocyte", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("medullary TEC myoid Type I$", "medullary TEC (myoid Type I)", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("medullary TEC myoid Type IIa", "medullary TEC (myoid Type IIa)", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("myocytes slow Type I", "myocyte (slow/Type I)", oneStudyPerTrait$Cell_type)
  oneStudyPerTrait$Cell_type = gsub("myocytes fast Type II", "myocyte (fast/Type II)", oneStudyPerTrait$Cell_type)
  
}

##### Heatmap ------------------------------------------------------------------
broad = oneStudyPerTrait
broad$Cell_type = paste0(broad$Organ, " ", broad$Cell_type)

# Each cell type should only appear once for the condensed figure
broad = broad %>%
  group_by(Cell_type) %>%
  mutate(total_z = sum(zscore)) %>%
  ungroup() %>%
  dplyr::arrange(desc(total_z)) %>%
  group_by(Cell_type) %>%
  slice_head(n = 1)


## Cluster traits
organs = sort(unique(broad$relevantOrgan))
column_order = c()
for(org in organs){
  print(org)
  column_ordering = broad[broad$relevantOrgan == org,][,c("meta_id", "Organ", "Cell_type", "zscore")]
  column_ordering = column_ordering %>%
    pivot_wider(names_from = "meta_id", values_from = "zscore", values_fill = 0)
  if(ncol(column_ordering)>3){
    dend = hclust(dist(t(column_ordering[,3:ncol(column_ordering)])))
    column_order = c(column_order, dend$label[dend$order])
  } else {
    column_order = c(column_order, names(column_ordering)[3])
  }
}

column_order = data.frame(meta_id = column_order, col_order = 1:length(column_order))
# Turn meta_id to trait labels
column_order = left_join(column_order, distinct(broad[,c("meta_id", "trait_label")]), by = "meta_id")

z_mat = broad %>% 
  pivot_wider(names_from = "trait_label", values_fill = 0,
              values_from = "zscore",
              id_cols = c("broad_cell_type", "Cell_type", "relevantOrgan")) %>%
  arrange(broad_cell_type)
z_mat = select(z_mat, "relevantOrgan", "broad_cell_type", "Cell_type", column_order$trait_label)
in_mat = as.matrix(z_mat[,-c(1:3)])
row.names(in_mat) = z_mat$Cell_type

DiseaseOrganAnnot = distinct(arrange(broad, relevantOrgan)[,c("trait_label", "relevantOrgan")])
top_ha = HeatmapAnnotation("Disease-relevant Organ" = as.factor(DiseaseOrganAnnot$relevantOrgan),
                           show_annotation_name = F, which = 'col',
                           col = list(
                             "Disease-relevant Organ" = c(
                               "Adrenal" = "#876941",
                               "Brain" = "#0C727C",
                               "Eye" = "#ff9f0f",
                               "Heart" = "#D51F26",
                               "Liver" = "#3b46a3",
                               "Lung" = "#f0643e",
                               "Muscle" = "#89C75F",
                               "Skin" = "#ad0773",
                               "Stomach" = "#208A42",
                               "Spleen" = "#3BBCA8",
                               "Thymus" = "#6E4B9E",
                               "Thyroid" = "#8A9FD1"
                             )))
organs = word(z_mat$Cell_type, 1, 1, sep = " ")
right_ha = HeatmapAnnotation("Organ" = as.factor(organs),
                             show_annotation_name = F, which = 'row',
                             col = list(
                               "Organ" = c(
                                 "Adrenal" = "#876941",
                                 "Brain" = "#0C727C",
                                 "Eye" = "#ff9f0f",
                                 "Heart" = "#D51F26",
                                 "Liver" = "#3b46a3",
                                 "Lung" = "#f0643e",
                                 "Muscle" = "#89C75F",
                                 "Skin" = "#ad0773",
                                 "Stomach" = "#208A42",
                                 "Spleen" = "#3BBCA8",
                                 "Thymus" = "#6E4B9E",
                                 "Thyroid" = "#8A9FD1"
                               )))

# Slight editing of broad cell type names
z_mat$broad_cell_type = gsub("Skeletal muscle", "Skeletal\nmuscle", z_mat$broad_cell_type)
z_mat$broad_cell_type = gsub("Myogenic progenitors", "Myogenic\nprogenitors", z_mat$broad_cell_type)
z_mat$broad_cell_type = gsub("Smooth muscle", "Smooth\nmuscle", z_mat$broad_cell_type)
z_mat$broad_cell_type = gsub("Neural crest derived", "Neural crest\nderived", z_mat$broad_cell_type)

svg(paste0(figout, "/heatmap_gchromvar_disease_only_fdr0.05_pipsum5_oneStudyPerTraitName_broad_celltypes.svg"),
    height = 11.5, width = 9.3)
Heatmap(in_mat,
        cluster_rows = T, cluster_columns = T,
        top_annotation = top_ha,
        right_annotation = right_ha,
        column_names_side = "bottom",
        row_names_side = "right",
        show_heatmap_legend = T,
        show_column_dend = F,
        show_row_dend = F,
        row_labels = row.names(in_mat),
        row_names_gp = gpar(fontsize = 9.5),
        row_title_gp = gpar(fontsize = 10), 
        row_split = z_mat$broad_cell_type,
        row_title_rot = 0,
        column_names_gp = gpar(fontsize = 9.5), 
        column_split = DiseaseOrganAnnot$relevantOrgan,
        column_title = "Disease-relevant Organ",
        col = c("#000075", "#FD619D", "#FFCD32"),
        name = "Z score")
dev.off()

###### Transposed version (final figure) ---------------------------------------
left_ha = HeatmapAnnotation("Disease-relevant Organ" = as.factor(DiseaseOrganAnnot$relevantOrgan),
                            show_annotation_name = F, which = 'row',
                            col = list(
                              "Disease-relevant Organ" = c(
                                "Adrenal" = "#876941",
                                "Brain" = "#0C727C",
                                "Eye" = "#ff9f0f",
                                "Heart" = "#D51F26",
                                "Liver" = "#3b46a3",
                                "Lung" = "#f0643e",
                                "Muscle" = "#89C75F",
                                "Skin" = "#ad0773",
                                "Stomach" = "#208A42",
                                "Spleen" = "#3BBCA8",
                                "Thymus" = "#6E4B9E",
                                "Thyroid" = "#8A9FD1"
                              )))

organs = word(z_mat$Cell_type, 1, 1, sep = " ")
bottom_ha = HeatmapAnnotation("Organ" = as.factor(organs),
                              show_annotation_name = F, which = 'column',
                              col = list(
                                "Organ" = c(
                                  "Adrenal" = "#876941",
                                  "Brain" = "#0C727C",
                                  "Eye" = "#ff9f0f",
                                  "Heart" = "#D51F26",
                                  "Liver" = "#3b46a3",
                                  "Lung" = "#f0643e",
                                  "Muscle" = "#89C75F",
                                  "Skin" = "#ad0773",
                                  "Stomach" = "#208A42",
                                  "Spleen" = "#3BBCA8",
                                  "Thymus" = "#6E4B9E",
                                  "Thyroid" = "#8A9FD1"
                                )))

svg(paste0(figout, "/heatmap_gchromvar_disease_only_fdr0.05_pipsum5_oneStudyPerTraitName_broad_celltypes_flipped.svg"),
    height = 7.1, width = 16)
ht = Heatmap(t(in_mat),
             cluster_rows = T, cluster_columns = T,
             left_annotation = left_ha,
             bottom_annotation = bottom_ha,
             column_names_side = "bottom",
             row_names_side = "right",
             show_heatmap_legend = T,
             show_column_dend = F,
             show_row_dend = F,
             border_gp = gpar(col = "black", lwd = 1),
             row_names_gp = gpar(fontsize = 9.5),
             row_split = DiseaseOrganAnnot$relevantOrgan,
             row_title = "Disease-relevant Organ",
             column_names_gp = gpar(fontsize = 9.5), 
             column_split = z_mat$broad_cell_type,
             column_title_gp = gpar(fontsize = 12), 
             row_title_gp = gpar(fontsize = 12), 
             col = c("#000075", "#FD619D", "#FFCD32"),
             name = "Z score")
draw(ht, padding = unit(c(5,1,5,1), "mm")) # b,l,t,r
dev.off()

pdf(paste0(figout, "/heatmap_gchromvar_disease_only_fdr0.05_pipsum5_oneStudyPerTraitName_broad_celltypes_flipped.pdf"),
    height = 7.1, width = 16)
draw(ht, padding = unit(c(5,1,5,1), "mm")) # b,l,t,r
dev.off()
