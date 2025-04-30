#!/usr/bin/env Rscript
# Preprocess scATAC using ArchR -----------------------------------------------------------------

# imports -----------------------------------------------------------------
#Load ArchR (and associated libraries)
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(mclust)
  library(Seurat)
  library(here)
})

# Get additional functions, etc.:
scriptPath <- here("code/utils")
source(paste0(scriptPath, "/plotting_config.R"))
source(paste0(scriptPath, "/archr_helpers.R"))
source(paste0(scriptPath, "/misc_helpers.R"))
source(paste0(scriptPath, "/matrix_helpers.R"))

# user inputs -------------------------------------------------------------
organ_code <- "##"
organ_batch <- "##" #e.g. "11" or c("12", "18")
wd <- here("output/01-preprocessing/02/##ORGAN##/atac_preprocess_output")
dir_input <- here("output/01-preprocessing/01/##ORGAN##/ATAC/")

dir_whitelist <- here("output/01-preprocessing/02/shared/whitelist_r1")
meta.path <- here("output/01-preprocessing/02/shared/SampleMetadata.csv")
cutoffs.path <- here("output/01-preprocessing/02/shared/sample_filter_cutoffs_metadata.csv")

# set working directory
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)

#Load Genome Annotations
pointSize <- 0.5
addArchRGenome("hg38")

# load meta data
meta <- read.csv(meta.path)
rownames(meta) <- meta$SampleName
samp.names <- meta$SampleName[meta$Batch %in% organ_batch]
samp.sex <- meta[samp.names, "Sex"] # this will be empty until we manually fill it after looking at rna
samp.age <- paste0("PCW",meta[samp.names,"PCW"])
names(samp.sex) <- samp.names
names(samp.age) <- samp.names

# main -------------------------------------------------------------------------------------

## Prepare data ----------------------------------------------------------------------------

#Get input files
inputFiles <- list.files(dir_input,
    #recursive = TRUE, 
    full.names = TRUE,
    pattern = "*gz$" #ignores tsv.gz.tbi using the $
    ) 

names(inputFiles) <- gsub(basename(inputFiles), pattern = ".fragments.tsv.gz", replacement = "")

# Create Arrow Files (~30 minutes)
# recommend you use as many threads as samples.
# This step will for each sample :
# 1. Read Accessible Fragments.
# 2. Identify Cells QC Information (TSS Enrichment, Nucleosome info).
# 3. Filter Cells based on QC parameters.
# 4. Create a TileMatrix 500-bp genome-wide.
# 5. Create a GeneScoreMatrix.

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 5,
  minFrags = 1000, # Default is 1000.
  addTileMat = FALSE, # Don't add tile or geneScore matrices yet. Will add them after we filter
  addGeneScoreMat = FALSE
)

#Create ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "unfiltered_output"
)

# Remove Arrow files after copying
unlink(paste0(wd, "/*.arrow"))

saveArchRProject(proj)
proj <- loadArchRProject(paste0(wd, "/unfiltered_output"))

# Add sample metadata
#proj$sex <- samp.sex[proj$Sample] %>% unlist() %>% as.factor()
proj$age <- samp.age[proj$Sample] %>% unlist()

# Now, add tile matrix and gene score matrix to ArchR project
proj <- addTileMatrix(proj, force=TRUE)
proj <- addGeneScoreMatrix(proj, force=TRUE)

# Visualize numeric metadata per grouping with a violin plot now that we have created an ArchR Project.
plotList <- list()

plotList[[1]] <- plotGroups(ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  #pal=samp_cmap,
  plotAs = "violin",
  addBoxPlot = TRUE,
  alpha = 0.4
  )

plotList[[2]] <- plotGroups(ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  #pal=samp_cmap,
  plotAs = "violin",
  addBoxPlot = TRUE,
  alpha = 0.4
  )

plotPDF(plotList = plotList, name = "Violin-TSS-log10(nFrag)", width = 4, height = 4,  ArchRProj = proj, addDOC = FALSE)

## Manually define cutoffs per sample and filter ArchRProject --------------------------------

### Manually enter into TSV the cutoff values for each sample after inspecting QC plots ###
# Load in the manual cutoff values
cutoffs <- read.csv(cutoffs.path, header=T, row.names=1)

# Filter archr project based on manual QC cutoffs
whitelist_cellnames <- c()

sample.list <- proj$Sample %>% unique
for (i in seq_along(sample.list)) {
    sample <- sample.list[i]
    nFragsCutoff <- cutoffs[sample,"atac_nFrags"]
    TSSCutoff <- cutoffs[sample,"atac_TSS"]

    idxPass <- which(
        proj$Sample %in% sample &
        proj$nFrags >= nFragsCutoff &
        proj$TSSEnrichment >= TSSCutoff
        )

    nTotal <- length(which(proj$Sample %in% sample))
    nFiltered <- nTotal - length(idxPass)
    percentFiltered <- (nFiltered/nTotal) * 100

    message("Removing ", round(percentFiltered, 2), "% of cells (", nFiltered, " cells) from ", sample, ", based on TSS cutoff = ", TSSCutoff, " and nFrags cutoff = ", nFragsCutoff)

    cellnames <- proj$cellNames[idxPass]
    whitelist_cellnames <- c(whitelist_cellnames, cellnames)
}

proj.filtered <- proj[whitelist_cellnames, ]
#plotTSSnFrag(proj, cutoffs.path) # replot ATAC QC
plotTSSnFrag(proj.filtered, cutoffs.path) # replot ATAC QC

# Save filtered ArchR Project
saveArchRProject(proj.filtered, outputDirectory = paste0(wd, "/filtered_r1"), dropCells = TRUE)

## Save cell name whitelist --------------------------------------------------------------

df.whitelist <- data.frame(atac_cb = proj.filtered$cellNames, atac_nFrags = proj.filtered$nFrags)
readr::write_tsv(df.whitelist, paste0(dir_whitelist, "/", organ_code, "_atac.txt"))
