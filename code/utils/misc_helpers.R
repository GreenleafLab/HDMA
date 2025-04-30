
#####################################
# Various helper functions in R
#####################################

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(tidyr)
  library(magrittr)
  library(Rcpp)
  library(RcppAlgos)
})


#####################################
# Miscellaneous
#####################################

`%ni%` <- Negate(`%in%`)

makeBins <- function(x, breaks=NULL, bins=10, bin.type="width"){
  # Return vector of labels for bins made from x
  ########################
  # x = numeric vector
  # breaks = manually set break values (signpost style, i.e. n+1 breaks for n bins)
  # bins = number of bins to return
  # bin.type = type of bins ('width' = equal width bins, 'size' = bins contain same number of elements)

  if(is.null(breaks)){
    stopifnot(bin.type %in% c('width', 'size'))
    if(bin.type == 'width') breaks <- seq(min(x), max(x), length.out=bins+1)
    if(bin.type == 'size') breaks <- quantile(x, probs=c(0,1/bins*c(1:bins)))
  }
  ids <- as.numeric(cut(x, breaks=breaks, include.lowest=TRUE))
  return(list(ids=ids, breaks=breaks))
}

invertList <- function(lst){
  # Swap names and values in list
  split(rep(names(lst), lengths(lst)), unlist(lst))
}

confusionMatrix <- function(i = NULL, j = NULL){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

getFreqs <- function(x){
  # Return a named vector of frequencies of x
  tab <- table(x) %>% as.data.frame.table()
  frqs <- tab$Freq
  names(frqs) <- tab[,1]
  frqs[order(frqs, decreasing=TRUE)]
}

jaccardIndex <- function(mat, i, j){
  # Calculate Jaccard Index between row i and column j in matrix mat
  # (Matrix is an intersection matrix of categories in rows i and columns j)
  AiB <- mat[i,j]
  AuB <- sum(mat[i,]) + sum(mat[,j]) - AiB
  AiB/AuB
}


#####################################
# Misc Genomics Tools
#####################################

convertMouseGeneList <- function(genes){
  # Function to convert mouse to human gene names
  # see: https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
  require("biomaRt")
  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
   
  genesV2 <- getLDS(attributes=c("mgi_symbol"), filters="mgi_symbol", values=genes, 
    mart=mouse, attributesL=c("hgnc_symbol"), martL=human, uniqueRows=TRUE)

  #humanx <- unique(genesV2[, 2])
  #return(humanx)
  return(genesV2)
}


gcContent <- function(gr, genome){
  # Calculate GC content for a provided genomic range and a corresponding BSgenome
  ##########################################################################
  # gr = genomic range
  # genome = BSgenome object
  freqs <- Biostrings::alphabetFrequency(Biostrings::getSeq(genome, gr))[,1:4]
  gc <- apply(freqs, 1, function(x){sum(x[c("C","G")])/sum(x)})
  gc
}


trim_N_seqs <- function(gr, genome, thresh=0){
  # Identify and remove ranges that contain any N values
  ##########################################################################
  # gr = genomic range
  # genome = BSgenome object
  # thresh = fraction of seq allowed to contain N's (e.g. 0.2 = 20% N's is acceptable)
  N_freq <- Biostrings::alphabetFrequency(Biostrings::getSeq(genome, gr))[,15] # 15th column contains NA counts
  N_frac <- N_freq / width(gr)
  gr[N_frac <= thresh]
}


trim_oob <- function(gr){
  # Identify and remove ranges that contain 'out of bounds' regions
  # (gr must contain seqinfo!)
  ##########################################################################
  # gr = genomic range
  idx <- GenomicRanges:::get_out_of_bound_index(gr)
  if(length(idx) != 0L){
    gr <- gr[-idx]
  }
  return(gr)
}


#####################################
# Sampling cells
#####################################

numCellsPerGroup <- function(groups, total.cells=NULL, cells.per.group=NULL, min.cells.per.group=50, warn.only=TRUE){

  # Number of cells to downsample from each group
  # Returns a table of number of cells to downsample from each group
  # Adapted from https://github.com/cssmillie/ulcerative_colitis/blob/master/downsample.r
  ######################
  # groups = vector of group membership of each cell barcode 
  # total.cells = total number of cells to downsample to
  # cells.per.group = number of cells per group to downsample

  if(is.null(total.cells) & is.null(cells.per.group)){
    stop("total.cells and cells.per.group cannot both be NULL!")
  }

  ncells <- table(groups) %>% sort() # Sorted low to high

  # Check if any groups are 'too small'
  if(any(ncells < min.cells.per.group)){
    smallGroups <- paste(names(ncells[ncells < min.cells.per.group]), collapse=';')
    if(warn.only){
      message(sprintf("Warning: groups %s had fewer than %s cells", smallGroups, min.cells.per.group))
    }else{
      stop(sprintf("Error: groups %s had fewer than %s cells. Exiting...", smallGroups, min.cells.per.group))
    }
  }

  # if sampling w/ cells.per.group, enforce thresholds
  if(!is.null(cells.per.group)){
    ncells[ncells > cells.per.group] <- cells.per.group
  }else{
    n <- table(groups) %>% sort()
    if(length(ncells) == 1){
      # In the case of one group, take total.cells from that group
      ncells <- total.cells
      names(ncells) <- names(n)
    }else{
      # For multiple groups, determine how many of each group to take
      # First: get cumulative totals of cells from 0 to the largest group
      u <- c(0, cumsum(n)[1:length(n)-1]) 
      # Second: Try to evenly split cells between groups.
      # If the smallest group isn't large enough, split 
      # subsequent groups into larger groups, appending the missing 
      # difference, and so on
      i <- (total.cells - u) / seq(length(n), 1, -1) < n 
      # If we have more cells than requested, split:
      if(sum(i) > 0){
        ncells[i] <- as.integer(ceiling((total.cells - sum(n[!i]))/sum(i)))
      }
    }
  }
  ncells
}


downsampleByGroup <- function(cells, groups, ngene=NULL, total.cells=NULL, cells.per.group=NULL, replace=FALSE){
  # Downsample cells evenly across groups
  # Adapted from https://github.com/cssmillie/ulcerative_colitis/blob/master/downsample.r
  #####################
  # cells = vector of cell names
  # groups = vector of cell groups
  # ngene = vector of number of genes detected per cell (if provided, will prioritize cells with more genes)
  # total.cells = total number of cells to downsample
  # cells.per.group = how many cells per group (cannot use with total.cells)
  # repalce = should cells be sampled with replacement

  if(is.null(total.cells) & is.null(cells.per.group)){
    stop("total.cells and cells.per.group cannot both be NULL!")
  }

  # Set ngene to ones if not provided
  if(is.null(ngene)){
    ngene <- structure(rep(1, length(cells)), names=cells)
  }
  if(is.null(names(ngene))){
    names(ngene) <- cells
  }

  # Calculate group sizes to sample
  groups <- as.factor(groups)
  cellsPerGroup <- numCellsPerGroup(groups, total.cells=total.cells, cells.per.group=cells.per.group)

  # Downsample cells within group
  downSampled <- sapply(levels(groups), function(x){

    # Shuffle cells within group
    gcells <- cells[groups == x]
    gcells <- sample(gcells, size=length(gcells), replace=FALSE)

    # Subsample, prioritizing high gene counts (random sample if ngene not provided)
    gcells[order(ngene[gcells], decreasing=TRUE)[1:cellsPerGroup[[x]]]]
    })

  # Combine and return
  downSampled <- unlist(downSampled) %>% na.omit() %>% as.character()
  return(downSampled)
}


selectCells <- function(covars, groupby=NULL, ngene=NULL, cells.use=NULL, max.cells=NULL, replace=FALSE){
  # Select cells from groups defined by covariates matrix
  # Adapted from https://github.com/cssmillie/ulcerative_colitis/blob/master/markers.r
  #######################
  # covars = data.frame of cells x covariates
  # groupby = covariate groups to sample by (i.e. vector of column names to sample evenly across)
  # cells.use = vector of cells allowed to be used
  # max.cells = maximum number of cells to subset from each group
  # repalce = should cells be sampled with replacement

  # Cells to use:
  cells <- rownames(covars)
  if(!is.null(cells.use)){
    cells <- intersect(cells, cells.use)
  }

  # Select cells evenly across groups (if we are subsampling cells)
  if(!is.null(max.cells)){

    # Prepare cell groups (i.e. define all unique combinations of groupby covariates)
    if(is.null(groupby)){
      # If groupby is not provided, will use all non-numeric covariates
      groupby <- sapply(covars, function(x) !is.numeric(x)) 
    }
    groups <- as.factor(apply(covars[cells, groupby, drop=FALSE], 1, function(x) paste(x, collapse='.')))
    names(groups) <- cells

    # Perform downsampling (prioritizes cells with more genes if ngene is provided)
    cells <- downsampleByGroup(cells, groups, ngene=ngene, total.cells=max.cells, replace=replace)

    # Print how many cells in each group:
    print(table(groups[cells]))
  }
  cells
}



#####################################
# GWAS and related tools
#####################################

parseVCF <- function(vcf_file, snpGR){
  # Read in a previously generated bcftools vcf file
  # The header of the VCF file has lines denoted by '##'. fread simply reads through these somehow.
  # The meat of the vcf format is a tab delimited text file with the following columns:
  # 1: Chromosome name
  # 2: 1-based position on the chromosome
  # 3. ID - identifier where available
  # 4: Reference base at this position
  # 5: Alternate base(s) at this position
  # 6: quality (phred-scaled quality score for assertion made in alt)
  # 7: filter status
  # 8: Info - additional information
  # 9: format - indicates the format of coming sample columns
  # 10+ : individual sample columns
  #######################################################################
  # Arguments provided:
  # vcf_file = previously generated vcf file
  # snpGR = snp genomic range
  require(vcfR)

  vcf <- read.vcfR(vcf_file, verbose=TRUE)
  vcf.fix <- as.data.frame(vcf@fix)
  chrs <- vcf.fix$CHROM
  positions <- as.integer(vcf.fix$POS)
  vcfGR <- GRanges(seqnames=chrs, IRanges(positions, end=positions))
  vcfGR$ref <- vcf.fix$REF
  vcfGR$alt <- vcf.fix$ALT
  vcfGR$linked_SNP <- snpGR$linked_SNP[findOverlaps(vcfGR, snpGR, maxgap=-1L, select="first", type="equal", ignore.strand=TRUE)]
  vcfGR$linked_refalt <- snpGR$linked_refalt[findOverlaps(vcfGR, snpGR, maxgap=-1L, select="first", type="equal", ignore.strand=TRUE)]

  # Add genotype data to GR
  gtdf <- as.data.table(vcf@gt)
  gtdf <- gtdf[,2:ncol(gtdf)] # drop format column that we don't care about

  # Rename samples to match expected sample names
  colnames(gtdf) <- sapply(colnames(gtdf), function(x){
    strsplit(x, split="_")[[1]] %>% head(.,-1) %>% paste(.,collapse="_")
    })

  # Remove other information we don't want
  gtdf <- apply(gtdf, c(1,2), function(x) strsplit(x, split=":")[[1]][1])
  # A '.' indicates where a sample was unable to be genotyped
  gtdf <- apply(gtdf, c(1,2), function(x) replace(x, grep("[.]", x), NA)) 
  gtdf

  values(vcfGR) <- cbind(values(vcfGR), gtdf)
  vcfGR
}



#####################################
# Milo differential accessibility
#####################################

runMilo <- function(obj, dim_reduction, k=30, prop=0.1, min_cells_per_samp=50, contrast="diseaseStatus", sample_name="Sample"){
  # Run all the necessary steps of Milo KNN graph-based differential abundance testing
  #####################################################
  # obj = A clustered seurat object
  # dim_reduction = the name of the dim_reduction to use for Milo
  # k = the number of neighbors for building KNN graph (30 is default)
  # prop = Sampling proportion for building neighborhoods (0.1 is recommended for dataselts < 30k cells)
  # min_cells_per_samp = minimum number of cells present in a given sample to use for Milo

  require(miloR)
  require(SingleCellExperiment)

  # Subset seurat object to remove samples that are too lowly represented:
  samp_freqs <- getFreqs(obj$Sample)
  valid_samps <- samp_freqs[samp_freqs > min_cells_per_samp] %>% names()
  Idents(rna_proj) <- "Sample"
  sub_obj <- subset(rna_proj, idents = c(valid_samps))

  # Create Milo object
  obj_milo <- Milo(as.SingleCellExperiment(sub_obj))

  # build KNN graph using LSI reduced dimensions
  d <- length(obj@reductions[[dim_reduction]])
  k <- 30 # 30 is default
  reduced.dims <- toupper(dim_reduction) # Milo capitalizes the dim reduc names
  obj_milo <- buildGraph(obj_milo, k=k, d=d, reduced.dim=reduced.dims)

  # Defining representative neighborhoods on the KNN graph
  prop <- prop # "We suggest using prop=0.1 for datasets of less than 30k cells. For bigger datasets using prop=0.05 should be sufficient (and makes computation faster)"
  obj_milo <- makeNhoods(obj_milo, prop=prop, k=k, d=d, refined=TRUE, reduced_dims=reduced.dims)

  # Counting cells in neighborhoods (to get how many cells from each sample are in each neighborhood)
  obj_milo <- countCells(obj_milo, meta.data=as.data.frame(colData(obj_milo)), sample=sample_name)

  # Defining experimental design
  # "We implement this hypothesis testing in a generalized linear model (GLM) framework, 
  # specifically using the Negative Binomial GLM implementation in edgeR.""
  milo_design <- data.frame(colData(obj_milo))[,c(sample_name, contrast)]
  colnames(milo_design) <- c("Sample", "contrast")
  milo_design <- distinct(milo_design)
  rownames(milo_design) <- milo_design$Sample

  # Computing neighborhood connectivity
  obj_milo <- calcNhoodDistance(obj_milo, d=d, reduced.dim=reduced.dims)

  # Testing differential abundance
  da_results <- testNhoods(obj_milo, design = ~ contrast, design.df=milo_design, reduced.dim=reduced.dims)

  # Visualize results:
  obj_milo <- buildNhoodGraph(obj_milo)

  return(list(milo_obj=obj_milo, da_results=da_results))
}


#################################################################################











