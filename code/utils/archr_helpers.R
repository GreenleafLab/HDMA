# Helper functions for ArchR

suppressPackageStartupMessages({
  library(ArchR)
  library(magrittr)
  library(SingleCellExperiment)
})


getMatrixValuesFromProj <- function(proj, matrixName="GeneScoreMatrix", names=NULL, imputeMatrix=FALSE){
  # Return the imputed matrix from an ArchR project
  # Must have already added imputedWeights, etc.
  # Names is a vector of feature names to return. If not provided, will use all available features
  # Warning though: imputing the matrix for all features may take a very long time

  # Returns a summarized experiment:
  se <- getMatrixFromProject(proj, useMatrix = matrixName, binarize = FALSE)

  # Get mat with cell names and row names
  mat <- assays(se)[[matrixName]]
  colnames(mat) <- rownames(colData(se))
  rownames(mat) <- rowData(se)$name # All matrix rowData has a name field

  # Subset by provided names
  if(!is.null(names)){
    # Check if any names are invalid
    validNames <- names[names %in% rownames(mat)]
    if(any(!names %in% validNames)){
      invalidNames <- names[!names %in% validNames]
      message(sprintf("Warning! name(s) %s are not present in matrix!", paste(invalidNames, collapse=',')))
    }
    mat <- mat[validNames,]
  }

  # Impute matrix values
  if(imputeMatrix){
    message("Imputing matrix...")
    imputeWeights <- getImputeWeights(proj)
    mat <- ArchR::imputeMatrix(mat = as.matrix(mat), imputeWeights = imputeWeights)
  }
  mat
}


getClusterPeaks <- function(proj, clusterNames, peakGR=NULL, replicateScoreQuantileCutoff=0, originalScore=FALSE, groupBy = NULL){
  # This function will return the subset of peaks from the full ArchR project that were 
  # initially called on the clusters provided in clusterNames.
  ######################################################################################
  # proj = ArchR project
  # clusterNames = name or names of clusters to pull peaks from. These cluster names must
  #   match the cluster names originally used to call peaks
  # peakGR = the ArchR peak genomic range obtained using 'getPeakSet'. Will return a subset of
  #   these peaks that overlap the peaks originally called using the provided clusters
  # replicateScoreQuantileCutoff = A numeric quantile cutoff for selecting peaks. 
  if(is.null(peakGR)){
    peakGR <- getPeakSet(proj)
  }
  if(!is.null(groupBy)){
    peakDir <- paste0(proj@projectMetadata$outputDirectory, "/PeakCalls/", groupBy)
  } else {peakDir <- paste0(proj@projectMetadata$outputDirectory, "/PeakCalls")}
  
  calledPeaks <- lapply(clusterNames, function(x){
      readRDS(paste0(peakDir, sprintf("/%s-reproduciblePeaks.gr.rds", x)))
    }) %>% as(., "GRangesList") %>% unlist()
  calledPeaks <- calledPeaks[calledPeaks$replicateScoreQuantile >= replicateScoreQuantileCutoff]
  peakGR <- peakGR[overlapsAny(peakGR, calledPeaks)]
  if(originalScore){
    message(sprintf("Getting original scores from clusters..."))
    # Replace 'score' column with the score of the original peak call in this cluster
    # (if multiple clusters, replaces with the maximum score)
    ol <- findOverlaps(peakGR, calledPeaks, type="any", maxgap=0, ignore.strand=TRUE)
    odf <- as.data.frame(ol)
    odf$og_score <- calledPeaks$score[odf$subjectHits]
    score_df <- odf %>% group_by(queryHits) %>% summarize(max_score=max(og_score)) %>% as.data.frame()
    peakGR$score[score_df$queryHits] <- score_df$max_score
  }
  peakGR
}


buildUMAPdfFromArchR <- function(proj, cellColData=NULL, embeddingName="UMAP", 
  useCells=NULL, dataMat=NULL, featureName=NULL, shuffle=TRUE, 
  lowerPctLim=NULL, upperPctLim=NULL){
  # Return a three column UMAP df from an ArchR project
  # If cellColData is not null, return the indicated column
  # dataMat is a pre-populated cell x feature matrix of values to plot. 
  # The featureName indicates which one

  # Get UMAP coordinates first:
  df <- proj@embeddings[[embeddingName]]$df
  if(is.null(useCells)){
    useCells <- rownames(df)
  }
  colnames(df) <- c("UMAP1", "UMAP2")
  df <- df[useCells,] %>% as.data.frame()
  if(!is.null(cellColData)){
    df[,3] <- proj@cellColData[useCells,cellColData] %>% as.vector()
    colnames(df) <- c("UMAP1", "UMAP2", cellColData)
  }
  if(!is.null(dataMat) & !is.null(featureName)){
    df <- merge(df, dataMat, by=0, all=TRUE)
    df <- df[,c("UMAP1", "UMAP2", featureName)] 
  }
  if(shuffle){
    df <- df[sample(nrow(df), replace=FALSE),]
  }
  # Force limits if indicated
  if(!is.null(lowerPctLim)){
    lowerLim <- quantile(df[,3], probs=c(lowerPctLim))
    df[,3][df[,3] <= lowerLim] <- lowerLim
  }
  if(!is.null(upperPctLim)){
    upperLim <- quantile(df[,3], probs=c(upperPctLim))
    df[,3][df[,3] >= upperLim] <- upperLim
  }
  df
}


scoreGeneSet <- function(expr, geneSet){
  # Generate scores for each cell in expr matrix (log2TP10K, genes x cells)
  # See: Smillie et al. Cell 2019

  # Subset expr matrix by genes in geneSet:
  validGenes <- geneSet[geneSet %in% rownames(expr)]
  subExpr <- expr[validGenes,]

  # Remove any genes that have no expression in any cells
  subExpr <- subExpr[rowSums(subExpr) > 0,]

  # Prevent highly expressed genes from dominating gene score signature by
  # scaling each gene by its root mean squared expression
  scaledSubExpr <- subExpr %>% t() %>% scale(., center=FALSE) %>% t()

  # Signature score is the mean scaled expression across all genes in signature
  scores <- colMeans(scaledSubExpr)
  return(scores)
}


# Functions for creating 'low-overlapping aggregates' of cells

computeKNN <- function(data=NULL, query=NULL, k=50, includeSelf=FALSE, ...){
  # Compute KNN for query points (usually a reduced dims matrix)
  # This returns a matrix of indices mapping query to neighbors in data
  # If query has n cells (rows) and k = 50, will be a n x 50 matrix
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}


getLowOverlapAggregates <- function(proj, target.agg=500, k=100, overlapCutoff=0.8, dimReduc="IterativeLSI", seed=1){
  # Generate low-overlapping aggregates of cells
  ##############################################
  # proj = ArchR project
  # target.agg = number of target aggregates (before filtering)
  # k = number of cells per aggreagate
  # overlapCutoff = Maximum allowable overlap between aggregates
  set.seed(seed)

  # Get reduced dims:
  rD <- getReducedDims(proj, reducedDims=dimReduc)

  # Subsample
  idx <- sample(seq_len(nrow(rD)), target.agg, replace = !nrow(rD) >= target.agg)

  # Get KNN Matrix:
  knnObj <- computeKNN(data=rD, query=rD[idx,], k=k)

  # Check whether aggregates pass overlap cutoff
  keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]

  # Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  # Name aggregates and return as a df of cell ids x aggs
  names(knnObj) <- paste0("agg", seq_len(length(knnObj)))
  knnDF <- data.frame(knnObj)[,c(3,2)]
  colnames(knnDF) <- c("cell_name", "group")
  knnDF$cell_name <- as.character(knnDF$cell_name)
  knnDF
}


# Cluster visualization helpers

relabelClusters <- function(proj, clusterName="Clusters"){
  # Relabel clusters to be ordered by cluster size

  ogClusters <- getCellColData(proj)[[clusterName]]
  tabDF <- base::table(ogClusters) %>% as.data.frame
  colnames(tabDF) <- c("Clusters", "count")
  tabDF["NewClusters"] <- rank(-tabDF$count)
  swapVec <- paste0("C", tabDF$NewClusters)
  names(swapVec) <- tabDF$Clusters

  # Now replace cluster names
  newClust <- sapply(ogClusters, function(x) swapVec[x]) %>% unname()
  proj <- addCellColData(proj, data=newClust, name=clusterName, cells=getCellNames(proj), force=TRUE)
  return(proj)
}


visualizeClustering <- function(proj, pointSize=0.75, prefix="", clusterName="Clusters", sampleName="Sample", embedding="UMAP", 
  sampleCmap=NULL, gest_age_cmap=NULL, barwidth=0.9, filename = "Plot-UMAP-Sample-Clusters.pdf"){
  # Plot various clustering results
  # Set colormap
  qualcmap <- cmaps_BOR$stallion
  quantcmap <- cmaps_BOR$solarExtra
  namedSampCmap <- TRUE
  namedAgeCmap <- TRUE

  if(is.null(sampleCmap)){
    sampleCmap <- qualcmap
    namedSampCmap <- FALSE
  }

  if(is.null(gest_age_cmap)){
    gest_age_cmap <- qualcmap
    namedAgeCmap <- FALSE
  }

  # Plot the UMAPs by Sample and Cluster:
  p1 <- plotEmbedding(proj, colorBy="cellColData", name=sampleName, embedding=embedding, plotAs="points", size=pointSize, pal=sampleCmap, labelMeans=FALSE)
  p2 <- plotEmbedding(proj, colorBy="cellColData", name=clusterName, embedding= embedding, plotAs="points", size=pointSize, labelMeans=FALSE)
  p3 <- plotEmbedding(proj, colorBy="cellColData", name="age", embedding=embedding, plotAs="points", size=pointSize, pal=gest_age_cmap, labelMeans=FALSE)

  proj@cellColData$log10nFrags <- log10(proj@cellColData$nFrags)
  p4 <- plotEmbedding(proj, colorBy="cellColData", name="log10nFrags", embedding=embedding, plotAs="points", size=pointSize, labelMeans=FALSE)
  p5 <- plotEmbedding(proj, colorBy="cellColData", name="TSSEnrichment", embedding=embedding, plotAs="points", size=pointSize, labelMeans=FALSE)
  ggAlignPlots(p1,p2,p3,p4,p5, type="h")
  plotPDF(p1,p2,p3,p4,p5, name = paste0(prefix, filename), ArchRProj=proj, addDOC=FALSE, width=5, height=5)

  return(proj)
}


# Functions for working with peak to gene linkages

getP2G_GR <- function(proj, corrCutoff=NULL, varCutoffATAC=0.25, varCutoffRNA=0.25, filtNA=TRUE){
  # Function to get peaks and genes involved in peak to gene links
  # (See: https://github.com/GreenleafLab/ArchR/issues/368)
  ############################################################
  # proj: ArchR project that alreayd has Peak2GeneLinks
  # corrCutoff: minimum numeric peak-to-gene correlation to return
  # varCutoffATAC: minimum variance quantile of the ATAC peak accessibility when selecting links
  # varCutoffRNA: minimum variance quantile of the RNA gene expression when selecting links
  p2gDF <- metadata(proj@peakSet)$Peak2GeneLinks
  p2gDF$symbol <- mcols(metadata(p2gDF)$geneSet)$name[p2gDF$idxRNA] %>% as.character()
  p2gDF$peakName <- (metadata(p2gDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2gDF$idxATAC]
  # Remove peaks with 'NA' correlation values
  if(filtNA){
    p2gDF <- p2gDF[!is.na(p2gDF$Correlation),]
    p2gDF <- p2gDF[!is.na(p2gDF$idxATAC),]
  }
  if(!is.null(corrCutoff)){
    p2gDF <- p2gDF[(p2gDF$Correlation > corrCutoff),]
  }
  # Filter by variance quantile
  p2gDF <- p2gDF[which(p2gDF$VarQATAC > varCutoffATAC & p2gDF$VarQRNA > varCutoffRNA),]
  # The genomic range contains just the peak ranges:
  p2gGR <- metadata(p2gDF)$peakSet[p2gDF$idxATAC]
  mcols(p2gGR) <- p2gDF
  p2gGR
}


grLims <- function(gr){
  # Get the minimum and maximum range from a GR
  if(length(gr) == 0){
    return(NA)
  }
  starts <- start(gr)
  ends <- end(gr)
  c(min(starts, ends), max(starts, ends))
}


getP2Gregions <- function(proj, genes, p2gGR=NULL, corrCutoff=0.4, buffer_space=0.05, min_width=25000, ...) {
  # Function to get regions containing entire peak to gene region,
  # i.e. a GR that contains all peak to gene links
  ###############################################################
  # p2gGR: genomic range containing all peak to gene links
  # genes: vector of genes to look up
  # buffer_space: fraction of total length to expand on each side of region

  # Get gene GR from ArchR project
  geneGR <- promoters(getGenes(proj)) # Promoters gets 2kb upstream and 200bp downstream
  geneGR <- geneGR[!is.na(geneGR$symbol)]

  # if p2gGR not provided, pull it from ArchR project
  if(is.null(p2gGR)){
    p2gGR <- getP2G_GR(proj, corrCutoff=corrCutoff, ...)
  }

  # Now for each gene, construct GR of all loops and gene TSS
  resultGR <- geneGR[match(genes, geneGR$symbol)]
  start(resultGR) <- sapply(resultGR$symbol, function(x){
      min(grLims(resultGR[resultGR$symbol == x]), grLims(p2gGR[p2gGR$symbol == x]), na.rm=TRUE)
    })
  end(resultGR) <- sapply(resultGR$symbol, function(x){
      max(grLims(resultGR[resultGR$symbol == x]), grLims(p2gGR[p2gGR$symbol == x]), na.rm=TRUE)
    })

  # Finally, resize by buffer space
  resultGR <- resize(resultGR, width=width(resultGR) + buffer_space*width(resultGR), fix="center")
  resultGR <- resize(resultGR, width=ifelse(width(resultGR) > min_width, width(resultGR), min_width), fix="center")
  resultGR
}

addSeuratGeneExpressionMatrix <- function(archr.proj, gtf.file, seu.file, assay="RNA"){
  # Function to add gene expression matrix from a Seurat object to an ArchR project
  # archr.proj: an ArchR object
  # gtf.file: path to the .gtf file used for RNA feature assignment during preprocessing
  # seu.file: path to the Seurat object .rds file with gene expression information, can also pass a Seurat object directly
  # assay: assay in the Seurat object to pull counts from, default to "RNA"
  
  # read RNA annotation
  message(sprintf("reading RNA annotations from %s", gtf.file))
  rna.annot <- rtracklayer::import(gtf.file)
  rna.annot <- rna.annot[rna.annot$type=="gene"]
  names(rna.annot) <- rna.annot$gene_name
  is.gene.related <- grep("gene_", colnames(mcols(rna.annot)))
  mcols(rna.annot) <- mcols(rna.annot)[,is.gene.related]
  
  # read genes from Seurat object
  if (typeof(seu.file)=="character"){ # if seu.file is a file path
    message(sprintf("reading Seurat object from %s", seu.file))
    seu <- readRDS(seu.file)
  } else if (typeof(seu.file)=="S4"){ # if seu.file is a Seurat object
    message("loading Seurat object from input")
    seu <- seu.file
  }
  sce <- seu %>% as.SingleCellExperiment(assay=assay)
  rm(seu) # save memory
  message(sprintf("# genes in Seurat object %s assay: %s", assay,dim(sce)[1]))
  
  # read genes from ArchR annotation
  atac.annot<- getArchRGenome(geneAnnotation = T)$gene
  names(atac.annot) <- atac.annot$symbol
  
  # clean up duplicate gene names with added suffix when imported into Seurat e.g. "TBCE" and "TBCE.1", remove the suffix
  dupmask <- !rownames(sce) %in% rna.annot$gene_name
  rownames(sce)[dupmask] <- (lapply(rownames(sce)[dupmask], function(n){strsplit(n,split="\\.")[[1]][1]}) %>% unlist)
  
  # double check that all Seurat gene names are found in rna.annot
  stopifnot("Seurat object has gene names not found in gtf file"=sum(!rownames(sce) %in% names(rna.annot))==0)
  
  # for Seurat object genes present in the ATAC annotation, use ATAC annotation genome ranges
  inatacmask <- rownames(sce) %in% names(atac.annot)
  message(sprintf("# genes in Seurat object %s assay and in ArchR annotation: %s", assay, sum(inatacmask)))
  ranges1 <- atac.annot[rownames(sce)[inatacmask]]
  
  # for Seurat object genes not present in the ATAC annotation, use RNA annotation genome ranges
  # e.g. these are mostly lncRNAs, antisense genes
  message(sprintf("# genes in Seurat object %s assay but not in ArchR annotation: %s", assay, sum(!inatacmask)))
  message(paste0("Examples of genes not in ArchR annotation: "))
  message(paste0(rownames(sce)[!inatacmask] %>% head(20), ","))
  ranges2 <- rna.annot[rownames(sce)[!inatacmask]]
  
  # rename/format/add some meta columns for consistency before merging
  genome(ranges2) <- "hg38"
  ranges2$symbol <- ranges2$gene_name
  ranges2$gene_name <- NULL
  ranges2$ensembl_gene_id <- ranges2$gene_id
  ranges2$gene_id <- NULL
  
  ranges1$gene_type <- rna.annot[names(ranges1)]$gene_type
  ranges1$ensembl_gene_id <- rna.annot[names(ranges1)]$gene_id
  
  # add combined row ranges
  range_combined <- c(ranges1, ranges2)
  rowRanges(sce) <- range_combined[rownames(sce)]
  
  message("adding gene expression matrix to ArchR project")
  archr.proj.gex <- addGeneExpressionMatrix(
    input = archr.proj,
    seRNA = sce,
    excludeChr = c("chrM", "chrY")
  )
  return(archr.proj.gex)
}

addSeuratGeneExpressionMatrixFixedFeat <- function(archr.proj, global.gexfeat.file, seu.file, assay="RNA"){
  # Function to add gene expression matrix from a Seurat object to an ArchR project with a fixed set of features
  #        output is a left join of the fixed set of features with Seurat object features, pad with zeros
  # archr.proj: an ArchR object
  # global.gexfeat.file: path to the GRanges object .rds file with a fixed set of features
  # seu.file: path to the Seurat object .rds file with gene expression information, can also pass a Seurat object directly
  # assay: assay in the Seurat object to pull counts from, default to "RNA"
  
  # read Seurat object
  if (typeof(seu.file)=="character"){ # if seu.file is a file path
    message(sprintf("reading Seurat object from %s", seu.file))
    seu <- readRDS(seu.file)
  } else if (typeof(seu.file)=="S4"){ # if seu.file is a Seurat object
    message("loading Seurat object from input")
    seu <- seu.file
  }
  sce <- seu %>% as.SingleCellExperiment(assay=assay)
  rm(seu) # save memory
  message(sprintf("# genes in Seurat object %s assay: %s", assay,dim(sce)[1]))
  
  # read genes from the fixed/global gene set
  global_gexfeat_ranges <- readRDS(global.gexfeat.file)
  
  # remove genes not in the fixed gene set
  filtermask <- rownames(sce) %in% global_gexfeat_ranges$symbol
  sce <- sce[filtermask]
  message(sprintf("# genes in Seurat object %s assay and the fixed gene set: %s", assay, dim(sce)[1]))
  
  # double check that all Seurat gene names are found in fixed gene set
  stopifnot("Seurat object has gene names not found in the given fixed gene set"=sum(!rownames(sce) %in% names(global_gexfeat_ranges))==0)
  
  # add genes in fixed gene set but not Seurat object, pad with zeros
  new_feat <- global_gexfeat_ranges[!global_gexfeat_ranges$symbol %in% rownames(sce)]
  new_counts <- matrix(0, nrow = length(new_feat), ncol = ncol(counts(sce)))
  rownames(new_counts) <- names(new_feat)
  
  mat <- rbind(counts(sce), new_counts)
  
  new_sce <- SingleCellExperiment(assays = list(counts = mat))
  new_sce <- new_sce[names(global_gexfeat_ranges), colnames(sce)]
  mainExpName(new_sce) <- mainExpName(sce)
  colData(new_sce) <- colData(sce)
  rowRanges(new_sce) <- global_gexfeat_ranges
  
  message("adding gene expression matrix to ArchR project")
  archr.proj.gex <- addGeneExpressionMatrix(
    input = archr.proj,
    seRNA = new_sce,
    excludeChr = c("chrM", "chrY"),
    force = T
  )
  return(archr.proj.gex)
}

plotTSSnFrag <- function(atac_proj, cutoffs.path, allplots=FALSE){
  # plots the TSS vs nFrag scatter QC plot, (optional) TSS Enrichment ratio and fragment size distribution for an ArchR project 
  # atac_proj: an ArchR object, must have "Sample" coldata
  # cutoffs.path: path to the cutoffs file with TSS and nFrag cutoffs for display in plot
  # allplots: default FALSE, if TRUE plot the TSS enrichment ratio and fragment size distribution too, which take ~10min
  # example cutoffs file format, must have columns "atac_TSS" and "atac_nFrags": 
  #           SampleName,atac_TSS,atac_nFrags
  #           Sample1,7,2000
  #           Sample2,6,2000
  
  sample.list <- atac_proj$Sample %>% unique
  cutoffs <- read.csv(cutoffs.path, header=T, row.names=1)
  
  for (i in seq_along(sample.list)) {
    sample <- sample.list[i]
    Metadata <- getCellColData(atac_proj) %>% as.data.frame %>% dplyr::filter(Sample==sample)
    ggtitle <- sprintf("%s\n%s\n%s",
                       paste0(sample, "\nnCells Pass Filter = ", dim(Metadata)[1]),
                       paste0("Median Frags = ", median(Metadata$nFrags)),
                       paste0("Median TSS Enrichment = ", median(Metadata$TSSEnrichment))
    )
    
    gg <- ggPoint(
      x = pmin(log10(Metadata$nFrags), 5) + rnorm(length(Metadata$nFrags), sd = 0.00001),
      y = Metadata$TSSEnrichment + rnorm(length(Metadata$nFrags), sd = 0.00001), 
      colorDensity = TRUE,
      xlim = c(2.5, 5),
      ylim = c(0, max(Metadata$TSSEnrichment) * 1.05),
      baseSize = 6,
      continuousSet = "sambaNight",
      xlabel = "Log 10 (Unique Fragments)",
      ylabel = "TSS Enrichment",
      title = ggtitle,
      rastr = TRUE, 
      size = 2) + 
      geom_hline(yintercept=cutoffs[sample,"atac_TSS"], lty = "dashed", size = 0.25) +
      geom_vline(xintercept=log10(cutoffs[sample,"atac_nFrags"]), lty = "dashed", size = 0.25)
    plotPDF(gg, name=paste0(sample,"-TSS_by_Unique_Frags.pdf"), width=4,height=4, ArchRProj = atac_proj,addDOC = FALSE)
  }
  
  # the following two plots are quite slow to run, can enable with the flag allplots
  if (allplots==TRUE){
    p1 <- plotFragmentSizes(ArchRProj = atac_proj)
    p2 <- plotTSSEnrichment(ArchRProj = atac_proj)
    plotPDF(p1, name="Frag_Size_Distribution.pdf", width=4, height=4, ArchRProj = atac_proj, addDOC=FALSE)
    plotPDF(p2, name="TSS_Enrichment.pdf", width=4, height=4, ArchRProj = atac_proj, addDOC=FALSE)
  }
  
}

addSeuratPCAtoArchRProj <- function(seurat_obj, reduction = "pca", atac_proj, convertCellName = TRUE) {
  # seurat_obj: Seurat object with the reduction to add to an ArchR Project
  # reduction: name of the reducedDims object in the seurat obj to add to the ArchR Project
  # atac_proj: an ArchR project to add the dimensionality reduction object
  # convertCellName: boolean to determine whether to convert Seurat cell barcode convention ("_") to ArchR cell barcode convention ("#"). This assumes 10x barcode structure.

  # Extract the reducedDims object from seurat object
  reducedDim <- Embeddings(seurat_obj, reduction = reduction)

  # convert cell names for ArchR
  if (convertCellName) {
    stringr::str_sub(rownames(reducedDim), -19, -19) <- "#"
  }

  # Add reducedDims matrix to ArchRProj
  atac_proj@reducedDims[[reduction]] <- SimpleList(
    matDR = reducedDim,
    date = Sys.time(),
    scaleDims = NA,
    corToDepth = NA
  )
  return(atac_proj)
}