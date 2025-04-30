# Helper functions for Seurat

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(magrittr)
  library(DoubletFinder)
  library(dplyr)
  library(tidyr)
  library(ggrastr)
})


scRNAdataPreProcessing <- function(
  obj, objPrefix, plotDir, # Seurat object, prefix, and plot dir for qc plots
  minFeatures=200, maxFeatures=Inf, minCounts=1000, maxCounts=Inf, maxPctMito=10, # Basic quality filters
  nfeatures=2500, dims=1:15, res=0.5, # Seurat clustering parameters
  runDoubletFinder=TRUE, # Should DoubletFinder be run?
  estDubRate=0.075, # DoubletFinder parameters
  runDecontX=TRUE, # Should DecontX be run?
  assays="RNA", # Which assays should be kept?
  ncores=1, use_logfile=TRUE
  ){

  # Perform some basic filtering on a seurat object.
  # Seurat object should be created from a single sample (e.g. one 10x run)
  # Optionally run other pre-processing tools:
  #
  # - DoubletFinder to estimate likely doublets 
  #   https://github.com/chris-mcginnis-ucsf/DoubletFinder
  #   https://www-cell-com.stanford.idm.oclc.org/cell-systems/fulltext/S2405-4712(19)30073-0
  #
  # - DecontX to reduce ambient RNA contamination
  #   https://github.com/campbio/celda
  #   https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6
  # 
  ################################
  # obj = seurat object
  # objPrefix = prefix for raq_qc plots
  # plotDir = directory for plotting
  # minFeatures = minimum number of features (genes) per cell
  # maxFeatures = maximum number of features (genes) per cell
  # minCounts = minimum number of UMIs per cell
  # maxCounts = maximum number of UMIs per cell
  # maxPctMito = maximum mitochondrial read percentage
  # nfeatures = number of variable features to be used in DoubletFinder
  # dims = which PCA dimensions will be used in DoubletFinder
  # estDubRate = estimated percent doublets (DF.classify will identify exactly this percent as doublets!)

  if(use_logfile){
    logfile <- paste0(plotDir, sprintf("/%s_preprocess_log_%s.txt", objPrefix, format(Sys.time(), "%Y%m%d-%H%M%S")))
    con <- file(logfile, open = "wt")
    sink(con, type="output")
    sink(con, type="message")
  }

  # Add percent.mt 
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  cellsBeforeFiltering <- dim(obj)[2]
  
  # Save some quick plots of raw qc
  histBreaks <- 100

  pdf(paste0(plotDir, sprintf("/%s_nCountHist.pdf", objPrefix)))
  df <- data.frame(cells=Cells(obj), log10nCount=log10(obj$nCount_RNA))
  p <- qcHistFilter(df, cmap = "blue", bins=histBreaks, border_color="black", lower_lim=log10(minCounts), upper_lim=log10(maxCounts))
  print(p)
  dev.off()

  pdf(paste0(plotDir, sprintf("/%s_nFeatureHist.pdf", objPrefix)))
  df <- data.frame(cells=Cells(obj), log10nFeatures=log10(obj$nFeature_RNA))
  p <- qcHistFilter(df, cmap = "blue", bins=histBreaks, border_color="black", lower_lim=log10(minFeatures), upper_lim=log10(maxFeatures))
  print(p)
  dev.off()

  pdf(paste0(plotDir, sprintf("/%s_pctMitoHist.pdf", objPrefix)))
  df <- data.frame(cells=Cells(obj), PctMito=obj$percent.mt)
  p <- qcHistFilter(df, cmap = "blue", bins=histBreaks, border_color="black", upper_lim=maxPctMito)
  print(p)
  dev.off()

  # Perform basic hard threshold filters
  message("Will filter based on:")
  message(sprintf("%s < unique genes < %s", minFeatures, maxFeatures))
  message(sprintf("%s < UMIs (counts) < %s", minCounts, maxCounts))
  message(sprintf("Percent mitochondrial < %s", maxPctMito))

  obj <- subset(obj, 
    subset = (
      nFeature_RNA > minFeatures & 
      nFeature_RNA < maxFeatures & 
      nCount_RNA > minCounts & 
      nCount_RNA < maxCounts & 
      percent.mt < maxPctMito
    )
  )
  cellsAfterFiltering <- dim(obj)[2]
  message(sprintf("%s filtered down to %s (%s%% remaining)", 
    cellsBeforeFiltering, cellsAfterFiltering, 
    round(100*(cellsAfterFiltering/cellsBeforeFiltering), 2)))

  # Perform standard Seurat pre-processing:
  obj <- seuratPreProcess(obj, vst.flavor="v2", dims=dims)

  # Close connections
  message("Finished preprocessing...")
  if(use_logfile){
    on.exit({ sink(type = "message"); sink(type = "output"); close(con) })
  }

  # Return processed Seurat object
  return(obj)
}


seuratPreProcess <- function(obj, vst.flavor="v2", dims=1:10, res=0.5, seed=1){
  # perform standard Seurat preprocessing 
  # (Normalization, Scaling, Dimensionality Reduction, Clustering)
  obj <- SCTransform(obj, vst.flavor = vst.flavor)
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims = dims)
  obj <- FindNeighbors(obj, dims = dims)
  obj <- FindClusters(obj, resolution=res, random.seed=1)
  return(obj)
}


runDoubletFinder <- function(obj, dims, estDubRate=0.075, ncores=1){
  # Run DoubletFinder on a provided (preprocessed) Seurat object
  # Return the seurat object with the selected pANN parameter and the 
  # DoubletFinder doublet classifications

  ### pK Identification (parameter-sweep) ###
  # "pK ~ This defines the PC neighborhood size used to compute pANN (proportion of artificial nearest neighbors), 
  # expressed as a proportion of the merged real-artificial data. 
  # No default is set, as pK should be adjusted for each scRNA-seq dataset"

  sweep.res.list <- paramSweep_v3(obj, PCs=dims, sct=FALSE, num.cores=ncores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  message(sprintf("Using pK = %s...", pK))

  # Get expected doublets (DF.classify will identify exactly this percent as doublets!)
  nExp_poi <- round(estDubRate * length(Cells(obj)))

  # DoubletFinder:
  obj <- doubletFinder_v3(obj, PCs = dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  # Rename results into more useful annotations
  pann <- grep(pattern="^pANN", x=names(obj@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(obj@meta.data), value=TRUE)
  obj$pANN <- obj[[pann]]
  obj$DF.classify <- obj[[classify]]
  obj[[pann]] <- NULL
  obj[[classify]] <- NULL

  return(obj)
}


runDecontX <- function(obj, seed=1){
  # Run DecontX on a provided Seurat object
  # From the DecontX vignette: 
  # "**Only the expression profile of *"real"* cells after cell calling are required to run DecontX. 
  # Empty cell droplet information (low expression cell barcodes before cell calling) are not needed.**"

  # DecontX can take either `SingleCellExperiment` object... or a single counts matrix as input. 
  # `decontX` will attempt to convert any input matrix to class `dgCMatrix` before beginning any analyses.
  counts <- GetAssayData(object = obj, slot = "counts")
  clusters <- Idents(obj) %>% as.numeric()

  # Run on only expressed genes
  x <- counts[rowSums(counts)>0,]
  message(sprintf("Running decontX on %s cells with %s non-zero genes...", dim(x)[2], dim(x)[1]))
  decon <- decontX(x, z=clusters, verbose=TRUE, seed=seed)

  # Save desired information back to Seurat Object
  # We will place the estimated 'decontaminated counts' in place of the original counts ('RNA')
  # and keep the original counts as a separate assay called 'origCounts'
  obj[["origCounts"]] <- CreateAssayObject(counts = counts)
  newCounts <- decon$decontXcounts
  # Add back unexpressed genes and sort according to original counts
  newCounts <- rbind(newCounts, counts[rowSums(counts)==0,])[rownames(counts),]
  obj[["RNA"]]@counts <- as(round(newCounts), "sparseMatrix")
  obj$estConp <- decon$contamination # Estimated 'contamination proportion, 0 to 1'

  return(obj)
}



plotClusterQC <- function(obj, subgroup, plotDir, pointSize=1.0, barwidth=0.9, sampleCmap=NULL, gest_age_cmap=NULL,
                          assay="RNA"){
  
  # Plot basic clustering plots
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

  ### Bar plot cluster counts ###
  tabDF <- base::table(obj$Clusters) %>% as.data.frame
  colnames(tabDF) <- c("Clusters", "count")

  pdf(paste0(plotDir, sprintf("/clusterBarPlot_%s.pdf", subgroup)))
  print(qcBarPlot(tabDF, cmap=qualcmap, barwidth=barwidth))
  dev.off()

  clustBySamp <- fractionXbyY(obj$Clusters, obj$Sample, add_total=TRUE, xname="Cluster", yname="Sample")

  pdf(paste0(plotDir, sprintf("/clustBySampleBarPlot_%s.pdf", subgroup)))
  print(stackedBarPlot(clustBySamp, cmap=sampleCmap, namedColors=namedSampCmap, barwidth=barwidth))
  dev.off()

  ### Stacked bar plot fraction gestational age in clusters ###
  clustByAge <- fractionXbyY(obj$Clusters, obj$age, add_total=TRUE, xname="Cluster", yname="Gestational Age")

  pdf(paste0(plotDir, sprintf("/clustByAgeBarPlot_%s.pdf", subgroup)))
  print(stackedBarPlot(clustByAge, cmap=gest_age_cmap, namedColors=namedAgeCmap, barwidth=barwidth))
  dev.off()

  ### Cluster UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$Clusters)
  # Randomize cells before plotting UMAP
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/cluster_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=qualcmap, point_size=pointSize))
  dev.off()

  ### Sample UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$Sample)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/sample_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=sampleCmap, namedColors=namedSampCmap, point_size=pointSize))
  dev.off()
  
  ### Cluster UMAP split by Sample ###
  pdf(paste0(plotDir, sprintf("/cluster_UMAP_splitsample_%s.pdf", subgroup)))
  for (samp in unique(obj$Sample)){
    obj.subset <- subset(obj, Sample==samp)
    umapDF <- data.frame(Embeddings(object=obj.subset, reduction="umap"), obj.subset$Clusters)
    # Randomize cells before plotting
    set.seed(1)
    umapDF <- umapDF[sample(nrow(umapDF)),]
    print(plotUMAP(umapDF, dataType="qualitative", cmap=qualcmap, point_size=pointSize, 
                   plotTitle=paste0(samp, " n = ", nrow(umapDF))))
  }
  dev.off()
  
  ### Age (PCW) UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$PCW)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]

  pdf(paste0(plotDir, sprintf("/PCW_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=gest_age_cmap, namedColors=namedAgeCmap, point_size=pointSize))
  dev.off()
  
  ### Age (PCD) UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$PCD)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]
  
  pdf(paste0(plotDir, sprintf("/PCD_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=gest_age_cmap, namedColors=namedAgeCmap, point_size=pointSize))
  dev.off()
  
  ### Sex UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$sex)
  # Randomize cells before plotting
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]
  
  pdf(paste0(plotDir, sprintf("/Sex_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=qualcmap, point_size=pointSize))
  dev.off()

  ### Percent Mitochondrial UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$percent.mt)
  pdf(paste0(plotDir, sprintf("/pctMito_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="quantitative", cmap=quantcmap, point_size=pointSize))
  dev.off()

  ### nCounts (UMIs) UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), log10(obj[[sprintf("nCount_%s", assay)]]))
  pdf(paste0(plotDir, sprintf("/log10nCount_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="quantitative", cmap=quantcmap, point_size=pointSize))
  dev.off()

  ### nFeatures (unique genes) UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj[[sprintf("nFeature_%s", assay)]])
  pdf(paste0(plotDir, sprintf("/nFeature_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="quantitative", cmap=quantcmap, point_size=pointSize))
  dev.off()

  ### Cell Cycle Phase UMAP ###
  umapDF <- data.frame(Embeddings(object=obj, reduction="umap"), obj$Phase)
  pdf(paste0(plotDir, sprintf("/ccPhase_UMAP_%s.pdf", subgroup)))
  print(plotUMAP(umapDF, dataType="qualitative", cmap=qualcmap, point_size=pointSize))
  dev.off()
  
    ### tree plot of clusters ###
  obj <- BuildClusterTree(obj, dims = 1:50)
  pdf(paste0(plotDir, sprintf("/cluster_tree_%s.pdf", subgroup)))
  PlotClusterTree(obj, direction = "downwards")
  dev.off()
  
  ### violin plot of UMI counts ###
  cmap <- getColorMap(cmaps_BOR$stallion, n=length(unique(obj$Clusters)))
  VlnPlot(obj, features=sprintf("nCount_%s", assay), group.by="Clusters", log = T, pt.size=0, 
          cols=cmap)
  ggsave(paste0(plotDir, sprintf("/vlnplot_umi_cluster_%s.pdf", subgroup)), height=5, width=10)
  
}

# Identifies consensus variable features from a list of sctransformed seurat objects
# Filters out blacklist (vector of gene names) genes
getConsensusVarFeatures <- function(objs, nfeatures = 3000, blacklist.genes = NULL){
  consFeatures <- SelectIntegrationFeatures(objs, nfeatures = nfeatures+1000) #Buffer to ensure enough features after filtering
  
  # Filter consensus features based on blacklist
  consFeatures <- consFeatures[which(!consFeatures %in% blacklist.genes)]

  return(consFeatures[1:nfeatures])
}

# SHAREseq functions--------------------
read_rna_data <- function(input.dir, use.sublib=FALSE, min.umi=100, min.cells=3, min.features=200){
  # read the matrix market files output from https://github.com/GreenleafLab/shareseq-pipeline 
  # input.dir: output folder from the preprocessing pipeline 
  # use.sublib: if true, read from RNA/sublibraries instead of RNA/samples
  # min.umi: filter out cell barcodes with less than min.umi reads
  # min.cells: filter out features detected in less than min.cells cells
  # min.features: filter out cells with less than min.features genes
  
  if (use.sublib){
    sample.list <- list.files(file.path(input.dir,"/RNA/sublibraries/"), pattern="matrix.mtx.gz", recursive=T) %>% 
      strsplit("/matrix.mtx.gz") %>% unlist
    
    prefix <- file.path(input.dir,"/RNA/sublibraries/")
    mtx.suffix <- "/matrix.mtx.gz"
    cells.suffix <- "/barcodes.tsv.gz"
    features.suffix <- "/features.tsv.gz"
  } else{
    sample.list <- list.files(file.path(input.dir,"/RNA/samples/"), pattern=".matrix.mtx.gz", recursive=F) %>% 
      strsplit(".matrix.mtx.gz") %>% unlist
    
    prefix <- file.path(input.dir,"/RNA/samples/")
    mtx.suffix <- ".matrix.mtx.gz"
    cells.suffix <- ".barcodes.tsv.gz"
    features.suffix <- ".features.tsv.gz"
  }
  
  obj.list <- c()
  for (i in 1:length(sample.list)) {
    n <- sample.list[i]
    message(n)
    data <- ReadMtx(mtx=paste0(prefix, n, mtx.suffix),
                    cells=paste0(prefix, n, cells.suffix),
                    features=paste0(prefix, n, features.suffix))
    obj <- CreateSeuratObject(counts=data, project=n, min.cells = min.cells, min.features = min.features)
    obj <- RenameCells(obj, new.names = paste0(n, "#", colnames(obj))) # rename to same as ArchR format
    
    # filter
    obj <- subset(obj, nCount_RNA>min.umi)
    
    # add percent mito
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    
    obj.list[[i]] <- obj
  }
  
  names(obj.list) <- sample.list
  return(obj.list)
}

read_rna_data_h5 <- function(input.dir, min.umi=100, min.cells=3, min.features=200){
  # read the h5 files output after cellbender filtering
  # input.dir: directory with the cellbender output files 
  # min.umi: filter out cell barcodes with less than min.umi reads
  # min.cells: filter out features detected in less than min.cells cells
  # min.features: filter out cells with less than min.features genes
  
  sample.list <- list.files(input.dir, pattern="_cellbender_filtered.h5", recursive=F) %>% 
    strsplit("_cellbender_filtered.h5") %>% unlist
  
  obj.list <- c()
  for (i in 1:length(sample.list)) {
    n <- sample.list[i]
    message(n)
    h5.file <- paste0(input.dir,"/", n, "_cellbender_filtered.h5")
    data <- Read10X_h5(h5.file)
    obj <- CreateSeuratObject(counts=data, project=n, min.cells = min.cells, min.features = min.features)
    obj <- RenameCells(obj, new.names = str_replace(colnames(obj), "_CL", "#CL")) # rename to same as ArchR format
    
    # filter
    obj <- subset(obj, nCount_RNA>min.umi)
    
    # add percent mito
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    
    obj.list[[i]] <- obj
  }
  
  names(obj.list) <- sample.list
  return(obj.list)
}

plot_rna_qc <- function(obj.list, plot.dir, umi.cutoffs=1e3, gene.cutoffs=500, pct.mt.cutoffs=20){
  # generate 6 plots per sample (cutoffs for display only, this function does not filter): 
  # - nUMI vs cell
  # - nGene vs cell
  # - nGene vs nUMI
  # - density plot of nUMI
  # - density plot of nGene
  # - density plot of pct.mt
  # if the cutoffs are passed as vectors (equal to the length of obj.list) use individual cutoffs
  umi.plots <- c()
  gene.plots <- c()
  gene.vs.umi.plots <- c()
  umi.density.plots <- c()
  gene.density.plots <- c()
  pctmt.density.plots <- c()
  
  for (i in seq_along(obj.list)){
    sample <- names(obj.list)[i]
    message(sample)
    obj <- obj.list[[i]]
    umis_per_cell <- obj$nCount_RNA
    genes_per_cell <- obj$nFeature_RNA
    
    # support different cutoffs per sample
    umi.cutoff <- umi.cutoffs
    gene.cutoff <- gene.cutoffs
    pct.mt.cutoff <- pct.mt.cutoffs
    if (length(umi.cutoffs)==length(obj.list)){
      umi.cutoff <- umi.cutoffs[i]} 
    if (length(gene.cutoffs)==length(obj.list)){
      gene.cutoff <- gene.cutoffs[i]} 
    if (length(pct.mt.cutoffs)==length(obj.list)){
      pct.mt.cutoff <- pct.mt.cutoffs[i]} 
    
    
    # knee plots
    p1 <- BPCells::plot_read_count_knee(umis_per_cell, cutoff = umi.cutoff) + ggtitle(sample) + ylab("nUMI") + geom_point(color="darkblue")
    p2 <- BPCells::plot_read_count_knee(genes_per_cell, cutoff = gene.cutoff) + ggtitle(sample) + ylab("nGene") + geom_point(color="darkblue")
    
    # gene vs umi
    downsample.idx <- plot_read_count_knee(umis_per_cell, return_data=TRUE)
    p3 <- ggplot(obj@meta.data[names(downsample.idx$data$reads),]) + aes(x=nCount_RNA,y=nFeature_RNA) + 
      geom_point(colour="darkblue") + theme_classic() + 
      xlab("nUMI") + ylab("nGene") + ggtitle(sample) + 
      geom_vline(xintercept=umi.cutoff, linetype="dashed") + geom_hline(yintercept=gene.cutoff, linetype="dashed") +
      xscale("log10") + yscale("log10") + annotation_logticks(sides="bl") 
    
    # gene density plot
    p4 <- ggplot(obj@meta.data, aes(x=nFeature_RNA)) + geom_density(stat="density") + 
      theme_classic() + xlab("nGene") + ylab("density") + xscale("log10") + annotation_logticks(sides="b") +
      geom_vline(xintercept=gene.cutoff, linetype="dashed") + ggtitle(sample) +
      annotate("text", x=gene.cutoff, y=Inf, 
               label=paste0(sum(genes_per_cell>gene.cutoff), " cells\nmedian: ", 
                            median(genes_per_cell[genes_per_cell>gene.cutoff])), hjust=-0.5, vjust=1)
    
    # umi density plot
    p5 <- ggplot(obj@meta.data, aes(x=nCount_RNA)) + geom_density(stat="density") + 
      theme_classic() + xlab("nUMI") + ylab("density") + xscale("log10") + annotation_logticks(sides="b") +
      geom_vline(xintercept=umi.cutoff, linetype="dashed") + ggtitle(sample) +
      annotate("text", x=umi.cutoff, y=Inf, 
               label=paste0(sum(umis_per_cell>umi.cutoff), " cells\nmedian: ", 
                            median(umis_per_cell[umis_per_cell>umi.cutoff])), hjust=-0.5, vjust=1)
    
    # percent mito density plot
    p6 <- ggplot(obj@meta.data, aes(x=percent.mt)) + geom_density(stat="density") + 
      theme_classic() + xlab("% mito") + ylab("density") +
      geom_vline(xintercept=pct.mt.cutoff, linetype="dashed") + ggtitle(sample) +
      annotate("text", x=pct.mt.cutoff, y=Inf, 
               label=paste0(sum(obj$percent.mt<pct.mt.cutoff), " cells"), hjust=1, vjust=1)
    
    umi.plots[[i]] <- p1
    gene.plots[[i]] <- p2
    gene.vs.umi.plots[[i]] <- p3
    gene.density.plots[[i]] <- p4
    umi.density.plots[[i]] <- p5
    pctmt.density.plots[[i]] <- p6
  }
  
  pdf(file.path(plot.dir, "rna_umi_knee.pdf"))
  for (p in umi.plots){
    print(p)
  }
  invisible(dev.off())
  
  pdf(file.path(plot.dir, "rna_gene_knee.pdf"))
  for (p in gene.plots){
    print(p)
  }
  invisible(dev.off())
  
  pdf(file.path(plot.dir, "rna_gene_vs_umi.pdf"))
  for (p in gene.vs.umi.plots){
    print(p)
  }
  invisible(dev.off())
  
  pdf(file.path(plot.dir, "rna_umi_density.pdf"))
  for (p in umi.density.plots){
    print(p)
  }
  invisible(dev.off())
  
  pdf(file.path(plot.dir, "rna_gene_density.pdf"))
  for (p in gene.density.plots){
    print(p)
  }
  invisible(dev.off())
  
  pdf(file.path(plot.dir, "rna_pctmt_density.pdf"))
  for (p in pctmt.density.plots){
    print(p)
  }
  invisible(dev.off())
}

plot_nfrag_numi <- function(df, plot.dir, ratio.cutoffs=0){
  # plot ATAC nFrag vs RNA nUMI for cells passing filter from both modalities
  # also plot the nFrag/nUMI ratio per cell as a density plot
  # cutoffs for display only, this function does not filter
  # df = a data frame input with the columns cb, rna_nUMIs, atac_nFrags
  # ratio.cutoffs = if the cutoffs are passed as vectors (equal to the length of obj.list) use individual cutoffs
  
  df$sample <- lapply(df$cb, function(n){strsplit(n, split="#")[[1]][1]})
  sample.list <- unique(df$sample) %>% unlist
  plot.list <- c()
  density.plot.list <- c()
  for (i in seq_along(sample.list)){
    sample <- sample.list[i]
    dfsub <- df[df$sample==sample,]
    umi.list <- dfsub$rna_nUMIs
    ratio.list <- dfsub$nFrags_nUMIs_ratio
    names(umi.list) <- rownames(dfsub)
    names(ratio.list) <- rownames(dfsub)
    
    # support different cutoffs per sample
    ratio.cutoff <- ratio.cutoffs
    if (length(ratio.cutoffs)==length(sample.list)){
      ratio.cutoff <- ratio.cutoffs[i]} 
    
    downsample.idx <- plot_read_count_knee(umi.list, return_data=TRUE)
    p1 <- ggplot(df[names(downsample.idx$data$reads),]) + aes(x=rna_nUMIs,y=atac_nFrags) + 
      geom_point(colour="darkblue") + theme_classic() + 
      xlab("nUMI") + ylab("nFrag") + ggtitle(sample) + 
      xscale("log10") + yscale("log10") + annotation_logticks(sides="bl") 
    
    p2 <- ggplot(df[names(downsample.idx$data$reads),]) + aes(x=nFrags_nUMIs_ratio)+ geom_density(stat="density") + 
      theme_classic() + xlab("nFrags/nUMI ratio") + ylab("density") +
      geom_vline(xintercept=ratio.cutoff, linetype="dashed") + ggtitle(sample) + 
      annotate("text", x=ratio.cutoff, y=Inf, 
               label=paste0(sum(ratio.list>ratio.cutoff), " cells\nmedian: ", 
                            median(ratio.list[ratio.list>ratio.cutoff])), hjust=-0.5, vjust=1)
    
    plot.list[[i]] <- p1
    density.plot.list[[i]] <- p2
  }
  
  pdf(file.path(plot.dir, "nfrag_vs_numi.pdf"))
  for (p in plot.list){
    print(p)
  }
  invisible(dev.off())
  
  pdf(file.path(plot.dir, "nfrag_vs_numi_ratio_density.pdf"))
  for (p in density.plot.list){
    print(p)
  }
  invisible(dev.off())
  
}


WNN <- function(atac.proj.dir, rna.proj.rds, atac.dimreduct.src="IterativeLSI", atac.dimreduct.dest="LSI",
                atac.umap.src = "UMAP", atac.umap.dest = "ATACUMAP", 
                rna.dimreduct.src="pca", rna.assay=NULL, atac.ndim=30, rna.ndim=30){
  ## calculate joint embeddings, from https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
  # atac.proj.dir: path to ATAC project directory (input to loadArchRProject)
  # rna.proj.rds: path to RNA project rds file (input to readRDS)
  # atac.dimreduct.src: name of the ATAC reduction in ArchR object to use, default to "IterativeLSI"
  # atac.dimreduct.dest: name of the ATAC reduction after adding to Seurat object, default to "LSI"
  # atac.umap.src: name of the ATAC UMAP embeddings in ArchR object to use (for visualization), default to "UMAP"
  # atac.umap.dest: name of the ATAC UMAP embeddings after adding to Seurat object (for visualization), default to "ATACUMAP"
  # rna.dimreduct.src: name of the RNA reduction in Seurat object to use, default to "pca"
  # rna.assay: which RNA assay to add the ATAC dim reduct to, if NULL default to DefaultAssay(rna.proj)
  # atac.ndim: number of ATAC dimensions to use, default to 30
  # rna.ndim: number of RNA dimensions to use, default to 30
  
  # read ArchR and Seurat objects
  message("reading ArchR project")
  atac <- loadArchRProject(atac.proj.dir)
  message("reading Seurat object")
  obj <- readRDS(rna.proj.rds)
  if (is.null(rna.assay)){
    rna.assay <- DefaultAssay(obj)
  }
  
  # get the LSI embeddings
  message("getting ATAC LSI embeddings")
  lsi <- atac@reducedDims[[atac.dimreduct.src]]$matSVD
  lsi <- lsi[colnames(obj),] # reorder based on seurat cell names
  obj[[atac.dimreduct.dest]] <- CreateDimReducObject(embeddings = lsi, key = paste0(atac.dimreduct.dest,"_"), assay = rna.assay)
  
  # also get the ATAC UMAP embeddings
  message("getting ATAC UMAP embeddings")
  aumap <- atac@embeddings[[atac.umap.src]]$df
  aumap <- aumap[colnames(obj),] %>% as.matrix
  obj[[atac.umap.dest]] <- CreateDimReducObject(embeddings = aumap, key = paste0(atac.umap.dest, "_"), global = T, assay = rna.assay)
  
  # Identify multimodal neighbors. These will be stored in the neighbors slot, 
  # and can be accessed using obj[['weighted.nn']]
  # The WNN graph can be accessed at obj[["wknn"]], 
  # and the SNN graph used for clustering at obj[["wsnn"]]
  # Cell-specific modality weights can be accessed at obj$RNA.weight and obj$ATAC.weight
  message("finding multimodal neighbors")
  obj <- FindMultiModalNeighbors(
    obj, reduction.list = list(rna.dimreduct.src, atac.dimreduct.dest), 
    dims.list = list(1:rna.ndim, 1:atac.ndim), modality.weight.name = c("RNA.weight", "ATAC.weight")
  )
  
  message("running UMAP on WNN output")
  obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  return(obj)
}

