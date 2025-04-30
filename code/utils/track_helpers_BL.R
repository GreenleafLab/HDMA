#' Get Peak2Gene links from ArchR project as a data frame passing filter 
#'        (Correlation>0.45 & FDR<1e-4 & VarQATAC>0.25 & VarQRNA>0.25)
#' 
#' @param atac_proj ArchR project.
#' @param filter Boolean value on whether to filter out any Peak2Gene links 
#'               where the RNA gene is not present in the ArchR genome annotation 
#'               (mostly antisense and lncRNA), default to T.
#' @returns A data frame with columns seqnames,start,end,idxATAC,idxRNA,Correlation,FDR,VarQATAC,VarQRNA,peakstart,genestart.
get_p2g_loops <- function(atac_proj, filter=T){
  p2gtmp <- metadata(atac_proj@peakSet)$Peak2GeneLinks
  if (filter) {
    # filter out RNA genes not present in the ArchR annotation (mostly antisense and lncRNA)
    rna_genes <- metadata(p2gtmp)$geneSet
    addArchRGenome("hg38")
    atac_genes <- getGenes()
    genespf <- which(rna_genes$name %in% atac_genes$symbol)
    p2gpf <- p2gtmp[p2gtmp$idxRNA %in% genespf,]
  } else{
    p2gpf <- p2gtmp
  } 
  
  
  # construct loop data frame
  loops <- p2gpf 
  peaks <- metadata(loops)$peakSet[loops$idxATAC]
  genes <- metadata(loops)$geneSet[loops$idxRNA]
  
  # double check all chromosomes match in each peak-gene pair
  stopifnot(identical(seqnames(peaks) %>% as.character, seqnames(genes) %>% as.character))
  
  loops$seqnames <- seqnames(peaks) %>% as.character
  loops$peakstart <- ranges(peaks)@start
  loops$genestart <- ranges(genes)@start
  loops$start <- pmin(loops$peakstart, loops$genestart)
  loops$end <- pmax(loops$peakstart, loops$genestart)
  
  # threshold
  loops <- loops %>% as.data.frame %>% dplyr::filter(Correlation>0.45 & FDR<1e-4 & VarQATAC>0.25 & VarQRNA>0.25)
  
  return(loops)
}


#' Trackplot helper using BPCells 
#' 
#' @param gene   the gene name to plot, default to NULL, must supply either gene or region, 
#'               if both are supplied, region is used for plotting window and gene is only used for expression
#' @param region a string or GRanges object with the region to plot, default to NULL, must supply either gene or region,
#'               if both are supplied, region is used for plotting window and gene is only used for expression
#' @param clusters an annotation vector (length of vector = # cells) used to group coverage plots
#' @param fragments a BPCells IterableFragments object with the raw ATAC fragment count data
#' @param cell_read_counts a numeric vector (length of vector = # cells) with total fragment counts per cell used to normalize ATAC data
#' @param transcripts Transcipt features given as GRanges, data.frame, or list. See BPCells::gene_region()
#' @param loops (optional) a GRanges object
#' @param loop_color (optional) the column name in mcols(loops) used to color the loops, only supports numeric columns for now
#' @param loop_label (optional) text label for loops track
#' @param annot (optional) a GRanges object for genomic annotations
#' @param annot_color (optional) Name of a metadata column in loci to use for coloring, or a data vector with same length as loci. Column must be numeric or convertible to a factor.
#' @param annot_labelby (optional) Name of a metadata column in loci to use for labeling, or a data vector with same length as loci. Column must hold string data.
#' @param annot_label (optional) text label for annot track
#' @param annot_strand (optional) boolean value, whether to show the strand as arrows
#' @param expression_matrix (optional) a BPCells IterableMatrix object with the normalized RNA count matrix (the plotting function do not normalize)
#' @param flank only used if gene!=NULL and region==NULL, extends the gene region to plot by this number of bp on both sides
#' @param plot_title title text
#' @param highlights a GRanges object including regions to highlight on coverage plot
#' @param cmap color map used for coverage plot for each cluster
#' 
#' @returns A ggplot object.
trackplot_helper_v2c <- function(gene=NULL, region=NULL, clusters, fragments, cell_read_counts,
                                 transcripts, loops=NULL, loop_color=NULL, loop_label="loop", 
                                 annot=NULL, annot_color=NULL, annot_labelby=NULL, annot_label="annot", annot_strand=FALSE,
                                 expression_matrix=NULL, flank=1e5, plot_title="", highlights=NULL, cmap=NULL) {
  
  #----- window to plot
  if (!is.null(region)){
    if (is.character(region)) {
      region <- BPCells:::normalize_ranges(region)
      region$end <- region$end + 1L
    } else {
      region <- BPCells:::normalize_ranges(region)
    }
  } 
  if (!is.null(gene)){
    if (is.null(region)) {
      region <- gene_region(transcripts, gene, extend_bp = flank)
    }
    plot_title <- gene
  }
  
  if (is.null(gene) && is.null(region)) {
    stop("Must supply either a gene or a region to plot!")
  }
  
  #----- colormap
  if (is.null(cmap)){
    pal <- BPCells::discrete_palette("stallion", n=length(unique(clusters)))
  } else{
    pal <- cmap
  }
  
  #----- window size scale bar
  base_size <- 11
  scale_plot <- trackplot_scalebar(region, font_pt=base_size)
  
  #----- coverage track plot
  bins = min(500, (region$end - region$start - 1))
  bulk_plot <- trackplot_coverage(
    fragments,
    region=region,
    groups=clusters,
    cell_read_counts=cell_read_counts,
    colors=pal,
    bins=bins
  )
  
  #----- (optional) peak highlights
  if (!is.null(highlights)){
    annotations <- list()
    for (i in seq_len(nrow(highlights))) {
      annotations <- c(annotations, list(
        ggplot2::annotate("rect", alpha=highlights$alpha[i],
                          xmin=highlights$xmin[i], xmax=highlights$xmax[i],
                          ymin=highlights$ymin[i], ymax=highlights$ymax[i],
                          fill=highlights$color[i])
      ))
    }
    bulk_plot <- Reduce(`+`, c(list(bulk_plot, annotations)))
  }
  
  #----- gene/transcript annotation
  gene_plot <- trackplot_gene(transcripts, region) + ggplot2::guides(color="none")
  plot_list <- list(scale_plot, bulk_plot, gene_plot)
  
  #----- (optional) loops
  if (!is.null(loops)){
    loop_plot <- trackplot_loop(loops, region, color_by=loop_color, track_label=loop_label)
    #plot_list <- list(scale_plot, bulk_plot, gene_plot, loop_plot)
    plot_list <- c(plot_list, list(loop_plot))
  }
  
  #----- (optional) track annotations 
  # e.g. peaks, motif annotations
  if (!is.null(annot)){
    annot_plot <- trackplot_genome_annotation(annot, region, color_by = annot_color, colors = NULL, label_by = annot_labelby,
                                              show_strand = annot_strand, annotation_size = 2.5, track_label = annot_label, return_data = FALSE)
    
    plot_list <- c(plot_list, list(annot_plot))
  }
  
  #----- (optional) gene expression
  if (!is.null(gene) & !is.null(expression_matrix)){
    expression <- collect_features(expression_matrix, gene)
    names(expression) <- "gene"
    
    expression_plot <- ggplot2::ggplot(expression, ggplot2::aes(clusters, gene, fill=clusters)) +
      ggplot2::geom_boxplot() +
      ggplot2::guides(y="none", fill="none") +
      ggplot2::labs(x=NULL, y="RNA normalized expression") +
      ggplot2::scale_fill_manual(values=pal, drop=FALSE) +
      trackplot_theme()
    
    trackplot_combine(plot_list, side_plot=expression_plot, title=gene)
  } else{
    trackplot_combine(plot_list, title=plot_title) &
      ggplot2::theme(legend.direction="vertical")
  }
  
}



trackplot_theme <- function(base_size=11) {
  ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing.y = ggplot2::unit(0, "pt"),
      strip.text.y.left = ggplot2::element_text(angle=0, hjust=1, size=ggplot2::rel(1.2)),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      axis.title.y.left = ggplot2::element_text(size=ggplot2::rel(1)),
      legend.direction="horizontal", 
      legend.title.position = "top"
    )
}

trackplot_contribs_BL <- function(bw, region, genome,
         track_label = "Contributions",
         facet_label = NULL,
         clip_quantile = 0.999,
         rel_height = 0.6,
         ymax = NULL,
         ymin = NULL) {
  
  region_gr <- str_to_gr(region)
  
  if (width(region_gr) > 500) warning("@ not recommended to plot basepair-level contribs for 500 bp")
  else message("@ plotting basepair-level contribs for width ", width(region_gr))
  
  contrib_filt <- bw[bw %over% region_gr]
  
  if (is.null(ymax)){
    ymax <- quantile(contrib_filt$score, clip_quantile)
  } 
  if (is.null(ymin)){
    ymin <- quantile(contrib_filt$score, 1-clip_quantile)
  }
  
  ymax_accuracy <- 10^as.integer(log10(0.01 * abs(ymax)))
  ymin_accuracy <- 10^as.integer(log10(0.01 * abs(ymin)))
  range_label <- sprintf("[%s-%s]", 
                         scales::label_comma(accuracy = ymin_accuracy, big.mark=" ")(ymin),
                         scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax)) 
  
  # clip values if needed
  contrib_filt$score <- pmin(contrib_filt$score, ymax)
  contrib_filt$score <- pmax(contrib_filt$score, ymin)
  
  # get DNA sequence in the region
  region_seq <- genome[[seqlevels(region_gr)]][(start(region_gr)):(end(region_gr))]
  
  # use the sequence & contribution scores to make a one-hot-encoding, where values are the contributions
  # adapted from Surag Nair:
  # https://github.dev/kundajelab/scATAC-reprog/blob/bcbe7ab4474fa2ddd0139e31f784b84af17c20b4/src/figures_factory/Fig5/Supp_Vignette_CTCF_footprint_R.ipynb
  mat_ohe <- matrix(0, length(region_seq), 4)
  dim(mat_ohe)
  colnames(mat_ohe) = c("A", "C", "G", "T")
  mat_ohe[cbind(seq(length(region_seq)), as.vector(matrix(region_seq)))] = contrib_filt$score
  contribs_ohe <- t(mat_ohe)
  
  if (!(length(region_seq) == width(region_gr))) stop("@ sequence and region length are not equal.")
  
  # make track using ggseqlogo
  plot <- ggseqlogo::ggseqlogo(contribs_ohe, method = "custom", seq_type = "dna") +
    ggplot2::geom_hline(yintercept = 0, color = "gray90") +
    # TODO: right now, the seqlogo has positions 1:N, not start:end, so we can't
    # rescale the x-axis and add the ticks easily. Need to map the tick labels to
    # genomic coordinates.
    # But, still must remove padding! Important for making sure things are aligned
    # to other plots, even if we can't scale the x axis
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y="none", fill="none") +
    BPCells:::trackplot_theme() +
    ggplot2::theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank()) +
    ggplot2::annotate("text", x=1, y=ymax, label=range_label, vjust=1.5, hjust=-0.1, size=11*.8/ggplot2::.pt)
    
  
  # make this plot a bit shorter
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(rel_height, "null"), takes_sideplot = FALSE)
  
  if (!is.null(facet_label)) trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  
  return(trackplot)
  
  
}





