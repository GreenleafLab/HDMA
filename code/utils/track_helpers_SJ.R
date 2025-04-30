# Helpers for plotting tracks with BPCells


# GRanges helpers --------------------------------------------------------------




#' Convert strings representing genomic regions to GRanges
#' 
#' @param regions character, vector of regions in the form chr1:10000-20000
#' @param sep character, two delimiters which separate chrom, start, end in the
#' strings provided to \code{regions}
#' 
#' @example 
#' str_to_gr("chr15:52785497-52791921")
str_to_gr <- function(regions, sep = c(":", "-")) {
  
  ranges.df <- data.frame(ranges = regions)
  separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  ) %>% 
    GRanges()
  
}



#' Convert GRanges to dataframes w/ the info to make
#' highlights on BPCells trackplots
#' 
#' @param regions_gr GRanges
#' @param color character, color to use for highlight. Default: red.
#' @param alpha numeric, value between [0, 1] specifying opacity for highlight
#' 
#' @example
#' gr_to_highlight("chr15:52785497-52791921")
gr_to_highlight <- function(regions_gr, color = "red", alpha = 0.2) {
  
  tibble(
    xmin = start(regions_gr),
    xmax = end(regions_gr),
    color = color,
    ymin = -Inf, ymax = Inf, alpha = alpha
  )
  
}


#' Convert strings representing genomic regions to dataframes w/ the info to make
#' highlights on BPCells trackplots
#' 
#' @param regions character, vector of regions in the form chr1:10000-20000
#' @param color character, color to use for highlight. Default: red.
#' @param alpha numeric, value between [0, 1] specifying opacity for highlight
#' 
#' @example
#' str_to_highlight("chr15:52785497-52791921")
str_to_highlight <- function(regions, sep = c(":", "-"), color = "red", alpha = 0.2) {
  
  regions_gr <- str_to_gr(regions, sep)
  
  gr_to_highlight(regions_gr, color = color, alpha = alpha)
  
}




#' Convert GRanges to dataframes w/ the info to make
#' peak/hit annotation on BPCells trackplots
#' 
#' @param regions_gr GRanges
#' @param color character, color to use for highlight. Default: red.
#' @param alpha numeric, value between [0, 1] specifying opacity for highlight
#' 
#' @example
#' str_to_highlight("chr15:52785497-52791921")
gr_to_anno <- function(regions_gr, color = "red") {
  
  tibble(
    xmin = start(regions_gr),
    xmax = end(regions_gr),
    color = color,
    ymin = -1, ymax = 1)
  
}


#' Make a prettier version of the coordinate string, using commas to separate 1kbs
#' so that it's more readable
#' 
#' @param region character, region in the form of chr1:10000-20000
#' 
#' @value "chr1:10,000-20,000"
str_to_pretty <- function(region) {
  
  region_gr <- str_to_gr(region)
  paste0(seqlevels(region_gr), ":", scales::comma(start(region_gr)), "-", scales::comma(end(region_gr)))
  
}





#' Given a gene name, get a GRanges object summarizing the region around it for plotting
#' 
#' @param transcripts Transcripts element of BPCells object
#' @param gene character, gene symbol
gene_to_gr <- function(transcripts, gene, flank) {
  
  gr <- BPCells::gene_region(genes = transcripts, gene_symbol = gene, extend_bp = flank) %>% 
    data.frame() %>% 
    GRanges()
  
  message("@ returning region around ", gene, " with width ", width(gr))
  
  return(gr)
  
}


#' Computed binned maximums for a numeric values associated with GRanges in provided bins,
#' similar to GenomicRanges::binnedAverage.
#' 
#' Adapted from 
#' https://divingintogeneticsandgenomics.com/post/compute-averages-sums-on-granges-or-equal-length-bins/
#' and https://stackoverflow.com/a/45683839
#' 
binned_max <- function(bins, numvar, mcolname) {
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewMaxs(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}



#' Wrapper function for \code{binned_max}, to calculate the binned maximums
#' given the region, the bigwig data, and the GRanges metadata column to use for
#' the aggregation.
#' 
#' @param region_gr GRanges
#' @param bw GRanges object obtained from \code{rtracklayer::import.bw(bigwig_path)}
#' @param mcolname character, metadata column name in \code{bw} to aggergate
#' @param tile_width numeric, size of bins (in bp)
#' 
#' @value
#' Returns a new GRanges, where each range is a bin and the metadata contains the
#' aggregated score.
calculate_binned_max <- function(region_gr, bw, mcolname = "score", tile_width) {
  
  bins <- GenomicRanges::tile(region_gr, width = tile_width)[[1]]
  signal_rle <- GenomicRanges::coverage(bw, weight = mcolname)
  seqlevels(bins, pruning.mode="coarse") <- names(signal_rle)
  binned_max <- binned_max(bins, numvar = signal_rle, mcolname = paste0("max_", mcolname))
  
  return(binned_max)
  
}





#' Given a GRanges, filter the metadata cols, optionally renaming them
gr_select <- function(gr, cols_keep, new_names = NULL) {
  
  mcols(gr) <- mcols(gr)[cols_keep]
  
  if (is.null(new_names)) new_names <- cols_keep
  
  colnames(mcols(gr)) <- new_names
  
  return(gr)
  
}

#' For tracks that are in relative coordinates (1, 2, 3,...) not genomic coordinates,
#' highlight the region with the direct start and end
highlight_relative_region <- function(start, end, alpha = 0.2, color = "red") {
  
  ggplot2::annotate("rect",
                    alpha = alpha,
                    xmin = start, xmax = end,
                    ymin = -Inf, ymax = Inf,
                    fill = color)
  
}






# Trackplots -------------------------------------------------------------------

#' BPCells::trackplot_coverage, with some modifications:
#' 1. Handle the case where there's only one group
#' 2. Optionally rasterize the area geom, containing the per-base signal via rasterize = TRUE
trackplot_coverage2 <- function(fragments, region, groups,
                                cell_read_counts,
                                group_order = NULL,
                                bins = 500,
                                clip_quantile = 0.999,
                                colors = discrete_palette("stallion"),
                                legend_label = "group",
                                rasterize = TRUE,
                                track_label = "Normalized Insertions (RPKM)",
                                zero_based_coords = !is(region, "GRanges"),
                                return_data = FALSE) {
  BPCells:::assert_is(fragments, "IterableFragments")
  BPCells:::assert_not_null(cellNames(fragments))
  region <- BPCells:::normalize_ranges(region)
  BPCells:::assert_true(length(region$chr) == 1)
  BPCells:::assert_is_wholenumber(bins)
  if (!is.null(clip_quantile)) {
    BPCells:::assert_is_numeric(clip_quantile)
    BPCells:::assert_len(clip_quantile, 1)
  }
  BPCells:::assert_is(colors, "character")
  BPCells:::assert_true(length(colors) >= length(unique(groups)))
  BPCells:::assert_is(legend_label, "character")
  
  groups <- as.factor(groups)
  BPCells:::assert_true(length(cellNames(fragments)) == length(groups))
  BPCells:::assert_true(length(cellNames(fragments)) == length(cell_read_counts))
  
  
  region$tile_width <- (region$end - region$start) %/% bins
  
  membership_matrix <- BPCells:::cluster_membership_matrix(groups, group_order)
  group_read_counts <- multiply_rows(membership_matrix, cell_read_counts) %>%
    colSums()
  group_norm_factors <- 1e9 / (group_read_counts * region$tile_width)
  
  if (is.null(names(colors))) {
    names(colors) <- colnames(membership_matrix)
  }
  colors <- colors[seq_len(ncol(membership_matrix))]
  
  bin_centers <- seq(region$start, region$end - 1, region$tile_width) + region$tile_width / 2
  bin_centers <- pmin(bin_centers, region$end - 1)
  
  mat <- (tile_matrix(fragments, region, explicit_tile_names = TRUE) %*% membership_matrix) %>%
    as("dgCMatrix") %>%
    as("matrix")
  # Discard any partial bins
  # SJ: use drop=FALSE to maintain matrix format, in the case of 1 column (1 group)
  mat <- mat[seq_along(bin_centers), , drop = FALSE]
  
  data <- tibble::tibble(
    pos = rep(bin_centers, ncol(mat)),
    group = factor(rep(colnames(mat), each = nrow(mat)), levels=levels(groups)),
    insertions = as.vector(mat),
    # Normalized insertions are "IPKM" - insertions per kilobase per million mapped reads
    normalized_insertions = as.vector(multiply_cols(mat, group_norm_factors))
  )
  
  if (return_data) {
    return(data)
  }
  
  ymax <- quantile(data$normalized_insertions, clip_quantile)
  # Set precision of y-axis range label to within 1% of the max value
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  print(ymax_accuracy)
  
  # Add a check to prevent extremely small accuracy values
  if (ymax_accuracy <= 0) {
    stop("Check values. Yamx < 0.")
  }
  
  range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))
  
  data$normalized_insertions <- pmin(data$normalized_insertions, ymax)
  
  plot <- ggplot2::ggplot(data) +
    ggplot2::geom_area(ggplot2::aes(pos, pmin(normalized_insertions, ymax), fill = group)) + 
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_comma(big.mark=" ")) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    ggplot2::annotate("text", x=region$start, y=ymax, label=range_label, vjust=1.5, hjust=-0.1, size=11*.8/ggplot2::.pt) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label, fill = legend_label) +
    ggplot2::guides(y="none", fill="none") +
    ggplot2::facet_wrap("group", ncol=1, strip.position="left") +
    BPCells:::trackplot_theme()
  
  # SJ: optionally rasterize the area geom
  if (rasterize) plot <- ggrastr::rasterize(plot, layer = "area", dpi = 400)
  
  BPCells:::wrap_trackplot(plot, ggplot2::unit(ncol(mat), "null"), takes_sideplot=TRUE)
  
}



#' Make a trackplot from an external bw or bedGraph file, that's compatible with BPCells plots
#'
#' This function calculates smoothed data (taking the maximum value within bins)
#' or basepair-resolution data and plots it along the genome.
#'
#' @param bw GRanges object obtained from \code{rtracklayer::import.bw(bigwig_path)} or
#' \code{rtracklayer::import.bedGraph(bg_path)}
#' @param gene character, gene symbol corresponding to the genomic region which should be plot 
#' @param region character, genomic coordinates in the form chr:start-end
#' @param track_label character, the y-axis label for the track
#' @param facet_label character, the facet label for the plot
#' @param transcripts transcripts from a BPCells object. Required if plotting based
#' on a gene.
#' @param flank numeric, specifying how much to extend on either side of the gene for plotting 
#' @param tile_width numeric, specifis bin size for aggregating data. Signal within
#' bins will be aggregated by taking the max in each bin
#' @param plot_as character, one of "area" or "bar", controlling whether signal should
#' be plot as an area or bar plot. Bar plot recommended for smaller regions.
#' @param clip_quantile numeric, quantile for clipping max values. Default: 0.999
#' @param color character, track color
#' @param return_data logical, whether to return input data before plotting.
#' 
#' @value
#' BPCells trackplot as returned by \code{BPCells:::wrap_trackplot}
trackplot_bw <- function(bw,
                         gene = NULL,
                         region = NULL,
                         track_label,
                         facet_label,
                         transcripts = NULL,
                         flank = NULL,
                         tile_width = 100,
                         plot_as = "area",
                         clip_quantile = 0.999,
                         ymin_zero = TRUE,
                         rasterize = TRUE,
                         return_data = FALSE,
                         color = "black") {
  
  message("@ preparing data...")
  
  # get a GRanges containing the region of interest, either from gene or region string
  if (!is.null(gene)) region_gr <- gene_to_gr(transcripts = transcripts, gene = gene, flank = flank)
  else if (!is.null(region)) {
    
    region_gr <- str_to_gr(region)
    
  } else if (is.null(gene) & is.null(region)) stop("Must provide either gene or region")
  
  message("@ plotting in region with width ", width(region_gr))
  
  # subset the bigwig
  bw_filt <- bw[bw %over% region_gr]
  
  # bw_filt <- BRGenomics::makeGRangesBRG(bw_filt)
  
  if (tile_width > 1) {
    
    # get binned values
    message("@ binning data with tile width ", tile_width)
    data_binned_max <- calculate_binned_max(region_gr, bw_filt, mcolname = "score", tile_width = tile_width)
    
    # calculate bin centers follwing BPCells::trackplot_coverage
    bin_centers <- seq(start(region_gr), end(region_gr) - 1, tile_width) + (tile_width/2)
    bin_centers <- pmin(bin_centers, end(region_gr) - 1)
    
    bw_data <- tibble::tibble(
      pos = bin_centers,
      signal = data_binned_max$max_score,
      facet_label = facet_label
    )
    
  } else {
    
    message("@ plotting without binning.")
    
    # get the coordinates of the basepair positions in the region
    positions <- seq(start(region_gr), end(region_gr)-1)
    
    # convert the data over the specified ranges to a dataframe
    bw_data <- tibble::tibble(pos = start(bw_filt), signal = bw_filt$score, facet_label = facet_label) %>%
      # complete missing positions for the regions with 0
      complete(pos = positions,
               fill = list("signal" = 0,
                           "facet_label" = facet_label)) %>%
      arrange(pos)
    
    # double check that lengths are the same
    nrow(bw_data) == length(positions)
    
  }
  
  if (return_data) return(bw_data)
  
  # get clipped ymax: note that we do this on the data in the bigwig *prior* to
  # completing the positions with 0s
  ymax <- quantile(bw_filt$score, clip_quantile)
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  
  if (!ymin_zero) {
    
    ymin <- quantile(bw_filt$score, 1-clip_quantile)
    ymin_accuracy <- 10^as.integer(log10(0.01 * ymin))
    
  } else {
    
    ymin = 0
    
  }
  
  range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))
  
  # clip values if needed
  bw_data$signal <- pmin(bw_data$signal, ymax)
  
  # plot track
  message("@ plotting...")
  plot <- ggplot2::ggplot(bw_data)
  
  if (plot_as == "area") {
    
    plot <- plot +
      ggplot2::geom_area(ggplot2::aes(x = pos, y = signal), fill = color)
    
  } else if (plot_as == "bar") {
    
    plot <- plot +
      ggplot2::geom_bar(ggplot2::aes(x = pos, y = signal), fill = color,
                        # don't leave spaces between bars
                        width = 1.1,
                        stat = "identity")
    
  }
  
  plot <- plot +
    ggplot2::scale_x_continuous(limits = c(start(region_gr), end(region_gr)), expand = c(0, 0), labels = scales::label_comma(big.mark=" ")) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
    ggplot2::annotate("text", x = start(region_gr), y = ymax, label = range_label, vjust=1.5, hjust=-0.1, size=11*.8/ggplot2::.pt) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y="none", fill="none") +
    # facetting is used to get a strip label on the left side, even if only one track plotted
    ggplot2::facet_wrap("facet_label", strip.position="left") +
    BPCells:::trackplot_theme()
  
  if (rasterize) plot <- ggrastr::rasterize(plot, layer = plot_as, dpi = 400)
  
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(1, "null"), takes_sideplot = FALSE) %>%
    BPCells:::set_trackplot_label(labels = facet_label)
  
  return(trackplot)
  
}


#' Make a trackplot from an external bw or bedGraph file, that's compatible with BPCells plots
#'
#' This function calculates smoothed data (taking the maximum value within bins)
#' or basepair-resolution data for two values and plots it along the genome, overlayed as
#' line plots.
#'
#' @param bw GRanges object obtained from \code{rtracklayer::import.bw(bigwig_path)} or
#' \code{rtracklayer::import.bedGraph(bg_path)}. Expecting score1 and score2 cols.
#' @param gene character, gene symbol corresponding to the genomic region which should be plot 
#' @param region character, genomic coordinates in the form chr:start-end
#' @param track_label character, the y-axis label for the track
#' @param facet_label character, the facet label for the plot
#' @param transcripts transcripts from a BPCells object. Required if plotting based
#' on a gene.
#' @param flank numeric, specifying how much to extend on either side of the gene for plotting 
#' @param tile_width numeric, specifis bin size for aggregating data. Signal within
#' bins will be aggregated by taking the max in each bin
#' @param plot_as character, one of "area" or "bar", controlling whether signal should
#' be plot as an area or bar plot. Bar plot recommended for smaller regions.
#' @param clip_quantile numeric, quantile for clipping max values. Default: 0.999
#' @param color character, track color
#' @param return_data logical, whether to return input data before plotting.
#' @param score_cmap character, length 2, names of the vector are the names of score1 and score2,
#' values are colors for score1 and score2. Default: ref = blue, alt = red.
#' 
#' @value
#' BPCells trackplot as returned by \code{BPCells:::wrap_trackplot}
trackplot_bw_paired <- function(bw,
                                gene = NULL,
                                region = NULL,
                                track_label,
                                facet_label,
                                transcripts = NULL,
                                flank = NULL,
                                alpha = 0.6,
                                ymin_zero = TRUE,
                                rasterize = FALSE,
                                score_cmap = c("ref" = "blue", "alt" = "red"),
                                return_data = FALSE) {
  
  message("@ preparing data...")
  
  # get a GRanges containing the region of interest, either from gene or region string
  if (!is.null(gene)) region_gr <- gene_to_gr(transcripts = transcripts, gene = gene, flank = flank)
  else if (!is.null(region)) {
    
    region_gr <- str_to_gr(region)
    
  } else if (is.null(gene) & is.null(region)) stop("Must provide either gene or region")
  
  message("@ plotting in region with width ", width(region_gr))
  
  # subset the bigwig
  bw_filt <- bw[bw %over% region_gr]
  
  message("@ plotting without binning.")
  
  # get the coordinates of the basepair positions in the region
  positions <- seq(start(region_gr), end(region_gr)-1)
  
  # convert the data over the specified ranges to a dataframe
  bw_data <- tibble::tibble(pos = start(bw_filt), signal1 = bw_filt$score1, signal2 = bw_filt$score2, facet_label = facet_label) %>%
    # complete missing positions for the regions with 0
    complete(pos = positions,
             fill = list("signal1" = 0,
                         "signal2" = 0,
                         "facet_label" = facet_label)) %>%
    arrange(pos)
  
  # transform to long format:
  bw_data <- bw_data %>% 
    pivot_longer(cols = c(signal1, signal2)) %>% 
    set_colnames(c("pos", "facet_label", "plot_group", "signal"))
  
  # rename groups
  bw_data[bw_data$plot_group == "signal1", ]$plot_group <- names(score_cmap)[1]
  bw_data[bw_data$plot_group == "signal2", ]$plot_group <- names(score_cmap)[2]
  
  # double check that lengths are the same
  nrow(bw_data) == length(positions)
  
  if (return_data) return(bw_data)
  
  # get ymax: note that we do this on the data in the bigwig *prior* to
  # completing the positions with 0s
  ymax <- max(bw_filt$score1, bw_filt$score2)
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  
  if (!ymin_zero) {
    
    ymin <- min(bw_filt$score1, bw_filt$score2)
    ymin_accuracy <- 10^as.integer(log10(0.01 * abs(ymin)))
    range_label <- glue("[{scales::label_comma(accuracy = ymin_accuracy, big.mark=' ')(ymin)}-{scales::label_comma(accuracy = ymax_accuracy, big.mark=' ')(ymax)}]")
    
  } else {
    
    ymin = 0
    range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))
    
  }
  
  # plot track
  message("@ plotting...")
  plot <- ggplot2::ggplot(bw_data, aes(group = plot_group))
  
  plot <- plot +
    ggplot2::geom_line(ggplot2::aes(x = pos, y = signal, color = plot_group), alpha = alpha) +
    scale_color_manual(values = score_cmap) +
    ggplot2::scale_x_continuous(limits = c(start(region_gr), end(region_gr)), expand = c(0, 0), labels = scales::label_comma(big.mark=" ")) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
    ggplot2::annotate("text", x = start(region_gr), y = ymax, label = range_label, vjust=1.5, hjust=-0.1, size=11*.8/ggplot2::.pt) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y="none", fill="none") +
    # facetting is used to get a strip label on the left side, even if only one track plotted
    ggplot2::facet_wrap("facet_label", strip.position="left") +
    BPCells:::trackplot_theme()
  
  if (!ymin_zero) plot <- plot + ggplot2::geom_hline(yintercept = 0)
  
  if (rasterize) plot <- ggrastr::rasterize(plot, layer = plot_as, dpi = 400)
  
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(1, "null"), takes_sideplot = FALSE) %>%
    BPCells:::set_trackplot_label(labels = facet_label)
  
  return(trackplot)
  
}




#' Make a trackplot for basepair-level neural network attributions, given
#' an external bw file, that's compatible with BPCells plots.
#' 
#' In these plots, at each position, the nucleotide in the reference genome
#' is plot at a height corresponding to its contribution score. This is implemented
#' using \code{ggseqlogo}.
#' 
#' TODO: we can't get genomic coordinates here because ggseqlogo treats x-axis
#' as positions starting at 1. So we need a way to map 1:N to the actual
#' genomic range for the labels.
#'
#' @param bw GRanges object obtained from \code{rtracklayer::import.bw(bigwig_path)},
#' corresponding to base-resolution contribution scores.
#' @param region character, genomic coordinates in the form chr:start-end
#' @param genome BSGenomes genome object e.g. BSgenome.Hsapiens.UCSC.hg38, used to
#' fetch sequence data.
#' @param track_label character, the y-axis label for the track
#' 
#' @value
#' BPCells trackplot as returned by \code{BPCells:::wrap_trackplot}
trackplot_contribs <- function(bw, region, genome,
                               track_label = "Contributions",
                               facet_label = NULL,
                               clip_quantile = 0.999,
                               rel_height = 0.6) {
  
  region_gr <- str_to_gr(region)
  
  if (width(region_gr) > 500) warning("@ not recommended to plot basepair-level contribs for 500 bp")
  else message("@ plotting basepair-level contribs for width ", width(region_gr))
  
  contrib_filt <- bw[bw %over% region_gr]
  
  ymax <- quantile(contrib_filt$score, clip_quantile)
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))
  
  # clip values if needed
  contrib_filt$score <- pmin(contrib_filt$score, ymax)
  
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
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y="none", fill="none") +
    BPCells:::trackplot_theme() +
    
    ggplot2::theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  
  # make this plot a bit shorter
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(rel_height, "null"), takes_sideplot = FALSE)
  
  if (!is.null(facet_label)) trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  
  return(trackplot)
  
  
}








#' Adapting trackplot_contribs2 to a second format of contribution scores as BED files as returned by 
#' variant scoring code.
#' 
#' @param gr GRanges object with at least four metadata columns: A, C, G, T which
#' hold the per nucleotide contribution scores per position
trackplot_contribs2 <- function(gr, region,
                                track_label = "Contributions",
                                facet_label = NULL,
                                clip_quantile = 0.999,
                                rel_height = 0.6) {
  
  region_gr <- str_to_gr(region)
  
  if (width(region_gr) > 500) warning("@ not recommended to plot basepair-level contribs for 500 bp")
  else message("@ plotting basepair-level contribs for width ", width(region_gr))
  
  contrib_filt <- gr[gr %over% region_gr]
  
  # rows are nucleotides and cols are positions
  contribs_ohe <- as.matrix(mcols(contrib_filt)[, c("A", "C", "G", "T")]) %>% 
    t()
  
  # make track using ggseqlogo
  plot <- ggseqlogo::ggseqlogo(contribs_ohe, method = "custom", seq_type = "dna") +
    ggplot2::geom_hline(yintercept = 0, color = "gray90") +
    # TODO: right now, the seqlogo has positions 1:N, not start:end, so we can't
    # rescale the x-axis and add the ticks easily. Need to map the tick labels to
    # genomic coordinates.
    # But, still must remove padding! Important for making sure things are aligned
    # to other plots, even if we can't scale the x axis
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y="none", fill="none") +
    BPCells:::trackplot_theme() +
    
    ggplot2::theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  
  # make this plot a bit shorter
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(rel_height, "null"), takes_sideplot = FALSE)
  
  if (!is.null(facet_label)) trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  
  return(trackplot)
  
  
}


#' Make a trackplot from an external bw or bedGraph file, that's compatible with BPCells plots
#'
#' This function calculates smoothed data (taking the maximum value within bins)
#' or basepair-resolution data for two values and plots it along the genome, overlayed as
#' line plots.
#'
#' @param bw GRanges object obtained from \code{rtracklayer::import.bw(bigwig_path)} or
#' \code{rtracklayer::import.bedGraph(bg_path)}. Expecting score1 and score2 cols.
#' @param gene character, gene symbol corresponding to the genomic region which should be plot 
#' @param region character, genomic coordinates in the form chr:start-end
#' @param track_label character, the y-axis label for the track
#' @param facet_label character, the facet label for the plot
#' @param transcripts transcripts from a BPCells object. Required if plotting based
#' on a gene.
#' @param flank numeric, specifying how much to extend on either side of the gene for plotting 
#' @param tile_width numeric, specifis bin size for aggregating data. Signal within
#' bins will be aggregated by taking the max in each bin
#' @param plot_as character, one of "area" or "bar", controlling whether signal should
#' be plot as an area or bar plot. Bar plot recommended for smaller regions.
#' @param clip_quantile numeric, quantile for clipping max values. Default: 0.999
#' @param color character, track color
#' @param return_data logical, whether to return input data before plotting.
#' @param score_cmap character, length 2, names of the vector are the names of score1 and score2,
#' values are colors for score1 and score2. Default: ref = blue, alt = red.
#' 
#' @value
#' BPCells trackplot as returned by \code{BPCells:::wrap_trackplot}
trackplot_bw_paired <- function(bw,
                                gene = NULL,
                                region = NULL,
                                track_label,
                                facet_label,
                                transcripts = NULL,
                                flank = NULL,
                                alpha = 0.6,
                                ymin_zero = TRUE,
                                rasterize = FALSE,
                                score_cmap = c("ref" = "blue", "alt" = "red"),
                                return_data = FALSE) {
  
  message("@ preparing data...")
  
  # get a GRanges containing the region of interest, either from gene or region string
  if (!is.null(gene)) region_gr <- gene_to_gr(transcripts = transcripts, gene = gene, flank = flank)
  else if (!is.null(region)) {
    
    region_gr <- str_to_gr(region)
    
  } else if (is.null(gene) & is.null(region)) stop("Must provide either gene or region")
  
  message("@ plotting in region with width ", width(region_gr))
  
  # subset the bigwig
  bw_filt <- bw[bw %over% region_gr]
  
  message("@ plotting without binning.")
  
  # get the coordinates of the basepair positions in the region
  positions <- seq(start(region_gr), end(region_gr)-1)
  
  # convert the data over the specified ranges to a dataframe
  bw_data <- tibble::tibble(pos = start(bw_filt), signal1 = bw_filt$score1, signal2 = bw_filt$score2, facet_label = facet_label) %>%
    # complete missing positions for the regions with 0
    complete(pos = positions,
             fill = list("signal1" = 0,
                         "signal2" = 0,
                         "facet_label" = facet_label)) %>%
    arrange(pos)
  
  # transform to long format:
  bw_data <- bw_data %>% 
    pivot_longer(cols = c(signal1, signal2)) %>% 
    set_colnames(c("pos", "facet_label", "plot_group", "signal"))
  
  # rename groups
  bw_data[bw_data$plot_group == "signal1", ]$plot_group <- names(score_cmap)[1]
  bw_data[bw_data$plot_group == "signal2", ]$plot_group <- names(score_cmap)[2]
  
  # double check that lengths are the same
  nrow(bw_data) == length(positions)
  
  if (return_data) return(bw_data)
  
  # get ymax: note that we do this on the data in the bigwig *prior* to
  # completing the positions with 0s
  ymax <- max(bw_filt$score1, bw_filt$score2)
  ymax_accuracy <- 10^as.integer(log10(0.01 * ymax))
  
  if (!ymin_zero) {
    
    ymin <- min(bw_filt$score1, bw_filt$score2)
    ymin_accuracy <- 10^as.integer(log10(0.01 * abs(ymin)))
    range_label <- glue("[{scales::label_comma(accuracy = ymin_accuracy, big.mark=' ')(ymin)}-{scales::label_comma(accuracy = ymax_accuracy, big.mark=' ')(ymax)}]")
    
  } else {
    
    ymin = 0
    range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))
    
  }
  
  # plot track
  message("@ plotting...")
  plot <- ggplot2::ggplot(bw_data, aes(group = plot_group))
  
  plot <- plot +
    ggplot2::geom_line(ggplot2::aes(x = pos, y = signal, color = plot_group), alpha = alpha) +
    scale_color_manual(values = score_cmap) +
    ggplot2::scale_x_continuous(limits = c(start(region_gr), end(region_gr)), expand = c(0, 0), labels = scales::label_comma(big.mark=" ")) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax), expand = c(0, 0)) +
    ggplot2::annotate("text", x = start(region_gr), y = ymax, label = range_label, vjust=1.5, hjust=-0.1, size=11*.8/ggplot2::.pt) +
    ggplot2::labs(x = "Genomic Position (bp)", y = track_label) +
    ggplot2::guides(y="none", fill="none") +
    # facetting is used to get a strip label on the left side, even if only one track plotted
    ggplot2::facet_wrap("facet_label", strip.position="left") +
    BPCells:::trackplot_theme()
  
  if (!ymin_zero) plot <- plot + ggplot2::geom_hline(yintercept = 0)
  
  if (rasterize) plot <- ggrastr::rasterize(plot, layer = plot_as, dpi = 400)
  
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(1, "null"), takes_sideplot = FALSE) %>%
    BPCells:::set_trackplot_label(labels = facet_label)
  
  return(trackplot)
  
}





#' Make a trackplot annotating motif hits, that's compatible with BPCells plots.
#' 
#' In these plots, each hit in the region is represented by a rectangle, and plot
#' alongside its name.
#'
#' @param hits GRanges object obtained from \code{rtracklayer::import.bed(hits_bed_path)},
#' corresponding to named motif hits.
#' @param region character, genomic coordinates in the form chr:start-end
#' @param track_label character, the y-axis label for the track
#' @param facet_label character, the facet label for the track
#' 
#' @value
#' BPCells trackplot as returned by \code{BPCells:::wrap_trackplot}
trackplot_hits <- function(hits, region, color, track_label = "Hits",
                           facet_label = NULL, rel_height = 0.6,
                           label_size = 4) {
  
  region_gr <- str_to_gr(region)
  
  hits_filt <- hits[hits %over% region_gr]
  
  # build the annotation table representing the coordinates of the rectangle for each hit
  hits_filt_anno <- tibble(
    xmin  = start(hits_filt),
    xmax  = end(hits_filt),
    alpha = hits_filt$score/1000,
    name  = hits_filt$name,
    ymin  = -0.5, ymax = 0.5) %>% 
    mutate(x = xmin + (xmax - xmin)/2)
  
  # plot hits
  plot <- ggplot2::ggplot(hits_filt_anno) +
    ggplot2::geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha), fill = color) +
    ggplot2::scale_y_continuous(limits = c(-3, 1)) +
    ggrepel::geom_text_repel(aes(label = name, x = x, y = ymin - 1), color = "black", fontface = "bold", size = label_size, hjust = 0,
                             force = 0.5,
                             min.segment.length = 1) +
    ggplot2::scale_alpha(range = c(0.5, 0.8), limits = c(0.7, 1)) +
    BPCells:::trackplot_theme() +
    ggplot2::scale_x_continuous(limits = c(start(region_gr), end(region_gr)),
                                expand = c(0, 0), labels = scales::label_comma(big.mark=" ")) +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL) +
    ggplot2::theme(axis.ticks.y = element_blank(),
                   axis.text.y = element_blank(),
                   legend.position = "none")
  
  # make this one a bit shorter
  trackplot <- BPCells:::wrap_trackplot(plot, ggplot2::unit(rel_height, "null"), takes_sideplot = FALSE)
  
  if (!is.null(facet_label)) trackplot <- BPCells:::set_trackplot_label(trackplot, labels = facet_label)
  
  return(trackplot)
  
}







# Annotation helpers -----------------------------------------------------------

#' Highlight a small genomic region on a trackplot
#' 
#' @param region character, specifying region to highlight
#'
#' @examples 
#' trackplot_coverage(
#'     region = region,
#'     fragments= bp_obj$frags,
#'     groups = bp_obj$cluster_metadata$Cluster) +
#'     highlight_region(region_small)
highlight_region <- function(region, alpha = 0.2, color = "red") {
  
  region_gr <- str_to_gr(region)
  
  ggplot2::annotate("rect",
                    alpha = alpha,
                    xmin = start(region_gr), xmax = end(region_gr),
                    ymin = -Inf, ymax = Inf,
                    fill = color)
  
}

highlight_relative_region <- function(start, end, alpha = 0.2, color = "red") {
  
  ggplot2::annotate("rect",
                    alpha = alpha,
                    xmin = start, xmax = end,
                    ymin = -Inf, ymax = Inf,
                    fill = color)
  
}





#' Combine track plots
#' 
#' Combines multiple track plots of the same region into a single grid.
#' Uses the `patchwork` package to perform the alignment.
#'
#' @param tracks List of tracks in order from top to bottom, generally ggplots as output from
#'    the other `trackplot_*()` functions.
#' @param side_plot Optional plot to align to the right (e.g. RNA expression per cluster). Will be aligned to the first
#'    `trackplot_coverage()` output if present, or else the first generic ggplot in the alignment. Should be in horizontal orientation and 
#'    in the same cluster ordering as the coverage plots.
#' @param title Text for overarching title of the plot
#' @param side_plot_width Fraction of width that should be used for the side plot relative to the main track area
#' @return A plot object with aligned genome plots. Each aligned row has
#'    the text label, y-axis, and plot body. The relative height of each row is given
#'    by heights. A shared title and x-axis are put at the top.
#' @seealso `trackplot_coverage()`, `trackplot_gene()`, `trackplot_loop()`, `trackplot_scalebar()`
#' @export
trackplot_combine2 <- function(tracks, side_plot = NULL, title = NULL, side_plot_width = 0.3) {
  for (plot in tracks) {
    BPCells:::assert_is(plot, "ggplot")
  }
  if (!is.null(side_plot)) {
    BPCells:::assert_is(side_plot, "ggplot")
  }
  
  # Calculate layout information on the plots
  heights <- list()
  collapse_upper_margin <- rep.int(TRUE, length(tracks))
  side_plot_row <- NULL
  areas <- NULL
  last_region <- NULL
  for (i in seq_along(tracks)) {
    if (is(tracks[[i]], "trackplot")) {
      if (tracks[[i]]$trackplot$takes_sideplot && is.null(side_plot_row)) {
        side_plot_row <- i
      }
      # If we switch regions, don't collapse margins into the track above
      if (!is.null(tracks[[i]]$trackplot$region)) {
        if (!is.null(last_region) && last_region != tracks[[i]]$trackplot$region && i > 1) collapse_upper_margin[i] <- FALSE 
        last_region <- tracks[[i]]$trackplot$region
      }
      
      # Preserve top and bottom margins if `keep_vertical_margin`
      if (tracks[[i]]$trackplot$keep_vertical_margin) {
        collapse_upper_margin[i] <- FALSE
        if (i < length(tracks)) collapse_upper_margin[i+1] <- FALSE
      }
    } else {
      if (is.null(side_plot_row)) side_plot_row <- i
    }
    heights <- c(heights, list(get_trackplot_height(tracks[[i]])))
    areas <- c(areas, list(patchwork::area(i, 1)))
  }
  heights <- do.call(grid::unit.c, heights)
  if (!is.null(side_plot) && is.null(side_plot_row)) {
    rlang::warn("Did not find a row to place the side_plot: no trackplot_coverage() or base ggplot tracks found. Defaulting to first row")
    side_plot_row <- 1L
  }
  
  # Collapse margins as needed among plots
  for (i in seq_along(tracks)) {
    plot.margin <- c(TRUE, TRUE, TRUE, TRUE) # Top, right, bottom, left
    if (!is.null(side_plot)) plot.margin[2] <- FALSE # Side plot should be flush on right side
    if (i < length(tracks) && collapse_upper_margin[i+1]) plot.margin[3] <- FALSE # Plot below should be flush
    if (collapse_upper_margin[i]) plot.margin[1] <- FALSE
    
    if (!plot.margin[3]) {
      tracks[[i]] <- tracks[[i]] +
        ggplot2::guides(x="none") +
        ggplot2::labs(x=NULL)
    }
    
    # Independent of showing the axis, we'll remove the bottom margin if the next row has the side_plot, since the
    # axis tick labels will add in some natural margin already
    # TODO: raise issue in BPCells
    if (!is.null(side_plot) && !is.null(side_plot_row) && (i+1 == side_plot_row)) plot.margin[3] <- FALSE  # ADDED CHECK HERE
    
    
    tracks[[i]] <- tracks[[i]] + ggplot2::theme(plot.margin=ggplot2::unit(5.5*plot.margin, "pt"))
  }
  
  # Reduce cut-off y-axis labels. Put plots with y axis labels later in the plot list, as they will take layer priority with patchwork
  has_y_axis <- vapply(tracks, function(t) is(t, "ggplot") && !is.null(t$labels$y), logical(1))
  tracks <- c(tracks[!has_y_axis], tracks[has_y_axis])
  areas <- c(areas[!has_y_axis], areas[has_y_axis])
  
  if (is.null(side_plot)) {
    widths <- c(1)
  } else {
    # Decide whether to put legends below/above side plot by adding up the height of all relatively-sized tracks
    height_above <- sum(as.vector(heights)[seq_along(heights) < side_plot_row & grid::unitType(heights) == "null"])
    height_below <- sum(as.vector(heights)[seq_along(heights) > side_plot_row & grid::unitType(heights) == "null"])
    if (height_above < height_below) {
      guide_position <- patchwork::area(side_plot_row+1L, 2, length(tracks))
    } else {
      guide_position <- patchwork::area(1L, 2, side_plot_row-1L)
    }
    
    widths <- c(1, side_plot_width)
    areas <- c(areas, list(patchwork::area(side_plot_row, 2), guide_position))
    # Make adjustments to the side plot style to fit in with tracks
    side_plot <- side_plot + 
      ggplot2::scale_x_discrete(limits=rev, position="top") +
      ggplot2::scale_y_continuous(position="right") +
      ggplot2::coord_flip() +
      ggplot2::labs(x=side_plot$labels$y, y=NULL) +
      ggplot2::theme(
        plot.margin=ggplot2::unit(c(0,0,0,0), "pt"),
        axis.ticks.length.x.top=grid::unit(-2.75, "pt")
      ) 
    guide_area <- patchwork::guide_area() + ggplot2::theme(plot.margin=ggplot2::unit(c(0,0,0,0), "pt"))
    tracks <- c(tracks, list(side_plot, guide_area))
  }
  
  patch <- patchwork::wrap_plots(tracks) +
    patchwork::plot_layout(ncol = 1, byrow = FALSE, heights = heights, widths=widths, guides = "collect", design=do.call(c, areas))
  
  if (!is.null(side_plot)) {
    # If a side plot is present, switch legend layout to use the horizontal space better
    patch <- patch * ggplot2::theme(
      legend.direction="horizontal", 
      legend.title.position = "top"
    )
  }
  if (!is.null(title)) {
    patch <- patch + patchwork::plot_annotation(title = title, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  return(patch)
}






