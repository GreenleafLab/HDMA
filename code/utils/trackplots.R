# Design: The aim is to make genome-browser style plots,
# where stacked plots are aligned to share a horizontal axis representing
# genome location.

# Each component plot function returns either a plain ggplot or a list of ggplot objects.
# These are then combined using draw_trackplot_grid, which handles alignment of
# all the plots with labels

#' Combine ggplot track plots into an aligned grid.
#' Uses patchwork to perform the alignment
#'
#' `r lifecycle::badge("deprecated")`
#'
#' This function has been renamed to `trackplot_combine()`
#'
#' @param ... Plots in order from top to bottom, generally plain ggplots.
#'    To better accomodate many bulk tracks, patchwork objects which contain multiple
#'    tracks are also accepted. In this case, plot labels will be drawn from the
#'    attribute `$patchwork$labels` if present, rather than the `labels` argument.
#' @param labels Text labels to display for each track
#' @param title Text for overarching title of the plot
#' @param heights Relative heights for each component plot. It is suggested to use 1 as standard height of a
#'    pseudobulk track.
#' @param label_width Fraction of width that should be used for labels relative to the main track area
#' @param label_style Arguments to pass to geom_text to adjust label text style
#' @return A plot object with aligned genome plots. Each aligned row has
#'    the text label, y-axis, and plot body. The relative height of each row is given
#'    by heights. A shared title and x-axis are put at the top.
#' @export
#' @keywords internal
draw_trackplot_grid <- function(..., labels, title = NULL,
                                heights = rep(1, length(plots)),
                                label_width = 0.2,
                                label_style = list(fontface = "bold", size = 4)) {
  lifecycle::deprecate_warn("0.2.0", "draw_trackplot_grid()", "trackplot_combine()")
  trackplot_combine(..., labels=labels, title=title, heights=heights, label_width=label_width, label_style=label_style)
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

# Internal helper function to return empty track plots if there's no data to be plotted
trackplot_empty <- function(region, label) {
  ggplot2::ggplot(tibble::tibble(label=label)) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_number()) +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL) +
    ggplot2::facet_wrap("label", strip.position="left") +
    trackplot_theme()
}

wrap_trackplot <- function(plot, height=NULL, takes_sideplot=FALSE) {
  if (!is.null(height)) {
    BPCells:::assert_is(height, "unit")
  }
  class(plot) <- c("trackplot", class(plot))
  plot$trackplot <- list(height=height, takes_sideplot=takes_sideplot)
  plot
}

get_patchwork_plots <- function(patchwork) {
  BPCells:::assert_is(patchwork, "patchwork")
  ret <- plot$patches
  plot$patches <- NULL
  class(plot) <- setdiff(class(plot), "patchwork")
  c(ret, plot)
}

#' Combine track plots
#' 
#' Combines multiple track plots of the same region into a single grid.
#' Uses the `patchwork` package to perform the alignment.
#'
#' @param tracks List of tracks in order from top to bottom, generally ggplots as output from
#'    the other `trackplot_*()` functions.
#' @param side_plot Optional plot to align to the right (e.g. RNA expression per cluster). Will be aligned to a
#'    `trackplot_coverage()` output if present, or else the first generic ggplot in the alignment. Should be in horizontal orientation and 
#'    in the same cluster ordering as the coverage plots.
#' @param title Text for overarching title of the plot
#' @param heights Relative heights for each component plot. It is suggested to use 1 as standard height of a
#'    pseudobulk track.
#' @param side_plot_width Fraction of width that should be used for the side plot relative to the main track area
#' @param label_width Fraction of width that should be used for labels relative to the main track area
#' @param label_style Arguments to pass to geom_text to adjust label text style
#' @return A plot object with aligned genome plots. Each aligned row has
#'    the text label, y-axis, and plot body. The relative height of each row is given
#'    by heights. A shared title and x-axis are put at the top.
#' @export
trackplot_combine <- function(tracks, side_plot = NULL, title = NULL, side_plot_width = 0.3) {
  for (plot in tracks) {
    BPCells:::assert_is(plot, "ggplot")
  }
  if (!is.null(side_plot)) {
    BPCells:::assert_is(side_plot, "ggplot")
  }
  
  heights <- list()
  side_plot_row <- 1L
  areas <- NULL
  for (i in seq_along(tracks)) {
    if (is(tracks[[i]], "trackplot")) {
      heights <- c(heights, list(tracks[[i]]$trackplot$height))
      if (tracks[[i]]$trackplot$takes_sideplot) {
        side_plot_row <- i
      } else if (side_plot_row == i) {
        side_plot_row <- side_plot_row + 1L
      }
    } else {
      heights <- c(heights, list(grid::unit(1L, "null")))
    }
    if (i != length(tracks)) {
      tracks[[i]] <- tracks[[i]] + 
        ggplot2::guides(x="none") +
        ggplot2::labs(x=NULL) +
        ggplot2::theme(plot.margin=ggplot2::unit(c(0,0,0,0), "pt"))
    } else {
      bottom_margin <- if(is.null(side_plot)) c(0,5.5,5.5,5.5) else c(0,0,5.5,5.5)
      tracks[[i]] <- tracks[[i]] + ggplot2::theme(plot.margin=ggplot2::unit(bottom_margin, "pt"))
    }
    areas <- c(areas, list(patchwork::area(i, 1)))
  }
  heights <- do.call(grid::unit.c, heights)
  
  if (is.null(side_plot)) {
    widths <- c(1)
  } else {
    widths <- c(1, side_plot_width)
    areas <- c(areas, list(patchwork::area(side_plot_row, 2), patchwork::area(side_plot_row+1L, 2, length(tracks))))
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
  if (!is.null(title)) {
    patch <- patch + patchwork::plot_annotation(title = title, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  return(patch)
}

#' Pseudobulk trackplot
#'
#' `r lifecycle::badge("deprecated")`
#'
#' This function has been renamed to `trackplot_coverage()`
#' 
#' Plot a pseudobulk genome track, showing the number of fragment insertions
#' across a region.
#'
#' @inheritParams cluster_membership_matrix
#' @inheritParams plot_embedding
#' @inheritParams convert_to_fragments
#' @param region GRanges of length 1 with region to plot, or list/data.frame with
#'    one entry each for chr, start, end. See `gene_region()` or [genomic-ranges] for details
#' @param fragments Fragments object
#' @param cell_read_counts Numeric vector of read counts for each cell (used for normalization)
#' @param bins Number of bins to plot across the region
#' @param clip_quantile (optional) Quantile of values for clipping y-axis limits. Default of 0.999 will crop out
#'    just the most extreme outliers across the region. NULL to disable clipping
#' @param colors Character vector of color values (optionally named by group)
#' @param legend_label Custom label to put on the legend
#'
#' @return Returns a combined plot of pseudobulk genome tracks. For compatability with
#' `draw_trackplot_grid()`, the extra attribute `$patches$labels` will be added to
#' specify the labels for each track. If `return_data` or `return_plot_list` is
#' `TRUE`, the return value will be modified accordingly.
#' @export
#' @keyword internal
trackplot_bulk <- function(fragments, region, groups,
                           cell_read_counts,
                           group_order = NULL,
                           bins = 200, clip_quantile = 0.999,
                           colors = discrete_palette("stallion"),
                           legend_label = "group",
                           zero_based_coords = !is(region, "GRanges"),
                           return_data = FALSE, return_plot_list = FALSE, apply_styling = TRUE) {
  lifecycle::deprecate_warn("0.2.0", "draw_trackplot_grid()", "trackplot_combine()")
  if (!isTRUE(apply_styling)) {
    lifecycle::deprecate_warn("0.2.0", "draw_trackplot_grid(apply_styling)", details="Argument is removed to simplify function. Styling can set manually on result")
  }
  if (!isFALSE(return_plot_list)) {
    lifecycle::deprecate_warn("0.2.0", "draw_trackplot_grid(return_plot_list)", details="Faceting is now used so returning a list of plots is unnecessary")
  }
  trackplot_coverage(fragments=fragments, region=region, groups=groups,
                     cell_read_counts=cell_read_counts,
                     group_order = group_order,
                     bins = bins, clip_quantile = clip_quantile,
                     colors = colors,
                     legend_label = legend_label,
                     zero_based_coords = zero_based_coords,
                     return_data = return_data, return_plot_list = return_plot_list)                        
}

#' Pseudobulk coverage trackplot
#'
#' Plot a pseudobulk genome track, showing the number of fragment insertions
#' across a region for each cell type or group. 
#'
#' @inheritParams cluster_membership_matrix
#' @inheritParams plot_embedding
#' @inheritParams convert_to_fragments
#' @param region GRanges of length 1 with region to plot, or list/data.frame with
#'    one entry each for chr, start, end. See `gene_region()` or [genomic-ranges] for details
#' @param fragments Fragments object
#' @param cell_read_counts Numeric vector of read counts for each cell (used for normalization)
#' @param scale_bar Whether to include a scale bar in the top track (`TRUE` or `FALSE`)
#' @param bins Number of bins to plot across the region
#' @param clip_quantile (optional) Quantile of values for clipping y-axis limits. Default of 0.999 will crop out
#'    just the most extreme outliers across the region. NULL to disable clipping
#' @param colors Character vector of color values (optionally named by group)
#' @param legend_label Custom label to put on the legend
#'
#' @return Returns a combined plot of pseudobulk genome tracks. For compatability with
#' `draw_trackplot_grid()`, the extra attribute `$patches$labels` will be added to
#' specify the labels for each track. If `return_data` or `return_plot_list` is
#' `TRUE`, the return value will be modified accordingly.
#' @export
trackplot_coverage <- function(fragments, region, groups,
                               cell_read_counts,
                               group_order = NULL,
                               bins = 500, clip_quantile = 0.999,
                               colors = discrete_palette("stallion"),
                               legend_label = "group",
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
  mat <- mat[seq_along(bin_centers), ]
  
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
  range_label <- sprintf("[0-%s]", scales::label_comma(accuracy = ymax_accuracy, big.mark=" ")(ymax))
  
  data$normalized_insertions <- pmin(data$normalized_insertions, ymax)
  
  plot <- ggplot2::ggplot(data) +
    ggplot2::geom_area(ggplot2::aes(pos, pmin(normalized_insertions, ymax), fill = group)) +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_comma(big.mark=" ")) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    ggplot2::annotate("text", x=region$start, y=ymax, label=range_label, vjust=1.5, hjust=-0.1, size=11*.8/ggplot2::.pt) +
    ggplot2::labs(x = "Genomic Position (bp)", y = "Normalized Insertions (RPKM)", fill = legend_label) +
    ggplot2::guides(y="none", fill="none") +
    ggplot2::facet_wrap("group", ncol=1, strip.position="left") +
    trackplot_theme() 
  
  wrap_trackplot(plot, ggplot2::unit(ncol(mat), "null"), takes_sideplot=TRUE)
}


#' Plot transcript models
#' @param transcripts GRanges, list, or data.frame of transcript features to plot.
#' Required attributes are:
#' - `chr`, `start`, `end`: genomic position
#' - `strand`: "+"/"-" or TRUE/FALSE for positive or negative strand
#' - `feature` (only entries marked as `"transcript"` or `"exon"` will be considered)
#' - `transcript_id`
#' - `gene_name`
#' See [genomic-ranges] for more details
#' @inheritParams trackplot_coverage
#' @param labels Character vector with labels for each item in transcripts. NA for items that should not be labeled
#' @param exon_size size for exon lines
#' @param transcript_size size for transcript lines
#' @param label_size size for transcript labels
#' @return Plot of gene locations
#' @export
trackplot_gene <- function(transcripts, region, exon_size = 2.5, gene_size = 0.5, label_size = 11*.8/ggplot2::.pt, track_label="Genes", return_data = FALSE) {
  region <- BPCells:::normalize_ranges(region)
  transcripts <- BPCells:::normalize_ranges(transcripts, metadata_cols = c("strand", "feature", "transcript_id", "gene_name"))
  
  # Adjust for the fact that exon_size and gene_size are now in units of linewidth = 0.75mm 
  # whereas size is given in units of 1mm (https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#linewidth)
  size_range <- c(min(exon_size, gene_size, label_size), max(exon_size, gene_size, label_size))
  linewidth_range <- size_range/.75
  exon_size <- exon_size/.75
  gene_size <- gene_size/.75
  
  data <- transcripts %>%
    tibble::as_tibble() %>%
    dplyr::filter(as.character(chr) == as.character(region$chr), end > region$start, start < region$end, feature %in% c("transcript", "exon")) %>%
    dplyr::mutate(size = dplyr::if_else(feature == "exon", exon_size, gene_size))
  
  if (nrow(data) == 0) {
    if (return_data) {
      data$y <- numeric(0)
      data$facet_label <- character(0)
      arrows <- tibble::tibble(start=numeric(0), end=numeric(0), y=integer(0), strand=logical(0), size=numeric(0))
      return(list(data=data, arrows=arrows))
    } else {
      return(trackplot_empty(region, label=track_label))
    }
  }
  
  ############################
  # Calculate y-position layout
  #############################
  # Steps:
  # 1. Calculate the maximum overlap depth of transcripts
  # 2. Iterate through transcript start/end in sorted order
  # 3. Randomly assign each transcript a y-coordinate between 1 and max overlap depth,
  #    with the restriction that a transcript can't have the same y-coordinate
  #    as a transcript it overlaps.
  tx_boundaries <- data %>% # Data frame of pos, tx_id, is_start
    dplyr::filter(feature == "transcript") %>%
    {
      dplyr::bind_rows(
        dplyr::transmute(., pos = start, transcript_id, start = TRUE),
        dplyr::transmute(., pos = end, transcript_id, start = FALSE)
      )
    } %>%
    dplyr::arrange(pos)
  total_positions <- max(cumsum(2 * tx_boundaries$start - 1))
  y_pos <- integer(0) # Names: tx_id, value: y position for transcript
  occupied_y <- rep(FALSE, total_positions) # Boolean vector, marking which y positions are open
  prev_seed <- BPCells:::get_seed()
  set.seed(12057235)
  for (i in seq_len(nrow(tx_boundaries))) {
    row <- tx_boundaries[i, ]
    if (row$start) {
      assigned_pos <- sample(which(!occupied_y), 1)
      y_pos[row$transcript_id] <- assigned_pos
      occupied_y[assigned_pos] <- TRUE
    } else {
      assigned_pos <- y_pos[row$transcript_id]
      occupied_y[assigned_pos] <- FALSE
    }
  }
  BPCells:::restore_seed(prev_seed)
  data <- dplyr::mutate(data, y=y_pos[transcript_id])
  
  #############################
  # Set up arrow coordinates
  #############################
  arrow_spacing <- (region$end - region$start) / 50
  arrow_list <- NULL
  transcript_coords <- dplyr::filter(data, feature == "transcript")
  for (i in seq_len(nrow(transcript_coords))) {
    if (transcript_coords$strand[i]) {
      endpoints <- seq(transcript_coords$start[i], transcript_coords$end[i], arrow_spacing)
    } else {
      endpoints <- rev(seq(transcript_coords$end[i], transcript_coords$start[i], -arrow_spacing))
    }
    arrow_list <- c(arrow_list, list(tibble::tibble(
      start = endpoints[-length(endpoints)],
      end = endpoints[-1],
      y = transcript_coords$y[i],
      strand = transcript_coords$strand[i],
      size = gene_size
    )))
  }
  arrows <- dplyr::bind_rows(arrow_list) %>%
    dplyr::filter(start >= region$start, start < region$end, end >= region$start, end < region$end)
  
  # Adjust the endpoints of any partially overlapping elements to fit within
  # the plot boundaries
  data <- dplyr::mutate(
    data,
    start = pmax(region$start, pmin(region$end, start)),
    end = pmax(region$start, pmin(region$end, end)),
    facet_label = track_label
  )
  
  if (return_data) {
    return(list(data=data, arrows=arrows))
  }
  
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = dplyr::if_else(strand, start, end), xend = dplyr::if_else(strand, end, start),
      y = y, yend = y,
      linewidth = size
    )
  ) +
    ggplot2::geom_segment(ggplot2::aes(color = dplyr::if_else(strand, "+", "-"))) +
    ggplot2::geom_segment(ggplot2::aes(color = dplyr::if_else(strand, "+", "-")), data=arrows, arrow=grid::arrow(length=grid::unit(.4*exon_size, "mm"))) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(data, feature == "transcript"),
      ggplot2::aes(label = gene_name),
      size = label_size,
      position = ggrepel::position_nudge_repel(y = 0.25)
    ) +
    ggplot2::scale_size(range = size_range, limits = size_range) +
    ggplot2::scale_linewidth(range = linewidth_range, limits = linewidth_range) +
    ggplot2::scale_color_manual(values = c("+" = "black", "-" = "darkgrey")) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_comma(big.mark=" ")) +
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL) +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL, color = "strand") +
    ggplot2::guides(size = "none", linewidth="none") +
    ggplot2::facet_wrap("facet_label", strip.position="left") +
    trackplot_theme()
  
  wrap_trackplot(plot, height=ggplot2::unit(1, "null"))
}

#' Plot loops
#'
#' @param loops `r document_granges()`
#' @param color_by Name of a metadata column in `loops` to use for coloring, or a numeric vector with same length as loops
#' @param colors Vector of hex color codes to use for the color gradient
#' @param allow_truncated If FALSE, remove any loops that are not fully contained within `region`
#' @param curvature Curvature value between 0 and 1. 1 is a 180-degree arc, and 0 is flat lines.
#' @inheritParams trackplot_coverage
#' 
#' @return Plot of arcs connecting genomic regions
#' @export
trackplot_loop <- function(loops, region, color_by=NULL, colors=c("#bfd3e6","#8c96c6","#88419d","#4d004b"), curvature=0.75, track_label="Links", return_data = FALSE) {
  region <- BPCells:::normalize_ranges(region)
  BPCells:::assert_true(is.null(color_by) || is.numeric(color_by) || is.character(color_by))
  BPCells:::assert_is_numeric(curvature)
  if (is.character(color_by)) {
    loops <- BPCells:::normalize_ranges(loops, metadata_cols=color_by) %>%
      dplyr::rename(color=dplyr::all_of(color_by))
  } else if (is.numeric(color_by)) {
    loops <- BPCells:::normalize_ranges(loops)
    loops[["color"]] <- color_by
    color_by <- argument_name(color_by, 2)
  }
  
  # Calculate curve points (segment of a circle) without scaling
  # Curve is scaled to start at x=0 and end at x=1
  resolution <- 100
  step <- (seq_len(resolution+1)-1)/resolution
  total_angle <- pi*curvature
  min_angle <- 1.5*pi - 0.5*total_angle
  max_angle <- min_angle + total_angle
  curve_x <- (cos(min_angle + step*total_angle) - cos(min_angle))/(cos(max_angle)-cos(min_angle))
  curve_y <- (sin(min_angle + step*total_angle) - sin(min_angle))/(cos(max_angle)-cos(min_angle))
  
  data <- loops %>%
    tibble::as_tibble() %>%
    dplyr::filter(as.character(chr) == as.character(region$chr), end > region$start | start < region$end) %>%
    dplyr::mutate(loop_id = dplyr::row_number()) %>%
    dplyr::cross_join(tibble::tibble(x=curve_x, y=curve_y)) %>%
    dplyr::mutate(
      x = x * (end - start) + start,
      y = y * (end - start),
      facet_label = track_label
    )

  if (return_data) {
    return(data)
  }
  
  if (nrow(data) == 0) {
    return(trackplot_empty(region, track_label))
  }
  
  ymin <- data %>%
    dplyr::filter(end <= region$end, start >= region$start) %>%
    dplyr::summarize(ymin=min(y)) %>%
    dplyr::pull(ymin)
  
  data <- data %>% 
    dplyr::mutate(
      y = pmax(y, 1.05*ymin),
      x = pmax(region$start, pmin(region$end, x))
    )

  if (is.null(color_by)) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x, y, group=loop_id))
  } else {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x, y, color=color, group=loop_id))
  }
  
  plot <- plot + 
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_comma(big.mark=" ")) +
    ggplot2::scale_y_continuous(limits = c(1.05*ymin, 0), labels=NULL, breaks=NULL, expand = c(0, 0)) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(x = "Genomic Position (bp)", y = NULL) +
    ggplot2::facet_wrap("facet_label", strip.position="left") + 
    trackplot_theme()
  
  if (!is.null(color_by)) {
    plot <- plot +
      ggplot2::scale_color_gradientn(limits=c(min(data$color), max(data$color)), colors=colors) +
      ggplot2::labs(color = color_by)
  } 
  
  wrap_trackplot(plot, ggplot2::unit(1, "null"))
}

#' Plot scale bar
#'
#' Plots a human-readable scale bar and coordinates of the region being plotted
#' 
#' @param font_pt Font size for scale bar labels in units of pt.
#' @inheritParams trackplot_coverage
#' 
#' @return Plot of loops
#' @export
trackplot_scalebar <- function(region, font_pt=11) {
  region <- BPCells:::normalize_ranges(region)
  breaks <- pretty(c(region$start, region$end))
  width <- diff(breaks)[1]
  
  width_text <- scales::label_number(scale_cut = scales::cut_si("b"))(width)
  bar_data <- tibble::tibble(
    right = region$end - 0.05*(region$end-region$start),
    left = right - width,
  )
  scale_label_data <- tibble::tibble(
    right = bar_data$left,
    text = sprintf("%s", width_text)
  )
  number_format <- scales::label_comma(big.mark=" ")
  region_data <- tibble::tibble(
    left = region$start,
    text = sprintf("%s: %s - %s", region$chr, number_format(region$start), number_format(region$end-1L))
  )
  
  plot <- ggplot2::ggplot() +
    ggplot2::geom_text(data=region_data, ggplot2::aes(x=left, y=0, label=text), size=10/ggplot2::.pt, hjust="left") +
    ggplot2::geom_text(data=scale_label_data, ggplot2::aes(x=right, y=0, label=text), size=10/ggplot2::.pt, hjust=1.1) +
    ggplot2::geom_errorbar(data=bar_data, ggplot2::aes(xmin=left, xmax=right, y=0), width=1) +
    ggplot2::scale_x_continuous(limits = c(region$start, region$end), expand = c(0, 0), labels=scales::label_comma(big.mark=" ")) +
    ggplot2::theme_void() 
  
  wrap_trackplot(plot, height=ggplot2::unit(font_pt*1.1, "pt"))
}

trackplot_helper_v1 <- function(gene, clusters, fragments, cell_read_counts, transcripts, loops, loop_color, expression_matrix, flank=1e5, region=NULL) {
  if (is.null(region)) {
    region <- gene_region(transcripts, gene, extend_bp = flank)
  } else {
    region <- BPCells:::normalize_ranges(region)
  }
  pal <- discrete_palette("stallion")
  
  base_size <- 11
  scale_plot <- trackplot_scalebar(region, font_pt=base_size)
  bulk_plot <- trackplot_coverage(
    fragments,
    region=region, 
    groups=clusters,
    cell_read_counts=cell_read_counts,
    colors=pal,
    bins=500
  )
  
  gene_plot <- trackplot_gene(transcripts, region) + ggplot2::guides(color="none")
  
  expression <- collect_features(expression_matrix, gene)
  names(expression) <- "gene"
  
  expression_plot <- ggplot2::ggplot(expression, ggplot2::aes(clusters, gene, fill=clusters)) +
    ggplot2::geom_boxplot() + 
    ggplot2::guides(y="none", fill="none") + 
    ggplot2::labs(x=NULL, y="RNA") +
    ggplot2::scale_fill_manual(values=pal, drop=FALSE) +
    trackplot_theme()
  
  loop_plot <- trackplot_loop(loops, region, color_by=loop_color, track_label="Co-Accessibility")
  
  trackplot_combine(list(scale_plot, bulk_plot, gene_plot, loop_plot), side_plot=expression_plot, title=gene)
}

gene_region <- function(genes, gene_symbol, extend_bp = 1e4, gene_mapping = human_gene_mapping) {
  genes <- BPCells:::normalize_ranges(genes, metadata_cols = c("strand", "gene_name"))
  idx <- BPCells::match_gene_symbol(gene_symbol, genes$gene_name)
  if (is.na(idx)) {
    rlang::stop("Could not locate gene")
  }
  if (length(extend_bp) == 1) {
    extend_bp <- c(extend_bp, extend_bp)
  }
  if (isFALSE(genes$strand[idx])) {
    extend_bp <- rev(extend_bp)
  }
  return(list(
    chr = as.character(genes$chr[idx]),
    start = genes$start[idx] - extend_bp[1],
    end = genes$end[idx] + extend_bp[2]
  ))
}

trackplot_helper_v2b <- function(gene=NULL, region=NULL, clusters, fragments, cell_read_counts,
                                 transcripts, loops, loop_color, expression_matrix=NULL,
                                 flank=1e5, loop_label="loop", plot_title="", highlights=NULL, cmap=NULL) {
  
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
  
  if (is.null(cmap)){
    pal <- BPCells::discrete_palette("stallion", n=length(unique(clusters)))
  } else{
    pal <- cmap
  }
  
  base_size <- 11
  scale_plot <- trackplot_scalebar(region, font_pt=base_size)
  bulk_plot <- trackplot_coverage(
    fragments,
    region=region,
    groups=clusters,
    cell_read_counts=cell_read_counts,
    colors=pal,
    bins=500
  )
  
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
  
  gene_plot <- trackplot_gene(transcripts, region) + ggplot2::guides(color="none")
  
  loop_plot <- trackplot_loop(loops, region, color_by=loop_color, track_label=loop_label)
  
  if (!is.null(gene) & !is.null(expression_matrix)){
    expression <- collect_features(expression_matrix, gene)
    names(expression) <- "gene"
    
    expression_plot <- ggplot2::ggplot(expression, ggplot2::aes(clusters, gene, fill=clusters)) +
      ggplot2::geom_boxplot() +
      ggplot2::guides(y="none", fill="none") +
      ggplot2::labs(x=NULL, y="RNA normalized expression") +
      ggplot2::scale_fill_manual(values=pal, drop=FALSE) +
      trackplot_theme()
    
    trackplot_combine(list(scale_plot, bulk_plot, gene_plot, loop_plot), side_plot=expression_plot, title=gene)
  } else{
    trackplot_combine(list(scale_plot, bulk_plot, gene_plot, loop_plot), title=plot_title) &
      ggplot2::theme(legend.direction="vertical")
  }
  
}
