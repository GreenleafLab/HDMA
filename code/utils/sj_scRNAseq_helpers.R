# Plotting & utility functions adapted from:
# https://github.com/fungenomics/HGG-oncohistones/blob/master/code/functions/scRNAseq.R
# Adapted to Greenleaf lab style

library(glue)
library(RColorBrewer)

# white -> red
whrd <- grDevices::colorRampPalette(c("white", "#FFC0C0", "red"))(n = 200)


# yellow -> red
ylrd    <- grDevices::colorRampPalette(brewer.pal(8, "OrRd"))(n = 100)

# purple -> blue -> green
purpgn    <- grDevices::colorRampPalette(brewer.pal(8, "PuBuGn")[2:6])(n = 100)

# blue -> red
rdbu    <- rev(grDevices::colorRampPalette(brewer.pal(8, "RdBu"))(n = 100))
rdbu2   <- colorRampPalette(c("navy", "white", "red"))(n = 100)
rdbu3   <- colorRampPalette(c("navy", "#2E2E97", "#4848A4", "white",
                              "#FF4E4E", "#FF3434", "red"))(n = 100)


# Plotting functions -----------------------------------------------------------

#' Colour cells in reduced dimension space by gene expression
#'
#' Plot a low-dimensional embedding of the cells,
#' coloured by expression of a gene, or mean expression of a group of marker
#' genes. Defaults to UMAP space, but see the \code{reduction} argument for
#' how to plot in t-SNE or PCA space instead. This function is based on 
#' \code{Seurat::FeaturePlot} and cytobox::tsneByMeanMarkerExpression.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA(), Seurat::RunTSNE(), or Seurat::RunUMAP() to the object)
#' @param genes String or character vector specifying gene(s) to use
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves UMAP by default.
#' Default: "umap"
#' @param label Logical, whether to label clusters on the plot. Default: TRUE.
#' @param label_repel Logical, if \code{label} is TRUE, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param palette String or character vector. If a string,
#' one of "viridis", "blues", or "redgrey", specifying which gradient
#' palette to use. Otherwise, a character vector of colours (from low to high)
#' to interpolate to create the scale. Default: redgrey.
#' @param point_size Numeric, size of points in scatterplot. Default: 1. (A smaller
#' value around 0.5 is better for plots which will be viewed at small scale.)
#' @param alpha Numeric, fixed alpha for points. Default: 0.6
#' @param legend Logical, whether or not to plot legend. Default: TRUE
#' @param hide_ticks Logical, whether to hide axis ticks. Default: FALSE
#' @param hide_axes Logical, whether to hide axis labels. Default: FALSE
#' @param title (Optional) String specifying the plot title
#' @param limits (Optional) A numeric vector of length two providing the limits to
#' use for the colour scale (documentation from \code{\link[ggplot2]{continous_scale}}. 
#' Default: 0 and max of the data.
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at levels(Ident(seurat))) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param dim1 Numeric, index of dimension from \code{reduction} to plot on
#' the x-axis. e.g. to plot the 3rd PC on the x-axis, pass 3. Default: 1.
#' @param dim2 Numeric, like \code{dim2}, but for the y-axis. Default: 2.
#' @param return_df Logical, whether to return the mean expression dataframe.
#'
#' @export
#' @return A ggplot object
#'
#' @author Selin Jessa, Samantha Worme
plot_mean_expression <- function(seurat, genes,
                                 assay = DefaultAssay(seurat),
                                 reduction = "umap",
                                 label = TRUE,
                                 label_repel = TRUE,
                                 label_size = 4,
                                 palette = "redgrey",
                                 point_size = 1,
                                 alpha = 0.6,
                                 legend = TRUE,
                                 hide_ticks = FALSE,
                                 hide_axes = FALSE,
                                 title = NULL,
                                 limits = NULL,
                                 label_short = FALSE,
                                 dim1 = 1,
                                 dim2 = 2,
                                 return_df = FALSE,
                                 theme = "BOR",
                                 rasterize = TRUE) {
  
  exp_df <- suppressWarnings(mean_gene_expression(seurat, genes, assay))
  
  # Get dimensionality reduction coordinates
  exp_df <- get_embedding(seurat, assay, exp_df, reduction, dim1, dim2) %>%
    # Order in which points will be plot, "front" points at the bottom
    dplyr::arrange(Mean_marker_expression)
  
  if (return_df) return(exp_df)
  
  # Get the variable names
  vars <- colnames(seurat[[reduction]]@cell.embeddings)[c(dim1, dim2)]
  
  # Set limits: if not provided, use default min/max
  if (is.null(limits)) limits <- c(NA, NA)
  
  # Plot
  gg <- exp_df %>%
    ggplot(aes(x = .data[[vars[1]]], y = .data[[vars[2]]]))
  
  if (rasterize) gg <- gg + geom_point_rast(aes(colour = Mean_marker_expression), size = point_size, alpha = alpha)
  else gg <- gg + geom_point(aes(colour = Mean_marker_expression), size = point_size, alpha = alpha)
  
  if (length(palette) == 1) {
    
    if (palette == "viridis") {
      
      gg <- gg + viridis::scale_color_viridis(limits = limits)
      
    } else if (palette == "blues") {
      
      gg <- gg + scale_colour_gradientn(
        colours = RColorBrewer::brewer.pal(n = 8, name = "Blues"),
        limits = limits)
      
    } else if (palette == "redgrey") {
      
      # NOTE: palette chosen is not the default gradient from gray -> red
      # but sets a midpoint at a lighter colour
      gg <- gg + scale_color_gradientn(
        colours = grDevices::colorRampPalette(c("gray93", "#D9B5B5", "red"))(n = 200),
        # colours = grDevices::colorRampPalette(c("gray83", "#E09797", "red"))(n = 200),
        limits = limits)
      
    } else {
      
      stop("Please pass the palette as a character vector ",
           "or specify one of: viridis, blues, redgrey")
      
    }
    
  } else if (length(palette) == 2) {
    
    gg <- gg + scale_color_gradient(low = palette[1], high = palette[2], limits = limits)
    
  } else {
    
    gg <- gg + scale_color_gradientn(colours = palette, limits = limits)
    
  }
  
  
  if (label) {
    
    if (reduction == "pca") {
      
      message("Plotting labels is currently only available for reduction = 'tsne';",
              " returning plot without labels.")
      
    } else {
      
      centers <- suppressMessages(compute_cluster_centers(seurat, reduction = reduction))
      gg <- gg + add_labels(centers, label_repel, label_size, label_short)
      
    }
  }
  
  axes <- gsub("_", " ", vars)
  
  gg <- gg +
    xlab(axes[1]) +
    ylab(axes[2]) +
    labs(colour = case_when(assay == "RNA"                            ~ "Expression",
                            assay == "SCENIC"                         ~ "Regulon AUC",
                            assay == "chromVAR"                       ~ "Motif activity",
                            assay %in% c("peaks", "ACTIVITY", "ATAC", "promoters") ~ "Accessibility",
                            TRUE                                      ~ "Score"))
  
  if (theme == "BOR") gg <- gg + theme_BOR()
  else gg <- gg + theme_min()
  
  if (!is.null(title)) gg <- gg + ggtitle(title)
  if (!legend) gg <- gg + hide_legend()
  if (hide_ticks) gg <- gg + hide_ticks()
  if (hide_axes) gg <- gg + xlab(NULL) + ylab(NULL)
  
  return(gg)
  
}



#' @describeIn tsneByMeanMarkerExpression Shortcut function for plotting mean expression
#' @export
feature <- function(seurat, genes,
                    assay = DefaultAssay(seurat),
                    per_gene = TRUE,
                    label = TRUE,
                    palette = "redgrey",
                    label_repel = FALSE,
                    label_size = 4,
                    label_short = FALSE,
                    legend = FALSE,
                    title = NULL,
                    reduction = "umap",
                    limits = c(NA, NA),
                    dim1 = 1,
                    dim2 = 2,
                    alpha = 0.6,
                    point_size = 0.3,
                    ncol = ifelse(length(genes) == 1, 1, ifelse(length(genes) %in% c(2, 4), 2, 3)),
                    hide_ticks = TRUE,
                    hide_axes = TRUE,
                    combine = TRUE,
                    theme = "BOR",
                    rasterize = FALSE) {
  
  if ((length(genes) >= 20) & per_gene) message("NOTE: you have input a lot of genes! ",
                                                "This function by default generates ",
                                                "one plot per gene. Set per_gene = FALSE ",
                                                "to plot a summary statistic of all genes.")
  
  
  if (per_gene) {
    
    genes_out <- find_genes(seurat, genes, assay = assay)
    if (length(genes_out$undetected > 0)) message(paste0("NOTE: [",
                                                         paste0(genes_out$undetected, collapse = ", "),
                                                         "] undetected in the data"))
    
    if(length(genes_out$detected) == 0) stop("No genes specified were ",
                                             "found in the data.")
    
    if ((ncol == 3) & (length(genes_out$detected) < 3)) ncol <- 2
    
    # this is a non-elegant way to ensure the supplied title is used in the one-gene case...
    if (!is.null(title) & length(genes) == 1) titles <- title
    else titles <- genes_out$detected
    
    plots <- lapply(seq_along(genes_out$detected),
                    function(i) plot_mean_expression(seurat,
                                                     genes_out$detected[i],
                                                     assay = assay,
                                                     reduction = reduction,
                                                     palette = palette,
                                                     title = titles[i],
                                                     legend = legend,
                                                     label = label,
                                                     label_short = label_short,
                                                     label_repel = label_repel,
                                                     label_size = label_size,
                                                     hide_ticks = hide_ticks,
                                                     hide_axes = hide_axes,
                                                     point_size = point_size,
                                                     limits = limits,
                                                     dim1 = dim1,
                                                     dim2 = dim2,
                                                     rasterize = rasterize))
    
    if (combine) cowplot::plot_grid(plotlist = plots, ncol = ncol)
    else return(plots)
    
  } else {
    
    plot_mean_expression(seurat,
                         genes,
                         assay = assay,
                         reduction = reduction,
                         palette = palette,
                         title = title,
                         legend = legend,
                         label = label,
                         label_repel = label_repel,
                         label_size = label_size,
                         label_short = label_short,
                         hide_ticks = hide_ticks,
                         hide_axes = hide_axes,
                         point_size = point_size,
                         limits = limits,
                         dim1 = dim1,
                         dim2 = dim2,
                         rasterize = rasterize)
  }
}


#' Plot a reduced dimensionality embedding for a Seurat object, coloring
#' by cluster or arbitrary variables.
#'
#' @param seurat Seurat object, where dimensionality reduction passed to \code{reduction} has been applied
#' @param reduction String, specifying a lot of \code{seurat@@reductions}, which
#' indicates which embedding to plot. Default: "umap". (Can also take "tsne" or "pca").
#' @param color_by (Optional) String, specifying the column in \code{seurat@@meta.data}
#' by which to color cells. Default: NULL, color cells by cluster (in \code{seurat@@ident}).
#' @param colors (Optional) Character vector of colors for points. If \code{color_by}
#' is NULL, cells will be colored by cluster; this should be a named character vector of colors for points. Names should
#' correspond to cluster names (e.g. \code{levels(Idents(seurat)}). If
#' specifying \code{color_by}, and the variable is discrete, this should be a character vector,
#' with either names or order corresponding to categorical
#' values in the column of \code{seurat@@meta.data} specified. Otherwise, if the variable
#' is continuous, pass the gradient to use, or a few colors (from low to high) from which a gradient
#' should be created, and specify \code{color_by_type = "continuous"}. The default is to
#' use ggplot2 colors.
#' @param color_by_type (Optional) String, one of "discrete" or "continuous".
#' If specifying \code{color_by} and providing colors to the \code{colors}
#' argument, specify whether the \code{color_by} variable is discrete or continuous.
#' Default: discrete. The function is pretty good at setting the palette correctly on
#' its own, but try modifying this in case of errors.
#' @param label Logical, whether to plot cluster labels. Default: TRUE
#' @param point_size Numeric, size of points in scatter plot. Default: 0.6 for datasets
#' with < 300 cells, and 1.3 for datasets with otherwise.
#' @param alpha Numeric, fixed alpha value for points: Default: 0.8
#' @param legend Logical, whether to plot legend. Default: FALSE if \code{color_by}
#' is NULL and \code{label} is TRUE, true otherwise.
#' @param label_repel Logical, if \code{label} is TRUE, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param title (Optional) String specifying title.
#' @param hide_ticks Logical, whether to hide axis ticks, i.e. both the text and the
#' small lines indicating the breaks along the x- and y-axes. Default: FALSE
#' @param hide_axes Logical, whether to hide axis labels. Default: TRUE
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes cluster
#' names (at seurat@@ident) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param cells Character vector of cell names if only a subset of cells should be
#' colored in the plot (these should correspond to seurat@@cell.names). Default: Plot all cells.
#' See the argument \code{clusters_to_label} for only labelling certain clusters.
#' See the \code{constrain_scale} argument for controlling the scales of the plot.
#' @param show_all_cells Logical. When passing cells to \code{highlight_cells},
#' if TRUE, plot the remaining cells in \code{na_color}, if FALSE,
#' only plot the highlighted cells. Default: TRUE.
#' @param order_by String, corresponding to a column in \code{seurat@@meta.data}, specifying
#' a variable to control the order in which cells are plot. (Thus, you can manually
#' specify the order, add it as a new column in \code{seurat@@meta.data}, and pass that).
#' If numeric, cells with high values are plot on top. If not, the column must
#' be a factor, and cells will be ordered according to the levels, with cells
#' in the first level plot on top. Default: if a numeric column is specified
#' to \code{color_by}, sort by that variable, otherwise, use the ordering of the cells
#' in the Seurat object.
#' @param clusters_to_label (Optional.) If \code{label} is TRUE,
#' clusters for which labels should be plot (if only a subset of clusters should be labelled).
#' Default: NULL (Label all clusters).
#' @param na_color String, specifying the color (built-in or hex code) to use to
#' plot points which have an NA value, for example
#' in the variable specified in \code{color_by}. Default: light gray ("gray80),
#' change to "white" to purposely hide those cells. If you do not want to plot
#' certain cells at all, pass names of cells to plot to the \code{cells} argument.
#' @param limits Numeric vector of length two providing the lower and upper limits of
#' the color scale, if coloring by a continuous variable. Default: min and max
#' of the values the variable takes on in the data.
#' @param constrain_scale Logical, if plotting a subset of cells, whether to
#' use the limits of the tSNE embedding computed on the whole dataset (useful
#' for constraining scales across plots while only plotting specific cells).
#' Default: TRUE
#' @param dim1 Numeric, dimension of embedding to use for x-axis. Default = 1.
#' @param dim2 Numeric, dimension of embedding to use for y-axis. Default = 2.
#'
#' @return A ggplot2 object
#' @export
#' @aliases plot_pca plot_tsne plot_umap
#'
#' @author Selin Jessa and Nisha Kabir
#'
#' @examples
#' tsne(pbmc)
#'
#' # Demonstrate label_short:
#' # Set cluster IDs to be too long
#' pbmc2 <- pbmc
#' levels(Seurat::Idents(pbmc2)) <- c("1-Cell type A", "2-Cell type B", "3-Cell type C", "4-Cell type D")
#' tsne(pbmc2)
#'
#' # Plot the prefixes only
#' tsne(pbmc2, label_short = TRUE)
#'
#' tsne(pbmc2, color_by = "nCount_RNA", color_by_type = "continuous")
plot_dr <- function (seurat,
                     reduction         = "umap",
                     color_by          = NULL,
                     colors            = NULL,
                     cmap              = NULL,
                     color_by_type     = "discrete",
                     label             = TRUE,
                     label_repel       = TRUE,
                     label_size        = 4,
                     # Decide default point size based on density / how many points are plot
                     point_size = ifelse(length(colnames(seurat)) > 300, 0.6, 1.3),
                     alpha             = 0.8,
                     # If coloring by clusters, hide legend by default, otherwise, show it
                     legend = ifelse((is.null(color_by)) && (label), FALSE, TRUE),
                     cells             = NULL,
                     show_all_cells    = TRUE,
                     order_by          = NULL,
                     clusters_to_label = NULL,
                     hide_ticks        = TRUE,
                     title             = NULL,
                     label_prefix_only = FALSE,
                     label_suffix_only = FALSE,
                     sep               = NULL,
                     na_color          = "gray80",
                     limits            = NULL,
                     constrain_scale   = TRUE,
                     hide_axes         = FALSE,
                     dim1              = 1,
                     dim2              = 2,
                     border_color      = NULL,
                     border_size       = NULL,
                     theme             = "BOR",
                     rasterize         = TRUE) {
  
  # Error if the selected reduction has not yet been computed
  if (!(reduction %in% names(seurat@reductions))) stop(reduction, " reduction has not been computed.")
  
  # Create a dataframe holding the embedding
  embedding <- data.frame(Cell = colnames(seurat),
                          dim1 = seurat[[reduction]]@cell.embeddings[, dim1],
                          dim2 = seurat[[reduction]]@cell.embeddings[, dim2],
                          Cluster = Idents(seurat), stringsAsFactors = FALSE)
  
  # If required, order the cells in the z-axis by the specified variable
  if (!is.null(order_by)) {
    
    # We can only order by a nuemric variable, or a discrete variable with specified order
    # of levels (i.e. a factor)
    if (!is.numeric(seurat@meta.data[[order_by]]) && !is.factor(seurat@meta.data[[order_by]])) {
      
      stop("The variable specified in 'order_by' is neither numeric ",
           "nor a factor. If the column is of type character, consider ",
           "converting it to a factor. Otherwise, pass the name of a numeric column.")
      
    }
    
    embedding[[order_by]] <- seurat@meta.data[[order_by]]
    
    
  } else if ((!is.null(color_by)) && is.numeric(seurat@meta.data[[color_by]])) {
    # If it's not specified, by default, order by the same variable used
    # to color cells
    
    order_by <- color_by
    
  }
  
  if (!is.null(color_by)) embedding[[color_by]] <- seurat@meta.data[[color_by]]
  
  # Highlight certain cells in the plot
  if (!is.null(cells)) {
    
    # Show all cells, but only color the selected ones, by setting the color_by
    # value to NA for all other cells
    if (show_all_cells) {
      
      if (is.null(color_by)) embedding[!(embedding$Cell %in% cells), ]$Cluster <- NA
      else embedding[!(embedding$Cell %in% cells), ][[color_by]] <- NA
      
    } else { # Otherwise, only display and color the selected ones
      
      embedding <- embedding %>% filter(Cell %in% cells)
      
    }
  }
  
  # We want to sort point such that any NAs will be plot first/underneath
  color_by2 <- ifelse(is.null(color_by), "Cluster", color_by)
  
  # If the ordering variable is specified, sort NA cells at the bottom, and then order by that variable
  if (!is.null(order_by)) embedding <- embedding %>% arrange_(paste0("!is.na(", color_by2, ")"), order_by)
  else embedding <- embedding %>% arrange_(paste0("!is.na(", color_by2, ")"))
  
  # Start the plot!
  gg <- ggplot(embedding, aes(x = dim1, y = dim2))
  
  # Make sure that there are idents available to use for labeling cells / clusters
  if (label && all(is.na(Idents(seurat)))) {
    
    label <- FALSE
    message("NOTE: identity of all cells is NA, setting 'label' to FALSE.")
    
  }
  
  # Deal with the palette
  if (is.null(color_by)) {
    
    if (is.null(colors)) {
      
      # Check if there's a cluster palette stored with the object
      if (!is.null(seurat@misc$colors)) colors <- seurat@misc$colors
      else if (!is.null(seurat@misc$colours)) colors <- seurat@misc$colours
      # Following the Greenleaf lab's style if a cmap is provided, choose a new color palette
      else if (!is.null(cmap)) {
        
        
        colors <- getColorMap(cmap, n = length(levels(Idents(seurat))))
        names(colors) <- levels(Idents(seurat))
        
      } else {
        
        # Assuming that the order of the levels is correct in the seurat object,
        # find the default colors for the clusters
        colors <- get_gg_colors(length(levels(Idents(seurat))))
        names(colors) <- levels(Idents(seurat))
        
      }
      
    }
    
    if (rasterize) gg <- gg + geom_point_rast(aes(color = Cluster), size = point_size, alpha = alpha)
    else gg <- gg + geom_point(aes(color = Cluster), size = point_size, alpha = alpha)
    
    # add colors
    gg <- gg +
      scale_color_manual(values = colors, na.value = na_color)
    
  } else {
    
    # If no limits are specified for the z-axis (variable used to color cells),
    # use the ggplot2 default, which is to fit range of the data
    if (is.null(limits)) lims <- c(NA, NA)
    else lims <- limits
    
    # Add the points to the plot object
    gg <- gg +
      geom_point(aes_string(color = color_by), size = point_size, alpha = alpha)
    
    # If a palette is provided, choose an appropriate scale based on whether
    # the coloring variable is discrete or continuous
    if (!is.null(colors)) {
      
      if (color_by_type == "discrete") gg <- gg + scale_color_manual(values = colors, na.value = na_color)
      else if (color_by_type == "continuous") {
        
        gg <- gg + scale_color_gradientn(colors  = colors,
                                         na.value = na_color,
                                         limits   = lims)
      }
      
    } else {
      
      if (color_by_type == "continuous") {
        
        gg <- gg + scale_color_gradientn(colors   = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "OrRd"))(n = 100),
                                         na.value = na_color,
                                         limits   = lims)
        
      }
    }
  }
  
  # Label clusters on the plot
  if (label) {
    
    centers <- compute_cluster_centers(seurat,
                                       reduction = reduction,
                                       dim1 = dim1, dim2 = dim2)
    
    gg <- gg + add_labels(centers           = centers,
                          label_repel       = label_repel,
                          label_size        = label_size,
                          label_prefix_only = label_prefix_only,
                          label_suffix_only = label_suffix_only,
                          sep               = sep,
                          clusters          = clusters_to_label)
  }
  
  # Add the default theme to the plot object
  if (theme == "BOR") gg <- gg + theme_BOR()
  else if (is.null(style)) gg <- gg + theme_min(border_size = border_size)
  
  # Set the right axes titles
  if (hide_axes) gg <- gg + xlab(NULL) + ylab(NULL)
  else if (reduction == "tsne") gg <- gg + xlab(glue("tSNE {dim1}")) + ylab(glue("tSNE {dim2}"))
  else if (reduction == "umap") gg <- gg + xlab(glue("UMAP {dim1}")) + ylab(glue("UMAP {dim2}"))
  else if (reduction == "phate") gg <- gg + xlab(glue("PHATE {dim1}")) + ylab(glue("PHATE {dim2}"))
  else if (reduction == "pca") {
    
    # If the reduction is PCA, automatically compute the variance explained
    var_exp <- getVarianceExplained(seurat)
    gg <- gg +
      xlab(glue("PC{dim1} ({round(var_exp$percent.var.explained[dim1], 1)}%)")) +
      ylab(glue("PC{dim2} ({round(var_exp$percent.var.explained[dim2], 1)}%)"))
    
  } else {
    gg <- gg + xlab(glue("{reduction} {dim1}")) + ylab(glue("{reduction} {dim2}"))
  }
  
  # Remove legend
  if (!legend) gg <- gg + hide_legend()
  else {
    
    gg <- gg + guides(color = guide_legend(override.aes = list(size=5)))
    
    if (!is.null(color_by)) {
      
      # This column typically contains the sample name in the seurat metadata;
      # use this more informative name for the legend title instead of "orig.ident"
      if (color_by == "orig.ident") gg <- gg + labs(color = "Sample")
      
    }
  }
  
  # Plot aesthetics
  if (!is.null(title)) gg <- gg + ggtitle(title)
  else {
    
    n_cells <- ncol(seurat)
    title <- glue("n = {n_cells}")
    gg <- gg + ggtitle(title)
    
  }
  
  if (hide_ticks) gg <- gg + hide_ticks()
  
  if (constrain_scale) gg <- gg + constrain_scale(seurat,
                                                  reduction = reduction,
                                                  dim1      = dim1,
                                                  dim2      = dim2)
  
  return(gg)
  
}


#' @describeIn plot_dr Plot a tSNE embedding
#' @export
tsne <- function(seurat, ...) {
  
  plot_dr(seurat, reduction = "tsne", ...)
  
}


#' @describeIn plot_dr Plot a PCA embedding
#' @export
pca <- function(seurat, ...) {
  
  plot_dr(seurat, reduction = "pca", ...)
  
}


#' @describeIn plot_dr Plot a UMAP embedding
#' @export
umap <- function(seurat, ...) {
  
  plot_dr(seurat, reduction = "umap", ...)
  
}




violin_grid <- function(seurat, genes,
                        group_col,
                        group_order = NULL,
                        order = "genes",
                        colours = NULL,
                        scale = "width",
                        title = NULL,
                        scales = "free_x") {
  
  expr <- FetchData(seurat, vars = c(group_col, genes)) %>%
    tidyr::gather(Marker, Expression, 2:ncol(.))
  expr$Cell <- colnames(seurat)
  expr$Cluster <- expr[[group_col]]
  expr[[group_col]] <- NULL
  
  # deal with all-0 expression
  maxes <- expr %>% group_by(Marker) %>% summarise(max = max(Expression)) %>% tibble::deframe()
  
  # convert to NA b/c of the jitter
  no_expr <- names(maxes)[maxes == 0]
  if (length(no_expr) > 0) expr[expr$Marker %in% no_expr, ]$Expression <- NA
  
  # complete for non-existent genes
  expr$Marker  <- factor(expr$Marker, levels = genes)
  expr$Cluster <- factor(expr$Cluster, levels = unique(expr$Cluster))
  expr <- expr %>% 
    complete(Cluster, nesting(Marker), fill = list(Expression = NA)) %>% 
    filter(!is.na(Cluster))
  
  genes_keep <- genes[genes %in% expr$Marker]
  
  if (length(order) == 1) {
    
    if (order == "sort") {
      
      sort_criteria <- glue::glue("desc({genes_keep})")
      
      cluster_order <- expr %>%
        filter(Expression != 0) %>%
        tidyr::spread(Marker, Expression) %>%
        group_by(Cluster) %>%
        dplyr::select(-Cell) %>% 
        summarise_all(funs(median), na.rm = TRUE) %>%
        arrange_(sort_criteria) %>%
        `[[`("Cluster")
      
      expr$Cluster <- factor(expr$Cluster, levels = rev(cluster_order))
      
    } else stop("Please set order = 'genes', or provide an ordering of clusters.")
    
  } else if (order == "genes") {
    
    
  } else if (!is.null(group_order)) {
    
    expr$Cluster <- factor(expr$Cluster, levels = rev(group_order))
    
  }
  
  expr$Marker <- factor(expr$Marker, levels = genes_keep)
  
  # Set colours to the levels of the clusters, so that they are preserved
  # even if an order was specified
  
  
  if (!is.null(colours) && length(colours) < length(unique(expr$Cluster))) {
    
    stop("Not enough colours! Please provide as many colours ",
         "as clusters in the dataset, or one per cluster specified in the ",
         "'subset_clusters' argument.")
    
    
  } else if (is.null(colours)) {
    
    colours <- ggColours(length(unique(expr$Cluster)))
    names(colours) <- unique(expr$Cluster)
    
  }
  
  gg <- expr %>%
    ggplot(aes(x = Cluster, y = Expression)) +
    geom_violin(aes(fill = Cluster), scale = scale, size = 0.5) +
    scale_fill_manual(values = colours) +
    facet_wrap(~ Marker, ncol = length(unique(expr$Marker)),
               scales = scales) +
    theme_min() +
    coord_flip() +
    ggplot2::theme(panel.border = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.line.x = element_blank(),
                   legend.position = "none",
                   strip.text.x = element_text(angle = 30, size = 8))
  
  if (!is.null(title)) gg <- gg + ggtitle(title)
  
  return(gg)
  
}


vln <- function(seurat, features, point_size = 0, cmap, ncol = NULL, assay = NULL, y_max = NULL) {
  
  colors <- getColorMap(cmap, n = length(levels(Idents(seurat))))
  names(colors) <- levels(Idents(seurat))
  
  if (is.null(assay)) assay <- DefaultAssay(seurat)
  
  VlnPlot(seurat, features, pt.size = point_size, cols = colors, ncol = ncol, assay = assay, y.max = y_max) +
    hide_legend()
  
}


# Summary stats ----
#' mean_gene_expression
#'
#' Calculate, for each cell, the mean expression of the given set of marker genes using the
#' normalized data.
#'
#' @param seurat Seurat object
#' @param genes String or character vector specifying gene(s) to use
#'
#' @return Data frame with two columns: "Cell" which specifies the cell ID
#' and "Mean_marker_expression" which is the expression value for the marker gene, 
#' or mean if multiple genes were provided.
#'
#' @export
#' @author Selin Jessa
#' @examples
#' mean_gene_expression(pbmc, genes = c("IL32", "MS4A1"))
mean_gene_expression <- function(seurat, genes, assay = "RNA") {
  
  # Get expression data from the Seurat object
  data.frame(Cell = colnames(GetAssayData(seurat, assay = assay)),
             Mean_marker_expression = rowMeans(fetch_data(seurat, genes, assay)),
             stringsAsFactors = FALSE)
  
}





# Utility functions ------------------------------------------------------------

#' find_genes
#'
#' Given a set of genes, find the ones which are detected in the sample,
#' and which are not.
#'
#' @param seurat Seurat object
#' @param genes Character vector of genes
#'
#' @return A named list with two elements: "detected" and "undetected"
#' each storing character vectors with the genes in each category
#' @export
#'
#' @author Selin Jessa
#' @examples
#' find_out <- find_genes(pbmc, c("IL32", "CD79B", "foo"))
#' find_out$detected
#' find_out$undetected
find_genes <- function(seurat, genes, assay = "RNA") {
  
  genes_detected <- genes[genes %in% rownames(GetAssayData(object = seurat, assay = assay))]
  genes_undetected <- setdiff(genes, genes_detected)
  
  list(detected = genes_detected, undetected = genes_undetected)
  
}



#' fetchData
#'
#' Subset the seurat@@data matrix by gene and cluster.
#' Similar to SeuratObject::FetchData except it doesn't thrown an error if a gene
#' is not found in the data, and accepts a specific assay argument.
#'
#' @param seurat Seurat object
#' @param genes Genes to filter
#' @param assay Assay to search for genes
#' @param clusters (Optional) Vector, include only cells with these identities
#' (e.g. cluster assignments). Searches in seurat@@active.ident.
#' @param return_cell Logical, whether or not to include a column with the cell ID.
#' Default: FALSE
#' @param return_cluster Logical, whether or not to include a column with the cluster.
#' Default: FALSE
#' @param scaled Logical, whether to fetch scaled data the assay.
#' Default: FALSE.
#'
#' @return Expression matrix for genes specified
#' @export
#'
#' @author Selin Jessa
#' @examples
#' fetchData(pbmc, c("IL32", "MS4A1"))
#' fetchData(pbmc, c("IL32"), c(1, 2))
#' fetchData(pbmc, c("IL32", "MS4A1"), c(1, 2), return_cluster = TRUE, return_cell = TRUE)
fetch_data <- function(seurat, features, assay = "RNA", clusters = NULL,
                       return_cell = FALSE, return_cluster = FALSE, scaled = FALSE) {
  
  message("@ Searching assay: ", assay)
  
  features_out <- find_genes(seurat, genes = features, assay)
  
  n_undetected <- length(features_out$undetected)
  
  if (n_undetected > 0) {
    
    if (n_undetected > 10) {
      
      message(paste0("@ NOTE: [",
                     paste0(head(features_out$undetected, 10), collapse = ", "),
                     "] and ", n_undetected - 10, " other features are undetected in ", seurat@project.name))
      
    } else {
      
      message(paste0("@ NOTE: [",
                     paste0(features_out$undetected, collapse = ", "),
                     "] undetected in ", seurat@project.name))
      
    }
  }
  
  if (length(features_out$detected) == 0) stop("No features specified were ",
                                               "found in the data.")
  
  if (scaled) exp <- as.matrix(seurat@assays[[assay]]@scale.data)
  else exp <- as.matrix(seurat@assays[[assay]]@data)
  
  exp_filt <- as.data.frame(t(exp[which(rownames(exp) %in% features_out$detected),]))
  
  # Keep all
  if(is.null(clusters)) clusters <- levels(seurat@active.ident)
  ident_idx <- which(seurat@active.ident %in% clusters)
  
  # Handle only one gene case, and properly return a data frame
  if(nrow(exp_filt) == 1) {
    
    if(!is.null(clusters)) exp_filt <- exp_filt[ident_idx]
    
    exp_filt <- as.data.frame(t(exp_filt))
    names(exp_filt) <- rownames(exp)[rownames(exp) %in% features]
    
  } else {
    if(!is.null(clusters)) exp_filt <- exp_filt[ident_idx,]
  }
  
  rownames(exp_filt) <- c() # Get rid of rownames
  
  if(return_cluster) exp_filt <- tibble::add_column(exp_filt, Cluster = seurat@active.ident[ident_idx], .before = 1)
  if(return_cell) exp_filt <- tibble::add_column(exp_filt, Cell = names(seurat@active.ident)[ident_idx], .before = 1)
  
  return(exp_filt)
  
}



#' get_embedding
#' 
#' Given a Seurat object and a data frame where the rows correspond to cells and a Cell column
#' is present, the dataframe is altered to add two columns giving coordinates in a dimensionality
#' reduced space.
#' 
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA(), Seurat::RunTSNE(), Seurat::RunUMAP() to the object).
#' @param df Data frame with at least one column, "Cell", giving the cell ID
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves UMAP by default.
#' Default: "umap"
#' @param dim1 Numeric, dimension of embedding to use for x-axis. Default = 1.
#' @param dim2 Numeric, dimension of embedding to use for y-axis. Default = 2.
#' 
#' @export
#' @author Selin Jessa, Samantha Worme
get_embedding <- function(seurat, assay = "RNA", df, reduction = "umap", dim1 = 1, dim2 = 2) {
  
  # Get the axes for the reduced space
  # See here: http://dplyr.tidyverse.org/articles/programming.html#setting-variable-names
  vars <- colnames(seurat[[reduction]]@cell.embeddings)[c(dim1, dim2)]
  
  df$Cell <- as.character(df$Cell)
  
  embedding <- data.frame(Cell = colnames(GetAssayData(seurat, assay = assay)), stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      !!vars[1] := seurat[[reduction]]@cell.embeddings[, dim1], 
      !!vars[2] := seurat[[reduction]]@cell.embeddings[, dim2])
  
  df <- dplyr::inner_join(df, embedding, by = "Cell")
  return(df)
  
}


#' Rotate the x axis labels in a ggplot
#'
#' @param angle Integer, value in degrees to rotate labels. Default: 90.
#'
#' @return A theme element to rotate labels
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(color = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + rotateX()
rotate_x <- function(angle = 90) {
  
  theme(axis.text.x = element_text(angle = angle, hjust = 1))
  
}



#' Remove the legend in a ggplot
#'
#' @return A theme element to hide legend
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(color = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + noLegend()
hide_legend <- function() {
  
  theme(legend.position = "none")
  
}

no_legend <- hide_legend


#' Remove axis ticks and tick labels from a ggplot
#'
#' @return A theme element to remove ticks
#' @export
#'
#' @author Selin Jessa
hide_ticks <- function() {
  
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())
  
}


#' compute_cluster_centers
#'
#' Get centers of clusters given a Seurat object, to use for labelling
#' in tSNE space. The cluster center is defined as the median X and Y coordinate
#' across cells in each cluster.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunTSNE() to the object).
#'
#' @return Data frame with three columns: Cluster, mean_tSNE_1, and mean_tSNE_2
#' @export
#'
#' @author Selin Jessa
#' @examples
#'
#' compute_cluster_centers(pbmc, reduction = "pca", dim1 = 1, dim2 = 3)
compute_cluster_centers <- function(seurat, reduction = "umap", dim1 = 1, dim2 = 2) {
  
  n_clusters <- length(unique(Idents(seurat)))
  
  # Attempts at tidyeval...
  # vars <- colnames(seurat[[reduction]]@cell.embeddings)[c(1, 2)]
  # col_names <- paste0("mean_", vars)
  
  # Get the embedding
  df <- as.data.frame(seurat[[reduction]]@cell.embeddings[, c(dim1, dim2)]) %>%
    mutate(Cell = colnames(seurat),
           Cluster = Idents(seurat))
  
  # Generalize these
  colnames(df)[c(1, 2)] <- c("Dim_1", "Dim_2")
  
  # Compute cluster centers
  centers <- df %>%
    group_by(Cluster) %>%
    dplyr::summarise(mean_x = median(Dim_1),
                     mean_y = median(Dim_2))
  
  return(centers)
  
}


#' Add cluster labels to a deminsionality reduction plot
#'
#' @param centers Data frame with at least three columns: "mean_x", "mean_y",
#' and "Cluster", as returned by \code{\link{clusterCenters}}
#' @param label_repel Logical, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at \code{seurat@@ident}) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param clusters (Optional) Clusters for which labels should be plot (if only
#' a subset of clusters should be labelled). Default: NULL (Label all clusters).
#'
#'
#' @author Selin Jessa and Nisha Kabir
#' @export
add_labels <- function(centers, label_repel = FALSE, label_size = 4, label_prefix_only = FALSE, label_suffix_only = FALSE, sep = NULL, clusters = NULL) {
  
  if (!is.null(clusters)) centers <- filter(centers, Cluster %in% clusters)
  
  if (label_prefix_only) centers <- suppressWarnings(tidyr::separate(centers, Cluster, into = c("Cluster", "Cluster_long"), extra = "drop", sep = sep))
  else if (label_suffix_only) centers <- suppressWarnings(tidyr::separate(centers, Cluster, into = c("Cluster_long", "Cluster"), extra = "drop", sep = sep))
  
  if (label_repel) {
    
    ggrepel::geom_label_repel(data = centers,
                              aes(x = mean_x, y = mean_y),
                              label = centers$Cluster,
                              size = label_size,
                              segment.color = 'grey50',
                              face = 'bold',
                              alpha = 0.8,
                              segment.alpha = 0.8,
                              label.size = NA,
                              force = 1,
                              # Leaving these unspecified for now, since it really depends on
                              # the dimensionality reduction
                              # nudge_x = 5, nudge_y = 5,
                              segment.size = 0.5,
                              arrow = arrow(length = unit(0.01, 'npc')))
    
  } else {
    
    geom_text(data = centers,
              aes(x = mean_x, y = mean_y, label = Cluster),
              size = label_size)
    
  }
  
}



#' Get the limits of a the first two dimensions in a dimensionality reduction
#'
#' When plotting an embedding, we may want to plot specific cells, but
#' constrain the scale to match plots of the whole dataset. Given a dim.
#' reduction, this function extracts the x and y limits to use for plotting.
#'
#' @param seurat Seurat object for which a dimensionality reduction has been
#' computed (e.g. PCA or tSNE)
#' @param reduction String, corresponding to the dimensionality reduction to use.
#' Default: "tsne".
#' @param dim1 Numeric, dimension of embedding to use for x-axis. Default = 1.
#' @param dim2 Numeric, dimension of embedding to use for y-axis. Default = 2.
#'
#' @return A list with two elements: "xlim", which is a character vector of
#' the limits for the x-axis, and "ylim", correspondingly for the y-axis
#' @export
#' @author Selin Jessa
#'
#' @examples
#' get_dr_lims(pbmc, reduction = "tsne")
get_dr_lims <- function(seurat, reduction, dim1 = 1, dim2 = 2) {
  
  dim1 <- seurat[[reduction]]@cell.embeddings[, dim1]
  dim2 <- seurat[[reduction]]@cell.embeddings[, dim2]
  
  return(list(xlim = c(min(dim1), max(dim1)),
              ylim = c(min(dim2), max(dim2))))
  
}


#' Constrain the scale of the plot to the dimensionality reduction limits
#'
#' @inheritParams get_dr_lims
#'
#' @export
#' @author Selin Jessa
constrain_scale <- function(seurat, reduction, dim1 = 1, dim2 = 2)  {
  
  limits <- get_dr_lims(seurat = seurat, reduction = reduction, dim1 = dim1, dim2 = dim2)
  lims(x = limits$xlim, y = limits$ylim)
  
}


#' Apply a clean theme to a ggplot2 object
#'
#' @references https://github.com/sjessa/ggmin
#'
#' @importFrom ggplot2 theme_light theme
#' @author Selin Jessa
#' @export
theme_min <- function(base_size = 11, base_family = "",
                      border_color = "grey90",
                      border_size = 1,
                      gridlines = FALSE) {
  
  my_theme <- theme_light(base_size = base_size, base_family = "") +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = border_color, size = border_size),
      axis.ticks = element_line(color = border_color),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text.x = element_text(color = "black", size = rel(1.2)),
      strip.text.y = element_text(color = "black", size = rel(1.2)),
      title = element_text(size = rel(0.9)),
      axis.text = element_text(color = "black", size = rel(0.8)),
      axis.title = element_text(color = "black", size = rel(0.9)),
      legend.title = element_text(color = "black", size = rel(0.9)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), color = "black"),
      legend.key = element_rect(color = NA, fill = NA),
      legend.background = element_rect(color = NA, fill = NA)
    )
  
  if (!gridlines) {
    
    my_theme + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
    
  } else {
    
    return(my_theme)
    
  }
  
}


#' theme_powerpoint
#'
#' A clean theme for ggplot2, with some handy defaults for Powerpoints:
#' essentially \code{\link{theme_min}} with larger text in certain elements,
#' and a little room between axis titles and the plot.
#'
#' @export
theme_powerpoint <- function(base_size = 13, base_family = "") {
  
  theme_light(base_size = base_size, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey90", size = 1),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black", size = rel(1.2), face = "bold"),
      strip.text.y = element_text(colour = "black", size = rel(1.2), face = "bold"),
      title = element_text(size = rel(1.52)),
      axis.text = element_text(colour = "black", size = rel(1.1)),
      axis.title = element_text(colour = "black", size = rel(1.2)),
      axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
      legend.title = element_text(colour = "black", size = rel(1.1)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(1.0), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}


# get top and bottom values
get_top_and_bottom <- function(df, wt, n = 10) {
  
  wt_var <- rlang::enquo(wt)
  
  df %>%
    {bind_rows(top_n(., n, !!wt_var),
               top_n(., -n, !!wt_var))}
  
}


# a distance function for Spearman rank correlation
# from https://davetang.org/muse/2010/11/26/hierarchical-clustering-with-p-values-using-spearman/
spearman <- function(x, ...) {
  x <- as.matrix(x)
  res <- as.dist(1 - cor(x, method = "spearman", use = "everything"))
  res <- as.dist(res)
  attr(res, "method") <- "spearman"
  return(res)
}



#' Get default ggplot2/Seurat colors
#'
#' Get evenly spaced colors from around the color wheel, which are the default
#' colors assigned to clusters by Seurat. The output of this function can be
#' passed to the \code{scale_color_manual()} and \code{scale_fill_manual()} functions
#' from ggplot2, as the \code{values} argument. (\code{\link{ggColors}} points
#' to this function.)
#'
#' @param n Number of colors to return
#'
#' @return Named character vector, where names are the names of clusters, from
#' 0 to n-1, and values are the hex codes for the colors.
#' @export
#'
#' @examples
#'
#' n_clust <- 5
#' get_gg_colors(n_clust)
#'
#' @references https://stackoverflow.com/a/8197703
#' @aliases get_gg_colours
#' @importFrom grDevices hcl
get_gg_colors <- function(n) {
  
  hues <- seq(15, 375, length = n + 1)
  colors <- hcl(h = hues, l = 65, c = 100)[1:n]
  names(colors) <- seq(0, n - 1) # Since the first cluster in Seurat is 0
  
  return(colors)
  
}


square <- function() {
  theme(aspect.ratio = 1)
}



# Analysis-specific ------------------------------------------------------------



#' @param df Data frame, containing column with the granular cluster label
#' that we want to summarize at a higher cell class level
#' @param cluster_col Character, name of column with granular label
summarize_cell_types <- function(df, cluster_col) {
  
  cc_quo <- rlang::sym(cluster_col)
  
  df %>%
    mutate(
      # Define some broader cell type classes
      Type = case_when(
        grepl("RG|[Rr]adial|NSC|prog", !!cc_quo) ~ "RGC",
        grepl("GLIP|OPAS|Glial", !!cc_quo) ~ "Glial progenitors",
        grepl("EXIP|INIP|NEURP|IP", !!cc_quo) & !grepl("VIP", !!cc_quo) ~ "Neuronal progenitors",
        grepl("MGE|CGE|LGE|SST|PV|INH|CEX|PEX|GABAN|NEU|SPN|NRGN|UBCN|MFN|SERN", !!cc_quo) ~ "Neurons",
        grepl("OPC-P", !!cc_quo) ~ "Proliferating OPC",
        grepl("OPC", !!cc_quo) ~ "OPC",
        grepl("NFOL|MOL", !!cc_quo) ~ "Oligodendrocytes",
        grepl("ASTR", !!cc_quo) ~ "Astrocytes",
        grepl("EPEN|ASEP", !!cc_quo) ~ "Ependymal",
        grepl("MGL|MAC", !!cc_quo) ~ "Immune",
        grepl("T|Fibr|Mur|Micro|Mixed|UNR|Other|ENDO|PERI|MNG", !!cc_quo) ~ "Vascular & other",
        TRUE ~ "Vascular & other")) %>%
    mutate(Type = factor(Type, levels = c("RGC",
                                          "Glial progenitors",
                                          "Proliferating OPC",
                                          "OPC",
                                          "Oligodendrocytes",
                                          "Astrocytes",
                                          "Ependymal",
                                          "Neuronal progenitors",
                                          "Neurons",
                                          "Immune",
                                          "Vascular & other",
                                          "Normal")))
  
}

#' @param df Data frame, containing column with the granular cluster label
#' that we want to summarize at a higher cell class level
#' @param cluster_col Character, name of column with granular label
summarize_cell_types.Seurat <- function(seurat, cluster_col) {
  
  cc_quo <- rlang::sym(cluster_col)
  
  types <- seurat@meta.data %>%
    mutate(
      # Define some broader cell type classes
      Type = case_when(
        grepl("RG|[Rr]adial|NSC|prog", !!cc_quo) ~ "RGC",
        grepl("GLIP|OPAS", !!cc_quo) ~ "Glial progenitors",
        grepl("EXIP|INIP|NEURP|IP", !!cc_quo) & !grepl("VIP", !!cc_quo) ~ "Neuronal progenitors",
        grepl("MGE|CGE|LGE|SST|PV|INH|CEX|PEX|GABAN|NEU|SPN|NRGN|UBCN|MFN|SERN", !!cc_quo) ~ "Neurons",
        grepl("OPC-P", !!cc_quo) ~ "Proliferating OPC",
        grepl("OPC", !!cc_quo) ~ "OPC",
        grepl("NFOL|MOL", !!cc_quo) ~ "Oligodendrocytes",
        grepl("ASTR", !!cc_quo) ~ "Astrocytes",
        grepl("EPEN|ASEP", !!cc_quo) ~ "Ependymal",
        grepl("MGL|MAC", !!cc_quo) ~ "Immune",
        grepl("T|Fibr|Mur|Micro|Mixed|UNR|Other|ENDO|PERI|MNG", !!cc_quo) ~ "Vascular & other",
        TRUE ~ "Vascular & other")) %>%
    mutate(Type = factor(Type, levels = c("RGC",
                                          "Glial progenitors",
                                          "Proliferating OPC",
                                          "OPC",
                                          "Oligodendrocytes",
                                          "Astrocytes",
                                          "Ependymal",
                                          "Neuronal progenitors",
                                          "Neurons",
                                          "Immune",
                                          "Vascular & other",
                                          "Normal"))) %>% 
    .$Type
  
  seurat$Type <- types
  return(seurat)
  
}


compute_pct1_per_sample <- function(seurat, genes) {
  
  expr_data <- FetchData(seurat,
                         vars = c(genes, "Sample",
                                  "COR_ref.joint_mouse_extended",
                                  "Malignant_normal_consensus",
                                  "Molecular", "Location", "GrowthFactorReceptor")) %>%
    summarize_cell_types("COR_ref.joint_mouse_extended") %>%
    select(Sample = Sample, Cell_type = Type, everything())
  
  # get # of cells per type and only keep ones which are more than 1% of the total # of cells
  celltype_freq <- (table(expr_data$Cell_type) / nrow(expr_data)) %>%
    data.frame() %>%
    filter(Freq > 0.01)
  
  pct1_joint_sample <- expr_data %>%
    group_by(Molecular, Sample, Cell_type) %>%
    # calculate the number of cells of each cell type per sample
    mutate(n_total = n()) %>%
    ungroup() %>%
    gather(Gene, Expression, all_of(genes)) %>%
    group_by(n_total, Molecular, GrowthFactorReceptor, Location, Sample, Cell_type, Gene) %>%
    # count the number of cells of each cell type where each gene is detected
    summarize(N_detected = sum(Expression > 0)) %>%
    # divide by the total # of cells in the cell tpye
    mutate(pct1 = round(N_detected / n_total, 2)) %>%
    ungroup() %>%
    select(-n_total) %>%
    filter(Cell_type != "Other" & Cell_type %in% celltype_freq$Var1)
  
  return(pct1_joint_sample)
  
}



#' boxplot, coloured by projected cell type
plot_pct_boxplot_per_sample <- function(seurat, genes, filter_quo = NULL, title = "", return_df = FALSE) {
  
  pct1_joint_sample <- compute_pct1_per_sample(seurat, genes)
  
  if (!is.null(filter_quo)) pct1_joint_sample <- pct1_joint_sample %>% filter(!! filter_quo)
  
  if (return_df) return(pct1_per_sample)
  
  pct1_joint_sample %>%
    mutate(Gene = factor(Gene, levels = genes)) %>% 
    mutate(Cell_type = factor(Cell_type, levels = rev(names(palette_type)))) %>% 
    # complete(Gene, nesting(Cell_type), fill = list(pct1 = 0)) %>%
    ggplot(aes(x = Cell_type, y = pct1)) +
    geom_boxplot(aes(fill = Cell_type, alpha = 0.5), lwd = 0.25) +
    geom_point(stat = "identity", aes(fill = Cell_type, size = N_detected), shape = 21) +
    scale_color_manual(values = palette_type, guide = FALSE) +
    scale_fill_manual(values = palette_type, guide = FALSE) +
    scale_alpha(guide = FALSE) +
    facet_wrap(~ Gene, ncol = 1) +
    coord_flip() +
    rotate_x() +
    # no_legend() +
    ylab("Proportion of cells in which the gene is detected") +
    ylim(c(0, 1)) +
    ggtitle(title) +
    theme(legend.position = "bottom")
  
}





#' Get genes highly correlated with a gene of interest
#'
#' @param seurat Seurat object
#' @param gene Gene, should be present (detected) in the Seurat object
#' @param method Character, what method to use to compute correlations, passed to the 
#' \code{method} argument of \code{cor()}
#'
#' @return Data frame with two columns, "Gene" and "Correlation". The first entry
#' will necessarily be the input gene, with correlation of 1.
#' 
#' @examples
#' get_highly_correlated_genes(seurat, "TOP2A")
get_highly_correlated_genes <- function(seurat, gene, method = "spearman") {
  
  # check for presence
  if (!(gene %in% rownames(seurat))) stop("Gene not detected")
  
  # get expression data
  expr2 <- seurat@assays$RNA@data
  m2 <- as.matrix(expr2)
  g2 <- m2[gene, ]
  
  # correlate input gene with all other genes
  correlations2 <- apply(m2, 1, function(x) cor(g2, x, method = method))
  
  out <- sort(correlations2, decreasing = TRUE) %>% 
    enframe("Gene", "Correlation")
  
  return(out)
  
}


prep_celltype_specificity <- function(genes) {
  
  feather::read_feather(pct1_feather, c("ID_20201028", unique(genes))) %>% 
    rowwise() %>% 
    dplyr::rename(Cluster = ID_20201028) %>% 
    filter(!grepl("EXCLUDE", Cluster) & Cluster %in% names(palette_joint_mouse_extended)) %>% 
    mutate(Timepoint = stringr::str_split_fixed(Cluster, "_", 2) %>% getElement(1)) %>% 
    relocate(Timepoint, .after = 1) %>% 
    gather(Gene, Pct1, 3:ncol(.))
  
}




#' Function for computing a celltype specificity score for individual genes
#' which is defined as Pct1 (detection rate in the highest-expressing cluster)
#' - Pct2 (detection rate in all other clusters)
#' 
#' @param df Data frame, with at least four columns: Timepoint, Cluster, Gene, Pct1.
#' Create a dummy variable called Timepoint if the data is not segregated by age.
#' @param n_cells_df Data frame with two columns, Cluster, matching values in 
#' \code{df$cluster}, and N_cells, containing the number of cells in each cluster. This
#' is used to compute the Pct2 value.
#' @param max_across_timepoints Logical, whether to return one score per gene
#' by taking the maximum specificity of the gene across all the timepoints. Default: TRUE.
#' If FALSE, there will be a specificity score for each gene for each timepoint/sample in the
#' dataset.
compute_celltype_specificity <- function(df, n_cells_df, max_across_timepoints = TRUE) {
  
  specificity_score <- function(df_per_timepoint) {
    
    if (sum(df_per_timepoint$Pct1, na.rm = TRUE) == 0 |
        all(is.na(df_per_timepoint$Pct1))) return(0)
    else {
      
      # find the cluster with the highest detection rate
      top_cluster <- df_per_timepoint %>% top_n(1, Pct1) %>% dplyr::slice(1)
      
      other_clusters <- df_per_timepoint %>%
        filter(Cluster != top_cluster$Cluster & !is.na(Pct1)) %>%
        left_join(n_cells_df, by = "Cluster") %>%
        # recover the # of cells where the gene is detected in each other cluster
        mutate(N_cells_exp = Pct1 * N_cells) %>%
        # aggregate that value to get the detection rate across all other clusters
        summarize(Pct2 = sum(N_cells_exp) / sum(N_cells))
      
      # return the difference
      return(round(top_cluster$Pct1 - other_clusters$Pct2, 2))    
      
    }
    
  }
  
  specificity_per_gene <- function(df, gene) {
    
    scores <- df %>% 
      filter(Gene == gene) %>% 
      group_by(Timepoint) %>% 
      summarize(Score = specificity_score(cur_data())) %>% 
      dplyr::select(Timepoint, Score)
    
    tops <- df %>% 
      filter(Gene == gene) %>% 
      group_by(Timepoint) %>% 
      arrange(desc(Pct1)) %>% 
      dplyr::slice(1) %>% 
      ungroup() %>% 
      dplyr::select(Top = Cluster)
    
    bind_cols(scores, tops) %>% 
      tibble::add_column(Gene = gene, .before = 1)
    
  }
  
  specificity_out <- map_dfr(unique(df$Gene), ~ specificity_per_gene(df, .x))
  
  # summarize cell types
  specificity_out <- summarize_cell_types(specificity_out, "Top")
  
  # taking max specificity across timepoints for each gene, to return one score per gene
  if (max_across_timepoints) specificity_out %>% 
    group_by(Gene) %>% 
    dplyr::slice(which.max(Score))
  else return(specificity_out)
  
}



calc_pct1 <- function(df, cluster_col) {
  
  require(data.table)
  
  clusters <- as.character(df[[cluster_col]])
  df[[cluster_col]] <- NULL
  
  x <- as.matrix(df)
  # Binarize
  x[x > 0] <- 1
  x <- as.data.table(x)
  # Add cluster info for cells
  x <- x[, Cluster := clusters]
  # Get prop of cells expressing a gene, within each cluster
  x <- x[, lapply(.SD, function(i) sum(i)/length(i)), by = Cluster]
  
  x <- as.data.frame(x)
  colnames(x)[which(colnames(x) == "Cluster")] <- cluster_col
  return(x)
  
}

plot_bubble <- function(genes, meanexp_feather, pct_feather, cluster_col, mean_exp, dendrogram_order) {
  
  rename_vars <- c(Cluster = cluster_col)
  
  x_mean <- feather::read_feather(meanexp_feather, c(cluster_col, genes)) %>%
    rename(!!rename_vars) %>% 
    tibble::column_to_rownames(var = "Cluster") %>% 
    apply(2, scales::rescale) %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "Cluster") %>%
    gather(Gene, Expression, 2:ncol(.))
  
  x_pct1 <- feather::read_feather(pct_feather, c(cluster_col, genes)) %>% 
    rename(!!rename_vars) %>% 
    gather(Gene, Pct1, 2:ncol(.))
  
  bubbleplot_data <- left_join(x_mean, x_pct1, by = c("Cluster", "Gene")) %>%
    filter(Cluster %in% dendrogram_order) %>%
    filter(Pct1 > 0) %>% 
    mutate(Cluster = factor(Cluster, levels = dendrogram_order)) %>%
    mutate(Gene = factor(Gene, levels = rev(genes))) %>% 
    replace_na(list(Expression = 0, Pct1 = 0))
  
  message("All clusters present: ", length(unique(bubbleplot_data$Cluster)) == length(dendrogram_order))
  
  bubbleplot_data %>%
    ggplot(aes(x = Cluster, y = Gene)) +
    geom_point(aes(size = Pct1, colour = Expression), alpha = 0.8) +
    scale_radius() +
    scale_color_gradientn(colours = tail(rdbu, 70)) +
    theme_min() +
    rotate_x() +
    theme(panel.grid.major.x = element_line(colour = "grey90"),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 13))
  
}