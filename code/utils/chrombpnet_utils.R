# Selin Jessa
# Helpers for working with motifs and ChromBPNet outputs

library(ggplot2)


#' Load hit calls from finemo-gpu
#' 
#' @param hits_path character, path to TSV containing hit calls with coordinates
#' @param one_index_peaks logical, whether to increment the peak ID in the hits
#' to change them from being 0-indexed to 1-indexed (the latter matching the indexing
#' used in R)
load_hits <- function(hits_path, one_index_peaks = TRUE) {
  
  hits <- data.table::fread(hits_path, data.table = FALSE)
  
  # convert peak id to 1-based indexing for working with peaks in R
  if (one_index_peaks) {
    message("@ converting peak IDs to 1-indexing")
    hits <- hits %>% mutate(peak_id = peak_id + 1)
  }
  
  return(hits)
  
}


#' Read a (possibly-zipped) BED file into a GRanges object.
read_bed <- function(bed, col_names) {
  
  df <- data.table::fread(bed, data.table = FALSE, header = FALSE)
  colnames(df) <- col_names
  
  # convert 0-based to 1-based (done automatically by rtracklayer)
  df$start <- df$start + 1
  gr <- GRanges(df)
  
  return(gr)
  
}


#' Write a GRanges object to a BED file
write_bed <- function(gr, out_bed)  {
  
  # convert to BED - flatten to GRanges and then write out
  # https://www.biostars.org/p/89341/#89351
  
  df <- data.frame(seqnames = seqnames(gr),
                   starts   = start(gr)-1, # switch to 0-based coords for BED
                   ends     = end(gr),
                   names    = mcols(gr)$motif_name,
                   scores   = gr$score,
                   strands  = strand(gr))
  
  write.table(df, file = out_bed, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
}



#' Helper function to sum up the number of hits per motif in a set of hits, with
#' some optional filtering options
#' 
#' @param pattern_class character, one of "pos_patterns" or "neg_patterns"
count_hits <- function(hits, pattern_class = NULL, peak_ids = NULL, motif_order = NULL) {
  
  # filter to a certain pattern class
  if (!is.null(pattern_class)) hits <- hits[hits$pattern_class == pattern_class, ]
  
  # filter to certain peak IDs
  if (!is.null(peak_ids)) hits <- hits %>% filter(peak_id %in% peak_ids)
  
  hits <- hits %>% 
    group_by(motif_name) %>% 
    count() %>% 
    ungroup() %>% 
    arrange(desc(n))
  
  # put the motifs in order and complete missing values (i.e. motifs with no hits) with 0s
  if (!is.null(motif_order)) {
    
    # filter motifs from hits not in the provided order
    motifs_dropped <- setdiff(hits$motif_name, motif_order)
    message("@ dropping ", glue_collapse(motifs_dropped, ", "))
    
    hits <- hits %>%
      filter(motif_name %in% motif_order) %>% 
      mutate(motif_name = factor(motif_name, levels = motif_order)) %>% 
      tidyr::complete(motif_name, fill = list("n" = 0)) %>% 
      arrange(motif_name)
    
    if (!all(hits$motif_name == motif_order)) stop("@ motif order does not match.")
    
  }
  
  return(hits)
  
}






#' Read CWMs from MEME format
#' 
#' Adapted from universalmotif::read_meme
#' https://github.com/bjmt/universalmotif/blob/master/R/read_meme.R
#' in order to handle reading in CWMs from MEME format, which don't adhere to 
#' the expectations of PPMs/PWMs.
#' 
#' @param memedb character, path to CWMs in meme db format (e.g. as extracted by
#' \code{modisco meme -t CWM})
read_cwm_meme <- function(memedb, skip=0) {
  
  raw_lines <- readLines(con <- file(memedb))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]
  raw_lines <- raw_lines[!grepl("^#", raw_lines)]
  raw_lines <- raw_lines[!grepl("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", raw_lines)]
  raw_lines <- raw_lines[!grepl("------------", raw_lines)]
  
  universalmotif:::check_meme_version(raw_lines)
  
  alph <- universalmotif:::get_meme_alph(raw_lines)
  alph.len <- universalmotif:::get_meme_alph_len(alph)
  alph.split <- switch(alph, DNA = Biostrings::DNA_BASES, RNA = Biostrings::RNA_BASES, AA = Biostrings::AA_STANDARD2,
                       safeExplode(alph))
  
  strands <- raw_lines[grepl("^strands:", raw_lines)]
  if (length(strands) > 0) {
    strands <- strsplit(strands, "\\s+")[[1]][-1]
  } else {
    message("Could not find strand info, assuming +.")
    strands <- "+"
  }
  if (all(c("+", "-") %in% strands)) {
    strands <- "+-"
  }
  bkg.start <- grep("^Background letter frequencies", raw_lines)
  if (length(bkg.start)) {
    bkg.offset <- 1
    bkg <- raw_lines[bkg.start + bkg.offset]
    bkg <- strsplit(bkg, "\\s+")[[1]]
    bkg <- bkg[bkg != ""]  # if the line if prepended with a space
    bkg <- as.numeric(bkg[seq_len(length(bkg)) %% 2 == 0])
    while (length(bkg) < alph.len) {
      bkg.offset <- bkg.offset + 1
      bkg.tmp <- raw_lines[bkg.start + bkg.offset]
      bkg.tmp <- strsplit(bkg.tmp, "\\s+")[[1]]
      bkg.tmp <- bkg.tmp[bkg.tmp != ""]
      bkg.tmp <- as.numeric(bkg.tmp[seq_along(bkg.tmp) %% 2 == 0])
      bkg <- c(bkg, bkg.tmp)
    }
    if (anyNA(bkg))
      stop("Could not parse background frequencies, check that they match alphabet")
  } else {
    message("Could not find background, assuming uniform frequencies.")
    bkg <- rep(1 / length(alph.split), length(alph.split))
  }
  
  motif_meta <- grep("^letter-probability matrix:", raw_lines)
  motif_names_i <- grep("^MOTIF ", raw_lines)
  motif_names <- lapply(raw_lines[motif_names_i], function(x) {
    x <- strsplit(x, "\\s+")[[1]]
    if (x[1] == "") x[3] else x[2]
  })
  motif_altnames <- lapply(raw_lines[motif_names_i], function(x) {
    x <- strsplit(x, "\\s+")[[1]]
    if (x[1] == "") x[4] else x[3]
  })
  motif_starts <- motif_meta + 1
  motif_stops <- sapply(raw_lines[motif_meta],
                        function(x) strsplit(x, "\\s+")[[1]][6])
  motif_stops <- motif_meta + as.numeric(motif_stops)
  
  motif_meta <- lapply(raw_lines[motif_meta],
                       function(x) {
                         x <- strsplit(x, "\\s+")[[1]]
                         c(nsites = x[8], eval = x[10])
                       })
  motif_list <- mapply(function(x, y) {
    z <- trimws(raw_lines[x:y])
    z <- sapply(z, function(x) strsplit(x, "\\s+")[[1]])
    if (nrow(z) != alph.len)
      stop("Alphabet length does not match motif length")
    z <- z[order(alph.split, method = "radix"), ]
    as.numeric(z)
  }, motif_starts, motif_stops, SIMPLIFY = FALSE)
  
  motifs <- purrr::map(motif_list, ~ t(matrix(.x, ncol = alph.len, byrow = TRUE))) %>% 
    set_names(unlist(motif_names))
  
  return(motifs)
  
}




#' Trim CWM based on a threshold
#' 
#' adapted from Jacob Schreiber in tf-modiscolite
#' https://github.com/jmschrei/tfmodisco-lite/blob/570535ee5ccf43d670e898d92d63af43d68c38c5/modiscolite/report.py#L213-L236
#' we use a slightly larger window than the actual trimmed CWM, as in:
#' https://github.com/jmschrei/tfmodisco-lite/blob/d98aeb17a9c79398ded3e21d0a30d004f39fdcd8/modiscolite/report.py#L233
#' 
#' @param cwm numeric matrix of 4 x N. 4 rows are bases and N columns are positions.
trim_cwm <- function(cwm, trim_threshold = 0.3, trim_min_length = 3, flank = 4) {
  
  score <- colSums(abs(cwm))
  trim_thresh <- max(score) * trim_threshold
  pass_inds <- which(score >= trim_thresh)
  # print(pass_inds)
  
  # add 4bp of context
  start <- max(min(pass_inds) - flank, 1)
  end <- min(max(pass_inds) + flank, length(score))
  
  trimmed <- cwm[, start:end]
  
  return(trimmed)
  
}



#' Function to get the highest affinity sequence from a CWM
#' by taking the base with highest absolute contributions at each position.
#'
#' @param cwm  numeric matrix of 4 x N. 4 rows are bases and N columns are positions.
get_high_affinity_seq <- function(cwm) {
  
  # define the bases
  bases <- c("A", "C", "G", "T")
  
  # get the highest affinity sequence
  paste0(bases[apply(cwm, MARGIN = 2, which.max)], collapse = "")
  
}



#' Reverse complement a motif represented as a matrix
revcomp <- function(motif) {
  
  # first reverse columns, then reverse rows, then transpose
  motif <- t(apply(apply(motif, 2, rev), 1, rev))
  
  # label columns and rows, necessary for some plotting & to avoid ambiguity
  rownames(motif) <- c("A", "C", "G", "T")
  colnames(motif) <- 1:ncol(motif)
  
  return(motif)
  
}




#' Plot a motif logo, optionally trimming
#' 
#' @param motif 4xN matrix, with rownames labelled as A, C, G, T
#' @param method character, passed one of "custom", "bits", "probability", passed
#' to ggseqlogo as plotting method. For CWMs, use "custom". Default: "custom".
#' @param revcomp logical, whether to reverse complement the mtoif before plotting
#' @param trim logical, whether to trim a CWM for a motif before plotting based
#' on a threshold, see \code{trim_cwm}
logo <- function(motif, trim = FALSE, revcomp = FALSE, method = "custom", ...) {
  
  if (trim) motif <- trim_cwm(motif)
  if (revcomp) motif <- revcomp(motif)
  
  ggseqlogo::ggseqlogo(motif, method = method, font = "roboto_bold", ...)
  
}




#' Plot several motifs in a vertical stack, with some helpful theme options
#' 
#' @param motifs list of 4xN matrices, with rownames labelled as A, C, G, T
#' @param method character, passed one of "custom", "bits", "probability", passed
#' to ggseqlogo as plotting method. For CWMs, use "custom". Default: "custom".
stack_logo <- function(motifs, method = "custom", title = NULL, hide_ticks = TRUE, ...) {
  
  gg <- ggplot() +
    geom_logo(motifs, method = method, font = "roboto_bold", ...) +
    theme_logo() + 
    facet_grid(seq_group ~ ., scales = "free_y", switch = "both") +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1),
          # reduce margins between facets, so that the logos can take up
          # a little bit more space
          panel.spacing=unit(0, "lines")
    )
  
  if (hide_ticks) gg <- gg + hide_ticks() + ylab(NULL)
  
  if (!is.null(title)) gg <- gg + ggtitle(title)
  
  return(gg)
  
}



#' Adapted from Signac::FindMotifs:
#' https://github.com/stuart-lab/signac/blob/8b98f63f89fb988fcfd5e1be6af6de316353008c/R/motifs.R#L316
#' 
#' Computes the number of hits containing the motif (observed) and
#' compares this to the total number of hits containing the
#' motif (background) using the hypergeometric test.
#' 
#' @param query_mat numeric matrix with N rows and 1 column. Rows are motifs and the column
#' contains hit counts in a region set of interest (e.g. for peaks from a given cell type).
#' containing counts of each motif in each cell type.
#' @param background_mat numeric matrix with N rows and M columns. Rows are motifs
#' and columns are different region sets (e.g. counts of motifs per cell type).
compute_motif_enrichment <- function(query_mat, background_mat, background_cols_drop = NULL) {
  
  print(dim(query_mat))
  print(dim(background_mat))
  
  if (!is.null(background_cols_drop)) {
    
    background_cols <- setdiff(colnames(background_mat), background_cols_drop)
    
  } else {
    
    background_cols <- colnames(background_mat)
    
  }
  
  # extract counts for the query set and background
  query_counts <- rowSums(query_mat)
  
  # add pseuodocount to avoid denominator of 0 in LFC
  background_counts <- rowSums(background_mat[, background_cols, drop = FALSE]) + 1
  
  # total number of sequences in cell type and background
  n_query <- sum(query_counts)
  n_background <- sum(background_counts)
  
  # compute percentages
  percent_query <- query_counts / n_query * 100
  percent_background <- background_counts / n_background * 100
  fold_enrichment <- percent_query / percent_background
  
  motifs_observed <- names(query_counts)[query_counts > 0]
  idx_motifs_observed <- which(names(query_counts) %in% motifs_observed)
  
  # run hypergeometric test
  p_list <- sapply(seq_along(motifs_observed), function(i) {
    motif <- motifs_observed[i]
    p <- phyper(
      q = query_counts[[motif]] - 1,
      m = background_counts[[motif]],
      n = sum(background_counts) - background_counts[[motif]],
      k = sum(query_counts),
      lower.tail = FALSE
    )
    return(p)
  })
  
  # adjust p-values for multiple testing
  p_adjusted <- p.adjust(p = p_list, method = "BH")
  
  # compile results
  results <- data.frame(
    motif_name = motifs_observed,
    observed = query_counts[motifs_observed],
    background = background_counts[motifs_observed],
    percent.observed = percent_query[motifs_observed],
    percent.background = percent_background[motifs_observed],
    fold.enrichment = fold_enrichment[motifs_observed],
    pvalue = p_list,
    p.adjust = p_adjusted,
    stringsAsFactors = FALSE
  ) %>%
    # sort by adjusted p-value and then descending fold enrichment
    arrange(p.adjust, desc(fold.enrichment))
  
  rownames(results) <- NULL
  
  return(results)
  
}



#' Import broadPeak BED file 
#' https://genome.ucsc.edu/FAQ/FAQformat.html#format13
import_broadpeak <- function(bed_path) {
  
  rtracklayer::import.bed(bed_path, extraCols = c("signalValue" = "numeric", "pValue" = "numeric", "qValue" = "numeric"))
  
}



#' Import narrowPeakBED file
#' https://genome.ucsc.edu/FAQ/FAQformat.html#format12
import_narrowpeak <- function(bed_path) {
  
  rtracklayer::import.bed(bed_path, extraCols = c("signalValue" = "numeric", "pValue" = "numeric", "qValue" = "numeric", "peak" = "numeric"))
  
}





#' Calculate the entropy given counts of hits per organ
#' 
#' We expect *lower* entropy when the uncertainty is low, meaning the
#' organ specificity is high and the distribution over organs is skewed.
#' 
#' We expect a *higher* entropy when the uncertainty is high, because
#' the distribution is near uniform.
#'
#' @param organ_counts numeric vector, containing counts of hits per organ.
calculate_entropy <- function(organ_counts) {
  
  # calculate total count
  total_count <- sum(organ_counts)
  
  # calculate probabilities
  probabilities <- organ_counts / total_count
  
  # calculate entropy
  entropy <- -sum(probabilities * log2(probabilities), na.rm = TRUE)
  
  return(entropy)
}




#' Plot a summary of multiple motifs and their hit call stats
#' 
#' - trimmed CWM of a representative motif per broad annotation
#' - known PWM for best matching known-motif
#' - total hits
#' - number of motif variants
#' - breakdown of hits in terms of organs
#' - breakdown of hits in terms of cell type compartments
#' 
#' @param motifs_filt data frame
#' @param hits_annotated data frame
#' @param to_revcomp character vector, specifying names of motifs to reverse
#' complement prior to plotting. Must match names in \code{motifs_filt$motif_name}
#' @param plot_known_motifs logical, whether to plot a stack of logos for known PWMs,
#' using the best match provided in \code{motifs_filt$motif0}
#' @param top_n numeric, whether to plot only the top N motifs (based on total hits),
#' to reduce plot size
#' @param rel_widths numeric vector, passed to cowplot::plot_grid to control 
#' relative widths of each plot to override defaults
#' @param order_by character, one of "hits" or "entropy"
plot_motif_summary <- function(motifs_filt,
                               hits_annotated,
                               hit_dyad_dist,
                               to_revcomp = NULL,
                               plot_known_motifs = TRUE,
                               top_n = NULL,
                               rel_widths = NULL,
                               known_motifs = vierstra_motifs,
                               known_motif_col = "motif0",
                               known_motif_name_col = "TF0",
                               order_by = "entropy") {
  
  # group by broad annotation and calculate summary stats
  motifs_to_plot <- motifs_filt %>% 
    group_by(annotation_broad) %>% 
    mutate(total_hits_broad = sum(total_hits),
           n_variants = n()) %>%
    # for visualization of the logo, use the motif w/ the most hits
    slice_max(n = 1, order_by = total_hits) %>% 
    arrange(desc(total_hits_broad))
  
  if (!all(motifs_to_plot$motif_name %in% names(cwm_list))) stop("Some motifs not found in CWM list")
  if (!all(to_revcomp %in% unique(motifs_to_plot$annotation_broad))) stop("Some motifs to revcomp not found in broad categories")
  
  if (!is.null(top_n)) motifs_to_plot <- motifs_to_plot[1:top_n, ]
  
  message("@ plotting ", nrow(motifs_to_plot), " motifs")
  
  # get an ordering
  if (order_by == "entropy") {
    
    # calculate number of hits per motif per organ
    hits_per_organ <- hits_annotated %>%
      filter(annotation_broad %in% motifs_to_plot$annotation_broad) %>% 
      dplyr::select(annotation_broad, organ, n) %>%
      group_by(annotation_broad, organ) %>%
      summarize(n_hits = sum(n))
    
    # sort by decreasing entropy
    motif_order <- split(hits_per_organ, hits_per_organ$annotation_broad) %>%
      map(~ .x %>% pull(n_hits)) %>%
      map_dbl(~ calculate_entropy(.x)) %>%
      sort(decreasing = TRUE) %>%
      names()
    
    # reorder data frame based on entropy
    motifs_to_plot <- motifs_to_plot %>% 
      mutate(annotation_broad = factor(annotation_broad, levels = motif_order)) %>% 
      arrange(annotation_broad)
    
  } else if (order_by == "distToTSS") {
    
    motif_order <- hits_annotated %>%
      filter(annotation_broad %in% motifs_to_plot$annotation_broad) %>%
      group_by(pattern_class, annotation_broad) %>%
      # sort by class then distToTSS
      mutate(pattern_class = factor(pattern_class, levels = c("pos_patterns", "neg_patterns"))) %>% 
      summarize(mean_dist = mean(median_distToTSS)) %>% 
      arrange(pattern_class, mean_dist) %>%
      pull(annotation_broad)
    
    # reorder data frame based on distance to TSS
    motifs_to_plot <- motifs_to_plot %>% 
      mutate(annotation_broad = factor(annotation_broad, levels = motif_order)) %>% 
      arrange(annotation_broad)
    
  } else if (order_by == "hits") {
    
    # use the order by total hits
    motif_order <- motifs_to_plot$annotation_broad
    
  }
  
  # get motifs in the new order
  cwm_list_subset <- cwm_list[motifs_to_plot$motif_name] %>% 
    map(trim_cwm, flank = 2)
  
  # use the annotation_broad name instead of the motif name
  # base_representative$motif_name_short <- gsub("^[^|]*\\|([^#]*)#.*$", "\\1", base_representative$motif_name)
  names(cwm_list_subset) <- plyr::mapvalues(names(cwm_list_subset),
                                            from = motifs_to_plot$motif_name,
                                            to = as.character(motifs_to_plot$annotation_broad))
  
  # reverse complement to match known motifs if needed
  if (!is.null(to_revcomp)) {
    cwm_list_subset[to_revcomp] <- map(cwm_list_subset[to_revcomp], revcomp)
  }
  
  # make the logos
  p1 <- stack_logo(cwm_list_subset, method = "custom", title = "De novo CWM")
  
  # plot how many total hits for this annotation
  p3 <- motifs_to_plot %>%
    distinct(annotation_broad, total_hits_broad) %>% 
    mutate(annotation_broad = factor(annotation_broad, levels = rev(motif_order))) %>%
    ggplot(aes(y = annotation_broad, x = "1")) +
    # color scale is on log10
    geom_tile(aes(fill = log10(total_hits_broad)), alpha = 1) +
    # labels are non log-transformed
    geom_text(aes(label = scales::comma(total_hits_broad)), color = "white", size = 3) +
    # don't include the very lightest colors; too hard to read
    scale_fill_gradientn(colors = viridis::plasma(100)[0:95], limits = c(0, 8)) +
    ylab(NULL) + xlab(NULL) + ggtitle("total # \nhits") +
    rotate_x() +
    hide_ticks() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = "bottom")
  
  # plot how many motifs in this annotation
  p4 <- motifs_to_plot %>%
    distinct(annotation_broad, n_variants) %>%
    mutate(annotation_broad = factor(annotation_broad, levels = rev(motif_order))) %>%
    ggplot(aes(y = annotation_broad, x = n_variants)) +
    # geom_point(aes(size = n_variants)) +
    # scale_radius() +
    geom_col(fill = "gray70") +
    geom_text(aes(label = n_variants), hjust = -0.5) +
    ylab(NULL) + xlab(NULL) + ggtitle("# motifs") +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "bottom") +
    xlim(c(0, 25))
  
  # plot breakdown of total hits based on cell type compartment
  p5 <- hits_annotated %>%
    filter(annotation_broad %in% motif_order) %>% 
    dplyr::select(annotation_broad, compartment, n) %>%
    group_by(annotation_broad, compartment) %>%
    summarize(n_hits = sum(n)) %>%
    mutate(annotation_broad = factor(annotation_broad, levels = rev(motif_order))) %>%
    ggplot(aes(y = annotation_broad, x = n_hits)) +
    geom_col(aes(fill = compartment), position = "fill") +
    scale_fill_manual(values = cmap_compartment2) +
    ylab(NULL) + xlab(NULL) + ggtitle("% hits \nby compartment") +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "bottom") +
    ggplot2::guides(fill=guide_legend(ncol=2))
  
  # plot breakdown of total hits based on organ
  p6 <- hits_annotated %>%
    filter(annotation_broad %in% motif_order) %>% 
    dplyr::select(annotation_broad, organ, n) %>%
    group_by(annotation_broad, organ) %>%
    summarize(n_hits = sum(n)) %>%
    mutate(annotation_broad = factor(annotation_broad, levels = rev(motif_order))) %>%
    ggplot(aes(y = annotation_broad, x = n_hits)) +
    geom_col(aes(fill = organ), position = "fill") +
    scale_fill_manual(values = cmap_organ) +
    ylab(NULL) + xlab(NULL) + ggtitle("% hits \nby organ") +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "bottom") +
    ggplot2::guides(fill=guide_legend(ncol=2))
  
  # plot breakdown of total hits based on peak/region type
  p7 <- hits_annotated %>%
    filter(annotation_broad %in% motif_order) %>% 
    dplyr::select(annotation_broad, organ, Distal, Exonic, Intronic, Promoter) %>%
    pivot_longer(cols = c(Distal, Exonic, Intronic, Promoter), names_to = "region_type", values_to = "n_hits") %>%
    group_by(annotation_broad, region_type) %>% 
    summarize(n_hits_total = sum(n_hits, na.rm = TRUE)) %>% 
    mutate(
      region_type = factor(region_type, levels = names(cmap_peaktype2)),
      annotation_broad = factor(annotation_broad, levels = rev(motif_order))) %>%
    ggplot(aes(y = annotation_broad, x = n_hits_total)) +
    geom_col(aes(fill = region_type), position = "fill", alpha = 0.7) +
    scale_fill_manual(values = cmap_peaktype2) +
    ylab(NULL) + xlab(NULL) + ggtitle("% hits \nby type") +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "bottom") +
    ggplot2::guides(fill=guide_legend(ncol=2))
  
  # plot distribution of median distance to TSS
  p8 <- hits_annotated %>%
    filter(annotation_broad %in% motif_order) %>% 
    dplyr::select(annotation_broad, organ, median_distToTSS) %>%
    mutate(annotation_broad = factor(annotation_broad, levels = rev(motif_order))) %>%
    ggplot(aes(y = annotation_broad, x = median_distToTSS)) +
    geom_boxplot(fill = "gray70", outliers = FALSE) +
    ylab(NULL) + xlab(NULL) + ggtitle("median distance \nto TSS") +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "bottom") +
    # chop x-limits
    coord_cartesian(xlim = c(0, 25000)) +
    # annotate thousands as K
    scale_x_continuous(labels = scales::label_number(scale_cut = cut_short_scale()))
  
  # plot distribution of median distance to summit
  p9 <- hits_annotated %>%
    filter(annotation_broad %in% motif_order) %>% 
    dplyr::select(annotation_broad, organ, median_distToPeakSummit) %>%
    mutate(annotation_broad = factor(annotation_broad, levels = rev(motif_order))) %>%
    ggplot(aes(y = annotation_broad, x = median_distToPeakSummit)) +
    geom_boxplot(fill = "gray70", outliers = FALSE) +
    ylab(NULL) + xlab(NULL) + ggtitle("median distance \nto summit") +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "bottom") +
    # set x-limits
    coord_cartesian(xlim = c(0, 300))
  
  # plot the distribution of binned distances to nucleosome dyads
  hit_dyad_dist_agg <- hit_dyad_dist %>% 
    filter(annotation_broad %in% motif_order) %>% 
    group_by(annotation_broad, bin_dist) %>% 
    summarize(sum_counts = sum(n)) %>% 
    # zscore counts per motif (row-wise)
    mutate(zscore = scale(sum_counts)) %>% 
    mutate(annotation_broad = factor(annotation_broad, levels = rev(motif_order)),
           bin_dist = as.numeric(bin_dist)) %>%
    ungroup()
  
  p10 <- hit_dyad_dist_agg %>% 
    ggplot(aes(x = bin_dist, y = annotation_broad)) +
    geom_tile(aes(fill = zscore)) +
    scale_fill_gradientn(colours = rdbu2, rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
    # geom_vline(xintercept = 75, linewidth = 1) +
    scale_x_continuous(breaks = seq(10, 250, by = 10)) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10),
          legend.position = "bottom") +
    xlab("distance bin (bp)") +
    ylab(NULL) +
    ggtitle("binned distances to \nnucleosome dyad")
  
  # pattern class
  p11 <- motifs_to_plot %>%
    distinct(annotation_broad, pattern_class) %>%
    mutate(annotation_broad = factor(annotation_broad, levels = rev(motif_order))) %>%
    mutate(class_short = ifelse(pattern_class == "pos_patterns", "+", "-")) %>%
    ggplot(aes(x = 1, y = annotation_broad)) +
    geom_tile(aes(fill = class_short), alpha = 0.1) +
    geom_text(aes(label = class_short, color = class_short), fontface = "bold", size = 8) +
    scale_fill_manual(values = c("+" = "darkgreen", "-" = "darkred"), aesthetics = c("color", "fill")) +
    ggtitle("class") +
    no_legend() +
    ylab(NULL) +
    xlab(NULL) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  if (is.null(rel_widths)) {
    
    if (plot_known_motifs) rel_widths = c(2.5, 2, 0.5, 1, 1, 1, 1, 1, 1, 1, 1)
    else if (!plot_known_motifs) rel_widths = c(3, 0.5, 1, 1, 1, 1, 1, 1, 1, 1)
    
  }
  
  # plot known motifs as a stack
  if (plot_known_motifs) {
    
    known_motifs_subset <- known_motifs[motifs_to_plot[[known_motif_col]]]
    
    names(known_motifs_subset) <- plyr::mapvalues(names(known_motifs_subset),
                                                  from = motifs_to_plot[[known_motif_col]],
                                                  to = motifs_to_plot[[known_motif_name_col]])
    
    # if non-unique, assign a prefix to allow plotting
    if (length(unique(names(known_motifs_subset))) < length(names(known_motifs_subset))) {
      
      names(known_motifs_subset) <- paste0(seq_along(known_motifs_subset), " ", names(known_motifs_subset))
      
    }
    
    known_motifs_subset <- map(known_motifs_subset, function(cwm) {
      
      rownames(cwm) <- c("A", "C", "G", "T")
      
      return(cwm)
      
    })
    
    # make the logos
    p2 <- stack_logo(known_motifs_subset, method = "bits", title = "Known PWM")
    
    plot_grid(p1, p2, p11, p3, p7, p8, p9, p10, p5, p6, p4, nrow = 1, align = "h", axis = "tb", rel_widths = rel_widths)
    
  } else {
    
    plot_grid(p1, p11, p3, p7, p8, p9, p10, p5, p6, p4, nrow = 1, align = "h", axis = "tb", rel_widths = rel_widths)
    
  }
  
}


