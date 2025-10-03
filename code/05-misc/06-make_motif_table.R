# Purpose: make an interactive table that also renders motif logos for easy
# exploration of the de novo motif compendium

# libraries
library(here)
library(DT)
library(readr)
library(dplyr)
library(knitr)
library(htmlwidgets)

script_path <- here("code/utils/")
source(file.path(script_path, "plotting_config.R"))
source(file.path(script_path, "hdma_palettes.R"))
source(file.path(script_path, "sj_scRNAseq_helpers.R"))

# Data preparation ----------------------------------------------------------

motif_anno <- read_tsv(here("output/03-chrombpnet/03-syntax/01/motifs_compiled_unique.tsv"))

motif_dir <- "motifs"

motif_anno2 <- motif_anno %>%
  mutate(
    total_log10_hits = log10(total_hits),
    fwd_path = if_else(!is.na(motif_name_safe), file.path(motif_dir, paste0(motif_name_safe, ".cwm.fwd.png")), NA_character_),
    rev_path = if_else(!is.na(motif_name_safe), file.path(motif_dir, paste0(motif_name_safe, ".cwm.rev.png")), NA_character_)
  ) %>%
  rowwise() %>%
  # encode images
  mutate(
    cwm_fwd = if (is.na(fwd_path) || !file.exists(fwd_path)) {
      "N/A"
    } else {
      paste0('<img src="', image_uri(fwd_path), '" height="40"></img>')
    },
    
    cwm_rev = if (is.na(rev_path) || !file.exists(rev_path)) {
      "N/A"
    } else {
      paste0('<img src="', image_uri(rev_path), '" height="40"></img>')
    }
  ) %>%
  dplyr::select(
    pattern_class,
    idx_uniq,
    motif_name,
    motif_name_safe,
    annotation,
    annotation_broad,
    category,
    query_consensus,
    cwm_fwd,
    cwm_rev,
    total_log10_hits,
    total_n_seqlets,
    n_component_celltypes,
    top_organ,
    cwm_entropy,
    entropy_ratio,
    pattern,
    best_match = TF0,
    best_match_TOMTOM_qval = qval0
  ) %>% 
  arrange(desc(n_component_celltypes))



# Create widget ----------------------------------------------------------------

dt_widget <- datatable(
  motif_anno2,
  # Tell DT not to escape the HTML in the image columns
  escape = FALSE,
  # Add filtering controls ('top', 'bottom', or 'none')
  filter = 'top',
  # Standard DT options
  options = list(
    pageLength = 600,
    autoWidth = TRUE,
    columnDefs = list(
      # Make the logo columns wider
      list(width = '300px', targets = c("cwm_fwd", "cwm_rev")),
      # Center the motif image
      list(className = 'dt-center', targets = c("cwm_fwd", "cwm_rev"))
    )
  ),
  # Give the table a caption
  caption = 'HDMA de novo motif compendium'
) %>%
  # Color numeric columns with a background bar
  formatStyle(
    'total_n_seqlets',
    background = styleColorBar(range(motif_anno2$total_n_seqlets, na.rm=TRUE), col2hex('olivedrab3')),
    backgroundSize = '98% 88%',
    backgroundRepeat = 'no-repeat',
    backgroundPosition = 'center'
  ) %>%
  # Color numeric columns with a background bar
  formatStyle(
    'total_log10_hits',
    background = styleColorBar(range(motif_anno2$total_log10_hits, na.rm=TRUE), col2hex('olivedrab3')),
    backgroundSize = '98% 88%',
    backgroundRepeat = 'no-repeat',
    backgroundPosition = 'center'
  ) %>%
  # Color numeric columns with a background bar
  formatStyle(
    'n_component_celltypes',
    background = styleColorBar(range(motif_anno2$n_component_celltypes, na.rm=TRUE), col2hex('olivedrab3')),
    backgroundSize = '98% 88%',
    backgroundRepeat = 'no-repeat',
    backgroundPosition = 'center'
  ) %>%
  # Color another numeric column by a color scale
  formatStyle(
    'cwm_entropy',
    backgroundColor = styleInterval(
      cut = seq(min(motif_anno2$cwm_entropy, na.rm=TRUE), max(motif_anno2$cwm_entropy, na.rm=TRUE), length.out=99),
      values = ylrd
    )
  ) %>%
  # Color another numeric column by a color scale
  formatStyle(
    'entropy_ratio',
    backgroundColor = styleInterval(
      cut = seq(min(motif_anno2$entropy_ratio, na.rm=TRUE), max(motif_anno2$entropy_ratio, na.rm=TRUE), length.out=99),
      values = ylrd
    )
  ) %>%
  # Color categorical columns by their value
  formatStyle(
    'category',
    backgroundColor = styleEqual(
      levels = names(cmap_category),
      values = as.character(col2hex(cmap_category))
    )
  ) %>%
  # Color categorical columns by their value
  formatStyle(
    'top_organ',
    backgroundColor = styleEqual(
      levels = names(cmap_organ),
      values = as.character(col2hex(cmap_organ))
    )
  ) %>%
  formatStyle(
    'pattern_class',
    color = styleEqual(c('pos_patterns', 'neg_patterns'), col2hex(c('tomato1', 'steelblue3'))),
    fontWeight = 'bold'
  )

saveWidget(dt_widget, here("MOTIFS.html"), selfcontained = TRUE)
