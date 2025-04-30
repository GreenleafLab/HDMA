# share color palettes for HDMA

# organ palette, modifying on cmaps_BOR$stallion
# also saved to hdma_shared/cmap_organ.rds
cmap_organ <- c("Heart"  = "#D51F26",
                "Brain"   = "#0C727C",
                "Eye"     = "#ff9f0f",
                "Adrenal" = "#876941",
                "Stomach" = "#208A42",
                "StomachEsophagus" = "#208A42",
                "Liver"   = "#3b46a3",
                "Lung"    = "#f0643e",
                "Muscle"  = "#89C75F", 
                "Skin"    = "#ad0773",
                "Spleen"  = "#3BBCA8",
                "Kidney"  = "#7E1416",
                "Thymus"  = "#6E4B9E",
                "Thyroid" = "#8A9FD1")

cmap_organcode <- c("HT" = "#D51F26",
                    "BR" = "#0C727C",
                    "EY" = "#ff9f0f",
                    "AG" = "#876941",
                    "ST" = "#208A42",
                    "LI" = "#3b46a3",
                    "LU" = "#f0643e",
                    "MU" = "#89C75F", 
                    "SK" = "#ad0773",
                    "SP" = "#3BBCA8",
                    "KI" = "#7E1416",
                    "TM" = "#6E4B9E",
                    "TR" = "#8A9FD1")

# age / PCW palette
cmap_pcw <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(n=12)
names(cmap_pcw) <- c("PCW10", "PCW12", "PCW13", "PCW14", "PCW15", "PCW17",
                     "PCW18", "PCW19", "PCW20", "PCW21", "PCW22", "PCW23")

cmap_pcw_discrete <- c("PCW10" = "#FFFFD9", # this is just from running the two lines above
                       "PCW12" = "#F1F9BB",
                       "PCW13" = "#DBF1B2",
                       "PCW14" = "#B9E3B5",
                       "PCW15" = "#85CFBA",
                       "PCW17" = "#57BEC0",
                       "PCW18" = "#33A8C2",
                       "PCW19" = "#1D8CBD",
                       "PCW20" = "#2167AC",
                       "PCW21" = "#23479D",
                       "PCW22" = "#1D2D83",
                       "PCW23" = "#081D58")          

# data types
cmap_genescore <- cmaps_BOR$horizonExtra # could try viridis::plasma() if needed
cmap_rna       <- cmaps_BOR$blueYellow
cmap_chromvar  <- cmaps_BOR$solarExtra


# cell types
cmap_compartment <- c(end = "#7F3C8D", epi = "#11A579", str = "#3969AC", imm = "#F2B701")

# more subtle colors
cmap_compartment2 <- c("end" = "gray80", "imm" = "black", "str" = "gray40", "epi" = "gray60")


cmap_category <- c("base"            = "firebrick1",
                   "base_with_flank" = "firebrick4",
                   "homocomposite"   = "darkorchid4",
                   "heterocomposite" = "royalblue3", # navy
                   "unresolved"      = "black",
                   "repeat"          = "gray30",
                   "partial"         = "gray50",
                   "exclude"         = "gray90")

# as in ArchR
cmap_peaktype <- c("Distal" = "#75B76C", "Exonic" = "#8DBFE6", "Intronic" = "#5B328C", "Promoter" = "#F5C768")

# blues
# cmap_peaktype2 <- RColorBrewer::brewer.pal(4, "Blues") 
cmap_peaktype2 <- c("Distal" = "darkgreen", "Intronic" = "darkolivegreen3", "Promoter" = "royalblue2", "Exonic" = "midnightblue")

# correlation (used for coloring ABC/P2G loops by correlation score)
cmap_cor <- c("#bfd3e6","#8c96c6","#88419d","#4d004b")
