#####################################
# Functions for GO term analysis
#####################################


suppressPackageStartupMessages({
  library(topGO)
  library(org.Hs.eg.db)
  library(stringr)
})

# topGO has a quirk where it can't report p-values less than 1e-30

calcTopGo <- function(
    allGenes, interestingGenes=NULL, pvals=NULL, geneSel=NULL,
    nodeSize=5, ontology="BP",
    alg="weight01", stat="fisher", topNodes=50
    ){
    # Calculate GO term enrichments using topGO on provided data
    # https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
    ############################################################
    # allGenes: vector of genenames to be used in GO term search. Expects gene 'symbol' 
    # interestingGenes: predefined list of 'instersting' genes. Incompatible with supplying pvalues.
    # geneSel: function for selecting 'interesting' genes. Can only really be a p-value cutoff...
    # pvals: vector of pvalues corresponding to geneList. If not provided, will assign everything to 1
    # nodeSize: will prune terms that have less than nodeSize number of genes
    # ontology: which GO ontology to use (MF, BP, CC)
    # alg: algorithm to be used for testing GO terms (topGO default is 'weight01')
    # stat: test statistic to use for significant GO terms
    # topNodes: how many GO terms to return in result table

    # Prepare geneList as expected for topGO (i.e. value vector with names of genes)
    if(!is.null(interestingGenes)){
        message(sprintf("Running GO enrichments with %s genes in universe of %s...", 
            length(interestingGenes), length(allGenes)))
        geneList <- factor(as.integer(allGenes %in% interestingGenes))
        names(geneList) <- allGenes
        # Create topGOdata object
        GOdata <- suppressMessages(new(
            "topGOdata",
            ontology = ontology,
            allGenes = geneList,
            annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
            nodeSize = nodeSize
            ))
    }else{
        geneList <- pvals
        names(geneList) <- allGenes
        message(sprintf("Running GO enrichments with %s genes in universe of %s...", 
            sum(geneSel(geneList)), length(allGenes)))
        GOdata <- suppressMessages(new(
            "topGOdata",
            ontology = ontology,
            allGenes = geneList,
            geneSel = geneSel,
            annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
            nodeSize = nodeSize
            ))
    }

    # Test for enrichment using Fisher's Exact Test
    GOresult <- suppressMessages(runTest(GOdata, algorithm=alg, statistic=stat))
    GenTable(GOdata, pvalue=GOresult, topNodes=topNodes, numChar=1000)
}



topGObarPlot <- function(goRes, cmap = NULL, nterms=10, border_color="black", 
    barwidth=0.85, title="", enrichLimits=c(0.0, 5.5), barLimits=NULL){
    # Plot GO results in bar plot form
    goRes$log2FoldEnrichment <- log2(goRes$Significant / goRes$Expected)
    goRes$log2FoldEnrichment <- ifelse(goRes$log2FoldEnrichment > enrichLimits[2], enrichLimits[2], goRes$log2FoldEnrichment)
    goRes$threshPval <- ifelse(goRes$pvalue == "< 1e-30", 1e-30, as.numeric(goRes$pvalue))
    goRes$log10pval <- -log10(goRes$threshPval)
    if(!is.null(barLimits)){
        goRes$log10pval <- ifelse(goRes$log10pval < barLimits[2], goRes$log10pval, barLimits[2])
    }

    # Only plot the top nterms (reverse order to plot most significant at top)
    goRes <- goRes[1:nterms,]
    goRes <- goRes[nrow(goRes):1,]

    if(is.null(cmap)){
        cmap <- cmaps_BOR$comet
    }
    p <- (
        ggplot(goRes, aes(x=Term, y=log10pval, fill=log2FoldEnrichment))
        + geom_bar(stat="identity", width=barwidth, color=border_color)
        + scale_x_discrete(
            limits=goRes$Term, # Required to prevent alphabetical sorting of terms
            labels= function(x) str_wrap(x, width=40) # wrap long GO term labels
            ) 
        + scale_fill_gradientn(colors=cmap, limits=enrichLimits)
        + xlab("")
        + ylab("-log10 pvalue")
        + ggtitle(title)
        + theme_BOR(border=FALSE)
        + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = 6/nterms, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1)) 
        + coord_flip()
    )
    if(!is.null(barLimits)){
        p <- p + scale_y_continuous(limits=barLimits, expand=c(0,0))
    }else{
        p <- p + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
    }
    p
}

