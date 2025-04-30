# Plotting settings

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(magrittr)
  library(ggplot2)
  library(ggrastr)
})


theme_BOR <- function(base_size=14, base_family="Helvetica", border = TRUE) {
  library(grid)
  library(ggthemes)
  # Should plots have a bounding border?
  if(border){
    panel.border <- element_rect(fill = NA, color = "black", size = 0.7)
    axis.line <- element_blank()
  }else{
    panel.border <- element_blank()
    axis.line <- element_line(color = "black", size = 0.5)
  }
  
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = panel.border,
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = axis.line,
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_blank(),
            strip.text = element_text()
    ))
  
}

scale_fill_BOR <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_BOR <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#---------------------------
# Colormaps
#---------------------------

cmaps_BOR <- list(
  # Many of these adapted from ArchR ColorPalettes.R by Jeff Granja or colors.R from BuenColors
  # https://github.com/GreenleafLab/ArchR/blob/master/R/ColorPalettes.R
  # https://github.com/caleblareau/BuenColors/blob/master/R/colors.R
  
  ## Sequential colormaps:
  solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', 
                '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'), #buencolors
  sunrise = c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
              "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
  horizon = c('#000075', '#2E00FF', '#9408F7', '#C729D6', '#FA4AB5', 
              '#FF6A95', '#FF8B74', '#FFAC53', '#FFCD32', '#FFFF60'),
  horizonExtra =c("#000436", "#021EA9", "#1632FB", "#6E34FC", "#C732D5",
                  "#FD619D", "#FF9965", "#FFD32B", "#FFFC5A"),
  blueYellow = c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                  "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D"),
  sambaNight = c('#1873CC','#1798E5','#00BFFF','#4AC596','#00CC00',
                  '#A2E700','#FFFF00','#FFD200','#FFA500'), #buencolors
  wolfgang_basic = c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", 
                    "#1D91C0", "#225EA8", "#253494", "#081D58"), #buencolors
  wolfgang_extra = c("#FFFFFF", "#FCFED3", "#E3F4B1", "#ABDEB6", "#60C1BF", 
                    "#2A9EC1", "#206AAD", "#243996", "#081D58"), #buencolors
  whitePurple = c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                  '#8c6bb1','#88419d','#810f7c','#4d004b'),
  whiteBlue = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf',
                '#3690c0','#0570b0','#045a8d','#023858'),
  whiteViolet = c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1', 
                  '#DD3497', '#AE017E', '#7A0177', '#49006A'),
  comet = c("#E6E7E8","#3A97FF","#8816A7","black"),

  flame_flame = c('#000033', '#0000A5', '#1E00FB', '#6F00FD', '#C628D6', 
    '#FE629D', '#FF9B64', '#FFD52C', '#FFFF5F'), # buencolors

  flame_short = c('#000033', '#0000A5', '#1E00FB', '#6F00FD', '#C628D6', 
    '#FE629D', '#FF9B64', '#FFD52C'), # Stop short of yellow (better for tracks, etc.)

  #7-colors
  greenBlue = c('#e0f3db','#ccebc5','#a8ddb5','#4eb3d3','#2b8cbe',
                '#0868ac','#084081'),
  
  #6-colors
  beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30","#F7962E","#FCEE2B"),
  
  #5-colors
  fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A"),
  greyMagma = c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"),
  fireworks2 = c("black", "#2488F0","#7F3F98","#E22929","#FCB31A"),
  purpleOrange = c("#581845", "#900C3F", "#C70039", "#FF5744", "#FFC30F"),
  beach = c("#87D2DB","#5BB1CB","#4F66AF","#F15F30","#F7962E","#FCEE2B"),
  zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"), #wesanderson
  darjeeling = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"), #wesanderson
  rushmore = c("#E1BD6D", "#EABE94", "#0B775E","#35274A" , "#F2300F"), #wesanderson
  FantasticFox1 = c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20"), #wesanderson
  BottleRocket2 = c("#FAD510", "#CB2314", "#273046", "#354823", "#1E1E1E"), #wesanderson
  Moonrise3 = c("#85D4E3", "#F4B5BD", "#9C964A", "#CDC08C", "#FAD77B"), #wesanderson
  fireworks = c("white","#2488F0","#7F3F98","#E22929","#FCB31A"),

  # Divergent sequential:
  coolwarm = c("#4858A7", "#788FC8", "#D6DAE1", "#F49B7C", "#B51F29"),
  brewer_yes = c("#053061", "#2971B1", "#6AACD0","#C1DDEB", "#F7F7F7", 
                  "#FACDB5", "#E58267", "#BB2933", "#67001F"), #buencolors
  brewer_celsius = c("#313695", "#5083BB", "#8FC3DD", "#D2ECF4", "#FFFFBF", 
                      "#FDD384", "#F88D51", "#DE3F2E", "#A50026"), #buencolors
  flame_blind = c("#0DB2AA", "#0AD7D3", "#00FFFF", "#B1FFFE", "#FFFFFF", 
                  "#FFA3EC", "#FF00D8", "#BD00EC", "#5F00FF"), #buencolors
  solar_flare = c('#3361A5', '#2884E7', '#1BA7FF', '#76CEFF', '#FFFFFF', 
                  '#FFE060', '#FA8E24', '#DA2828', '#A31D1D'), #buencolors
  brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', 
                '#FACDB5', '#E58267', '#BB2933', '#67001F'), #buencolors
  
  ## Qualitative colormaps:
  
  # see: https://carto.com/carto-colors/
  cartoPrism = c('#7F3C8D', '#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310', 
                  '#008695', '#CF1C90', '#F97B72', '#4B4B8F'),
  cartoSafe = c('#88CCEE', '#CC6677', '#DDCC77', '#117733', '#332288', '#AA4499', '#44AA99',
                 '#999933', '#882255', '#661100', '#6699CC'),
  cartoBold = c('#7F3C8D' ,'#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310',
                 '#008695', '#CF1C90', '#f97b72', '#4b4b8f'),
  cartoAntique = c('#855C75', '#D9AF6B', '#AF6458', '#736F4C', '#526A83', '#625377', '#68855C',
                    '#9C9C5E', '#A06177', '#8C785D', '#467378'),
  cartoPastel = c('#66C5CC', '#F6CF71', '#F89C74', '#DCB0F2', '#87C55F', '#9EB9F3', '#FE88B1',
                   '#C9DB74', '#8BE0A4', '#B497E7', '#D3B484'),
  cartoVivid = c('#E58606', '#5D69B1', '#52BCA3', '#99C945', '#CC61B0', '#24796C', '#DAA51B',
                  '#2F8AC4', '#764E9F', '#ED645A', '#CC3A8E'),
  # 15 color
  circus = c("#D52126", "#88CCEE", "#FEE52C", "#117733", "#CC61B0", "#99C945", "#2F8AC4", "#332288", 
              "#E68316", "#661101", "#F97B72", "#DDCC77", "#11A579", "#89288F", "#E73F74"),
  iron_man = c('#371377','#7700FF','#9E0142','#FF0080', '#DC494C',"#F88D51","#FAD510","#FFFF5F",'#88CFA4',
               '#238B45',"#02401B","#0AD7D3","#046C9A", "#A2A475", 'grey35'),
  # The following 3 were designed by Ryan Corces.
  stallion = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B", "#FEE500","#8A9FD1","#C06CAB", "#D8A767",
               "#90D5E4", "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8", "#6E4B9E","#0C727C", "#7E1416", "#E6C2DC"),
  calm = c("#7DD06F", "#844081", "#688EC1", "#C17E73", "#484125", "#6CD3A7", "#597873","#7B6FD0", "#CF4A31", "#D0CD47",
           "#722A2D", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D", "#403552", "#76D73C", "#96CED5", "#CE54D1", "#C48736"),
  kelly = c("#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
            "#FF7A5C", "#53377A", "#FF8E00","#B32851", "#F4C800", "#7F180D", "#93AA00", "#593315", "#F13A13")
)

#--------------------------
# Colormap helper functions
#--------------------------

mostDifferentColors <- function(cols, n=20, colorspace="Lab", startingCols=NULL){
  stopifnot(length(cols) > n)
  rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue=255)
  
  # Convert sRGB to another colorspace (more 'perceptually uniform' colorspace, e.g. "Lab")
  rgbCols <- t(col2rgb(cols))
  conv <- grDevices::convertColor(rgbCols, from="sRGB", to=colorspace, scale.in=255)
  
  # Now select n 'furthest neighbors' colors
  # This performs an iterative procedure for picking colors that maximize
  # 'distance' to already selected colors. The first color is picked randomly.
  # If starting cols provided, add these to the list of picked cols
  if(!is.null(startingCols)){
    stConv <- grDevices::convertColor(t(col2rgb(startingCols)), from="sRGB", to=colorspace, scale.in=255)
    pickedColors <- list()
    for(i in seq_len(nrow(stConv))){
      pickedColors[[i]] <- stConv[i,]
    }
    remainingColors <- conv
  }else{
    idx <- sample(1:nrow(conv), 1)
    pickedColors <- list(conv[idx,])
    remainingColors <- conv[-idx,]
  }
  pickedLen <- length(pickedColors)
  
  # Iteratively add the furthest color from the selected colors
  for(i in seq(pickedLen, n - 1)){
    distList <- list()
    for(j in seq_along(pickedColors)){
      colJ <- pickedColors[[j]]
      distMat <- dist(rbind(colJ, remainingColors), method="euclidean") %>% as.matrix
      distList[[j]] <- distMat[2:nrow(distMat),1]
    }
    # Maximize the minimum distance between each color
    distMat <- do.call(cbind, distList)
    distMins <- apply(distMat, 1, FUN = min)
    idx <- which(max(distMins) == distMins)
    pickedColors[[i + 1]] <- remainingColors[idx,]
    remainingColors <- remainingColors[-idx,]
  }
  pickedLab <- do.call(rbind, pickedColors)
  pickedRgb <- round(grDevices::convertColor(pickedLab, from = colorspace, to = "sRGB", scale.out = 255),0)
  hex <- apply(pickedRgb, 1, rgb2hex)
  hex
}


mostSimilarColors <- function(color, colorOptions, n=5, colorspace="Lab"){
  stopifnot(length(colorOptions) > n)
  rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)

  colorOptions <- colorOptions[colorOptions != color]
  
  # Convert sRGB to another colorspace (more 'perceptually uniform' colorspace)
  rgb <- t(col2rgb(color))
  rgbCols <- t(col2rgb(colorOptions))
  fullMat <- rbind(rgb, rgbCols)
  rownames(fullMat) <- 1:nrow(fullMat)
  conv <- grDevices::convertColor(fullMat, from = "sRGB", to = colorspace, scale.in = 255)

  # Calcualte distances and pick n most similar to starting color
  distMat <- dist(conv, method = "euclidean") %>% as.matrix()
  pickedIdx <- distMat[1,2:ncol(distMat)] %>% sort() %>% head(.,n=n) %>% names() %>% as.integer()
  colorOptions[pickedIdx-1]
}


pairwiseColorInterpolations <- function(cols, colorspace = "Lab"){
  # Get all pairwise interpolations between a vector of colors
  rgb2hex <- function(rgb) rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
  interpolate <- function(c1, c2, colorspace){
    rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
  }
  paired <- sapply(cols, function(x) sapply(cols, function(y) interpolate(x, y, colorspace)))
  unique(as.vector(paired))
}


getColorMap <- function(cmap, n, type='qualitative'){
  stopifnot(n >= 1)
  # Return a character vector of n colors based on
  # the provided colormap. If n > length(cmap), do
  # some smart interpolation to get enough colors
  names(cmap) <- NULL # Having names on colors causes problems for some plotting routines

  if(type == 'qualitative'){
    # If qualitative colormap, do 'mostDifferent' interpolation
    if(length(cmap) < n){
      cmap <- mostDifferentColors(
        pairwiseColorInterpolations(cmap), 
        colorspace = "Apple RGB", n = n, startingCols = cmap
      )
    }

  }else{
    # Otherwise, return sequential colors based on provided palette
    colfunc <- colorRampPalette(cmap)
    cmap <- colfunc(n)
  }
  cmap[1:n]
}

plotColorMap <- function(cols){
  # Plot each of the colors in a colormap
  cols <- base::unname(cols)
  n <- length(cols)
  df <- data.frame(
    x = seq_len(n),
    y = rep(1, n),
    z = factor(seq_len(n))
  )
  p <- (
    ggplot(df, aes(x=x,y=y,color=z)) 
    + geom_tile(aes(fill=z))
    + theme_BOR()
    + scale_color_manual(values = cols)
    + scale_fill_manual(values = cols)
  )
  p
}

#-------------------
# Plotting functions
#-------------------

plotUMAP <- function(df, dataType = "qualitative", cmap = NULL, covarLabel = "", point_size=0.5, 
  namedColors=FALSE, plotTitle=NULL, colorLims=NULL, na.value="grey35", useRaster=TRUE){
  # Given a df containing the UMAP x and y coords and a third column, 
  # plot the UMAP

  if(useRaster){
    p <- (
      ggplot(df, aes(x = df[,1], y = df[,2], color = df[,3]))
      + geom_point_rast(size = point_size)
      + theme_BOR()
      + theme(
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        aspect.ratio=1
        )
      + xlab("UMAP1")
      + ylab("UMAP2")
      )
  }else{
    p <- (
      ggplot(df, aes(x = df[,1], y = df[,2], color = df[,3]))
      + geom_point(size = point_size)
      + theme_BOR()
      + theme(
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        aspect.ratio=1
        )
      + xlab("UMAP1")
      + ylab("UMAP2")
      )
  }

  # Set plot title
  if(!is.null(plotTitle)){
    p <- p + ggtitle(plotTitle)
  }else{
    p <- p + ggtitle(sprintf("n = %s", nrow(df)))
  }
  # If colormap provided, update colors
  if(!is.null(cmap)){
    if(namedColors){
      # Named colormap corresponding to discrete values in third column
      p <- p + scale_color_manual(values=cmap, limits=names(cmap), name=covarLabel, na.value=na.value)
      p <- p + guides(fill = guide_legend(title=covarLabel), 
                      colour = guide_legend(override.aes = list(size=5)))
    }else{
      # Remove names
      names(cmap) <- NULL
      if(dataType == "qualitative"){
        # check to make sure you have enough colors for qualitative mapping
        nvals <- length(unique(df[,3]))
        cmap <- getColorMap(cmap, n=nvals)
        p <- p + scale_color_manual(values=cmap, name=covarLabel, na.value=na.value)
        p <- p + guides(fill = guide_legend(title=covarLabel), 
                          colour = guide_legend(override.aes = list(size=5)))
      }else{
        if(!is.null(colorLims)){
          p <- p + scale_color_gradientn(colors=cmap, name=covarLabel, limits=colorLims, na.value=na.value)
        }else{
          p <- p + scale_color_gradientn(colors=cmap, name=covarLabel, na.value=na.value)
        }
      }
    }
  }
  p
}

# added by Betty
quickumap <- function(obj, reduct, groupby, title){
  umapDF <- data.frame(Embeddings(object=obj, reduction=reduct), obj[[groupby]])
  set.seed(1)
  umapDF <- umapDF[sample(nrow(umapDF)),]
  p <- plotUMAP(umapDF, plotTitle = title, dataType = "qualitative", cmap = cmaps_BOR$stallion, point_size=0.1)
  return(p)
}


plotEachQualitative <- function(umapDF, colors=NULL, defaultColor="red", bgColor="grey", cmap = camps_BOR$stallion, pointSize=0.5){
  # Plot separate UMAPs for each variable in the third column
  # Return plots in a named list
  plotList <- list("all" = plotUMAP(umapDF, dataType = "qualitative", cmap = cmap, point_size=pointSize))
  unq <- unique(umapDF[,3]) %>% as.character()
  if(is.null(colors)){
    colors <- rep(defaultColor, length(unq))
  }
  for(i in seq_along(unq)){
    s <- unq[i]
    sc <- colors[i]
    message(sprintf("Plotting %s...", s))
    pDF <- umapDF
    pDF[,3] <- ifelse(pDF[,3] == s, s, "other") %>% factor(., order = TRUE, levels = c(s, "other"))
    # Sort sample to front
    pDF <- pDF[order(pDF[,3], decreasing=TRUE),]
    nCells <- sum(pDF[,3] == s)
    pTitle <- sprintf("%s; n = %s", s, nCells)
    plotList[[s]] <- plotUMAP(pDF, dataType = "qualitative", cmap = c(sc, bgColor), plotTitle=pTitle, point_size=pointSize)
  }
  plotList
}


qcHistFilter <- function(df, cmap = NULL, bins=100, border_color="black", lower_lim=NULL, upper_lim=NULL){
  # Histogram, with cutoffs if included

  # Fix colormap if provided
  if(is.null(cmap)){
    cmap <- "blue"
  }
  
  p <- (
    ggplot(df, aes(x=df[,2]))
    + geom_histogram(bins=bins, fill=cmap, color=border_color)
    + xlab(colnames(df)[2])
    + ylab("Frequency")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 0.6,
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
  )

  x <- df[,2]

  if(!is.null(lower_lim)){
    if(is.finite(lower_lim)){
      p <- p + geom_vline(aes(xintercept=lower_lim), color="red", linetype="dashed")
      thresh <- round((sum(x < lower_lim) / length(x)) * 100, 2)
      p <- p + annotate("text",  x=-Inf, y = Inf, label = paste0(thresh, "%"), vjust=1, hjust=-1)
    }
  }
  if(!is.null(upper_lim)){
    if(is.finite(upper_lim)){
      p <- p + geom_vline(aes(xintercept=upper_lim), color="red", linetype="dashed")
      thresh <- round((sum(x > upper_lim) / length(x)) * 100, 2)
      p <- p + annotate("text",  x=Inf, y = Inf, label = paste0(thresh, "%"), vjust=1, hjust=1)
    }
  }
  p
}


qcBarPlot <- function(df, cmap = NULL, border_color="black", barwidth=0.5){
  # Plot a bar plot (df is a 2+ column dataframe with column 1 = x and column 2 = y)
  nsamp <- nrow(df)
  # Fix colormap if provided
  if(!is.null(cmap)){
    if(length(cmap) > 1){
      cmap <- getColorMap(cmap, n = nsamp)
    }
  }else{
    cmap <- "blue"
  }
  
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2]))
    + geom_bar(stat = "identity", fill = cmap, width=barwidth, color=border_color)
    + scale_fill_manual(values = cmap)
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
  )
  p
}


qcViolinPlot <- function(df, cmap = NULL, makeLog = FALSE){
  # Plot a violin plot
  nsamp <- length(unique(df[,1]))
  aspectRatio <- 6/nsamp
  # Assume that the first column is the sample and the second column is the variable of interest
  if(makeLog){
    df[,2] <- log10(df[,2])
    colnames(df)[2] <- paste0("log10 ", colnames(df)[2]) 
  }
  
  # Plot a violin / box plot
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2], color = df[,1]))
    + geom_violin(aes(fill = df[,1]))
    + geom_boxplot(width = 0.8, alpha = 0)
    + scale_color_manual(values = cmap)
    + scale_fill_manual(values = alpha(cmap, 0.2))
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1)) 
  )
  p
  
  # Adjust colors if necessary:
  if(!is.null(cmap)){
    cmap <- getColorMap(cmap, n = nsamp)
  }else{
    cmap <- rep("blue", times = nsamp)
  }
  p <- suppressMessages(p + scale_color_manual(values = cmap))
  p <- suppressMessages(p + scale_fill_manual(values = alpha(cmap, 0.3)))
  p
}


stackedBarPlot <- function(df, xcol = 1, fillcol = 2, ycol = 3, cmap = NULL, border_color="black", covarLabel = "", namedColors=FALSE, barwidth=0.5){
  # Plot a stacked bar plot
  # Expects a 'melted' dataframe/matrix as input
  nsamp <- length(unique((df[,xcol])))
  # Assume that we want to show all xaxis labels
  xID <- unique((df[,xcol]))

  # Fix colormap if provided
  if(!namedColors){
    if(!is.null(cmap)){
      cmap <- getColorMap(cmap, n = length(unique((df[,fillcol]))))
    }else{
      cmap <- getColorMap(cmaps_BOR$stallion, n = length(unique((df[,fillcol]))))
    }
  }
  
  p <- (
    ggplot(df, aes(x=df[,xcol], y=df[,ycol], fill=df[,fillcol]))
    + geom_bar(stat = "identity", position="fill", width=barwidth, color=border_color)
    + xlab(xcol)
    + ylab(ycol)
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio = 6/nsamp, # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
  )

  # If colormap provided, update colors
  if(namedColors){
    # Named colormap corresponding to discrete values in third column
    p <- p + scale_color_manual(values = cmap, limits = names(cmap), name = covarLabel)
    p <- p + scale_fill_manual(values = cmap, limits = names(cmap), name = covarLabel)
  }else{
    p <- p + scale_color_manual(values = cmap, name = covarLabel)
    p <- p + scale_fill_manual(values = cmap, name = covarLabel)
  }
  p 
}


groupedBarPlot <- function(df, xcol=1, ycol=2, fillcol=3, cmap = NULL, border_color="black", barwidth=0.5){
  # Plot a bar plot
  nsamp <- nrow(df)
  ngroups <- length(unique(df[,fillcol]))
  # Fix colormap if provided
  if(!is.null(cmap)){
    cmap <- getColorMap(cmap, n = ngroups, type='qualitative')
  }else{
    cmap <- getColorMap(cmaps_BOR$stallion, n=ngroups, type='qualitative')
  }
  
  p <- (
    ggplot(df, aes(x=df[,xcol], y=df[,ycol], fill=df[,fillcol]))
    + geom_bar(
      stat = "identity", 
      position=position_dodge2(width=barwidth + barwidth/(ngroups*2), preserve="single"), 
      width=barwidth, color=border_color
      )
    + scale_fill_manual(values = cmap)
    + xlab(xcol)
    + ylab(ycol)
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = 8/(nsamp + nsamp/ngroups), # What is the best aspect ratio for a bar chart?
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + scale_y_continuous(expand = c(0, 0)) # Make bars start at the axis
    + guides(
      fill = guide_legend(title=fillcol)
      )
  )
  p
}


dotPlot <- function(df, xcol, ycol, color_col, size_col, xorder=NULL, yorder=NULL, cmap=NULL, 
  color_label=NULL, size_label=NULL, aspectRatio=NULL, sizeLims=NULL, colorLims=NULL){
  # Plot rectangular dot plot where color and size map to some values in df
  # (Assumes xcol, ycol, color_col and size_col are named columns)

  # If neither x or y col order is provided, make something up
  # Sort df:
  if(is.null(xorder)){
    xorder <- unique(df[,xcol]) %>% sort()
  }
  if(is.null(yorder)){
    yorder <- unique(df[,ycol]) %>% sort()
  }
  if(is.null(aspectRatio)){
    aspectRatio <- length(yorder)/length(xorder) # What is the best aspect ratio for this chart?
  }
  df[,xcol] <- factor(df[,xcol], levels=xorder)
  df[,ycol] <- factor(df[,ycol], levels=yorder)
  df <- df[order(df[,xcol], df[,ycol]),]

  # Make plot:
  p <- (
    ggplot(df, aes(x=df[,xcol], y=df[,ycol], color=df[,color_col], size=ifelse(df[,size_col] > 0, df[,size_col], NA)))
    + geom_point()
    + xlab(xcol)
    + ylab(ycol)
    + theme_BOR(border=TRUE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,0,0.25,1), "cm"), 
            aspect.ratio = aspectRatio,
            axis.text.x = element_text(angle = 90, hjust = 1)) 
    + guides(
      fill = guide_legend(title=""), 
      colour = guide_legend(title=color_label, override.aes = list(size=5)),
      size = guide_legend(title=size_label)
      )
  )
  if(!is.null(cmap)){
    if(!is.null(colorLims)){
      p <- p + scale_color_gradientn(colors=cmap, limits=colorLims, oob=scales::squish, name = "")
    }else{
      p <- p + scale_color_gradientn(colors=cmap, name = "")
    }
  }
  if(!is.null(sizeLims)){
    p <- p + scale_size_continuous(limits=sizeLims)
  }
  p
}


# This is used primarily for making colormaps for ComplexHeatmap
makeColFun <- function(start, end, cmap, midpoint = NULL){
  # Make a color ramp function from provided start and end breaks,
  # and optionally a midpoint
  cmapLen <- length(cmap)
  if(!is.null(midpoint)){
    interpolate <- function(c1, c2, colorspace = "Lab"){
      rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
    }
    if(length(cmap) %% 2 == 0){
      # Interpolate middle colors if necessary to get midpoint
      preMidIdx <- floor(cmapLen / 2)
      midCol <- interpolate(cmap[preMidIdx], cmap[preMidIdx + 1])
      cmap <- c(cmap[1:preMidIdx], midCol, cmap[(preMidIdx + 1):cmapLen])
      cmapLen <- length(cmap)
    }
    midIdx <- ceiling(cmapLen / 2)
    breaks <- c(seq(start, midpoint, length.out = midIdx), seq(midpoint, end, length.out = midIdx)[2:midIdx])
  } else {
    breaks <- seq(start, end, length.out = cmapLen)
  }
  colorRamp2(breaks, cmap)
}


# Heatmap wrapper:
BORHeatmap <- function(
  mat, # Data to plot (matrix or dataframe)
  limits = NULL, # Enforced limits for colormap (2 dimensional array)
  clusterCols = TRUE, # Should columns be clustered
  clusterRows = TRUE, # Should rows be clustered
  labelCols = FALSE, # Should columns be labeled
  labelRows = FALSE, # Should rows be labeled
  dataColors = NULL, # Colormap for plotting data
  dataColorMidPoint = NULL, # The data value to be the middle of the color map
  customRowLabel = NULL,
  customRowLabelIDs = NULL,
  customColLabel = NULL,
  customColLabelIDs = NULL,
  customLabelWidth = 0.15,
  useRaster = TRUE, # Should heatmap be rasterized
  rasterDevice = "CairoPNG",
  rasterQuality = 5, # Raster quality. Higher is {better?}
  fontSize = 6, # Font size for labels
  showColDendrogram = FALSE, # Should the column dendrogram be shown
  showRowDendrogram = FALSE, # Should the row dendrogram be shown
  borderColor = NA, # Color for lines between cells
  mapname = " ", # 'Name' to give heatmap
  legendTitle = " ", # Name of legend
  ...
){
  
  #Packages
  suppressPackageStartupMessages(require(ComplexHeatmap))
  suppressPackageStartupMessages(require(circlize))
  
  # Make sure mat is actually a matrix
  if(!is.matrix(mat)){
    message("'mat' needs to be a matrix. Converting...")
    mat <- as.matrix(mat)
  }
  
  # Prepare color function
  if(!is.null(limits)){
    ll <- limits[1]
    ul <- limits[2]
  }else{
    ll <- min(mat, na.rm=TRUE)
    ul <- max(mat, na.rm=TRUE)
  }
  # If no colormap provided, use solarExtra
  if(is.null(dataColors)){
    dataColors <- c("1"='#3361A5', "2"='#248AF3', "3"='#14B3FF', 
                    "4"='#88CEEF', "5"='#C1D5DC', "6"='#EAD397', 
                    "7"='#FDB31A', "8"= '#E42A2A', "9"='#A31D1D')
  }
  dataColFun <- makeColFun(ll, ul, dataColors, midpoint = dataColorMidPoint)
  
  message("Preparing Heatmap...")
  hm <- Heatmap(
    # Main components:
    matrix = mat,
    name = mapname,
    col = dataColFun,
    
    # Legend options:
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      legend_width = unit(1, "cm"),
      title = legendTitle
    ),
    rect_gp = gpar(col = borderColor), 
    
    # Column options:
    show_column_names = labelCols,
    cluster_columns = clusterCols,
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    #column_names_gp = gpar(fontsize = fontSize), 
    
    # Row options:
    show_row_names = labelRows,
    cluster_rows = clusterRows,
    show_row_dend = showRowDendrogram,
    clustering_method_rows = "ward.D2",
    #row_names_gp = gpar(fontsize = fontSize), 
    
    # Raster info:
    use_raster = useRaster,
    raster_device = rasterDevice,
    raster_quality = rasterQuality,

    # Other
    ...
  )
  
  # Add row labels if provided:
  if(!is.null(customRowLabel)){
    if(is.null(customRowLabelIDs)){
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    hm <- hm + rowAnnotation(
      link = anno_mark(at = customRowLabel, labels = customRowLabelIDs, labels_gp = gpar(fontsize = fontSize)),
      width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs)
    )
  }
  
  return(hm)
}



################################################################################
# Volcano / MA plots
################################################################################

volcanoPlot <- function(df, cmap=NULL, cmap_style='qualitative', title=NULL, covarLabel="", 
  namedColors=FALSE, colorColName="color", minxmax=NULL, minymax=NULL, point_size=1){
  # Plot a volcano plot of differential genes
  # df is a 2+ column df:
  # col 1 = x axis (e.g. fold change)
  # col 2 = significance (e.g. adj log10 pval)
  # col 3 = labels (should points be labeled)
  # col 4 = point color
  # If camp is a named vector, will match values in column 4
  # min(x/y)max indicates the lowest allowable value of (x/y)max
  n <- nrow(df)
  na_col <- "grey88"

  # Convert data.table back to df if necessary (data.tables cause problems)
  df <- as.data.frame(df)

  # Get color col
  color_col <- match(colorColName, colnames(df))

  # Set colors
  if(ncol(df) < 4){
    df[,4] <- NA
  }
  nc <- length(unique(df[,color_col]))
  if(is.null(cmap)){
    cmap <- getColorMap(cmaps_BOR$stallion, n=nc)
  }
  
  # Plot a volcano plot
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2], group=df[,color_col]))
    + geom_point_rast(aes(color=df[,color_col]), size=point_size)
    #+ geom_point(aes(color=df[,color_col]), size=point_size)
    + geom_text_repel(aes(label=df[,3]), size=3, max.overlaps=Inf)
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + ggtitle(title)
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            plot.margin=unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio=1,
            #legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(expand=c(0, 0.75)) # Make bars start at the axis
  )
  if(namedColors){
    # Named colormap corresponding to discrete values in third column
    p <- p + scale_color_manual(values=cmap, limits=names(cmap), 
      name=covarLabel, na.value=na_col, drop=FALSE)
    p <- p + guides(fill = guide_legend(title=covarLabel), 
                    colour = guide_legend(override.aes=list(size=5)))
  }else{
    p <- p + scale_color_manual(values=colors, na.value=na_col)
  }
  # Enforce x and y lims if indicated
  if(!is.null(minxmax)){
    xrng <- layer_scales(p)$x$get_limits()
    xmin <- xrng[1]
    xmax <- xrng[2]
    xmin <- ifelse(abs(xmin) < minxmax, -minxmax, xmin)
    xmax <- ifelse(xmax < minxmax, minxmax, xmax)
    p <- p + xlim(xmin, xmax)
  }
  if(!is.null(minymax)){
    yrng <- layer_scales(p)$y$get_limits()
    ymin <- yrng[1]
    ymax <- yrng[2]
    ymax <- ifelse(ymax < minymax, minymax, ymax)
    suppressMessages(
      p <- p + scale_y_continuous(expand = c(0, 0), limits=c(0, ymax*1.05))
    ) 
  }
  p
}


MAPlot <- function(df, cmap=NULL, cmap_style='qualitative', title=NULL, covarLabel="", 
  namedColors=FALSE, colorColName="color", minxmax=NULL, minymax=NULL, point_size=1){
  # Plot a MA plot of differential genes
  # df is a 2+ column df:
  # col 1 = x axis (e.g. base mean expression)
  # col 2 = y axis (e.g. fold change)
  # col 3 = labels (should points be labeled)
  # col 4 = point color
  # If camp is a named vector, will match values in column 5
  # min(x/y)max indicates the lowest allowable value of (x/y)max
  # sig_cutoff
  n <- nrow(df)
  na_col <- "grey88"

  # Convert data.table back to df if necessary (data.tables cause problems)
  df <- as.data.frame(df)

  # Get color col
  color_col <- match(colorColName, colnames(df))

  # Set colors
  if(ncol(df) < 4){
    df[,4] <- NA
  }
  nc <- length(unique(df[,color_col]))
  if(is.null(cmap)){
    cmap <- getColorMap(cmaps_BOR$stallion, n=nc)
  }
  
  # Plot a MA plot
  p <- (
    ggplot(df, aes(x=df[,1], y=df[,2], group=df[,color_col]))
    + geom_point_rast(aes(color=df[,color_col]), size=point_size)
    #+ geom_point(aes(color=df[,color_col]), size=point_size)
    + geom_text_repel(aes(label=df[,3]), 
      size=3, max.overlaps=Inf,
      box.padding=0.75, force=0.75
      )
    + geom_hline(yintercept=0.0, linetype="dashed")
    + xlab(colnames(df)[1])
    + ylab(colnames(df)[2])
    + ggtitle(title)
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            plot.margin=unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio=1,
            #legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    #+ scale_y_continuous(expand=c(0, 0.75)) # Make bars start at the axis
  )
  if(namedColors){
    # Named colormap corresponding to discrete values in third column
    p <- p + scale_color_manual(values=cmap, limits=names(cmap), 
      name=covarLabel, na.value=na_col, drop=FALSE)
    p <- p + guides(fill = guide_legend(title=covarLabel), 
                    colour = guide_legend(override.aes=list(size=5)))
  }else{
    p <- p + scale_color_manual(values=colors, na.value=na_col)
  }
  # Enforce x and y lims if indicated
  if(!is.null(minxmax)){
    xrng <- layer_scales(p)$x$get_limits()
    xmin <- xrng[1]
    xmax <- xrng[2]
    xmin <- ifelse(xmin > -minxmax, -minxmax, xmin)
    xmax <- ifelse(xmax < minxmax, minxmax, xmax)
    p <- p + xlim(xmin, xmax)
  }
  if(!is.null(minymax)){
    yrng <- layer_scales(p)$y$get_limits()
    ymin <- yrng[1]
    ymax <- yrng[2]
    ymin <- ifelse(ymin > -minymax, -minymax, ymin)
    ymax <- ifelse(ymax < minymax, minymax, ymax)
    suppressMessages(
      p <- p + scale_y_continuous(expand = c(0, 0), limits=c(ymin*1.05, ymax*1.05))
    ) 
  }
  p
}



