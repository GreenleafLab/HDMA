
#####################################
# Matrix helper functions in R
#####################################


suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(tidyr)
  library(magrittr)
  library(Rcpp)
  library(RcppAlgos)
})


#####################################
# Matrix helper functions:
#####################################

# Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}


# Helper function for getting the mean of matrix groups
groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}


# Helper function for applying a custom function to matrix groups
groupFun <- function (mat, fun, groups = NULL, ...){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        # Apply custom function to columns matching group
        fun(mat[, which(groups == x), drop = F], ...)
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}


unmelt <- function(long_df, row_col, col_col, val_col){
  # 'unmelt' a long form data.frame
  # Requires columns of long_df to be named
  ##################################
  # long_df = long-format data.frame
  # row_col = name of column with id's that will become rows
  # col_col = name of column with id's that will become columns
  # val_col = name of column with values that will fill in matrix
  wide_df <- pivot_wider(long_df, id_cols=all_of(c(row_col, col_col)), names_from=all_of(col_col), values_from=all_of(val_col)) %>% as.data.frame()
  rownames(wide_df) <- wide_df[,1]
  wide_df <- wide_df[,2:ncol(wide_df)]
  return(wide_df)
}


prettyOrderMat <- function(mat, scale=TRUE, cutOff=1, lmat=NULL, clusterCols=TRUE){
  # Reorder mat in a prettier way for plotting
  # Adapted from Jeff's ArchR .binarySort
  ###################################
  # mat = matrix (like) object to sort
  # scale = should mat be scaled before building logical mat
  # cutOff = cutoff for lmat
  # lmat = logical matrix for ordering rows (binary sorting)
  # clusterCols = should columns be clustered?
  mat <- as.matrix(mat)

  if(is.null(lmat)){
    # Compute row Z-scores
    if(scale){
      lmat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    }else{
      lmat <- mat
    }
    # Logical matrix of values passing cutoff 
    lmat <- lmat >= cutOff
  }

  # Transpose:
  mat <- t(mat)
  lmat <- t(lmat)

  # Identify column ordering:
  if(clusterCols){
    hc <- hclust(dist(mat))
    colIdx <- hc$order
    mat <- t(mat[colIdx,])
    lmat <- t(lmat[colIdx,])
  }else{
    mat <- t(mat)
    lmat <- t(lmat)
    hc <- NULL
  }

  # Identify row ordering:
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  mat <- mat[rowIdx,]

  return(list(mat=mat, hclust=hc))
}


fractionXbyY <- function(x, y, add_total=FALSE, xname="x", yname="y", ylab="proportion"){
  # Returns a melted dataframe with the proportion of x by y
  # e.g. if x is cluster and y is samples, will return the proportion of each sample 
  # making up the total amount of each cluster (i.e. the proportional contribution of each y to each x)
  ###################################
  # x, y: paired vectors
  # add_total: bool indicating whether you would like to add a 'total' entry (i.e. proportion of y 
  # across all x)
  # xname, yname: char with x and y labels
  # ylab: char with proportion label
  XbyYdf <- data.frame("xgroup" = x, "ygroup" = y) %>% 
  group_by(xgroup, ygroup) %>% # Group by both x and y
  dplyr::summarize(n = n()) %>% ungroup() %>% # summarize (i.e. tally for each group(s))
  pivot_wider(names_from=xgroup, values_from=n, values_fill= list(n=0)) %>% # Expand into multiple columns
  as.data.frame()

  rownames(XbyYdf) <- XbyYdf[,1]
  XbyYmat <- as.matrix(XbyYdf[,-1]) %>% t()

  if(add_total){
    rnms <- rownames(XbyYmat)
    XbyYmat <- rbind(XbyYmat, colSums(XbyYmat))
    rownames(XbyYmat) <- c(rnms, "total")
  }
  XbyYmat <- (XbyYmat/rowSums(XbyYmat)) %>% reshape2::melt()
  colnames(XbyYmat) <- c(xname, yname, ylab) 
  # Force cluster to be qualitative
  XbyYmat[,xname] <- as.factor(XbyYmat[,xname])
  return(XbyYmat)
}


sparseLogX <- function(spmat, logtype="log2", scale=FALSE, scaleFactor=10^4){
  # Adapted from https://rdrr.io/github/YosefLab/FastProjectR/src/R/Utilities.R
  # Takes the log2(x + 1) of a potentially sparse matrix without creating an 
  # intermediate dense matrix
  ###################################
  # spmat = sparse matrix
  # logtype = which log to use (log, log2, log10)
  # scale = should matrix be scaled first?
  # scaleFactor = the amount to depth-normalize to
  stopifnot(any(logtype == c("log", "log2", "log10")))

  if(scale == TRUE){
      spmat <- t(t(spmat)/Matrix::colSums(spmat)) * scaleFactor
  }

  if(is(spmat, "sparseMatrix")){
    matsum <- summary(spmat) # Get the sparse matrix summary
    if(logtype == "log"){
      logx <- log(matsum$x + 1) 
    }else if(logtype == "log2"){
      logx <- log2(matsum$x + 1) 
    }else{
      logx <- log10(matsum$x + 1)
    }
    logmat <- sparseMatrix(i = matsum$i, j = matsum$j, # convert back to sparse matrix
                           x = logx, dims = dim(spmat),
                           dimnames = dimnames(spmat))
  }else{
    if(logtype == "log"){
      logmat <- log(spmat + 1) 
    }else if(logtype == "log2"){
      logmat <- log2(spmat + 1) 
    }else{
      logmat <- log10(spmat + 1) 
    }
  }
  return(logmat)
}

#####################################
# Expression matrix functions
#####################################


averageExpr <- function(counts_mat, group_vec, log=TRUE){
  # Can't use Seurat's average Expression function since it expects certain types of data
  # This function works on the raw counts matrix.
  # If log==TRUE, returns log2(group means + 1)

  # First, get counts per 10k
  cp10k <- sparseLogX(counts_mat, logtype="log2", scale=TRUE, scaleFactor=10000)

  # Get group means
  group_means <- groupMeans(cp10k, group_vec)

  # Log transform, if necessary
  if(log){
    group_means <- log2(group_means + 1)
  }

  return(group_means)
}


pctExpr <- function(mat){
  # Report the 'percent expressed' or fraction of cells w/ non-zero expression
  # mat is genes x cells
  nexpr <- apply(mat, 1, function(x) sum(x > 0))
  nexpr / ncol(mat)
}


avgAndPctExpressed <- function(count_mat, groups, feature_normalize=FALSE, min_pct=NULL){
  # Create dataframe of average and percent of cells expressing each feature in matrix
  # Takes the raw counts matrix as input (features x cells)

  # First, get average expression in log2CP10k space
  message("Calculating group average expression...")
  avgExpr <- averageExpr(count_mat, groups)

  # If indicated, normalize mean expression to maximum detected
  if(feature_normalize){
    message("Normalizing to maximum...")
    # First, need to drop features with zero expression
    avgExpr <- avgExpr[rowSums(avgExpr) > 0,]
    avgExpr <- avgExpr / apply(avgExpr,1,max)
  }

  # Next, get percent of cells expressing each feature
  message("Calculating percent of cells per group expressing feature...")
  pctExprMat <- groupFun(count_mat, pctExpr, groups=groups) * 100

  # Filter genes that are too lowly expressed, if indicated
  if(!is.null(min_pct)){
    message(sprintf("Filtering features expressed by less than %s%% of cells in any cluster...", min_pct))
    pctExprMat <- pctExprMat[apply(pctExprMat, 1, function(x) max(x) > min_pct),]
  }

  # Melt each matrix and then merge
  message("Reshaping and returning...")
  expMelt <- reshape2::melt(avgExpr)
  colnames(expMelt) <- c("feature", "grp", "avgExpr")
  pctMelt <- reshape2::melt(pctExprMat)
  colnames(pctMelt) <- c("feature", "grp", "pctExpr")

  # Merge dfs:
  mergedDF <- base::merge(expMelt, pctMelt, by=c("feature", "grp"))
  return(mergedDF)
}


#####################################
# Matrix Correlation tools
#####################################

# Rcpp script for performing fast row-wise pearson correlations:

sourceCpp(code='
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// Borrowed from https://github.com/GreenleafLab/ArchR/blob/master/src/Correlation.cpp
// who in turn adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
// [[Rcpp::export]]
Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
  
  if(X.ncol() != Y.ncol()){
    stop("Columns of Matrix X and Y must be equal length!");
  }

  if(max(idxX) > X.nrow()){
    stop("Idx X greater than nrow of Matrix X");
  }

  if(max(idxY) > Y.nrow()){
    stop("Idx Y greater than nrow of Matrix Y");
  }
    
  // Transpose Matrices
  X = transpose(X);
  Y = transpose(Y);
  
  const int nx = X.ncol();
  const int ny = Y.ncol();

  // Centering the matrices
  for (int j = 0; j < nx; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }

  for (int j = 0; j < ny; ++j) {
    Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
  }

  // Compute 1 over the sample standard deviation
  Rcpp::NumericVector inv_sqrt_ss_X(nx);
  for (int i = 0; i < nx; ++i) {
    inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
  }

  Rcpp::NumericVector inv_sqrt_ss_Y(ny);
  for (int i = 0; i < ny; ++i) {
    inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
  }

  //Calculate Correlations
  const int n = idxX.size();
  Rcpp::NumericVector cor(n);
  for(int k = 0; k < n; k++){
    cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X( idxX[k] - 1) * inv_sqrt_ss_Y( idxY[k] - 1);
  } 

  return(cor);

}'
)


corMatrix <- function(mat, nonredundant=TRUE, subsetx=NULL, subsety=NULL, nThreads=1){
    # Calculate correlation matrix and associated statistics using Rcpp
    ###########################
    # mat: matrix with rows to be correlated
    # nonredunant: if true, will not calculate redundant correlations (requires)
    # subsetx: vector of rownames or indices to use for calculating correlations
    # subsety: vector or rownames or indices to use for calculating correlations

    # First, make sure no rows are zero
    if(any(Matrix::rowSums(mat) == 0)){
      stop("Error: matrix contains rows of all 0's!")
    }
    if(is.null(rownames(mat))){
        rownames(mat) <- 1:nrow(mat)
    }
    mat <- as.matrix(mat)

    if(is.null(subsetx)){
      subsetx <- rownames(mat)
    }
    if(is.null(subsety)){
      subsety <- rownames(mat)
    }
    xmat <- mat[subsetx,]
    ymat <- mat[subsety,]

    # Get indices to correlate
    message("Determining indices to correlate...")
    maxRows <- max(nrow(xmat), nrow(ymat))
    # idx <- combn(1:maxRows, 2) %>% t() # Native combn is actually very slow...
    idx <- RcppAlgos::comboGeneral(1:maxRows, 2, nThreads=nThreads)
    idx <- idx[idx[,1] <= nrow(xmat),]
    idx <- idx[idx[,2] <= nrow(ymat),]
    xidx <- idx[,1]
    yidx <- idx[,2]

    df <- data.frame(
        "x" = rownames(xmat)[xidx],
        "y" = rownames(ymat)[yidx]
        )
    message(sprintf("Calculating %s correlations...", nrow(df)))
    df$Correlation <- rowCorCpp(xidx, yidx, xmat, ymat)
    message("Finished. Calculating statistics...")
    df$TStat <- (df$Correlation / sqrt((pmax(1-df$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(mat)-2))) #T-statistic P-value
    df$Pval <- 2*pt(-abs(df$TStat), ncol(mat) - 2)
    df$FDR <- p.adjust(df$Pval, method = "fdr")
    df <- df[, c("x", "y", "Correlation", "FDR")]
    return(df)
}

cor2Matrices <- function(mat1, mat2, subset1=NULL, subset2=NULL){
    # Calculate row correlations and associated statistics of two distinct matrices using Rcpp
    ###########################
    # mat1: first matrix with rows to be correlated
    # mat2: second matrix with rows to be correlated
    # subset1: vector of rownames or indices to use for calculating correlations
    # subset2: vector or rownames or indices to use for calculating correlations

    # First, make sure no rows are zero
    if(any(Matrix::rowSums(mat1) == 0)){
      stop("Error: matrix 1 contains rows of all 0's!")
    }
    if(any(Matrix::rowSums(mat2) == 0)){
      stop("Error: matrix 2 contains rows of all 0's!")
    }
    if(is.null(rownames(mat1))){
        rownames(mat1) <- 1:nrow(mat1)
    }
    mat1 <- as.matrix(mat1)
    if(is.null(rownames(mat2))){
        rownames(mat2) <- 1:nrow(mat2)
    }
    mat2 <- as.matrix(mat2)

    if(is.null(subset1)){
      subset1 <- rownames(mat1)
    }
    if(is.null(subset2)){
      subset2 <- rownames(mat2)
    }
    mat1 <- mat1[subset1,]
    mat2 <- mat2[subset2,]

    # Get indices to correlate
    message("Determining indices to correlate...")
    idx <- expand.grid(rownames(mat1), rownames(mat2))
    idx1 <- match(idx[,1], rownames(mat1))
    idx2 <- match(idx[,2], rownames(mat2))

    df <- data.frame(
        "x" = rownames(mat1)[idx1],
        "y" = rownames(mat2)[idx2]
        )
    message(sprintf("Calculating %s correlations...", nrow(df)))
    df$Correlation <- rowCorCpp(idx1, idx2, mat1, mat2)
    message("Finished. Calculating statistics...")
    df$TStat <- (df$Correlation / sqrt((pmax(1-df$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(mat1)-2))) #T-statistic P-value
    df$Pval <- 2*pt(-abs(df$TStat), ncol(mat1) - 2)
    df$FDR <- p.adjust(df$Pval, method = "fdr")
    df <- df[, c("x", "y", "Correlation", "FDR")]
    return(df)
}


corPairwise <- function(mat1, mat2){
    # Calculate row correlations and associated statistics of two paired matrices using Rcpp
    # Assumes that rows in mat1 are already paired to those in mat2 (i.e. row 1 of mat 1 will be 
    # correlated to row 1 of mat 2, and so on)
    ###########################
    # mat1: first matrix with rows to be correlated
    # mat2: second matrix with rows to be correlated

    # First, make sure no rows are zero
    if(any(Matrix::rowSums(mat1) == 0)){
      stop("Error: matrix 1 contains rows of all 0's!")
    }
    if(any(Matrix::rowSums(mat2) == 0)){
      stop("Error: matrix 2 contains rows of all 0's!")
    }
    if(is.null(rownames(mat1))){
        rownames(mat1) <- 1:nrow(mat1)
    }
    mat1 <- as.matrix(mat1)
    if(is.null(rownames(mat2))){
        rownames(mat2) <- 1:nrow(mat2)
    }
    mat2 <- as.matrix(mat2)

    df <- data.frame(
        "ix" = rownames(mat1)
        )
    message(sprintf("Calculating %s correlations...", nrow(df)))
    df$Correlation <- rowCorCpp(1:nrow(mat1), 1:nrow(mat2), mat1, mat2)
    message("Finished. Calculating statistics...")
    df$TStat <- (df$Correlation / sqrt((pmax(1-df$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(mat1)-2))) #T-statistic P-value
    df$Pval <- 2*pt(-abs(df$TStat), ncol(mat1) - 2)
    df$FDR <- p.adjust(df$Pval, method = "fdr")
    df <- df[, c("ix","Correlation", "FDR")]
    return(df)
}
