#' Impute NAs with per-gene mean (simple)
#' 
#' @param expr Numeric matrix (genes x samples), can contain NAs
#' @return Numeric matrix with all NA values replaced by the mean expression of 
#'         the corresponding gene (row-wise mean).
impute_gene_mean <- function(expr) {
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  gene_means <- rowMeans(expr, na.rm = TRUE)
  for (i in seq_len(nrow(expr))) {
    nas <- is.na(expr[i, ])
    if (any(nas)) expr[i, nas] <- gene_means[i]
  }
  expr
}

#' Run PCA on expression matrix (genes x samples)
#' 
#' @param expr Numeric matrix (genes x samples)
#' @param center Logical, whether to center variables
#' @param scale. Logical, whether to scale variables
#' @param n_pcs Number of principal components to retain
#' @return prcomp object with truncated PC scores
run_pca <- function(expr, center = TRUE, scale. = TRUE, n_pcs = 10) {
  expr_imp <- impute_gene_mean(expr)
  
  # PCA on samples: transpose so samples are rows
  x <- t(expr_imp)
  
  p <- stats::prcomp(x, center = center, scale. = scale.)
  # keep first n_pcs scores
  p$x <- p$x[, seq_len(min(n_pcs, ncol(p$x))), drop = FALSE]
  p
}

#' Basic PCA plot (PC1 vs PC2)
#' 
#' @param pca_res prcomp object returned by run_pca()
#' @param labels Optional character vector of sample labels
#' @return Invisibly returns NULL (produces a plot)
plot_pca <- function(pca_res, labels = NULL) {
  scores <- pca_res$x
  x <- scores[, 1]
  y <- scores[, 2]
  
  plot(x, y,
       xlab = "PC1",
       ylab = "PC2",
       main = "PCA of Samples",
       pch = 19)
  
  if (!is.null(labels)) {
    text(x, y, labels = labels, pos = 3, cex = 0.8)
  }
}

#' Run k-means clustering on PCA scores
#'
#' @param pca_res A prcomp object returned by run_pca().
#' @param k Integer; number of clusters to form.
#' @param pcs Integer vector indicating which principal components
#'        to use for clustering (e.g., 1:2 for PC1 and PC2).
#' @param seed Integer; random seed for reproducibility.
#'
#' @return A kmeans object containing:
#'         \itemize{
#'           \item cluster: cluster assignment for each sample
#'           \item centers: coordinates of cluster centers
#'           \item tot.withinss: total within-cluster sum of squares
#'         }
#'
#' @details k-means clustering is performed on selected principal
#'          component scores to group samples based on similarity
#'          in reduced-dimensional space.
run_kmeans <- function(pca_res, k = 3, pcs = 1:2, seed = 123) {
  set.seed(seed)
  x <- pca_res$x[, pcs, drop = FALSE]
  stats::kmeans(x, centers = k, nstart = 25)
}

#' Plot k-means clusters on PCA scatter plot
#'
#' @param pca_res A prcomp object returned by run_pca().
#' @param km_res A kmeans object returned by run_kmeans().
#'
#' @return Invisibly returns NULL; produces a scatter plot of
#'         PC1 vs PC2 with samples labeled or grouped by cluster.
#'
#' @details This function visualizes clustering structure by
#'          projecting samples onto the first two principal
#'          components and overlaying k-means cluster assignments.
plot_kmeans_on_pca <- function(pca_res, km_res) {
  scores <- pca_res$x
  x <- scores[, 1]
  y <- scores[, 2]
  cl <- km_res$cluster
  
  plot(x, y,
       xlab = "PC1",
       ylab = "PC2",
       main = "k-means clusters (on PCA)",
       pch = 19)
  
  # Overplot by cluster (no custom colors needed)
  for (g in sort(unique(cl))) {
    idx <- which(cl == g)
    points(x[idx], y[idx], pch = 19)
    text(mean(x[idx]), mean(y[idx]), labels = paste0("C", g), pos = 3)
  }
}
