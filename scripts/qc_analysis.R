#' Compute QC metrics for expression matrix
#' 
#' @param expr Numeric matrix (genes x samples), can contain NAs
#' @param detected_threshold Values > threshold count as "detected" after transform
#' @param corr_to_mean Whether to compute correlation of each sample to the mean
#' @return List with sample_qc (data.frame) and gene_qc (data.frame)
compute_qc_metrics <- function(expr,
                               detected_threshold = 0,
                               corr_to_mean = TRUE) {
  
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  
  # ---- Sample-level metrics (per column) ----
  n_genes <- nrow(expr)
  
  sample_missing_frac <- colMeans(is.na(expr))
  sample_mean <- colMeans(expr, na.rm = TRUE)
  sample_median <- apply(expr, 2, median, na.rm = TRUE)
  sample_sd <- apply(expr, 2, sd, na.rm = TRUE)
  
  # Detected genes per sample: count of (non-NA) values over threshold
  sample_detected_genes <- colSums(!is.na(expr) & expr > detected_threshold)
  
  # Generate QC summary in a data frame
  sample_qc <- data.frame(
    sample_id = colnames(expr),
    n_genes = n_genes,
    missing_frac = sample_missing_frac,
    mean = sample_mean,
    median = sample_median,
    sd = sample_sd,
    detected_genes = sample_detected_genes,
    detected_frac = sample_detected_genes / n_genes,
    stringsAsFactors = FALSE
  )
  
  # ---- Gene-level metrics (per row) ----
  n_samples <- ncol(expr)
  
  gene_missing_frac <- rowMeans(is.na(expr))
  gene_mean <- rowMeans(expr, na.rm = TRUE)
  gene_median <- apply(expr, 1, median, na.rm = TRUE)
  gene_sd <- apply(expr, 1, sd, na.rm = TRUE)
  
  # Number of samples that detect a gene: count of (non-NA) values over threshold
  gene_detected_samples <- rowSums(!is.na(expr) & expr > detected_threshold)
  
  gene_qc <- data.frame(
    gene_id = rownames(expr),
    n_samples = n_samples,
    missing_frac = gene_missing_frac,
    mean = gene_mean,
    median = gene_median,
    sd = gene_sd,
    detected_samples = gene_detected_samples,
    detected_frac = gene_detected_samples / n_samples,
    stringsAsFactors = FALSE
  )
  
  return(list(sample_qc = sample_qc, gene_qc = gene_qc))
}