get_project_root <- function() {
  # Works when sourced from file (including when knitting)
  this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NULL)
  
  if (!is.null(this_file)) {
    return(dirname(this_file))  # directory containing run_pipeline.R
  }
  
  # Fallback: current working directory
  getwd()
}

#' Run full simulation -> missingness -> normalization -> QC -> PCA/kmeans pipeline
#'
#' @param num_genes Number of genes
#' @param num_samples Number of samples
#' @param missing_frac Fraction missing to introduce
#' @param seed Random seed
#' @param norm_method Normalization method (passed to normalize_expression)
#' @param detected_threshold Threshold used for QC detection metrics
#' @param do_pca Whether to run PCA
#' @param do_kmeans Whether to run k-means (usually on PCA scores)
#' @param k Number of clusters for k-means
#' @param out_dir Optional output directory to save plots/tables
#' @return A list of all intermediate and final results
run_pipeline <- function(num_genes = 1000,
                         num_samples = 30,
                         missing_frac = 0.1,
                         seed = 123,
                         norm_method = "cpm_log1p",
                         detected_threshold = 0,
                         do_pca = TRUE,
                         do_kmeans = TRUE,
                         k = 3,
                         out_dir = NULL) {
  
  # ---- 1) Simulate ----
  expr_raw <- simulate_expression(num_genes = num_genes,
                                  num_samples = num_samples,
                                  seed = seed)
  
  # ---- 2) Missingness ----
  expr_miss <- introduce_missingness(expr_raw,
                                     missing_frac = missing_frac,
                                     seed = seed + 1)
  
  # ---- 3) Normalize ----
  expr_norm <- normalize_expression(expr_miss,
                                    method = norm_method)
  
  # ---- 4) QC metrics ----
  qc <- compute_qc_metrics(expr_norm,
                           detected_threshold = detected_threshold,
                           corr_to_mean = TRUE)
  
  # Optional: save QC tables
  if (!is.null(out_dir)) {
    write.csv(qc$sample_qc, file.path(out_dir, "tables", "sample_qc.csv"), row.names = FALSE)
    write.csv(qc$gene_qc, file.path(out_dir, "tables", "gene_qc.csv"), row.names = FALSE)
  }
  
  # ---- 5) PCA (recommended on normalized data) ----
  pca_res <- NULL
  if (do_pca) {
    pca_res <- run_pca(expr_norm)  # defined below (or in dim_reduction.R)
    
    if (!is.null(out_dir)) {
      png(file.path(out_dir, "figures", "pca_scatter.png"), width = 900, height = 700)
      plot_pca(pca_res)            # defined below
      dev.off()
    }
  }
  
  # ---- 6) k-means (usually on PCA scores) ----
  km_res <- NULL
  if (do_kmeans) {
    if (is.null(pca_res)) {
      stop("do_kmeans=TRUE requires do_pca=TRUE (k-means uses PCA scores).")
    }
    km_res <- run_kmeans(pca_res, k = k)  # defined below
    
    if (!is.null(out_dir)) {
      png(file.path(out_dir, "figures", "kmeans_pca.png"), width = 900, height = 700)
      plot_kmeans_on_pca(pca_res, km_res) # defined below
      dev.off()
    }
  }
  
  # Return everything for your Rmd to use
  list(
    expr_raw = expr_raw,
    expr_miss = expr_miss,
    expr_norm = expr_norm,
    sample_qc = qc$sample_qc,
    gene_qc = qc$gene_qc,
    pca = pca_res,
    kmeans = km_res
  )
}
