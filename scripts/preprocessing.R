#' Anything that transforms the matrix

#' Normalize expression matrix
#' 
#' @param expr Numeric matrix (genes x samples)
#' @param method One of: "log1p", "zscore_gene", "cpm_log1p"
#' @param pseudocount Value added before log transform
#' @return Normalized numeric matrix
normalize_expression <- function(expr,
                                 method = c("log1p", "zscore_gene", "cpm_log1p"),
                                 pseudocount = 1) {
  
  method = match.arg(method)
  
  # Convert expr to a matrix if it is not already
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  
  # Log(1 + x)
  if (method == "log1p") {
    return(log(expr + pseudocount))
  }
  
  # Counts Per Million + log(1 + x)
  if (method == "cpm_log1p") {
    # Column sums ignoring NAs
    lib_sizes <- colSums(expr, na.rm = TRUE)
    # Avoid divide-by-zero
    lib_sizes[lib_sizes == 0] <- NA
    
    cpm <- sweep(expr, 2, lib_sizes, FUN = "/") * 1e6
    return(log(cpm + pseudocount))
  }
  
  # Z Score
  if (method == "zscore_gene") {
    gene_means <- rowMeans(expr, na.rm = TRUE)
    gene_sds <- apply(expr, 1, sd, na.rm = TRUE)
    gene_sds[gene_sds ==0] <- NA
    
    centered <- sweep(expr, 1, gene_means, FUN = "-")
    return(sweep(centered, 1, gene_sds, FUN = "/"))
  }
}


#' optional: filter_low_detected_genes(), impute_smple()