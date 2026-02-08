#' Simulate gene expression
#'
#' @param num_genes Number of genes
#' @param num_samples Number of samples
#' @param seed Random seed for reproducibility
#' @return Numeric matrix of expression values
simulate_expression <- function(num_genes = 1000,
                                num_samples = 20,
                                seed = 123) {
  
  set.seed(seed)
  
  # Simulate expression values
  expr <- matrix(
    rlnorm(num_genes * num_samples, meanlog = 1, sdlog = 1),
    nrow = num_genes,
    ncol = num_samples
  )
  
  rownames(expr) <- paste0("Gene_", seq_len(num_genes))
  colnames(expr) <- paste0("Sample_", seq_len(num_samples))
  
  return(expr)
}


#' Inject missing residues into sequence matrix using a 'Missing Completely At Random' (MCAR) approach
#' @param expr_matrix Numeric expression matrix (genes x samples)
#' @param missing_frac Proportion of values to be set as NA
#' @param seed Random seed for reproducibility
#' @return Expression matrix with introduced missing values
introduce_missingness <- function(expr,
                                  missing_frac = 0.1,
                                  seed = 123) {
  
  set.seed
  
  # Copy matrix to avoid overwriting the original
  expr_miss <- expr
  
  total_values <-  length(expr_miss)
  n_missing <- floor(total_values * missing_frac)
  
  # Generate random linear indices
  missing_indices <- sample(seq_len(total_values), size = n_missing)
  
  expr_miss[missing_indices] <- NA
  
  return(expr_miss)
}