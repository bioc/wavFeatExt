#' Principal Component Analysis (PCA) for a List of Matrices
#'
#' Applies principal component analysis (PCA) to each matrix in a
#' collection of data sets and returns cumulative principal components
#' for use in downstream classification or feature extraction.
#'
#' @param x Either a numeric matrix, a list of numeric matrices (with rows
#'   corresponding to observations and columns to variables), or the output
#'   from \code{\link{simulateCNA}}. Each matrix will undergo PCA separately.
#' @param k Integer specifying the number of principal components to extract
#'   from each matrix. Must be at least 2.
#' @param ... Additional arguments passed to \code{\link[stats]{prcomp}}.
#'
#' @details
#' For each matrix in \code{x}, the function calls
#' \code{\link[stats]{prcomp}} to compute principal component scores.
#' Let \eqn{Z} denote the resulting score matrix (rows = observations,
#' columns = principal components).
#'
#' The function then constructs a list of \code{k - 1} matrices per data
#' set, where the \eqn{j}-th element (for \eqn{j = 2, \dots, k}) contains
#' the first \eqn{j} principal components, i.e. \code{Z[, 1:j]}. This
#' cumulative structure mirrors the behaviour of \code{\link{getIca}} and
#' provides feature sets with increasing dimensionality.
#'
#' @return A list of length equal to the number of matrices in \code{x}.
#'   For each input matrix, an element is returned which is itself a list
#'   of length \code{k - 1}. The \eqn{j}-th element (for \eqn{j = 2, \dots, k})
#'   of the inner list is a numeric matrix containing the first \eqn{j}
#'   principal components for that data set.
#'
#' @note
#' The cumulative component representation is designed to be used in
#' combination with \code{\link{classifyPcaIca}} for comparing
#' classification performance across different feature dimensions.
#'
#' @references
#' Jolliffe, I. T. (2002).
#' \emph{Principal Component Analysis}. 2nd edition. Springer.
#'
#' @seealso
#' \code{\link{getIca}},
#' \code{\link{classifyPcaIca}},
#' \code{\link{simulateCNA}},
#' \code{\link[stats]{prcomp}}
#'
#' @examples
#' set.seed(10)
#'
#' ## Small simulated data for fast execution
#' sim.dat <- simulateCNA(
#'   n.obs = 20,
#'   p = 32,
#'   n.sim = 1,
#'   n.block = 8,
#'   verbose = FALSE
#' )
#'
#' ## Obtain PCA-based features
#' pca <- getPca(sim.dat, k = 5)
#'
#' length(pca)
#' length(pca[[1]])
#' sapply(pca[[1]], dim)
#'
#' @author Maharani Ahsani Ummi and Arief Gusnanto
#'
#' @export
#' @importFrom stats prcomp
getPca <- function(x, k, ...) {
  if (missing(k) || length(k) != 1L || k < 2L) {
    stop("'k' must be a single integer >= 2.")
  }
  k <- as.integer(k)
  
  ## Normalise input to list of matrices
  if (is.matrix(x) || is.data.frame(x)) {
    x_list <- list(as.matrix(x))
  } else if (is.list(x)) {
    x_list <- lapply(x, function(z) {
      if (!is.matrix(z)) {
        z <- as.matrix(z)
      }
      if (!is.numeric(z)) {
        stop("All matrices in 'x' must be numeric.")
      }
      z
    })
  } else {
    stop("'x' must be a matrix, data frame, or list of matrices.")
  }
  
  all.res <- vector("list", length(x_list))
  
  for (i in seq_along(x_list)) {
    Xi <- x_list[[i]]
    
    ## PCA
    pc_fit <- prcomp(Xi, ...)
    scores <- pc_fit$x  # scores are safer than Xi %*% rotation
    
    ## cek tersedia cukup komponen
    if (ncol(scores) < k) {
      stop(
        sprintf(
          "Matrix %d has only %d PCs available; 'k' was %d.",
          i, ncol(scores), k
        )
      )
    }
    
    ## Buat cumulative PC feature sets: PC1–2, PC1–3, ..., PC1–k
    comp_list <- lapply(
      seq.int(2L, k),
      function(j) scores[, seq_len(j), drop = FALSE]
    )
    all.res[[i]] <- comp_list
  }
  
  all.res
}