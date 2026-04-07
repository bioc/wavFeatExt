#' Independent Component Analysis (ICA) for a List of Matrices
#'
#' Applies independent component analysis (ICA) to each matrix in a
#' collection of data sets and returns cumulative independent components
#' for use in downstream classification or feature extraction.
#'
#' @param x Either a numeric matrix, a list of numeric matrices (with rows
#'   corresponding to observations and columns to variables), or the output
#'   from \code{\link{simulateCNA}}. Each matrix will undergo ICA separately.
#' @param k Integer specifying the number of independent components to extract
#'   from each matrix. Must be at least 2.
#' @param ... Additional arguments passed to \code{\link[ica]{icafast}} from
#'   the \pkg{ica} package.
#'
#' @details
#' For each matrix in \code{x}, the function calls
#' \code{\link[ica]{icafast}} to perform independent component analysis,
#' requesting \code{k} independent components. Let \eqn{S} denote the
#' resulting matrix of independent components (rows = observations,
#' columns = components).
#'
#' The function then constructs a list of \code{k - 1} matrices per data set,
#' where the \eqn{j}-th element (for \eqn{j = 2, \dots, k}) contains the
#' first \eqn{j} independent components, i.e. \code{S[, 1:j]}. This
#' cumulative structure mirrors the behaviour of \code{\link{getPca}} and
#' provides feature sets with increasing dimensionality.
#'
#' When \code{x} is a single matrix or data frame, it is internally wrapped
#' into a list of length one, and the result is still a list of lists for
#' consistency.
#'
#' @return A list of length equal to the number of matrices in \code{x}.
#'   For each input matrix, an element is returned which is itself a list
#'   of length \code{k - 1}. The \eqn{j}-th element (for \eqn{j = 2, \dots, k})
#'   of the inner list is a numeric matrix containing the first \eqn{j}
#'   independent components for that data set.
#'
#' @note
#' The cumulative component representation (2, 3, \dots, \code{k}
#' components) is designed to be used in combination with
#' \code{\link{classifyPcaIca}} for comparing classification performance
#' across different feature dimensions.
#'
#' @references
#' Hyv\"arinen, A. and Oja, E. (2000).
#' Independent component analysis: algorithms and applications.
#' \emph{Neural Networks} \strong{13}(4--5), 411--430.
#'
#' @seealso
#' \code{\link{getPca}},
#' \code{\link{classifyPcaIca}},
#' \code{\link{simulateCNA}},
#' \code{\link[ica]{icafast}}
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
#' ## Obtain ICA-based features
#' ica <- getIca(sim.dat, k = 5)
#'
#' length(ica)
#' length(ica[[1]])
#' sapply(ica[[1]], dim)
#'
#' @author Maharani Ahsani Ummi and Arief Gusnanto
#'
#' @export
#' @importFrom ica icafast
getIca <- function(x, k, ...) {
  
  if (missing(k) || length(k) != 1L || k < 2L) {
    stop("'k' must be a single integer >= 2.")
  }
  k <- as.integer(k)
  
  ## Normalise x to a list of numeric matrices (rows = obs, cols = vars)
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
    stop("'x' must be a matrix, data frame, or a list of such objects.")
  }
  
  n.set <- length(x_list)
  if (n.set == 0L) {
    stop("'x' must contain at least one matrix.")
  }
  
  all.res <- vector("list", n.set)
  
  for (i in seq_len(n.set)) {
    Xi <- x_list[[i]]
    
    if (ncol(Xi) < k) {
      stop("For all matrices in 'x', the number of columns must be >= 'k'.")
    }
    
    ## ICA using the 'ica' package
    ica_fit <- icafast(Xi, nc = k, ...)
    S <- ica_fit$S  # source signals: rows = obs, cols = components
    
    ## Build cumulative component matrices: [, 1:2], [, 1:3], ..., [, 1:k]
    comp_list <- vector("list", k - 1L)
    for (j in seq.int(2L, k)) {
      comp_list[[j - 1L]] <- as.matrix(S[, seq_len(j), drop = FALSE])
    }
    
    all.res[[i]] <- comp_list
  }
  
  all.res
}