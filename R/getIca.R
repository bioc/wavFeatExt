#' @export
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