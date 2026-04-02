#' @export
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