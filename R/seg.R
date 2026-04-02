#' @export
seg <- function(test.noise, denoise = "CBS", verbose = FALSE) {
  if (!is.matrix(test.noise)) {
    test.noise <- as.matrix(test.noise)
  }
  
  nrow.x <- nrow(test.noise)
  n <- ncol(test.noise)
  
  if (nrow.x == 0L || n == 0L) {
    stop("'test.noise' must have positive number of rows and columns.")
  }
  
  ## Currently only CBS is implemented
  if (!identical(denoise, "CBS")) {
    stop("Currently only 'CBS' is supported for 'denoise'.")
  }
  
  ## Segmentation
  test.seg <- matrix(NA_real_, nrow.x, n)
  
  for (i in seq_len(nrow.x)) {
    test.seg[i, ] <- CBS(test.noise[i, ], chr = rep(1L, n),verbose = verbose)
  }
  
  test.seg
}
