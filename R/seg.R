#' Segmentation Per Observation Using CBS
#'
#' Applies circular binary segmentation (CBS) to each row of a data matrix.
#' This function is a convenience wrapper around \code{\link{CBS}} for
#' performing segmentation per observation.
#'
#' @param test.noise Numeric matrix or data frame containing the noisy data
#'   to be segmented. Rows correspond to observations (e.g. samples or
#'   patients) and columns to genomic locations or predictor variables.
#' @param denoise Character string specifying the segmentation method.
#'   Currently only \code{"CBS"} is supported.
#' @param verbose Logical. If \code{TRUE}, progress messages from the
#'   segmentation procedure are shown. If \code{FALSE}, they are suppressed.
#'
#' @details
#' For each row of \code{test.noise}, the function calls \code{\link{CBS}}
#' to perform circular binary segmentation using the \pkg{DNAcopy}
#' infrastructure. The result is a matrix of segmented values with the
#' same dimensions as the input, where each row represents the segmented
#' copy-number (or intensity) profile for that observation.
#'
#' The \code{denoise} argument is retained for future extensions, where
#' alternative segmentation methods might be implemented.
#'
#' @return A numeric matrix of the same dimension as \code{test.noise},
#'   containing the segmented profiles for all observations. Each row
#'   corresponds to one observation and each column to one position.
#'
#' @references
#' Olshen, A. B., Venkatraman, E. S., Lucito, R. and Wigler, M. (2004).
#' Circular binary segmentation for the analysis of array-based DNA copy
#' number data. \emph{Biostatistics} \strong{5}(4), 557--572.
#'
#' @seealso
#' \code{\link{CBS}},
#' \code{\link{wavFeatExt}}
#'
#' @examples
#' set.seed(10)
#'
#' cn <- replicate(
#'   10,
#'   rnorm(100, mean = c(rep(2, 50), rep(1, 50)), sd = 0.4)
#' )
#' cn <- t(cn)
#'
#' cn.seg <- seg(cn, denoise = "CBS", verbose = FALSE)
#'
#' dim(cn.seg)
#'
#' par(mfrow = c(1, 2))
#' plot(cn[1, ], type = "p", main = "Observation 1",
#'      xlab = "Position", ylab = "Copy-number")
#' lines(cn.seg[1, ], col = 2, lwd = 2)
#'
#' plot(cn[2, ], type = "p", main = "Observation 2",
#'      xlab = "Position", ylab = "Copy-number")
#' lines(cn.seg[2, ], col = 2, lwd = 2)
#'
#' @author Maharani Ummi \email{maharaniahsani@@gmail.com}
#'
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
