#' Feature Extraction Using Non-decimated Haar Wavelet Transform
#'
#' Performs feature extraction from copy-number (or other one-dimensional)
#' profiles using the non-decimated (stationary) Haar wavelet transform.
#' The function computes detail or scaling coefficients at multiple scales
#' for each observation and organises them into a structured list.
#'
#' @param data Either a numeric matrix, a list of numeric matrices (with rows
#'   corresponding to observations and columns to variables or genomic
#'   locations), or the output from \code{\link{simulateCNA}}. Each matrix
#'   represents a set of profiles for which wavelet-based features are to be
#'   extracted.
#' @param type Character string indicating which type of coefficients to
#'   extract from the non-decimated wavelet transform. One of \code{"detail"}
#'   or \code{"scaling"}.
#'
#' @details
#' For each matrix in \code{data}, the function applies a non-decimated
#' Haar wavelet transform (via \code{\link{nhwt}}) to each row. This yields,
#' for each scale, a matrix of wavelet coefficients of the same dimension
#' as the original data matrix.
#'
#' The input \code{data} is first normalised to a list of matrices. If a
#' single matrix is provided, it is internally wrapped into a list of length
#' one. When the input is the output of \code{\link{simulateCNA}}, each
#' element corresponds to one simulated data set.
#'
#' The resulting object is a list of length equal to the number of matrices
#' in the input. Each element of this list is itself a list of length equal
#' to the number of available scales. At each scale, a numeric matrix
#' contains the wavelet coefficients (detail or scaling) for all observations.
#'
#' @return A list of length equal to the number of matrices in \code{data}.
#'   For each data matrix, an element is returned which is a list of length
#'   \eqn{J}, where \eqn{J} is the number of scales determined by the
#'   wavelet transform. Each sub-list element is a numeric matrix of the
#'   same dimension as the original data matrix.
#'
#' @references
#' Nason, G. P. (2008).
#' \emph{Wavelet Methods in Statistics with R}. Springer.
#'
#' @seealso
#' \code{\link{nhwt}},
#' \code{\link{simulateCNA}},
#' \code{\link{classifyWavFeatExt}},
#' \code{\link{classifyPcaIca}}
#'
#' @examples
#' set.seed(10)
#'
#' sim.dat <- simulateCNA(
#'   n.obs = 20,
#'   p = 32,
#'   n.sim = 1,
#'   n.block = 8,
#'   verbose = FALSE
#' )
#'
#' det.coef <- wavFeatExt(sim.dat, type = "detail")
#' sca.coef <- wavFeatExt(sim.dat, type = "scaling")
#'
#' length(det.coef)
#' length(sca.coef)
#'
#' @author Maharani Ahsani Ummi and Arief Gusnanto
#'
#' @export
wavFeatExt <- function(data, type = c("detail", "scaling")) {
  # Feature extraction using non-decimated Haar wavelet transform (NHWT)
  # data: matrix, list of matrices, or output from sim.CNA()
  # type: "detail" or "scaling"
  
  type <- match.arg(type)
  
  ## Normalise input to a list of matrices
  if (is.matrix(data)) {
    data_list <- list(data)
  } else if (is.data.frame(data)) {
    data_list <- list(as.matrix(data))
  } else if (is.list(data)) {
    data_list <- lapply(data, function(x) {
      if (!is.matrix(x)) as.matrix(x) else x
    })
  } else {
    stop("'data' must be a matrix, data frame, or a list of matrices.")
  }
  
  n.sim <- length(data_list)
  if (n.sim == 0L) {
    stop("'data' must contain at least one data matrix.")
  }
  
  ## Apply NHWT to each data matrix
  res <- vector("list", n.sim)
  for (j in seq_len(n.sim)) {
    res[[j]] <- nhwt(data_list[[j]], type = type)
  }
  
  res
}
