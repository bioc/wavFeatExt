#' Non-decimated Haar Wavelet Transform (Legacy)
#'
#' Performs a non-decimated Haar wavelet transform (stationary wavelet
#' transform) on one or more one-dimensional signals. This function is
#' kept for compatibility with older code and is not used in the main
#' workflow of the package.
#'
#' @param data Numeric vector, matrix, or data frame containing the data
#'   to be decomposed. If a matrix or data frame is supplied, rows are
#'   interpreted as observations and columns as ordered locations.
#' @param type Character string indicating which coefficients to extract
#'   from the transform. Use \code{"detail"} for detail coefficients or
#'   \code{"scaling"} for scaling coefficients.
#'
#' @details
#' The function applies a non-decimated (stationary) Haar wavelet transform
#' to each row of the input. If the number of columns is not a power of two,
#' the series is extended to the next power of two by constant padding with
#' value 1 on both sides before computing the transform. After the transform,
#' the resulting coefficients are re-aligned and cropped back to the original
#' length.
#'
#' The underlying transform is computed using \code{\link[wavethresh]{wd}}
#' with \code{family = "DaubExPhase"}, \code{filter.number = 1}, and
#' \code{type = "station"}. For each scale, either detail or scaling
#' coefficients can be extracted via \code{\link[wavethresh]{accessD}} or
#' \code{\link[wavethresh]{accessC}}, respectively.
#'
#' @return A list of length \eqn{J}, where \eqn{J = \lfloor \log_2(n) \rfloor}
#'   and \eqn{n} is the number of columns in \code{data}. Each element is a
#'   numeric matrix of dimension \code{nrow(data) x n} containing the detail
#'   or scaling coefficients at a given scale.
#'
#' @note
#' This function is considered legacy and is not used by the main feature
#' extraction and classification functions in the package. It is provided
#' for backward compatibility with older analysis pipelines.
#'
#' @references
#' Nason, G. P. (2008).
#' \emph{Wavelet Methods in Statistics with R}. Springer.
#'
#' @seealso
#' \code{\link{wavFeatExt}},
#' \code{\link{plot.nhwt}},
#' \code{\link[wavethresh]{wd}},
#' \code{\link[wavethresh]{accessD}},
#' \code{\link[wavethresh]{accessC}}
#'
#' @examples
#' ## Simple example: non-decimated Haar wavelet coefficients
#' ## for a piecewise constant signal
#' obj <- c(rep(1, 10), rep(3, 20))
#'
#' ## Detail coefficients at all scales
#' obj.nhwt.det <- nhwt(obj, type = "detail")
#'
#' length(obj.nhwt.det)
#' sapply(obj.nhwt.det, dim)
#'
#' @author Maharani Ahsani Ummi and Arief Gusnanto
#'
#' @export
#' @importFrom wavethresh wd accessD accessC
nhwt <- function(data, type = c("detail", "scaling")) {
  # Non-decimated Haar wavelet transform (stationary transform)
  # Legacy function, kept for compatibility with older code.
  
  type <- match.arg(type)
  
  ## Coerce to matrix like original version:
  ## - if vector: 1 x n (1 sample, n locations)
  ## - if matrix/data.frame: rows = samples, cols = locations
  if (is.null(nrow(data))) {
    data <- t(as.matrix(data))  # vector -> 1 x n
  } else if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  n.sampl <- nrow(data)
  n <- ncol(data)
  
  if (n.sampl <= 0L || n <= 0L) {
    stop("'data' must have positive number of rows and columns.")
  }
  if (n < 2L) {
    stop("'data' must have at least 2 columns/locations.")
  }
  
  ## Extend to next power of two by constant padding with value 1
  n.up <- ceiling(log2(n))
  n.extend <- 2^n.up
  n.res <- n.extend - n
  
  if (n.res == 0L) {
    data.extend <- data
    left <- 0L
  } else {
    data.extend <- matrix(NA_real_, nrow = n.sampl, ncol = n.extend)
    left <- floor(n.res / 2)
    right <- n.res - left
    for (i in seq_len(n.sampl)) {
      data.extend[i, ] <- c(rep(1, left), data[i, ], rep(1, right))
    }
  }
  
  if (n.up <= 1L) {
    stop("Series is too short for a multi-scale non-decimated transform.")
  }
  
  coef.ndwt <- vector("list", n.up - 1L)
  data.ndwt <- vector("list", n.up - 1L)
  
  for (i in seq_len(n.up - 1L)) {
    coef.ndwt[[i]] <- matrix(NA_real_, nrow = n.sampl, ncol = n.extend)
    data.ndwt[[i]] <- matrix(NA_real_, nrow = n.sampl, ncol = n)
  }
  
  scale <- 1L
  for (k in rev(seq_len(n.up - 1L))) {
    for (i in seq_len(n.sampl)) {
      wt <- wd(
        data.extend[i, ],
        filter.number = 1,
        family = "DaubExPhase",
        type = "station"
      )
      
      temp <- if (type == "detail") {
        accessD(wt, level = k)
      } else {
        accessC(wt, level = k)
      }
      
      shift <- 2^(scale - 1L) - 1L
      lead <- 2^(scale - 1L)
      
      if (n.res == 0L) {
        coef.ndwt[[scale]][i, (1L + shift):n] <-
          temp[seq_len(n - shift)]
        coef.ndwt[[scale]][i, seq_len(lead)] <-
          temp[(n - shift):n]
        data.ndwt[[scale]][i, ] <- coef.ndwt[[scale]][i, ]
      } else {
        coef.ndwt[[scale]][i, (1L + shift):n.extend] <-
          temp[seq_len(n.extend - shift)]
        coef.ndwt[[scale]][i, seq_len(lead)] <-
          temp[(n.extend - shift):n.extend]
        
        start.idx <- left + 1L
        end.idx <- left + n
        data.ndwt[[scale]][i, ] <-
          coef.ndwt[[scale]][i, start.idx:end.idx]
      }
    }
    scale <- scale + 1L
  }
  
  class(data.ndwt) <- "nhwt"
  data.ndwt
}