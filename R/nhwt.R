#' @export
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