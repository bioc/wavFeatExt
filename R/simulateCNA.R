#' Simulate Copy-number Alteration (CNA) Data in a Two-group Setting
#'
#' Generates simulated copy-number alteration (CNA) data for a two-group
#' setting, under a multivariate normal model with block correlation
#' structure. Mean differences between the two groups are introduced in
#' selected blocks, and each simulated profile is subsequently segmented
#' using circular binary segmentation (CBS).
#'
#' @param n.obs Number of observations in each simulated data set. Must be
#'   even, so that the first half belong to group 1 and the second half to
#'   group 2.
#' @param p Number of variables (genomic locations) in each simulated data
#'   set. It is recommended to use a power of two for compatibility with
#'   subsequent wavelet-based analyses.
#' @param n.sim Number of data sets to generate. Each data set has
#'   \code{n.obs} rows and \code{p} columns.
#' @param effect.diff Magnitude of the mean difference between the two groups
#'   in the blocks where differences are present.
#' @param n.block Number of correlation blocks along the genome. Must be a
#'   divisor of \code{p}, so that the block size \code{p / n.block} is an
#'   integer.
#' @param true.rho Correlation between genomic regions within a block.
#' @param block.cor Correlation between adjacent blocks. A value of zero
#'   implies that blocks are independent from each other.
#' @param block.diff Pattern of blocks that carry mean differences between
#'   the two groups. Currently one of \code{"A"} or \code{"B"}.
#' @param true.mu Optional numeric vector of length \code{p} giving the
#'   baseline mean level of CNA across genomic locations, common to both
#'   groups. If \code{NULL}, a simple four-block mean pattern is used and
#'   repeated along the genome.
#' @param verbose Logical. If \code{TRUE}, progress messages are shown during
#'   simulation and segmentation. If \code{FALSE}, they are suppressed.
#'
#' @details
#' The function first constructs a block-structured covariance matrix
#' \eqn{\Sigma} of dimension \code{p x p}, with within-block correlations
#' given by \code{true.rho} and between-block correlations given by
#' \code{block.cor}. Copy-number data are then simulated from a multivariate
#' normal distribution with mean vector \code{true.mu} and covariance
#' \eqn{\Sigma}.
#'
#' To create a two-group setting, a mean shift vector is added to the first
#' half of the observations (group 1), while the second half (group 2)
#' remains at the baseline mean. The pattern of mean shifts depends on
#' \code{block.diff}.
#'
#' Finally, for each simulated data set, circular binary segmentation (CBS)
#' is applied to every observation (row) via \code{\link{seg}} and
#' \code{\link{CBS}}. The returned data therefore correspond to segmented
#' CNA profiles rather than raw multivariate normal values.
#'
#' @return A list of length \code{n.sim}. Each element is a numeric matrix of
#'   dimension \code{n.obs x p}, containing the segmented CNA profiles for
#'   one simulated data set.
#'
#' @seealso
#' \code{\link{seg}},
#' \code{\link{CBS}},
#' \code{\link{wavFeatExt}},
#' \code{\link{classifyPcaIca}},
#' \code{\link{classifyWavFeatExt}}
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
#' length(sim.dat)
#' dim(sim.dat[[1]])
#'
#' plot(sim.dat[[1]][1, ], type = "l",
#'      xlab = "Genomic location", ylab = "Segmented CNA",
#'      main = "Simulated segmented CNA profile")
#'
#' @author Arief Gusnanto
#'
#' @export
#' @importFrom MASS mvrnorm
simulateCNA <- function(n.obs      = 100,
                    p          = 1024,
                    n.sim      = 1,
                    effect.diff = 1,
                    n.block    = 32,
                    true.rho   = 0.9,
                    block.cor  = 0.4,
                    block.diff = "A",
                    true.mu    = NULL,
                    verbose = FALSE) {
  # Function to generate CNA data
  # Author: Arief Gusnanto (a.gusnanto@leeds.ac.uk)
  
  ## Basic checks
  n.obs <- as.integer(n.obs)
  p     <- as.integer(p)
  n.sim <- as.integer(n.sim)
  n.block <- as.integer(n.block)
  
  if (n.obs <= 0L || p <= 0L) {
    stop("'n.obs' and 'p' must be positive integers.")
  }
  if (n.obs %% 2L != 0L) {
    stop("Number of observations 'n.obs' must be even (two-group setting).")
  }
  if (n.sim <= 0L) {
    stop("'n.sim' must be a positive integer.")
  }
  if (!block.diff %in% c("A", "B")) {
    stop("Argument 'block.diff' must be either 'A' or 'B'.")
  }
  if (p %% n.block != 0L) {
    stop("'n.block' must divide 'p' exactly (so that block size is integer).")
  }
  
  size.block <- p / n.block
  
  ## true.mu: mean levels across genomic locations
  if (is.null(true.mu)) {
    ## default mean profile: 4-block pattern repeated
    base.mu <- c(rep(1,   size.block),
                 rep(1.5, size.block),
                 rep(1,   size.block),
                 rep(0.5, size.block))
    true.mu <- rep(base.mu, length.out = p)
  }
  if (length(true.mu) != p) {
    stop("The length of 'true.mu' must be equal to 'p'.")
  }
  
  ## Covariance matrix of the simulated (noisy) data
  true.Sigma <- matrix(0, p, p)
  
  ## Within-block correlations
  temp.within <- matrix(true.rho, size.block, size.block)
  for (k in seq_len(n.block)) {
    idx <- ((k - 1L) * size.block + 1L):(k * size.block)
    true.Sigma[idx, idx] <- temp.within
  }
  
  ## Between-block correlations
  if (n.block > 1L) {
    temp.between <- matrix(block.cor, size.block, size.block)
    for (k in seq_len(n.block - 1L)) {
      idx1 <- ((k - 1L) * size.block + 1L):(k * size.block)
      idx2 <- (k * size.block + 1L):((k + 1L) * size.block)
      true.Sigma[idx2, idx1] <- temp.between
      true.Sigma[idx1, idx2] <- temp.between
    }
  }
  diag(true.Sigma) <- 1
  
  ## True differences between groups
  if (block.diff == "B") {
    true.d1 <- c(
      rep(effect.diff,  size.block),
      rep(-effect.diff, size.block),
      rep(-effect.diff, size.block),
      rep(effect.diff,  size.block),
      rep(0, p - 4L * size.block)
    )
  } else {  # block.diff == "A"
    true.d1 <- c(
      rep(effect.diff,  size.block),
      rep(-effect.diff, size.block),
      rep(0, p - 2L * size.block)
    )
  }
  
  ## Simulate data sets
  sim.dat <- vector("list", n.sim)
  
  for (j in seq_len(n.sim)) {
    if (isTRUE(verbose)) {
      message("Simulating data set ", j)
    }
    Z.sim <- mvrnorm(n.obs, mu = true.mu, Sigma = true.Sigma)
    
    ## Add mean shift to first group
    Z.sim[seq_len(n.obs / 2L), ] <- sweep(Z.sim[seq_len(n.obs / 2L), , drop = FALSE],
                                          2L, true.d1, FUN = "+")
    
    ## Segment each profile using CBS
    sample.CBS <- seg(Z.sim, denoise = "CBS")
    sim.dat[[j]] <- sample.CBS
  }
  
  sim.dat
}
