#' @export
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
