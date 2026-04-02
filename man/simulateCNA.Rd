\name{simulateCNA}
\alias{simulateCNA}
\title{
Simulate Copy-number Alteration (CNA) Data in a Two-group Setting
}
\description{
Generates simulated copy-number alteration (CNA) data for a two-group
setting, under a multivariate normal model with block correlation
structure. Mean differences between the two groups are introduced in
selected blocks, and each simulated profile is subsequently segmented
using circular binary segmentation (CBS).
}
\usage{
simulateCNA(n.obs = 100, p = 1024, n.sim = 1, effect.diff = 1,
            n.block = 32, true.rho = 0.9, block.cor = 0.4,
            block.diff = "A", true.mu = NULL, verbose = FALSE)
}
\arguments{
  \item{n.obs}{
    Number of observations in each simulated data set. Must be even,
    so that the first half belong to group 1 and the second half to
    group 2.
  }
  \item{p}{
    Number of variables (genomic locations) in each simulated data set.
    It is recommended to use a power of two (e.g. 512, 1024, 2048)
    for compatibility with subsequent wavelet-based analyses.
  }
  \item{n.sim}{
    Number of data sets to generate. Each data set has \code{n.obs}
    rows and \code{p} columns.
  }
  \item{effect.diff}{
    Magnitude of the mean difference between the two groups in the
    blocks where differences are present (see \code{block.diff}).
  }
  \item{n.block}{
    Number of correlation blocks along the genome. Must be a divisor
    of \code{p}, so that the block size \code{p / n.block} is an
    integer.
  }
  \item{true.rho}{
    Correlation between genomic regions within a block.
  }
  \item{block.cor}{
    Correlation between adjacent blocks. A value of zero implies that
    blocks are independent from each other.
  }
  \item{block.diff}{
    Pattern of blocks that carry mean differences between the two
    groups. Currently one of \code{"A"} or \code{"B"}.
  }
  \item{true.mu}{
    Optional numeric vector of length \code{p} giving the baseline
    mean level of CNA across genomic locations, common to both groups.
    If \code{NULL} (default), a simple four-block mean pattern is
    used and repeated along the genome.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, progress messages are shown during
    simulation and segmentation. If \code{FALSE}, they are suppressed.
  }
}
\details{
The function first constructs a block-structured covariance matrix
\eqn{\Sigma} of dimension \code{p x p}, with within-block correlations
given by \code{true.rho} and between-block correlations given by
\code{block.cor}. Copy-number data are then simulated from a
multivariate normal distribution with mean vector \code{true.mu} and
covariance \eqn{\Sigma}.

To create a two-group setting, a mean shift vector is added to the
first half of the observations (group 1), while the second half
(group 2) remains at the baseline mean. The pattern of mean shifts
depends on \code{block.diff}.

Finally, for each simulated data set, circular binary segmentation
(CBS) is applied to every observation (row) via \code{\link{seg}}
and \code{\link{CBS}}. The returned data therefore correspond to
segmented CNA profiles rather than raw multivariate normal values.
}
\value{
A list of length \code{n.sim}. Each element is a numeric matrix of
dimension \code{n.obs x p}, containing the segmented CNA profiles for
one simulated data set.
}
\author{
Arief Gusnanto
}
\seealso{
\code{\link{seg}},
\code{\link{CBS}},
\code{\link{wavFeatExt}},
\code{\link{classifyPcaIca}},
\code{\link{classifyWavFeatExt}}
}
\examples{
set.seed(10)

sim.dat <- simulateCNA(
  n.obs = 20,
  p = 32,
  n.sim = 1,
  n.block = 8,
  verbose = FALSE
)

length(sim.dat)
dim(sim.dat[[1]])

plot(sim.dat[[1]][1, ], type = "l",
     xlab = "Genomic location", ylab = "Segmented CNA",
     main = "Simulated segmented CNA profile")
}
\keyword{datagen}