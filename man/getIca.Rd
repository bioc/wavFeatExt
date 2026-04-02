\name{getIca}
\alias{getIca}
\title{
Independent Component Analysis (ICA) for a List of Matrices
}
\description{
Applies independent component analysis (ICA) to each matrix in a
collection of data sets and returns cumulative independent components
for use in downstream classification or feature extraction.
}
\usage{
getIca(x, k, ...)
}
\arguments{
  \item{x}{
    Either a numeric matrix, a list of numeric matrices (with rows
    corresponding to observations and columns to variables), or the
    output from \code{\link{simulateCNA}}. Each matrix will undergo ICA
    separately.
  }
  \item{k}{
    Integer specifying the number of independent components to extract
    from each matrix. Must be at least 2.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[ica]{icafast}} from the
    \pkg{ica} package.
  }
}
\details{
For each matrix in \code{x}, the function calls \code{\link[ica]{icafast}}
to perform independent component analysis, requesting \code{k}
independent components. Let \eqn{S} denote the resulting matrix of
independent components (rows = observations, columns = components).

The function then constructs a list of \code{k - 1} matrices per data
set, where the \eqn{j}-th element (for \eqn{j = 2, \dots, k}) contains
the first \eqn{j} independent components, i.e. \code{S[, 1:j]}. This
cumulative structure mirrors the behaviour of \code{\link{getPca}} and
provides feature sets with increasing dimensionality (2 components,
3 components, and so on).

When \code{x} is a single matrix or data frame, it is internally
wrapped into a list of length one, and the result is still a list of
lists for consistency.
}
\value{
A list of length equal to the number of matrices in \code{x}. For each
input matrix, an element is returned which is itself a list of length
\code{k - 1}. The \eqn{j}-th element (for \eqn{j = 2, \dots, k}) of
the inner list is a numeric matrix of dimension \code{nrow(x) x j}
containing the first \eqn{j} independent components for that data set.
}
\references{
Hyv\"arinen, A. and Oja, E. (2000).
Independent component analysis: algorithms and applications.
\emph{Neural Networks} \bold{13}(4--5), 411--430.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\note{
The cumulative component representation (2, 3, \dots, \code{k}
components) is designed to be used in combination with
\code{\link{classifyPcaIca}} for comparing classification performance
across different feature dimensions.
}
\seealso{
\code{\link{getPca}},
\code{\link{classifyPcaIca}},
\code{\link{simulateCNA}},
\code{\link[ica]{icafast}}
}
\examples{
set.seed(10)

## Small simulated data for fast execution
sim.dat <- simulateCNA(
  n.obs = 20,
  p = 32,
  n.sim = 1,
  n.block = 8,
  verbose = FALSE
)

## Obtain ICA-based features
ica <- getIca(sim.dat, k = 5)

length(ica)
length(ica[[1]])
sapply(ica[[1]], dim)
}
\keyword{multivariate}
\keyword{dimension reduction}
