\name{seg}
\alias{seg}
\title{Segmentation Per Observation Using CBS}
\description{
Applies circular binary segmentation (CBS) to each row of a data matrix.
This function is a convenience wrapper around \code{\link{CBS}} for
performing segmentation per observation.
}
\usage{
seg(test.noise, denoise = "CBS", verbose = FALSE)
}
\arguments{
  \item{test.noise}{
    Numeric matrix or data frame containing the noisy data to be
    segmented. Rows correspond to observations (e.g. samples or
    patients) and columns to genomic locations or predictor variables.
  }
  \item{denoise}{
    Character string specifying the segmentation method. Currently
    only \code{"CBS"} is supported and used to indicate circular
    binary segmentation.
  }
  \item{verbose}{
    Logical. If \code{TRUE}, progress messages from the segmentation
    procedure are shown. If \code{FALSE}, they are suppressed.
  }
}
\details{
For each row of \code{test.noise}, the function calls \code{\link{CBS}}
to perform circular binary segmentation using the \pkg{DNAcopy}
infrastructure. The result is a matrix of segmented values with the
same dimensions as the input, where each row represents the segmented
copy-number (or intensity) profile for that observation.

The \code{denoise} argument is retained for future extensions, where
alternative segmentation methods might be implemented.
}
\value{
A numeric matrix of the same dimension as \code{test.noise}, containing
the segmented profiles for all observations. Each row corresponds to one
observation and each column to one position, as in the input.
}
\references{
Olshen, A. B., Venkatraman, E. S., Lucito, R. and Wigler, M. (2004).
Circular binary segmentation for the analysis of array-based DNA copy
number data.
\emph{Biostatistics} \bold{5}(4), 557--572.
}
\author{
Maharani Ummi \email{maharaniahsani@gmail.com}
}
\seealso{
\code{\link{CBS}},
\code{\link{wavFeatExt}}
}
\examples{
set.seed(10)

cn <- replicate(
  10,
  rnorm(100, mean = c(rep(2, 50), rep(1, 50)), sd = 0.4)
)
cn <- t(cn)

cn.seg <- seg(cn, denoise = "CBS", verbose = FALSE)

dim(cn.seg)

par(mfrow = c(1, 2))
plot(cn[1, ], type = "p", main = "Observation 1",
     xlab = "Position", ylab = "Copy-number")
lines(cn.seg[1, ], col = 2, lwd = 2)

plot(cn[2, ], type = "p", main = "Observation 2",
     xlab = "Position", ylab = "Copy-number")
lines(cn.seg[2, ], col = 2, lwd = 2)
}
\keyword{segmentation}