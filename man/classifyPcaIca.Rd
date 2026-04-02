\name{classifyPcaIca}
\alias{classifyPcaIca}
\title{
Classification on PCA- and ICA-based Feature Sets
}
\description{
Performs classification using several machine learning methods on
feature sets derived from principal component analysis (PCA) and
independent component analysis (ICA), as well as on the original
(or segmented) data matrix. Supported methods include Lasso,
elastic net, random forest, neural network, partial least squares,
and k-nearest neighbours.
}
\usage{
classifyPcaIca(data, y, pca, ica,
               method = c("lasso", "elnet", "RF", "NN", "PLS", "KNN"),
               k = 5,
               ite = length(data))
}
\arguments{
  \item{data}{
    A list of numeric matrices, each with observations in rows and
    variables in columns, or an object with equivalent structure
    (e.g. the output from \code{\link{simulateCNA}}). Each element of
    the list corresponds to one simulated data set or replicate.
  }
  \item{y}{
    Response vector for the classification task, of length equal to
    the number of rows in \code{data[[1]]}. It must be either a
    binary numeric vector coded as \code{0} and \code{1} or a binary
    factor with two levels.
  }
  \item{pca}{
    Result of applying \code{\link{getPca}} to \code{data}, i.e. a
    list of PCA-based feature sets for each data set. For each
    data set, the corresponding element is a list of cumulative
    principal components obtained from \code{\link{getPca}}.
  }
  \item{ica}{
    Result of applying \code{\link{getIca}} to \code{data}, i.e. a
    list of ICA-based feature sets for each data set. For each
    data set, the corresponding element is a list of cumulative
    independent components obtained from \code{\link{getIca}}.
  }
  \item{method}{
    Character string specifying the classification method to use.
    One of \code{"lasso"} for Lasso logistic regression,
    \code{"elnet"} for elastic net, \code{"RF"} for random forest,
    \code{"NN"} for neural network, \code{"PLS"} for partial least
    squares, or \code{"KNN"} for k-nearest neighbours.
  }
  \item{k}{
    Number of folds for k-fold cross-validation. Default is 5.
    The number of observations must be divisible by \code{k}.
  }
  \item{ite}{
    Number of cross-validation replications. The default is
    \code{length(data)}, i.e. one replication per data set in the
    list.
  }
}
\details{
For each replication, the function constructs a k-fold
cross-validation partition of the observations in \code{data[[1]]}.
For every feature set considered (each cumulative PCA feature set,
each cumulative ICA feature set, and the original/segmented data
matrix), the chosen classification method is trained on the training
folds and evaluated on the held-out fold.

The cross-validated misclassification error (CE) and area under the
ROC curve (AUC) are computed for each fold and then averaged across
folds. This procedure is repeated for \code{ite} replications, and
the resulting averages are stored in matrices for CE and AUC.

The columns of the output matrices correspond to the feature sets
used: the first \code{length(pca[[1]])} columns to PCA-based feature
sets, the next \code{length(ica[[1]])} columns to ICA-based feature
sets, and the final column to the original (or segmented) data
matrix (labelled \code{"seg"}).
}
\value{
A list with the following components:
\item{CE}{
  Numeric matrix of cross-validated misclassification errors.
  Rows correspond to replications (of length \code{ite}) and
  columns to feature sets (PCA-based, ICA-based, and the
  original/segmented data).
}
\item{AUC}{
  Numeric matrix of cross-validated areas under the ROC curve.
  The dimensions and column names match those of \code{CE}.
}
\item{method}{
  Character string giving the classification method used.
}
}
\references{
Tibshirani, R. (1996).
Regression shrinkage and selection via the Lasso.
\emph{Journal of the Royal Statistical Society, Series B}
\bold{58}, 267--288.

Breiman, L. (2001).
Random forests.
\emph{Machine Learning} \bold{45}(1), 5--32.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\seealso{
\code{\link{simulateCNA}},
\code{\link{getPca}},
\code{\link{getIca}},
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

pca <- getPca(sim.dat, k = 4)
ica <- getIca(sim.dat, k = 4)

y <- factor(rep(c("Group1", "Group2"), each = 10))

res <- classifyPcaIca(
  sim.dat,
  y,
  pca,
  ica,
  method = "KNN",
  k = 5,
  ite = 1
)

res
}
\keyword{classification}