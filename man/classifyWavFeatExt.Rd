\name{classifyWavFeatExt}
\alias{classifyWavFeatExt}
\title{
Classification on Wavelet-based Features
}
\description{
Performs classification using several machine learning methods on
features extracted via non-decimated Haar wavelet transformation,
namely detail and scaling coefficients, as well as on the original
(or segmented) data matrix. Optionally, all wavelet coefficients can be
combined into a single feature set.
}
\usage{
classifyWavFeatExt(data, y, det, sca,
                   method = c("lasso", "elnet", "RF", "NN", "PLS", "KNN"),
                   k = 5,
                   ite = length(data),
                   all = FALSE, verbose = FALSE)
}
\arguments{
  \item{data}{
    A non-empty list of numeric matrices with equal number of rows
    (rows = observations, columns = variables). Each element of the list
    corresponds to one data set / replicate (simulation).
  }
  \item{y}{
    Response vector for the classification task, of length equal to
    the number of rows in \code{data[[1]]}. It must be either a binary
    numeric vector coded as \code{0} and \code{1}, or a binary factor
    with exactly two levels. Internally, both a factor and a 0/1 numeric
    version are created.
  }
  \item{det}{
    Detail coefficients of \code{data}, typically obtained by
    \code{\link{wavFeatExt}(data, type = "detail")}. Must be a list of
    the same length as \code{data}; for each data set, the corresponding
    element is a list of matrices, one per wavelet scale.
  }
  \item{sca}{
    Scaling coefficients of \code{data}, typically obtained by
    \code{\link{wavFeatExt}(data, type = "scaling")}. Must be a list of
    the same length as \code{data}; for each data set, the corresponding
    element is a list of matrices, one per wavelet scale.
  }
  \item{method}{
    Character string specifying the classification method to use.
    One of \code{"lasso"} (Lasso logistic regression via \pkg{glmnet}),
    \code{"elnet"} (elastic net via \pkg{glmnet}, \code{alpha = 0.5}),
    \code{"RF"} (random forest via \pkg{randomForest}),
    \code{"NN"} (feed-forward neural network via \pkg{neuralnet}),
    \code{"PLS"} (PLS-DA via \pkg{caret} \code{plsda}),
    or \code{"KNN"} (k-nearest neighbours via \pkg{caret} \code{knn3}).
  }
  \item{k}{
    Number of folds for k-fold cross-validation (default \code{5}).
    The number of observations \code{nrow(data[[1]])} must be
    divisible by \code{k} because the implementation uses equal-sized
    folds.
  }
  \item{ite}{
    Number of cross-validation replications (default \code{length(data)}).
    For replication \code{i}, the function uses data set
    \code{((i - 1) \%\% length(data)) + 1} when \code{length(data) > 1}.
    If \code{length(data) == 1}, the single data matrix is reused for
    all replications.
  }
  \item{all}{
    Logical; if \code{FALSE} (default), each wavelet scale is evaluated
    separately (all detail scales, all scaling scales) plus the
    original/segmented matrix (\code{"seg"}).
    If \code{TRUE}, all wavelet coefficients (all detail + all scaling)
    are concatenated into a single feature set (\code{"ALL"}), and the
    original/segmented matrix (\code{"seg"}) is evaluated as the second
    feature set.
  }
  \item{verbose}{
    Logical; if TRUE, progress messages are printed. Default is FALSE.
  }
}
\details{
For each replication, the function randomly assigns observations to
\code{k} folds of equal size. For every feature set considered, the
chosen classification method is trained on the training folds and
evaluated on the held-out fold.

The misclassification error (CE) is computed as the proportion of
incorrect predicted classes on the test fold. In addition, the area
under the ROC curve (AUC) is computed from the predicted class
probabilities using \pkg{pROC}.

Fold-wise CE and AUC are averaged across the \code{k} folds to obtain
one CE and one AUC per feature set and replication. This is repeated
for \code{ite} replications, producing two matrices (CE and AUC).

\strong{Feature sets and column names:}
\itemize{
  \item If \code{all = FALSE}: columns correspond to
    \code{D1, D2, ..., Dq} (detail scales),
    then \code{S1, S2, ..., Sp} (scaling scales),
    and the final column \code{"seg"} (the original/segmented matrix).
    Here \code{q = length(det[[1]])} and \code{p = length(sca[[1]])}.
  \item If \code{all = TRUE}: columns are \code{"ALL"} (all wavelet
    coefficients combined) and \code{"seg"}.
}
}
\value{
An object of class \code{"wavFeatExtClassifier"}, a list with components:
\item{CE}{
  Numeric matrix of cross-validated misclassification errors.
  Rows correspond to replications (\code{ite}) and columns to feature
  sets (see Details for column naming).
}
\item{AUC}{
  Numeric matrix of cross-validated AUC values.
  Dimensions and column names match those of \code{CE}.
}
\item{method}{
  Character string giving the classification method used.
}
\item{all}{
  Logical indicating whether \code{all = TRUE} was used.
}
}
\references{
Nason, G. P. (2008).
\emph{Wavelet Methods in Statistics with R}. Springer.

Tibshirani, R. (1996).
Regression shrinkage and selection via the Lasso.
\emph{Journal of the Royal Statistical Society, Series B}
\bold{58}, 267--288.
}
\author{
Maharani Ahsani Ummi and Arief Gusnanto
}
\seealso{
\code{\link{simulateCNA}},
\code{\link{wavFeatExt}},
\code{\link{classifyPcaIca}},
\code{\link{nhwt}}
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

## Wavelet coefficients
det.coef <- wavFeatExt(sim.dat, type = "detail")
sca.coef <- wavFeatExt(sim.dat, type = "scaling")

## Binary response
y <- factor(rep(c("Group1", "Group2"), each = 10))

## Classification
res <- classifyWavFeatExt(
  sim.dat,
  y,
  det.coef,
  sca.coef,
  method = "KNN",
  k = 5,
  ite = 1
)

res
}
\keyword{classification}
\keyword{wavelets}
