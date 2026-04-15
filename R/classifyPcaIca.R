#' Classification on PCA- and ICA-based Feature Sets
#'
#' Performs classification using several machine learning methods on
#' feature sets derived from principal component analysis (PCA) and
#' independent component analysis (ICA), as well as on the original
#' (or segmented) data matrix. Supported methods include Lasso,
#' elastic net, random forest, neural network, partial least squares,
#' and k-nearest neighbours.
#'
#' @param data A list of numeric matrices, each with observations in rows
#'   and variables in columns, or an object with equivalent structure
#'   (e.g. the output from \code{\link{simulateCNA}}). Each element of
#'   the list corresponds to one simulated data set or replicate.
#' @param y Response vector for the classification task, of length equal
#'   to the number of rows in \code{data[[1]]}. It must be either a binary
#'   numeric vector coded as \code{0} and \code{1} or a binary factor with
#'   two levels.
#' @param pca Result of applying \code{\link{getPca}} to \code{data}, i.e.
#'   a list of PCA-based feature sets for each data set. For each data set,
#'   the corresponding element is a list of cumulative principal components.
#' @param ica Result of applying \code{\link{getIca}} to \code{data}, i.e.
#'   a list of ICA-based feature sets for each data set. For each data set,
#'   the corresponding element is a list of cumulative independent components.
#' @param method Character string specifying the classification method to use.
#'   One of \code{"lasso"}, \code{"elnet"}, \code{"RF"}, \code{"NN"},
#'   \code{"PLS"}, or \code{"KNN"}.
#' @param k Number of folds for k-fold cross-validation. Default is 5.
#'   The number of observations must be divisible by \code{k}.
#' @param ite Number of cross-validation replications. Default is
#'   \code{length(data)}, i.e. one replication per data set in the list.
#'
#' @details
#' For each replication, the function constructs a k-fold cross-validation
#' partition of the observations in \code{data[[1]]}. For every feature set
#' considered (each cumulative PCA feature set, each cumulative ICA feature
#' set, and the original/segmented data matrix), the chosen classification
#' method is trained on the training folds and evaluated on the held-out fold.
#'
#' The cross-validated misclassification error (CE) and area under the ROC
#' curve (AUC) are computed for each fold and then averaged across folds.
#' This procedure is repeated for \code{ite} replications, and the resulting
#' averages are stored in matrices for CE and AUC.
#'
#' The columns of the output matrices correspond to the feature sets used:
#' the first \code{length(pca[[1]])} columns to PCA-based feature sets,
#' the next \code{length(ica[[1]])} columns to ICA-based feature sets,
#' and the final column to the original (or segmented) data matrix
#' (labelled \code{"seg"}).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{CE}{Numeric matrix of cross-validated misclassification errors.
#'   Rows correspond to replications (of length \code{ite}) and columns
#'   to feature sets.}
#'   \item{AUC}{Numeric matrix of cross-validated areas under the ROC curve.
#'   The dimensions and column names match those of \code{CE}.}
#'   \item{method}{Character string giving the classification method used.}
#' }
#'
#' @references
#' Tibshirani, R. (1996).
#' Regression shrinkage and selection via the Lasso.
#' \emph{Journal of the Royal Statistical Society, Series B}
#' \strong{58}, 267--288.
#'
#' Breiman, L. (2001).
#' Random forests.
#' \emph{Machine Learning} \strong{45}(1), 5--32.
#'
#' @seealso
#' \code{\link{simulateCNA}},
#' \code{\link{getPca}},
#' \code{\link{getIca}},
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
#' pca <- getPca(sim.dat, k = 4)
#' ica <- getIca(sim.dat, k = 4)
#'
#' y <- factor(rep(c("Group1", "Group2"), each = 10))
#'
#' res <- classifyPcaIca(
#'   sim.dat,
#'   y,
#'   pca,
#'   ica,
#'   method = "KNN",
#'   k = 5,
#'   ite = 1
#' )
#'
#' res
#'
#' @author Maharani Ahsani Ummi and Arief Gusnanto
#'
#' @export
classifyPcaIca <- function(data, y, pca, ica,
                           method = c("lasso", "elnet", "RF", "NN", "PLS", "KNN"),
                           k = 5,
                           ite = length(data)) {
  # Classification on PCA/ICA-based features using multiple ML methods
  # data : list of matrices (rows = obs, cols = vars)
  # y    : binary response (0/1 numeric or factor with 2 levels)
  # pca  : list as returned by get.pca()
  # ica  : list as returned by get.ica()
  # k    : number of CV folds
  # ite  : number of CV replications
  
  method <- match.arg(method)
  
  ## --- Basic structure checks ---
  if (!is.list(data) || length(data) < 1L) {
    stop("'data' must be a non-empty list of matrices.")
  }

  
  n.data <- nrow(data[[1]])
  if (n.data < k) {
    stop("Number of observations must be at least equal to 'k'.")
  }
  
  if (length(y) != n.data) {
    stop("Length of 'y' must match the number of rows in 'data[[1]]'.")
  }
  
  ## --- y: binary response (both factor & numeric versions) ---
  if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      stop("'y' must be a binary factor (2 levels).")
    }
    y_fac <- y
    y_bin <- as.numeric(y) - 1L  # 0 / 1
  } else {
    y <- as.numeric(y)
    vals <- sort(unique(y))
    if (!all(vals %in% c(0, 1))) {
      stop("'y' must be either a 0/1 numeric vector or a binary factor.")
    }
    y_bin <- y
    y_fac <- factor(y_bin, levels = c(0, 1))
  }
  
  ## --- PCA & ICA structures ---
  if (!is.list(pca) || !is.list(ica)) {
    stop("'pca' and 'ica' must be lists as returned by get.pca() and get.ica().")
  }
  
  n.sim <- length(data)
  if (length(pca) != n.sim || length(ica) != n.sim) {
    stop("Lengths of 'data', 'pca', and 'ica' must be the same.")
  }
  
  l.pca <- length(pca[[1]])
  l.ica <- length(ica[[1]])
  n.features <- l.pca + l.ica + 1L  # PCA, ICA, seg/original
  
  if (!is.numeric(k) || k < 2) {
    stop("'k' must be an integer >= 2.")
  }
  k <- as.integer(k)
  
  if (n.data %% k != 0L) {
    stop("'nrow(data[[1]])' must be divisible by 'k' (for equal-sized folds).")
  }
  fold_size <- n.data / k
  
  ite <- as.integer(ite)
  if (ite < 1L) {
    stop("'ite' must be a positive integer.")
  }
  
  ## --- Result containers ---
  all.mce  <- matrix(NA_real_, nrow = ite, ncol = n.features)
  all.mauc <- matrix(NA_real_, nrow = ite, ncol = n.features)
  
  for (i in seq_len(ite)) {
    ## choose which data set to use in this replication
    idx.sim <- if (n.sim == 1L) 1L else ((i - 1L) %% n.sim) + 1L
    
    ## random k-fold assignment for this iteration
    fold_id <- integer(n.data)
    
    for (lev in levels(y_fac)) {
      idx <- which(y_fac == lev)
      idx <- sample(idx)
      fold_id[idx] <- rep(seq_len(k), length.out = length(idx))
    }
    
    mce  <- numeric(n.features)
    mauc <- numeric(n.features)
    
    for (s in seq_len(n.features)) {
      ## --- Select feature set ---
      if (s <= l.pca) {
        x <- as.matrix(pca[[idx.sim]][[s]])
      } else if (s > l.pca && s <= l.pca + l.ica) {
        x <- as.matrix(ica[[idx.sim]][[s - l.pca]])
      } else {
        x <- if (n.sim == 1L) as.matrix(unname(data[[1]])) else as.matrix(data[[idx.sim]])
      }
      
      ce_vec  <- numeric(k)
      auc_vec <- numeric(k)
      
      for (f in seq_len(k)) {
        test.idx  <- which(fold_id == f)
        train.idx <- which(fold_id != f)
        
        x.train <- x[train.idx, , drop = FALSE]
        x.test  <- x[test.idx, , drop = FALSE]
        
        y.train_fac <- y_fac[train.idx]
        y.test_fac  <- y_fac[test.idx]
        
        y.train_bin <- y_bin[train.idx]
        y.test_bin  <- y_bin[test.idx]
        
        data.train <- data.frame(y = y.train_bin, x.train)
        
        ## --- Fit and predict ---
        if (method == "lasso") {
          
          fit <- glmnet::cv.glmnet(x.train, y.train_bin,
                                   family = "binomial",
                                   type.measure = "class",
                                   alpha = 1)
          pred_prob <- as.numeric(
            predict(fit, s = fit$lambda.1se,
                    newx = x.test, type = "response")
          )
          pred_class <- ifelse(pred_prob > 0.5, 1, 0)
          ce_vec[f] <- mean(pred_class != y.test_bin)
          
        } else if (method == "elnet") {
          
          fit <- glmnet::cv.glmnet(x.train, y.train_bin,
                                   family = "binomial",
                                   type.measure = "class",
                                   alpha = 0.5)
          pred_prob <- as.numeric(
            predict(fit, s = fit$lambda.1se,
                    newx = x.test, type = "response")
          )
          pred_class <- ifelse(pred_prob > 0.5, 1, 0)
          ce_vec[f] <- mean(pred_class != y.test_bin)
          
        } else if (method == "RF") {
          
          fit <- randomForest::randomForest(x = x.train, y = y.train_fac)
          rf.class <- predict(fit, newdata = x.test)
          rf.vote  <- predict(fit, newdata = x.test, type = "vote")
          ## probability of the second level (class "1")
          pred_prob <- rf.vote[, 2]
          ce_vec[f] <- mean(rf.class != y.test_fac)
          
        } else if (method == "NN") {
          
          data.train.nn <- data.frame(y = y.train_bin, x.train)
          fit <- neuralnet::neuralnet(y ~ ., data.train.nn)
          nn.pred <- neuralnet::compute(fit, x.test)
          pred_prob <- as.numeric(nn.pred$net.result[, 1])
          pred_class <- ifelse(pred_prob > 0.5, 1, 0)
          ce_vec[f] <- mean(pred_class != y.test_bin)
          
        } else if (method == "KNN") {
          
          fit <- caret::knn3(x.train, y.train_fac, k = 3)
          knn.class <- predict(fit, x.test, type = "class")
          knn.prob  <- predict(fit, x.test, type = "prob")
          ## probability of second level (class "1")
          pred_prob <- knn.prob[, 2]
          ce_vec[f] <- mean(knn.class != y.test_fac)
          
        } else if (method == "PLS") {
          
          fit <- caret::plsda(x.train, y.train_fac, ncomp = 1)
          pls.class <- predict(fit, x.test, type = "class")
          pls.prob  <- predict(fit, x.test, type = "prob")
          ## probability of second level (class "1")
          pred_prob <- pls.prob[, 2, 1]
          ce_vec[f] <- mean(pls.class != y.test_fac)
          
        }
        
        ## --- AUC ---
        roc.obj <- pROC::roc(y.test_bin, pred_prob, quiet = TRUE)
        auc_vec[f] <- as.numeric(pROC::auc(roc.obj))
      }
      
      mce[s]  <- mean(ce_vec)
      mauc[s] <- mean(auc_vec)
    }
    
    all.mce[i, ]  <- mce
    all.mauc[i, ] <- mauc
  }
  
  ## Column names (note: s = 1..l.pca are cumulative PC feature sets)
  pca_names <- paste0("PC", seq_len(l.pca))
  ica_names <- paste0("I",  seq_len(l.ica) + 1L)
  
  colnames(all.mce)  <- c(pca_names, ica_names, "seg")
  colnames(all.mauc) <- c(pca_names, ica_names, "seg")
  
  res <- list(
    CE     = all.mce,
    AUC    = all.mauc,
    method = method
  )

  class(res) <- "pcaIcaClassifier"
  res
}
