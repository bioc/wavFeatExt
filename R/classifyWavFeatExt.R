#' @export
classifyWavFeatExt <- function(data, y, det, sca,
                               method = c("lasso", "elnet", "RF", "NN", "PLS", "KNN"),
                               k = 5,
                               ite = length(data),
                               all = FALSE, verbose = FALSE) {
  
  method <- match.arg(method)
  
  ## --- Basic structure checks ---
  if (!is.list(data) || length(data) < 1L) stop("'data' must be a non-empty list of matrices.")
  if (!is.matrix(data[[1]])) stop("'data[[1]]' must be a matrix (rows = observations, cols = variables).")
  
  n.data <- nrow(data[[1]])
  if (n.data < k) stop("Number of observations must be at least equal to 'k'.")
  if (length(y) != n.data) stop("Length of 'y' must match the number of rows in 'data[[1]]'.")
  
  ## --- y: binary response (factor & numeric versions) ---
  if (is.factor(y)) {
    if (nlevels(y) != 2L) stop("'y' must be a binary factor (2 levels).")
    y_fac <- y
    y_bin <- as.numeric(y) - 1L
  } else {
    y <- as.numeric(y)
    vals <- sort(unique(y))
    if (!all(vals %in% c(0, 1))) stop("'y' must be either a 0/1 numeric vector or a binary factor.")
    y_bin <- y
    y_fac <- factor(y_bin, levels = c(0, 1))
  }
  
  ## --- det & sca structures ---
  if (!is.list(det) || !is.list(sca)) stop("'det' and 'sca' must be lists as returned by wavFeatExt().")
  
  n.sim <- length(data)
  if (length(det) != n.sim || length(sca) != n.sim) stop("Lengths of 'data', 'det', and 'sca' must be the same.")
  
  l.det <- length(det[[1]])
  l.sca <- length(sca[[1]])
  
  if (!is.numeric(k) || k < 2) stop("'k' must be an integer >= 2.")
  k <- as.integer(k)
  if (n.data %% k != 0L) stop("'nrow(data[[1]])' must be divisible by 'k' (for equal-sized folds).")
  fold_size <- n.data / k
  
  ite <- as.integer(ite)
  if (ite < 1L) stop("'ite' must be a positive integer.")
  
  combine_all_coef <- function(det_one, sca_one) {
    mats <- c(det_one, sca_one)
    mats <- lapply(mats, function(m) as.matrix(m))
    nr <- vapply(mats, nrow, integer(1))
    if (length(unique(nr)) != 1L) stop("Row counts differ across coefficient matrices.")
    do.call(cbind, mats)
  }
  
  if (!isTRUE(all)) {
    n.features <- l.det + l.sca + 1L
    det_names <- paste0("D", seq_len(l.det))
    sca_names <- paste0("S", seq_len(l.sca))
    feat_names <- c(det_names, sca_names, "seg")
  } else {
    n.features <- 2L
    feat_names <- c("ALL", "seg")
  }
  
  ## --- Result containers ---
  all.mce  <- matrix(NA_real_, nrow = ite, ncol = n.features)
  all.mauc <- matrix(NA_real_, nrow = ite, ncol = n.features)
  
  for (i in seq_len(ite)) {
    
    idx.sim <- if (n.sim == 1L) 1L else ((i - 1L) %% n.sim) + 1L
    fold_id <- sample(rep(seq_len(k), each = fold_size))
    
    mce  <- numeric(n.features)
    mauc <- numeric(n.features)
    
    for (s in seq_len(n.features)) {
      
      ## --- Select feature set ---
      if (!isTRUE(all)) {
        if (s <= l.det) {
          x <- as.matrix(det[[idx.sim]][[s]])
        } else if (s <= l.det + l.sca) {
          x <- as.matrix(sca[[idx.sim]][[s - l.det]])
        } else {
          x <- if (n.sim == 1L) as.matrix(unname(data[[1]])) else as.matrix(data[[idx.sim]])
        }
      } else {
        if (s == 1L) {
          x <- combine_all_coef(det[[idx.sim]], sca[[idx.sim]])
        } else {
          x <- if (n.sim == 1L) as.matrix(unname(data[[1]])) else as.matrix(data[[idx.sim]])
        }
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
        
        if (method == "lasso") {
          
          fit <- glmnet::cv.glmnet(x.train, y.train_bin,
                                   family = "binomial",
                                   type.measure = "class",
                                   alpha = 1)
          pred_prob <- as.numeric(predict(fit, s = fit$lambda.1se, newx = x.test, type = "response"))
          pred_class <- ifelse(pred_prob > 0.5, 1, 0)
          ce_vec[f] <- mean(pred_class != y.test_bin)
          
        } else if (method == "elnet") {
          
          fit <- glmnet::cv.glmnet(x.train, y.train_bin,
                                   family = "binomial",
                                   type.measure = "class",
                                   alpha = 0.5)
          pred_prob <- as.numeric(predict(fit, s = fit$lambda.1se, newx = x.test, type = "response"))
          pred_class <- ifelse(pred_prob > 0.5, 1, 0)
          ce_vec[f] <- mean(pred_class != y.test_bin)
          
        } else if (method == "RF") {
          
          fit <- randomForest::randomForest(x = x.train, y = y.train_fac)
          rf.class <- predict(fit, newdata = x.test)
          rf.vote  <- predict(fit, newdata = x.test, type = "vote")
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
          pred_prob <- knn.prob[, 2]
          ce_vec[f] <- mean(knn.class != y.test_fac)
          
        } else if (method == "PLS") {
          
          fit <- caret::plsda(x.train, y.train_fac, ncomp = 1)
          
          pls.class <- predict(fit, x.test, type = "class")
          pls.prob  <- predict(fit, x.test, type = "prob")
          
          pred_prob <- pls.prob[,2,1]
          
          ce_vec[f] <- mean(pls.class != y.test_fac)
          
        } else {
          stop("Unknown 'method'.")
        }
        
        roc.obj <- pROC::roc(y.test_bin, pred_prob, quiet = TRUE)
        auc_vec[f] <- as.numeric(pROC::auc(roc.obj))
      }
      
      mce[s]  <- mean(ce_vec)
      mauc[s] <- mean(auc_vec)
    }
    
    all.mce[i, ]  <- mce
    all.mauc[i, ] <- mauc
  }
  
  colnames(all.mce)  <- feat_names
  colnames(all.mauc) <- feat_names
  
  res <- list(
    CE     = all.mce,
    AUC    = all.mauc,
    method = method,
    all    = isTRUE(all)
  )
  class(res) <- "wavFeatExtClassifier"
  res
}