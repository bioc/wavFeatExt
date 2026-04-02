#' Plot Classification Results from \code{classifyPcaIca}
#'
#' Create boxplots of cross-validated performance measures
#' (classification error or AUC) for the feature sets used in
#' \code{\link{classifyPcaIca}}, namely PCA-based, ICA-based,
#' and original/segmented data features.
#'
#' @param x The result object returned by \code{\link{classifyPcaIca}},
#'   containing the components \code{CE}, \code{AUC}, and \code{method}.
#' @param type Character string specifying which performance measure to plot.
#'   One of \code{"CE"} or \code{"AUC"}.
#' @param d.type Character string indicating which feature types to include
#'   in the plot: \code{"both"}, \code{"PCA"}, or \code{"ICA"}.
#' @param adjust Logical flag reserved for future extensions of the plotting
#'   function. Currently not used.
#' @param ... Additional graphical parameters passed to
#'   \code{\link[graphics:boxplot]{graphics::boxplot}}.
#'
#' @details
#' This function assumes that \code{\link{classifyPcaIca}} assigns column
#' names to its \code{CE} and \code{AUC} matrices, with prefixes
#' \code{"PC"} for PCA-based features, \code{"I"} for ICA-based features,
#' and the column name \code{"seg"} for the original (or segmented) data matrix.
#'
#' @return
#' This function is used for its side effect of producing a plot and returns
#' \code{invisible(NULL)}.
#'
#' @seealso
#' \code{\link{classifyPcaIca}},
#' \code{\link{classifyWavFeatExt}},
#' \code{\link{simulateCNA}},
#' \code{\link{getPca}},
#' \code{\link{getIca}}
#'
#' @examples
#' set.seed(10)
#' sim.dat <- simulateCNA(
#'   n.obs = 20,
#'   p = 32,
#'   n.sim = 1,
#'   n.block = 8,
#'   verbose = FALSE
#' )
#' pca <- getPca(sim.dat, k = 4)
#' ica <- getIca(sim.dat, k = 4)
#' y <- factor(rep(c("Group1", "Group2"), each = 10))
#' res <- classifyPcaIca(
#'   sim.dat,
#'   y,
#'   pca,
#'   ica,
#'   method = "KNN",
#'   k = 5,
#'   ite = 1
#' )
#' plot(res, type = "CE")
#' plot(res, type = "AUC", d.type = "PCA")
#'
#' @method plot pcaIcaClassifier
#' @export
plot.pcaIcaClassifier <- function(x,
                                type  = c("CE", "AUC"),
                                d.type = c("both", "PCA", "ICA"),
                                adjust = TRUE,  # reserved, not used currently
                                ...) {
  # x: result from classif.pcaica()
  type  <- match.arg(type)
  d.type <- match.arg(d.type)
  
  ## --- Basic checks ---
  if (!is.list(x) || is.null(x$CE) || is.null(x$AUC)) {
    stop("'x' must be the result object returned by 'classif.pcaica()'.")
  }
  
  CE  <- x$CE
  AUC <- x$AUC
  
  if (is.vector(CE))  CE  <- matrix(CE,  nrow = 1L)
  if (is.vector(AUC)) AUC <- matrix(AUC, nrow = 1L)
  
  if (ncol(CE) != ncol(AUC)) {
    stop("'CE' and 'AUC' in 'x' must have the same number of columns.")
  }
  
  cn <- colnames(AUC)
  if (is.null(cn)) {
    stop("Column names of 'AUC' (and 'CE') must be set in 'classif.pcaica()'.")
  }
  
  ## --- Identify PCA / ICA / seg columns by name ---
  pc_idx  <- grep("^PC[0-9]+$", cn)
  ic_idx  <- grep("^I[0-9]+$",  cn)
  seg_idx <- which(cn == "seg")
  
  if (length(seg_idx) != 1L) {
    stop("Could not uniquely identify the 'seg' column in the results.")
  }
  
  if (length(pc_idx) == 0L && d.type %in% c("both", "PCA")) {
    warning("No 'PC*' columns found; skipping PCA curves.")
  }
  if (length(ic_idx) == 0L && d.type %in% c("both", "ICA")) {
    warning("No 'I*' columns found; skipping ICA curves.")
  }
  
  idx <- switch(d.type,
                "both" = c(pc_idx, ic_idx, seg_idx),
                "PCA"  = c(pc_idx, seg_idx),
                "ICA"  = c(ic_idx, seg_idx)
  )
  idx <- idx[idx > 0L]  # buang indeks kosong kalau ada
  
  if (length(idx) == 0L) {
    stop("No columns selected to plot for the chosen 'd.type'.")
  }
  
  subset_CE  <- CE[,  idx, drop = FALSE]
  subset_AUC <- AUC[, idx, drop = FALSE]
  labels     <- cn[idx]
  
  ## --- Pilih matriks mana yang dipakai ---
  if (type == "AUC") {
    mat <- subset_AUC
    main_title <- paste("AUC:", x$method)
  } else {
    mat <- subset_CE
    main_title <- paste("Classification Error:", x$method)
  }
  
  ## --- Boxplot ---
  boxplot(mat, names = labels, main = main_title, ...)
  
  ## --- Garis referensi horizontal ---
  # Seg column is always the last of subset (karena selalu kita masukkan terakhir)
  seg_col <- ncol(mat)
  seg_vals <- mat[, seg_col]
  
  seg_med <- stats::median(seg_vals, na.rm = TRUE)
  feat_meds <- matrixStats::colMedians(mat, na.rm = TRUE)
  
  if (type == "AUC") {
    # merah: median AUC seg
    # biru: AUC terbaik (median terbesar)
    best_val <- max(feat_meds, na.rm = TRUE)
    abline(h = seg_med,  lty = 2, col = "red")
    abline(h = best_val, lty = 3, col = "blue")
  } else {
    # merah: median CE seg
    # biru: CE terbaik (median terkecil)
    best_val <- min(feat_meds, na.rm = TRUE)
    abline(h = seg_med,  lty = 2, col = "red")
    abline(h = best_val, lty = 3, col = "blue")
  }
  
  invisible(NULL)
}
