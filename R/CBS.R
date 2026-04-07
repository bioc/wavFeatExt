#' Circular Binary Segmentation of Copy-number Profiles
#'
#' Performs circular binary segmentation (CBS) using the
#' \pkg{DNAcopy} package on a numeric vector representing
#' a single copy-number profile.
#'
#' @param obj Numeric vector containing the copy-number
#'   (or log-ratio) values along the genome for a single sample.
#' @param chr Vector of the same length as \code{obj} giving the
#'   chromosome index for each position. By default all positions
#'   are assumed to belong to chromosome 1.
#' @param verbose Logical. If \code{TRUE}, progress messages from
#'   the segmentation procedure are shown. If \code{FALSE}, they
#'   are suppressed.
#'
#' @details
#' This function is a thin wrapper around \code{\link[DNAcopy]{CNA}}
#' and \code{\link[DNAcopy]{segment}} from the \pkg{DNAcopy} package.
#' The input vector is first converted to a \code{"CNA"} object using
#' artificial map locations \code{1:length(obj)}. Circular binary
#' segmentation is then applied, and the estimated segment means are
#' expanded back to the original resolution to form a segmented
#' profile of the same length as the input.
#'
#' @return A numeric vector of length \code{length(obj)} containing
#'   the segmented copy-number (or log-ratio) values, where each
#'   element equals the estimated segment mean for the corresponding
#'   genomic position.
#'
#' @references
#' Olshen, A. B., Venkatraman, E. S., Lucito, R. and Wigler, M. (2004).
#' Circular binary segmentation for the analysis of array-based DNA
#' copy number data. \emph{Biostatistics}, \strong{5}(4), 557--572.
#'
#' @seealso \code{\link{wavFeatExt}}
#'
#' @examples
#' ## Generate noisy copy-number data under normal error
#' set.seed(10)
#' cn <- rnorm(100, mean = c(rep(2, 50), rep(1, 50)), sd = 0.4)
#'
#' ## Segmentation using CBS
#' cn.cbs <- CBS(cn)
#'
#' ## Plot raw and segmented profiles
#' plot(cn, type = "p", xlab = "Genomic position", ylab = "Copy-number")
#' lines(cn.cbs, col = 2, lwd = 2)
#'
#' @author Maharani Ummi \email{maharaniahsani@@gmail.com}
#'
#' @export
#' @importFrom DNAcopy CNA segment
CBS <- function(obj, chr = rep(1L, length(obj)), verbose = FALSE) {
  CNA.obj <- DNAcopy::CNA(obj, chrom = chr, maploc = seq_along(obj), data.type = "binary")
  seg.obj <- DNAcopy::segment(CNA.obj, verbose = if (isTRUE(verbose)) 1L else 0L,)
  num.seg <- length(seg.obj$segRows$startRow)
  cbs.seg <- numeric(length(obj))
  for (k in seq_len(num.seg)) {
    cbs.seg[seg.obj$segRows$startRow[k]:seg.obj$segRows$endRow[k]] <- seg.obj$output$seg.mean[k]
  }
  cbs.seg
}
