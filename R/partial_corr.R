#' Partial Correlation Network
#'
#' Build a partial correlation network with scRNAseq data.
#'
#' @param mtx single-cell expression matrix (log2 normalized data), genes by cells. Typically genes are transcription factors.
#' @param threshold parameter passed to \code{threshold} in \code{qgraph}, could be "sig" without multiple testing correction, or "BH" for fdr.
#'        See \code{qgraph} for detail.
#' @param alpha threshold to select significant edges.
#' @details
#' Partial correlation network based on spearman correlation.
#' @return
#' A matrix.
#' @export
skel_net <- function(mtx, threshold = "BH", alpha = 0.05) {
  if (!is.matrix(mtx)) stop("(EE) Input should be a matrix.")

  #- Values in a row are identicial, might happen in subsampling.
  if (any(apply(mtx, 1, sd) %==% 0)) {
    message("(II) NOT run: zero standard deviation found.")
  } else {
    cor_mtx <- cor(t(mtx), method = "spearman")

    if (all(eigen(cor_mtx)$values > 1E-6)) {
      qres <- qgraph::qgraph(cor_mtx,
                             graph      = "pcor",
                             sampleSize = ncol(mtx),
                             labels     = rownames(cor_mtx),
                             threshold  = threshold,
                             alpha      = alpha,
                             DoNotPlot  = TRUE)

      #- Return a symm weight matrix.
      w_mtx <- qgraph::getWmat(qres, directed = FALSE) %>% extract(rownames(mtx), rownames(mtx))
      return(w_mtx)
    } else {
      message("(II) The correlation matrix is not positive definite.")
    }
  }
}

#' PCN on large datasets
#'
#' Build a consensus netork based on partial correlation networks on subsamples of large scRNAseq data.
#'
#' @param mtx single-cell expression matrix (log2 normalized data), genes by cells. Typically genes are transcription factors.
#' @param n_cell number of cells to sample. If NULL, boostrapping the data and keep unique sets of cells.
#' @param n_draw times of sampling.
#' @param ... additional parameters to \code{skel_net}.
#' @return
#' A matrix.
#' @export
pcn <- function(mtx, n_cell = NULL, n_draw = 50, ...) {
  if (!is.matrix(mtx)) stop("(EE) Input should be a matrix.")

  suppressPackageStartupMessages(
    res <- foreach(i = seq_len(n_draw)) %dorng% {
      if (is.null(n_cell) || n_cell > ncol(mtx)) {
        subcell <- sample(colnames(mtx), ncol(mtx), replace = TRUE) %>% unique
      } else {
        subcell <- sample(colnames(mtx), n_cell)
      }

      skel_net(mtx[, subcell], ...)
    }
  )

  #- Remove empty elements.
  res[lengths(res) == 0] <- NULL
  consensus_net <- Reduce("+", res) / length(res)
  return(consensus_net)
}
