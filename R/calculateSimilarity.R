#' Calculate similarity between two GRNs
#'
#' @param grn_before Adjacency matrix or GRN object.
#' @param grn_after Adjacency matrix or GRN object.
#'
#' @return A list with Pearson and cosine similarity plus Frobenius change.
calculateSimilarity <- function(grn_before, grn_after) {
  a <- .extract_adjacency(grn_before)
  b <- .extract_adjacency(grn_after)

  common <- intersect(rownames(a), rownames(b))
  if (length(common) < 2) {
    stop("Networks must have at least two overlapping genes.")
  }

  a <- a[common, common, drop = FALSE]
  b <- b[common, common, drop = FALSE]

  idx <- upper.tri(a, diag = FALSE)
  va <- as.numeric(a[idx])
  vb <- as.numeric(b[idx])

  pearson <- stats::cor(va, vb)
  cosine <- sum(va * vb) / (sqrt(sum(va^2)) * sqrt(sum(vb^2)) + 1e-12)
  frobenius_change <- sqrt(sum((a - b)^2))

  list(
    pearson = as.numeric(pearson),
    cosine = as.numeric(cosine),
    frobenius_change = as.numeric(frobenius_change),
    overlap_genes = common
  )
}
