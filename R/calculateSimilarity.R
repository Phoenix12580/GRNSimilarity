#' Calculate similarity between two GRNs with multiple metrics
#'
#' @param grn_before Adjacency matrix or GRN object.
#' @param grn_after Adjacency matrix or GRN object.
#' @param metrics Character vector of metrics. Supported values include
#'   `"pearson"`, `"cosine"`, `"jaccard"`, `"spearman"`, `"frobenius"`.
#' @param edge_threshold Numeric threshold used to binarize weighted edges for
#'   Jaccard similarity.
#'
#' @return A list containing selected metrics and aligned gene names.
#' @export
calculateSimilarity <- function(grn_before,
                                grn_after,
                                metrics = c("pearson", "cosine", "jaccard"),
                                edge_threshold = 0) {
  a <- .extract_adjacency(grn_before)
  b <- .extract_adjacency(grn_after)

  common <- intersect(rownames(a), rownames(b))
  if (length(common) < 2) {
    stop("Networks must have at least two overlapping genes.")
  }

  if (!is.numeric(edge_threshold) || length(edge_threshold) != 1 || is.na(edge_threshold)) {
    stop("`edge_threshold` must be a single numeric value.")
  }

  allowed_metrics <- c("pearson", "cosine", "jaccard", "spearman", "frobenius")
  if (!is.character(metrics) || length(metrics) == 0) {
    stop("`metrics` must be a non-empty character vector.")
  }
  metrics <- unique(tolower(metrics))
  invalid <- setdiff(metrics, allowed_metrics)
  if (length(invalid) > 0) {
    stop(sprintf("Unsupported metrics: %s. Allowed: %s",
                 paste(invalid, collapse = ", "),
                 paste(allowed_metrics, collapse = ", ")))
  }

  a <- a[common, common, drop = FALSE]
  b <- b[common, common, drop = FALSE]

  idx <- upper.tri(a, diag = FALSE)
  va <- as.numeric(a[idx])
  vb <- as.numeric(b[idx])

  out <- list(overlap_genes = common)

  if ("pearson" %in% metrics) {
    out$pearson <- suppressWarnings(stats::cor(va, vb, method = "pearson", use = "pairwise.complete.obs"))
  }
  if ("spearman" %in% metrics) {
    out$spearman <- suppressWarnings(stats::cor(va, vb, method = "spearman", use = "pairwise.complete.obs"))
  }
  if ("cosine" %in% metrics) {
    out$cosine <- sum(va * vb) / (sqrt(sum(va^2)) * sqrt(sum(vb^2)) + 1e-12)
  }
  if ("jaccard" %in% metrics) {
    ea <- va > edge_threshold
    eb <- vb > edge_threshold
    union_count <- sum(ea | eb)
    out$jaccard <- if (union_count == 0) 1 else sum(ea & eb) / union_count
  }
  if ("frobenius" %in% metrics) {
    out$frobenius <- sqrt(sum((a - b)^2))
  }

  out
}
