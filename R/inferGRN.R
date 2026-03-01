#' Infer a single-cell gene regulatory network
#'
#' @param expr_data Numeric matrix (genes x cells).
#' @param method One of "hdWGCNA", "ARACNE", or "SCENIC".
#' @param power Soft-thresholding power used by the hdWGCNA-like method.
#' @param k Number of strongest outgoing edges retained per gene.
#' @param regulon_db Optional regulon specification for SCENIC-like mode.
#'   Can be a named list (`tf -> character vector of targets`) or a
#'   two-column data.frame with columns `tf` and `target`.
#'
#' @return A list with adjacency matrix and edge table.
inferGRN <- function(expr_data,
                     method = c("hdWGCNA", "ARACNE", "SCENIC"),
                     power = 6,
                     k = 20,
                     regulon_db = NULL) {
  method <- match.arg(method)

  if (!is.matrix(expr_data) || !is.numeric(expr_data)) {
    stop("`expr_data` must be a numeric matrix with genes in rows and cells in columns.")
  }
  if (nrow(expr_data) < 2 || ncol(expr_data) < 3) {
    stop("`expr_data` must contain at least 2 genes and 3 cells.")
  }

  if (is.null(rownames(expr_data))) {
    rownames(expr_data) <- paste0("Gene", seq_len(nrow(expr_data)))
  }

  cor_mat <- stats::cor(t(expr_data), method = "pearson", use = "pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  diag(cor_mat) <- 0

  adjacency <- switch(
    method,
    hdWGCNA = abs(cor_mat)^power,
    ARACNE = .aracne_like_adjacency(cor_mat),
    SCENIC = .scenic_like_adjacency(cor_mat, regulon_db)
  )

  adjacency <- .keep_top_k(adjacency, k = k)

  edges <- as.data.frame(as.table(adjacency), stringsAsFactors = FALSE)
  colnames(edges) <- c("regulator", "target", "weight")
  edges <- edges[edges$weight > 0, , drop = FALSE]
  edges <- edges[order(edges$weight, decreasing = TRUE), , drop = FALSE]

  list(
    method = method,
    adjacency = adjacency,
    edges = edges,
    parameters = list(power = power, k = k)
  )
}

.aracne_like_adjacency <- function(cor_mat) {
  mi_proxy <- -0.5 * log(pmax(1 - cor_mat^2, .Machine$double.eps))
  mi_proxy[!is.finite(mi_proxy)] <- max(mi_proxy[is.finite(mi_proxy)], na.rm = TRUE)
  mi_proxy[is.na(mi_proxy)] <- 0

  thresh <- stats::quantile(mi_proxy[upper.tri(mi_proxy)], probs = 0.7, na.rm = TRUE)
  mi_proxy[mi_proxy < thresh] <- 0
  mi_proxy
}

.scenic_like_adjacency <- function(cor_mat, regulon_db = NULL) {
  genes <- rownames(cor_mat)
  adj <- matrix(0, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))

  if (is.null(regulon_db)) {
    tf_idx <- grepl("^TF", genes, ignore.case = TRUE)
    if (!any(tf_idx)) {
      tf_idx <- seq_len(min(10, length(genes)))
    }
    adj[tf_idx, ] <- pmax(cor_mat[tf_idx, , drop = FALSE], 0)
    diag(adj) <- 0
    return(adj)
  }

  if (is.data.frame(regulon_db)) {
    if (!all(c("tf", "target") %in% colnames(regulon_db))) {
      stop("When `regulon_db` is data.frame, columns `tf` and `target` are required.")
    }
    tf_target_pairs <- split(as.character(regulon_db$target), as.character(regulon_db$tf))
  } else if (is.list(regulon_db)) {
    tf_target_pairs <- regulon_db
  } else {
    stop("`regulon_db` must be NULL, a named list, or a data.frame with columns tf and target.")
  }

  for (tf in names(tf_target_pairs)) {
    if (!tf %in% genes) next
    targets <- intersect(tf_target_pairs[[tf]], genes)
    if (length(targets) == 0) next
    adj[tf, targets] <- pmax(cor_mat[tf, targets], 0)
  }

  diag(adj) <- 0
  adj
}

.keep_top_k <- function(mat, k = 20) {
  if (!is.numeric(k) || length(k) != 1 || k <= 0) {
    stop("`k` must be a positive integer.")
  }
  k <- as.integer(k)

  out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  for (i in seq_len(nrow(mat))) {
    vals <- mat[i, ]
    if (all(vals <= 0)) next
    ord <- order(vals, decreasing = TRUE)
    keep <- ord[seq_len(min(k, length(ord)))]
    keep <- keep[vals[keep] > 0]
    out[i, keep] <- vals[keep]
  }
  diag(out) <- 0
  out
}
