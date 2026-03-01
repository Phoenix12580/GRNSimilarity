#' Rank TF impact by subgroup-specific network changes
#'
#' @param expr_data Numeric matrix (genes x cells).
#' @param cells_col Cell subgroup annotation vector, length equals number of cells.
#' @param cells_A Label value(s) in `cells_col` OR numeric indices for subgroup A.
#' @param cells_B Label value(s) in `cells_col` OR numeric indices for subgroup B.
#' @param tf_list Character vector of TFs for virtual knockout.
#' @param method One of "hdWGCNA", "ARACNE", or "SCENIC".
#' @param power Soft-thresholding power (used by hdWGCNA-like inference).
#' @param k Number of retained edges per regulator.
#' @param regulon_db Optional regulon mapping for SCENIC-like inference.
#'
#' @return A list with baseline similarity, per-TF results, and rankings.
networkRanking <- function(expr_data,
                           cells_col,
                           cells_A,
                           cells_B,
                           tf_list,
                           method = c("hdWGCNA", "ARACNE", "SCENIC"),
                           power = 6,
                           k = 20,
                           regulon_db = NULL) {
  method <- match.arg(method)
  .validate_inputs(expr_data, cells_col, cells_A, cells_B, tf_list)

  idx_A <- .resolve_cell_index(cells_col, cells_A)
  idx_B <- .resolve_cell_index(cells_col, cells_B)

  expr_A <- expr_data[, idx_A, drop = FALSE]
  expr_B <- expr_data[, idx_B, drop = FALSE]

  grn_A <- inferGRN(expr_A, method = method, power = power, k = k, regulon_db = regulon_db)
  grn_B <- inferGRN(expr_B, method = method, power = power, k = k, regulon_db = regulon_db)

  baseline <- calculateSimilarity(grn_A, grn_B)

  tf_results <- lapply(tf_list, function(tf) {
    ko_A <- simulateKO(grn_A, tf)
    ko_B <- simulateKO(grn_B, tf)
    sim_after <- calculateSimilarity(ko_A$grn_after, ko_B$grn_after)

    data.frame(
      tf = tf,
      baseline_cosine = baseline$cosine,
      postKO_cosine = sim_after$cosine,
      delta_cosine = baseline$cosine - sim_after$cosine,
      baseline_frobenius = baseline$frobenius_change,
      postKO_frobenius = sim_after$frobenius_change,
      delta_frobenius = sim_after$frobenius_change - baseline$frobenius_change,
      stringsAsFactors = FALSE
    )
  })

  ranking <- do.call(rbind, tf_results)
  ranking <- ranking[order(ranking$delta_cosine, decreasing = TRUE), , drop = FALSE]
  rownames(ranking) <- NULL

  per_tf_networks <- lapply(tf_list, function(tf) {
    ko_A <- simulateKO(grn_A, tf)
    ko_B <- simulateKO(grn_B, tf)
    list(tf = tf, A_after = ko_A$grn_after, B_after = ko_B$grn_after)
  })
  names(per_tf_networks) <- tf_list

  gene_scores <- rowMeans(do.call(cbind, lapply(tf_list, function(tf) {
    ko_A <- simulateKO(grn_A, tf)
    ko_B <- simulateKO(grn_B, tf)
    common <- intersect(names(ko_A$impact_score), names(ko_B$impact_score))
    score <- (ko_A$impact_score[common] + ko_B$impact_score[common]) / 2
    full <- rep(0, nrow(expr_data))
    names(full) <- rownames(expr_data)
    full[common] <- score
    full
  })))

  largest_change_tf <- ranking$tf[which.max(ranking$delta_cosine)]

  list(
    baseline_similarity = baseline,
    ranking = ranking,
    largest_change = list(
      tf = largest_change_tf,
      networks = per_tf_networks[[largest_change_tf]]
    ),
    gene_scores = sort(gene_scores, decreasing = TRUE),
    grn_A = grn_A,
    grn_B = grn_B
  )
}

.validate_inputs <- function(expr_data, cells_col, cells_A, cells_B, tf_list) {
  if (!is.matrix(expr_data) || !is.numeric(expr_data)) {
    stop("`expr_data` must be a numeric matrix (genes x cells).")
  }
  if (length(cells_col) != ncol(expr_data)) {
    stop("`cells_col` length must equal number of columns (cells) in `expr_data`.")
  }
  if (missing(cells_A) || missing(cells_B)) {
    stop("Both `cells_A` and `cells_B` must be provided.")
  }
  if (missing(tf_list) || length(tf_list) == 0) {
    stop("`tf_list` must contain at least one TF.")
  }
}

.resolve_cell_index <- function(cells_col, selector) {
  if (is.numeric(selector)) {
    idx <- unique(as.integer(selector))
  } else {
    idx <- which(cells_col %in% selector)
  }

  if (length(idx) < 3) {
    stop("Each subgroup must include at least 3 cells.")
  }
  idx
}
