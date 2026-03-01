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
#' @param use_external Logical, whether to use package backends (`hdWGCNA`, `corto`, `SCENIC`).
#' @param method_params Named list of method-specific parameters passed to [inferGRN()].
#' @param similarity_metrics Similarity metrics passed to [calculateSimilarity()].
#' @param primary_metric Metric used for ranking TF effects.
#' @param edge_threshold Threshold for Jaccard binarization.
#'
#' @return A list with baseline similarity, per-TF rankings, ggplot data, and GRNs.
#' @export
networkRanking <- function(expr_data,
                           cells_col,
                           cells_A,
                           cells_B,
                           tf_list,
                           method = c("hdWGCNA", "ARACNE", "SCENIC"),
                           power = 6,
                           k = 20,
                           regulon_db = NULL,
                           use_external = TRUE,
                           method_params = list(),
                           similarity_metrics = c("cosine", "jaccard", "pearson"),
                           primary_metric = "cosine",
                           edge_threshold = 0) {
  method <- match.arg(method)
  .validate_inputs(expr_data, cells_col, cells_A, cells_B, tf_list)

  similarity_metrics <- unique(tolower(similarity_metrics))
  primary_metric <- tolower(primary_metric)
  if (!primary_metric %in% similarity_metrics) {
    similarity_metrics <- unique(c(similarity_metrics, primary_metric))
  }

  idx_A <- .resolve_cell_index(cells_col, cells_A)
  idx_B <- .resolve_cell_index(cells_col, cells_B)

  expr_A <- expr_data[, idx_A, drop = FALSE]
  expr_B <- expr_data[, idx_B, drop = FALSE]

  grn_A <- inferGRN(expr_A, method = method, power = power, k = k, regulon_db = regulon_db,
                    use_external = use_external, method_params = method_params)
  grn_B <- inferGRN(expr_B, method = method, power = power, k = k, regulon_db = regulon_db,
                    use_external = use_external, method_params = method_params)

  baseline <- calculateSimilarity(grn_A, grn_B, metrics = similarity_metrics, edge_threshold = edge_threshold)

  tf_list <- unique(as.character(tf_list))
  tf_list <- tf_list[!is.na(tf_list) & nzchar(tf_list)]
  if (length(tf_list) == 0) {
    stop("`tf_list` contains only missing/empty values.")
  }

  tf_rows <- vector("list", length(tf_list))
  per_tf_networks <- vector("list", length(tf_list))
  gene_impacts <- vector("list", length(tf_list))

  for (i in seq_along(tf_list)) {
    tf <- tf_list[i]
    ko_A <- simulateKO(grn_A, tf)
    ko_B <- simulateKO(grn_B, tf)
    sim_after <- calculateSimilarity(ko_A$grn_after, ko_B$grn_after,
                                     metrics = similarity_metrics,
                                     edge_threshold = edge_threshold)

    metric_rows <- do.call(rbind, lapply(similarity_metrics, function(metric_name) {
      base_v <- if (!is.null(baseline[[metric_name]])) baseline[[metric_name]] else NA_real_
      post_v <- if (!is.null(sim_after[[metric_name]])) sim_after[[metric_name]] else NA_real_
      direction <- if (metric_name == "frobenius") 1 else -1
      delta_v <- direction * (post_v - base_v)
      data.frame(
        tf = tf,
        metric = metric_name,
        baseline = base_v,
        post_ko = post_v,
        delta = delta_v,
        stringsAsFactors = FALSE
      )
    }))

    tf_rows[[i]] <- metric_rows
    per_tf_networks[[i]] <- list(tf = tf, A_after = ko_A$grn_after, B_after = ko_B$grn_after)

    common <- intersect(names(ko_A$impact_score), names(ko_B$impact_score))
    score <- (ko_A$impact_score[common] + ko_B$impact_score[common]) / 2
    full <- stats::setNames(rep(0, nrow(expr_data)), rownames(expr_data))
    full[common] <- score
    gene_impacts[[i]] <- full
  }

  names(per_tf_networks) <- tf_list
  ranking_long <- do.call(rbind, tf_rows)
  ranking <- ranking_long[ranking_long$metric == primary_metric, , drop = FALSE]
  ranking <- ranking[order(ranking$delta, decreasing = TRUE), , drop = FALSE]
  rownames(ranking) <- NULL

  gene_scores <- rowMeans(do.call(cbind, gene_impacts))
  gene_scores <- sort(gene_scores, decreasing = TRUE)
  gene_scores_df <- data.frame(gene = names(gene_scores), score = as.numeric(gene_scores), stringsAsFactors = FALSE)

  largest_change_tf <- if (nrow(ranking) > 0) ranking$tf[1] else NA_character_

  list(
    baseline_similarity = baseline,
    ranking = ranking,
    ranking_long = ranking_long,
    largest_change = list(
      tf = largest_change_tf,
      networks = per_tf_networks[[largest_change_tf]]
    ),
    gene_scores = gene_scores,
    plot_data = list(
      tf_ranking = ranking_long,
      gene_scores = gene_scores_df,
      edges_A = grn_A$plot_data,
      edges_B = grn_B$plot_data
    ),
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
    if (any(idx < 1 | idx > length(cells_col))) {
      stop("Numeric subgroup indices are out of valid cell range.")
    }
  } else {
    idx <- which(cells_col %in% selector)
  }

  if (length(idx) < 3) {
    stop("Each subgroup must include at least 3 cells.")
  }
  idx
}
