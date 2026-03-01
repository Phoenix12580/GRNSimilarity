#' Simulate virtual knockout of transcription factors
#'
#' @param grn GRN object returned by [inferGRN()] or a square adjacency matrix.
#' @param tf_list Character vector of TFs to knockout.
#'
#' @return A list with knockout network and per-gene impact score.
#' @export
simulateKO <- function(grn, tf_list) {
  adjacency <- .extract_adjacency(grn)

  if (missing(tf_list) || length(tf_list) == 0) {
    stop("`tf_list` must contain at least one transcription factor.")
  }
  tf_list <- unique(as.character(tf_list))
  tf_list <- tf_list[!is.na(tf_list) & nzchar(tf_list)]
  if (length(tf_list) == 0) {
    stop("`tf_list` contains only missing/empty values.")
  }

  tf_hit <- intersect(tf_list, rownames(adjacency))
  if (length(tf_hit) == 0) {
    warning("None of `tf_list` found in GRN row names. Returning unchanged network.")
    zero <- stats::setNames(rep(0, nrow(adjacency)), rownames(adjacency))
    return(list(
      knocked_tfs = character(0),
      grn_before = adjacency,
      grn_after = adjacency,
      impact_score = zero
    ))
  }

  grn_after <- adjacency
  grn_after[tf_hit, ] <- 0
  grn_after[, tf_hit] <- 0

  impact_score <- rowSums(abs(adjacency - grn_after)) + colSums(abs(adjacency - grn_after))

  list(
    knocked_tfs = tf_hit,
    grn_before = adjacency,
    grn_after = grn_after,
    impact_score = impact_score
  )
}

.extract_adjacency <- function(grn) {
  if (is.matrix(grn)) {
    adjacency <- grn
  } else if (is.list(grn) && !is.null(grn$adjacency) && is.matrix(grn$adjacency)) {
    adjacency <- grn$adjacency
  } else {
    stop("`grn` must be an adjacency matrix or an object returned by `inferGRN`.")
  }

  if (!is.numeric(adjacency)) {
    stop("Adjacency matrix must be numeric.")
  }
  if (nrow(adjacency) != ncol(adjacency)) {
    stop("Adjacency matrix must be square.")
  }
  if (anyNA(adjacency)) {
    warning("Adjacency contains NA; replacing NA with 0.")
    adjacency[is.na(adjacency)] <- 0
  }

  if (is.null(rownames(adjacency))) {
    genes <- paste0("Gene", seq_len(nrow(adjacency)))
    rownames(adjacency) <- genes
    colnames(adjacency) <- genes
  }
  if (is.null(colnames(adjacency))) {
    colnames(adjacency) <- rownames(adjacency)
  }

  adjacency
}
