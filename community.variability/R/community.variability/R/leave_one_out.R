## ================================================================
## Sensitivity analyses for community arrays
## ================================================================
## Leave-one-out operations are direct array slices. Metric computation
## remains separate and receives X[time, site, taxon].
leave_one_out <- function(X,
                          margin,
                          metric_fn = calc_metacommunity_metrics,
                          show_progress = FALSE) {
  margin_idx <- if (is.character(margin)) {
    match(margin, names(dimnames(X)))
  } else {
    margin
  }
  removed_col <- paste0(names(dimnames(X))[margin_idx], "_removed")
  labels <- dimnames(X)[[margin_idx]]
  baseline <- metric_fn(X) %>%
    dplyr::mutate(!!removed_col := "Baseline", .before = 1)
  ## Optionally show lightweight progress during leave-one-out slicing.
  ##
  ## The progress bar tracks omitted components:
  ##   i = 1, ..., length(labels)
  progress_bar <- NULL
  if (show_progress) {
    progress_bar <- utils::txtProgressBar(min = 0, max = length(labels), style = 3)
    on.exit(close(progress_bar), add = TRUE)
  }
  omitted <- purrr::map_dfr(seq_along(labels), function(i) {
    idx <- rep(list(TRUE), length(dim(X)))
    idx[[margin_idx]] <- -i
    omitted_metrics <- metric_fn(do.call(`[`, c(list(X), idx, list(drop = FALSE)))) %>%
      dplyr::mutate(!!removed_col := labels[i], .before = 1)
    if (show_progress) {
      utils::setTxtProgressBar(progress_bar, i)
    }
    omitted_metrics
  })
  dplyr::bind_rows(baseline, omitted)
}
add_baseline_delta <- function(sensitivity_results,
                               removed_col,
                               baseline_label = "Baseline",
                               group_cols = NULL) {
  join_cols <- c(group_cols, "varname")
  baseline <- sensitivity_results %>%
    dplyr::filter(.data[[removed_col]] == baseline_label) %>%
    dplyr::select(dplyr::all_of(join_cols), baseline = estimate)
  sensitivity_results %>%
    dplyr::left_join(baseline, by = join_cols) %>%
    dplyr::mutate(
      delta = .data$estimate - .data$baseline,
      abs_delta = abs(.data$delta)
    )
}
