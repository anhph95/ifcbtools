## ================================================================
## Bootstrap resampling for community arrays
## ================================================================
## This file only resamples. Metric computation is delegated to
## metric_fn, which receives X[time, site, taxon].
bootstrap_by_dimension <- function(X,
                                   margin,
                                   metric_fn = calc_metacommunity_metrics,
                                   n_boot = 1000,
                                   seed = 123,
                                   baseline_in_boot = TRUE,
                                   show_progress = FALSE) {
  set.seed(seed)
  margin_idx <- if (is.character(margin)) {
    match(margin, names(dimnames(X)))
  } else {
    margin
  }
  n_margin <- dim(X)[[margin_idx]]
  baseline <- metric_fn(X)
  ## Optionally show lightweight progress during bootstrap resampling.
  ##
  ## The progress bar tracks completed bootstrap replicates:
  ##   i = 1, ..., n_boot
  progress_bar <- NULL
  if (show_progress) {
    progress_bar <- utils::txtProgressBar(min = 0, max = n_boot, style = 3)
    on.exit(close(progress_bar), add = TRUE)
  }
  boot_replicates <- purrr::map_dfr(seq_len(n_boot), function(i) {
    sampled_idx <- sample.int(n_margin, size = n_margin, replace = TRUE)
    idx <- rep(list(TRUE), length(dim(X)))
    idx[[margin_idx]] <- sampled_idx
    boot_metrics <- metric_fn(do.call(`[`, c(list(X), idx, list(drop = FALSE)))) %>%
      dplyr::select(varname, estimate) %>%
      tidyr::pivot_wider(names_from = varname, values_from = estimate) %>%
      dplyr::mutate(
        boot_id = i,
        sample_type = "Bootstrap",
        .before = 1
      )
    if (show_progress) {
      utils::setTxtProgressBar(progress_bar, i)
    }
    boot_metrics
  })
  baseline_wide <- baseline %>%
    dplyr::select(varname, estimate) %>%
    tidyr::pivot_wider(names_from = varname, values_from = estimate) %>%
    dplyr::mutate(
      boot_id = 0L,
      sample_type = "Baseline",
      .before = 1
    )
  boot_out <- if (baseline_in_boot) {
    dplyr::bind_rows(baseline_wide, boot_replicates)
  } else {
    boot_replicates
  }
  summary <- baseline %>%
    dplyr::left_join(
      tibble::tibble(
        varname = baseline$varname,
        lwr = purrr::map_dbl(
          baseline$varname,
          ~ stats::quantile(boot_replicates[[.x]], 0.025, na.rm = TRUE)
        ),
        upr = purrr::map_dbl(
          baseline$varname,
          ~ stats::quantile(boot_replicates[[.x]], 0.975, na.rm = TRUE)
        )
      ),
      by = "varname"
    )
  list(
    baseline = baseline,
    boot = boot_out,
    boot_replicates = boot_replicates,
    summary = summary
  )
}
