## ================================================================
## Generic NMDS helpers
## ================================================================
run_nmds_ordination <- function(community_matrix,
                                distance = "euclidean",
                                k = 2,
                                trymax = 1000,
                                seed = 123,
                                trace = TRUE) {
  set.seed(seed)
  scaled <- vegan::decostand(
    community_matrix,
    method = "hellinger",
    na.rm = TRUE
  )
  vegan::metaMDS(
    scaled,
    distance = distance,
    k = k,
    trymax = trymax,
    trace = trace
  )
}
extract_nmds_sites <- function(model,
                               metadata,
                               site_col = "site",
                               time_col = "sample_time") {
  as.data.frame(vegan::scores(model, display = "sites")) %>%
    dplyr::rename(Axis1 = 1, Axis2 = 2) %>%
    dplyr::bind_cols(metadata) %>%
    dplyr::mutate(
      "{site_col}" := factor(.data[[site_col]]),
      "{time_col}" := .data[[time_col]]
    )
}
extract_nmds_path_segments <- function(ord_sites,
                                       site_col = "site",
                                       time_col = "sample_time") {
  all_segments <- ord_sites %>%
    dplyr::arrange(.data[[site_col]], .data[[time_col]]) %>%
    dplyr::group_by(.data[[site_col]]) %>%
    dplyr::mutate(
      x = .data$Axis1,
      y = .data$Axis2,
      xend = dplyr::lead(.data$Axis1),
      yend = dplyr::lead(.data$Axis2)
    ) %>%
    dplyr::filter(!is.na(.data$xend)) %>%
    dplyr::ungroup()
  list(
    all_segments = all_segments,
    start_points = ord_sites %>%
      dplyr::arrange(.data[[site_col]], .data[[time_col]]) %>%
      dplyr::group_by(.data[[site_col]]) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup(),
    end_segments = all_segments %>%
      dplyr::group_by(.data[[site_col]]) %>%
      dplyr::slice_tail(n = 1) %>%
      dplyr::ungroup()
  )
}
fit_env_vectors <- function(model,
                            env_data,
                            subset = NULL,
                            permutations = 999,
                            choices = c(1, 2),
                            p_threshold = 0.05) {
  fit <- vegan::envfit(
    model,
    env_data,
    permutations = permutations,
    choices = choices,
    subsets = subset,
    na.rm = TRUE
  )
  env_df <- as.data.frame(fit$vectors$arrows)
  env_df$var <- rownames(env_df)
  env_df$pval <- fit$vectors$pvals
  list(
    fit = fit,
    vectors = env_df %>%
      dplyr::filter(.data$pval < p_threshold)
  )
}
