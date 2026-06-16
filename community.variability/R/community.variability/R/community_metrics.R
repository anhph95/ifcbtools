## ================================================================
## Community array construction and metric output helpers
## ================================================================
##
## Data cleaning can remain tidy/data-frame oriented. Once the
## observation table is prepared, these helpers convert it to the
## numerical metacommunity array used by the variability functions:
##
##   X[time, site, taxon]
##
## where:
##   time  = repeated surveys t
##   site  = local communities i
##   taxon = species/taxa j
##
## Metric output order follows the alpha-gamma-phi partitioning:
##   gamma = regional/metacommunity variability
##   alpha = local/community variability
##   phi   = spatial synchrony component
##
## BD_spatial_weighted is reported as BD_beta here only if you want
## to keep the previous output name. It is a biomass-weighted summary
## of spatial compositional differences among sites through time, not
## the Lamy et al. BD_phi ratio.
community_metric_order <- c(
  "CV_gamma", "CV_alpha", "CV_phi",
  "BD_gamma", "BD_alpha", "BD_phi", "BD_beta"
)
make_community_array <- function(data_wide,
                                 taxa_cols,
                                 time_step_col_name = "timestep",
                                 site_id_col_name = "site") {
  ## Identify the time and site levels that define the first two
  ## dimensions of X. Sorting gives reproducible array ordering.
  time_ids <- sort(unique(data_wide[[time_step_col_name]]))
  site_ids <- sort(unique(data_wide[[site_id_col_name]]))
  ## Initialize the biomass array:
  ##   X[t, i, j] = biomass of taxon j at site i and time t
  X <- array(
    0,
    dim = c(length(time_ids), length(site_ids), length(taxa_cols)),
    dimnames = list(
      timestep = as.character(time_ids),
      site = as.character(site_ids),
      taxon = taxa_cols
    )
  )
  ## Map each row in the prepared wide table to its array position.
  time_index <- match(data_wide[[time_step_col_name]], time_ids)
  site_index <- match(data_wide[[site_id_col_name]], site_ids)
  ## Extract taxon biomass columns as a numeric matrix.
  ## Non-finite values are treated as zero biomass.
  biomass <- as.matrix(data_wide[, taxa_cols, drop = FALSE])
  storage.mode(biomass) <- "double"
  biomass[!is.finite(biomass)] <- 0
  ## Fill X row by row.
  ##
  ## If duplicate site x time rows occur, their biomass values are
  ## added. This preserves total biomass and avoids silently dropping
  ## replicate records.
  for (row_idx in seq_len(nrow(data_wide))) {
    X[time_index[row_idx], site_index[row_idx], ] <-
      X[time_index[row_idx], site_index[row_idx], ] + biomass[row_idx, ]
  }
  X
}
calc_metacommunity_metrics <- function(X) {
  ## Calculate alpha and gamma components first.
  ##
  ## Aggregate variability:
  ##   CV_gamma = regional variability in total metacommunity biomass
  ##   CV_alpha = local variability in total site biomass
  ##
  ## Compositional variability:
  ##   BD_gamma = regional variability in metacommunity composition
  ##   BD_alpha = biomass-weighted local compositional variability
  metrics <- c(
    CV_gamma = CV_gamma(X),
    CV_alpha = CV_alpha(X),
    BD_gamma = BD_gamma(X),
    BD_alpha = BD_alpha(X),
    BD_beta = BD_spatial_weighted(X)
  )
  ## Calculate spatial synchrony components using the multiplicative
  ## partitions:
  ##
  ##   CV_gamma = CV_alpha * CV_phi
  ##   BD_gamma = BD_alpha * BD_phi
  metrics <- c(
    metrics,
    CV_phi = metrics[["CV_gamma"]] / metrics[["CV_alpha"]],
    BD_phi = metrics[["BD_gamma"]] / metrics[["BD_alpha"]]
  )
  ## Return a tidy metric table in a fixed order for downstream joins,
  ## plotting, and comparison across datasets.
  tibble::tibble(
    varname = community_metric_order,
    estimate = unname(metrics[community_metric_order])
  )
}
wide_metric_table <- function(metric_table,
                              value_col = "estimate") {
  ## Convert the long metric table:
  ##   varname | estimate
  ##
  ## to a one-row wide table:
  ##   CV_gamma | CV_alpha | CV_phi | ...
  metric_table %>%
    dplyr::select(varname, dplyr::all_of(value_col)) %>%
    tidyr::pivot_wider(
      names_from = varname,
      values_from = dplyr::all_of(value_col)
    )
}
calc_spatial_bd_by_time <- function(X) {
  ## Return the time-resolved spatial compositional variability
  ## calculated inside spatial_bd_by_time().
  ##
  ## BD:
  ##   spatial compositional variability among sites at each timestep
  ##
  ## total_metacommunity_biomass:
  ##   X_.t., used to weight each timestep
  ##
  ## weights:
  ##   w_t = X_.t. / sum_t X_.t.
  ##
  ## BD_x_wt:
  ##   weighted contribution of each timestep to BD_spatial_weighted
  site_biomass_by_time <- apply(X, c(1, 2), sum, na.rm = TRUE)
  site_relative_biomass <- sweep(
    X,
    c(1, 2),
    site_biomass_by_time,
    "/"
  )
  site_relative_biomass[!is.finite(site_relative_biomass)] <- 0
  site_hellinger <- sqrt(site_relative_biomass)
  BD <- vapply(seq_len(dim(site_hellinger)[1]), function(time_idx) {
    site_by_taxon <- site_hellinger[time_idx, , , drop = FALSE][1, , ]
    taxon_spatial_var <- apply(site_by_taxon, 2, stats::var, na.rm = TRUE)
    sum(taxon_spatial_var, na.rm = TRUE)
  }, numeric(1))
  total_metacommunity_biomass <- rowSums(site_biomass_by_time, na.rm = TRUE)
  weights <- total_metacommunity_biomass /
    sum(total_metacommunity_biomass, na.rm = TRUE)
  tibble::tibble(
    timestep = dimnames(X)$timestep,
    BD = BD,
    total_metacommunity_biomass = total_metacommunity_biomass,
    weights = weights,
    BD_x_wt = BD * weights
  )
}
