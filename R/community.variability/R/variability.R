## ================================================================
## Metacommunity variability metrics following Lamy et al. 2021
## ================================================================
##
## Input:
##   X[time, site, taxon]
##
## Ecological meaning:
##   time  = repeated surveys t
##   site  = local communities i
##   taxon = species/taxa j
##
## Aggregate variability metrics use total biomass.
## Compositional variability metrics use Hellinger-transformed
## relative biomass so that compositional change is separated from
## variation in total biomass.
##
## Metrics calculated here:
##   CV_gamma  = regional aggregate variability
##   CV_alpha  = local aggregate variability
##   CV_phi    = spatial aggregate synchrony = CV_gamma / CV_alpha
##
##   BD_gamma  = regional compositional variability
##   BD_alpha  = biomass-weighted local compositional variability
##   BD_phi    = spatial compositional synchrony = BD_gamma / BD_alpha
##
## Additional helper:
##   spatial_bd_by_time = per-time spatial compositional variability
##                        among local communities
CV_gamma <- function(X) {
  ## Regional aggregate variability, CV_gamma^2
  ##
  ## First collapse the full metacommunity array to one regional
  ## biomass time series:
  ##   X_.t. = sum_i sum_j X_itj
  total_metacommunity_biomass <- apply(X, 1, sum, na.rm = TRUE)
  ## CV_gamma^2 is the squared temporal coefficient of variation
  ## of total metacommunity biomass:
  ##   CV_gamma^2 = (sd_t(X_.t.) / mean_t(X_.t.))^2
  mu_TT <- mean(total_metacommunity_biomass, na.rm = TRUE)
  sigma_TT <- stats::sd(total_metacommunity_biomass, na.rm = TRUE)
  (sigma_TT / mu_TT)^2
}
CV_alpha <- function(X) {
  ## Local aggregate variability, CV_alpha^2
  ##
  ## Collapse taxa within each local community to obtain local
  ## total biomass through time:
  ##   X_it. = sum_j X_itj
  site_biomass_by_time <- apply(X, c(1, 2), sum, na.rm = TRUE)
  ## Regional mean biomass provides the common denominator used
  ## to scale the summed local standard deviations.
  total_metacommunity_biomass <- rowSums(site_biomass_by_time, na.rm = TRUE)
  mu_TT <- mean(total_metacommunity_biomass, na.rm = TRUE)
  ## Temporal standard deviation of total biomass within each site:
  ##   sigma_Ti = sd_t(X_it.)
  site_sd <- apply(site_biomass_by_time, 2, stats::sd, na.rm = TRUE)
  ## CV_alpha^2 is the squared sum of local biomass standard
  ## deviations scaled by mean regional biomass:
  ##   CV_alpha^2 = (sum_i sigma_Ti / mu_TT)^2
  (sum(site_sd, na.rm = TRUE) / mu_TT)^2
}
CV_phi <- function(X) {
  ## Spatial aggregate synchrony, phi
  ##
  ## Measures how strongly fluctuations in total biomass are
  ## synchronized among local communities.
  ##
  ## Multiplicative partition:
  ##   CV_gamma^2 = CV_alpha^2 * phi
  ##
  ## Therefore:
  ##   phi = CV_gamma^2 / CV_alpha^2
  CV_gamma(X) / CV_alpha(X)
}
BD_gamma <- function(X) {
  ## Regional compositional variability, BD_gamma^h
  ##
  ## Collapse sites into the metacommunity to obtain regional
  ## biomass of each taxon at each time:
  ##   X_.tj = sum_i X_itj
  regional_taxon_biomass <- apply(X, c(1, 3), sum, na.rm = TRUE)
  ## Regional total biomass at each time:
  ##   X_.t. = sum_j X_.tj
  total_metacommunity_biomass <- rowSums(regional_taxon_biomass, na.rm = TRUE)
  ## Convert regional taxon biomass to regional relative biomass,
  ## then apply the Hellinger transformation:
  ##   z_.tj = sqrt(X_.tj / X_.t.)
  ##
  ## This removes the effect of total biomass and focuses on
  ## changes in relative taxon composition through time.
  regional_relative_biomass <- sweep(
    regional_taxon_biomass,
    1,
    total_metacommunity_biomass,
    "/"
  )
  regional_relative_biomass[!is.finite(regional_relative_biomass)] <- 0
  regional_hellinger <- sqrt(regional_relative_biomass)
  ## BD_gamma^h is the sum across taxa of temporal variances in
  ## Hellinger-transformed regional composition:
  ##   BD_gamma^h = sum_j Var_t(z_.tj)
  taxon_temporal_var <- apply(regional_hellinger, 2, stats::var, na.rm = TRUE)
  sum(taxon_temporal_var, na.rm = TRUE)
}
BD_alpha <- function(X) {
  ## Local compositional variability, BD_alpha^h
  ##
  ## First compute total biomass within each local community and time:
  ##   X_it. = sum_j X_itj
  site_biomass_by_time <- apply(X, c(1, 2), sum, na.rm = TRUE)
  ## Convert each site-by-time sample to taxon relative biomass,
  ## then apply the Hellinger transformation:
  ##   z_itj = sqrt(X_itj / X_it.)
  site_relative_biomass <- sweep(
    X,
    c(1, 2),
    site_biomass_by_time,
    "/"
  )
  site_relative_biomass[!is.finite(site_relative_biomass)] <- 0
  site_hellinger <- sqrt(site_relative_biomass)
  ## For each local community i, calculate compositional variability:
  ##   BD_i^h = sum_j Var_t(z_itj)
  site_bd <- vapply(seq_len(dim(site_hellinger)[2]), function(site_idx) {
    taxon_by_time <- site_hellinger[, site_idx, , drop = FALSE][, 1, ]
    taxon_temporal_var <- apply(taxon_by_time, 2, stats::var, na.rm = TRUE)
    sum(taxon_temporal_var, na.rm = TRUE)
  }, numeric(1))
  ## Weight each local community by its contribution to total
  ## metacommunity biomass:
  ##   w_i = mean_t(X_it.) / sum_i mean_t(X_it.)
  mean_site_biomass <- colMeans(site_biomass_by_time, na.rm = TRUE)
  site_weights <- mean_site_biomass / sum(mean_site_biomass, na.rm = TRUE)
  ## BD_alpha^h is the biomass-weighted mean of local
  ## compositional variability:
  ##   BD_alpha^h = sum_i w_i BD_i^h
  sum(site_bd * site_weights, na.rm = TRUE)
}
BD_phi <- function(X) {
  ## Spatial compositional synchrony, BD_phi^h
  ##
  ## Measures how strongly local communities follow similar
  ## compositional trajectories through time.
  ##
  ## Multiplicative partition:
  ##   BD_gamma^h = BD_alpha^h * BD_phi^h
  ##
  ## Therefore:
  ##   BD_phi^h = BD_gamma^h / BD_alpha^h
  BD_gamma(X) / BD_alpha(X)
}
spatial_bd_by_time <- function(X) {
  ## Per-timestep spatial compositional variability among sites
  ##
  ## This is not the same as BD_phi^h above.
  ## Instead, this calculates spatial compositional differences
  ## among local communities at each time step, then returns the
  ## biomass weights needed to average those values through time.
  ## Total biomass within each site and time:
  ##   X_it. = sum_j X_itj
  site_biomass_by_time <- apply(X, c(1, 2), sum, na.rm = TRUE)
  ## Hellinger-transform local taxon composition:
  ##   z_itj = sqrt(X_itj / X_it.)
  site_relative_biomass <- sweep(
    X,
    c(1, 2),
    site_biomass_by_time,
    "/"
  )
  site_relative_biomass[!is.finite(site_relative_biomass)] <- 0
  site_hellinger <- sqrt(site_relative_biomass)
  ## For each time step t, calculate spatial compositional
  ## variability among sites:
  ##   BD_t^h = sum_j Var_i(z_itj)
  BD <- vapply(seq_len(dim(site_hellinger)[1]), function(time_idx) {
    site_by_taxon <- site_hellinger[time_idx, , , drop = FALSE][1, , ]
    taxon_spatial_var <- apply(site_by_taxon, 2, stats::var, na.rm = TRUE)
    sum(taxon_spatial_var, na.rm = TRUE)
  }, numeric(1))
  ## Weight each time step by its contribution to total
  ## metacommunity biomass:
  ##   w_t = X_.t. / sum_t X_.t.
  total_metacommunity_biomass <- rowSums(site_biomass_by_time, na.rm = TRUE)
  weights <- total_metacommunity_biomass /
    sum(total_metacommunity_biomass, na.rm = TRUE)
  list(
    BD = BD,
    total_metacommunity_biomass = total_metacommunity_biomass,
    weights = weights,
    BD_x_wt = BD * weights
  )
}
BD_spatial_weighted <- function(X) {
  ## Biomass-weighted mean spatial compositional variability
  ## among local communities through time.
  ##
  ## This summarizes the output of spatial_bd_by_time().
  out <- spatial_bd_by_time(X)
  sum(out$BD_x_wt, na.rm = TRUE)
}