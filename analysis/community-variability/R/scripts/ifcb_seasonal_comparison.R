## ================================================================
## IFCB seasonal comparison plotting workflow
## ================================================================
##
## This script reads saved bootstrap, sensitivity, and spatial
## compositional variability outputs, then creates cross-season
## comparison plots.
##
## Metrics follow the same notation used in variability.R:
##   CV_gamma, BD_gamma = regional/metacommunity variability
##   CV_alpha, BD_alpha = local/community variability
##   CV_phi,   BD_phi   = spatial synchrony partitions
##   BD_beta           = biomass-weighted spatial BD through time
rm(list = ls())
library(dplyr)
library(tidyr)
library(tidyverse)
library(readr)
library(ggplot2)
library(cowplot)
library(grid)
library(community.variability)
source("R/scripts/ifcb_common.R")
logger <- setup_workflow_logging("ifcb_seasonal_comparison_R")
logger$info("Starting R seasonal-comparison workflow")
options(error = function() {
  error <- geterrmessage()
  logger$error(simpleError(error))
})
## ------------------------------------------------
## 1. Project settings
## ------------------------------------------------
## Define season order and plotting colors. The factor order is used
## throughout so seasonal comparisons are displayed as:
##   JFM -> AMJ -> JAS -> OND
seasons <- c("JFM", "AMJ", "JAS", "OND")
show_bootstrap <- T
point_alpha <- 0.01
plot_top_n_taxa <- TRUE
n_top_taxa <- 10
season_colors <- c(
  JFM = "#0072B2",
  AMJ = "#E69F00",
  JAS = "#009E73",
  OND = "#CC79A7"
)
season_labels <- c(
  JFM = "Winter",
  AMJ = "Spring",
  JAS = "Summer",
  OND = "Fall"
)
taxon_effect_colors <- c(
  Positive = "red",
  Negative = "blue",
  `Not significant` = "gray70"
)
logger$config(list(
  command = logger$command(commandArgs()),
  working_directory = getwd(),
  results_dir = normalizePath(results_dir, winslash = "/", mustWork = FALSE),
  seasons = seasons,
  show_bootstrap = show_bootstrap,
  point_alpha = point_alpha,
  plot_top_n_taxa = plot_top_n_taxa,
  n_top_taxa = n_top_taxa
))
## ------------------------------------------------
## 2. Read bootstrap outputs
## ------------------------------------------------
## Read seasonal bootstrap summaries and replicate-level metric tables.
##
## Each bootstrap row represents metrics recalculated from a resampled
## metacommunity array:
##   X*[t, i, j] = X[sampled t, i, j]
summary_all <- purrr::map_dfr(seasons, function(season) {
  readr::read_csv(
    file.path(results_dir, paste0("summary_", season, ".csv")),
    show_col_types = FALSE
  ) %>%
    mutate(season = season)
})
boot_all <- purrr::map_dfr(seasons, function(season) {
  readr::read_csv(
    file.path(results_dir, paste0("boot_", season, ".csv")),
    show_col_types = FALSE
  ) %>%
    mutate(season = season)
}) %>%
  filter(.data$sample_type == "Bootstrap") %>%
  mutate(season = factor(.data$season, levels = seasons))
summary_wide <- summary_all %>%
  select(season, varname, estimate, lwr, upr) %>%
  mutate(varname = factor(.data$varname, levels = community_metric_order)) %>%
  pivot_wider(
    names_from = varname,
    values_from = c(estimate, lwr, upr),
    names_glue = "{varname}_{.value}"
  ) %>%
  mutate(season = factor(.data$season, levels = seasons))
## ------------------------------------------------
## 3. Read and plot all-season sensitivity analyses
## ------------------------------------------------
## Read leave-one-year and leave-one-taxon outputs, then calculate each
## perturbation's departure from its seasonal baseline:
##   delta = metric_without_component - metric_baseline
##
## Large absolute deltas identify years t or taxa j that strongly
## influence variability or synchrony estimates.
loo_year_all <- purrr::map_dfr(seasons, function(season) {
  readr::read_csv(
    file.path(results_dir, paste0("leave_one_year_out_", season, ".csv")),
    show_col_types = FALSE
  )
}) %>%
  mutate(
    season = factor(.data$season, levels = seasons),
    varname = factor(.data$varname, levels = community_metric_order)
  ) %>%
  add_baseline_delta(removed_col = "year_removed", group_cols = "season")
loo_taxon_all <- purrr::map_dfr(seasons, function(season) {
  readr::read_csv(
    file.path(results_dir, paste0("leave_one_taxon_out_", season, ".csv")),
    show_col_types = FALSE
  )
}) %>%
  mutate(
    season = factor(.data$season, levels = seasons),
    varname = factor(.data$varname, levels = community_metric_order)
  ) %>%
  add_baseline_delta(removed_col = "taxon_removed", group_cols = "season")
plot_vars_sensitivity <- c("BD_gamma", "CV_gamma", "BD_phi", "CV_phi")
## Plot year sensitivity for regional variability and synchrony:
##   X_without_year = X[-t, i, j]
year_levels <- sort(unique(
  loo_year_all$year_removed[loo_year_all$year_removed != "Baseline"]
))
loo_year_plot_data <- loo_year_all %>%
  filter(.data$year_removed != "Baseline", .data$varname %in% plot_vars_sensitivity) %>%
  mutate(
    varname = factor(.data$varname, levels = plot_vars_sensitivity),
    year_removed = factor(.data$year_removed, levels = year_levels),
    ymin = if_else(.data$delta >= 0, .data$baseline, .data$estimate),
    ymax = if_else(.data$delta >= 0, .data$estimate, .data$baseline)
  )
year_dodge_width <- 0.8
year_bar_width <- year_dodge_width / length(year_levels)
year_breaks_global <- loo_year_plot_data %>%
  mutate(
    season_num = as.numeric(.data$season),
    year_num = as.numeric(.data$year_removed),
    year_offset = (
      .data$year_num - (length(year_levels) + 1) / 2
    ) * year_bar_width,
    xmid = .data$season_num + .data$year_offset
  ) %>%
  distinct(.data$xmid, .data$year_removed) %>%
  arrange(.data$xmid)
make_year_delta_plot <- function(dat,
                                 variable_name,
                                 show_year_labels = FALSE,
                                 show_season_top_axis = FALSE) {
  dat_var <- dat %>%
    filter(.data$varname == variable_name) %>%
    mutate(
      season_num = as.numeric(.data$season),
      year_num = as.numeric(.data$year_removed),
      year_offset = (
        .data$year_num - (length(year_levels) + 1) / 2
      ) * year_bar_width,
      xmin = .data$season_num + .data$year_offset - year_bar_width * 0.45,
      xmax = .data$season_num + .data$year_offset + year_bar_width * 0.45
    )
  baseline_dat <- dat_var %>%
    distinct(.data$season, .data$varname, .data$baseline)
  y_min_var <- min(c(dat_var$estimate, dat_var$baseline), na.rm = TRUE)
  y_max_var <- max(c(dat_var$estimate, dat_var$baseline), na.rm = TRUE)
  y_pad_var <- (y_max_var - y_min_var) * 0.18
  if (y_pad_var == 0) y_pad_var <- abs(y_max_var) * 0.18
  if (y_pad_var == 0) y_pad_var <- 1
  ggplot(dat_var) +
    geom_vline(
      xintercept = seq(1.5, length(seasons) - 0.5, by = 1),
      linewidth = 0.35,
      color = "gray70"
    ) +
    geom_segment(
      data = baseline_dat,
      aes(
        x = as.numeric(.data$season) - 0.42,
        xend = as.numeric(.data$season) + 0.42,
        y = .data$baseline,
        yend = .data$baseline
      ),
      linewidth = 0.5,
      color = "black"
    ) +
    geom_rect(
      aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      fill = "gray70",
      color = "black",
      linewidth = 0.25,
      alpha = 0.85
    ) +
    scale_x_continuous(
      breaks = year_breaks_global$xmid,
      labels = if (show_year_labels) {
        as.character(year_breaks_global$year_removed)
      } else {
        rep("", nrow(year_breaks_global))
      },
      sec.axis = dup_axis(
        breaks = seq_along(seasons),
        labels = if (show_season_top_axis) {
          unname(season_labels[seasons])
        } else {
          rep("", length(seasons))
        },
        name = NULL
      ),
      limits = c(0.5, length(seasons) + 0.5),
      expand = c(0, 0)
    ) +
    coord_cartesian(
      ylim = c(
        y_min_var - y_pad_var * 0.15,
        y_max_var + y_pad_var
      ),
      clip = "off"
    ) +
    labs(x = NULL, y = metric_plot_labels[[variable_name]]) +
    box_theme(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks.x.bottom = element_line(),
      axis.text.x.top = element_text(angle = 0, margin = margin(b = 4)),
      axis.ticks.x.top = element_line(),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      axis.ticks.y.left = element_line(),
      axis.ticks.length.y = unit(3, "pt"),
      plot.margin = margin(t = 4, r = 4, b = 4, l = 4)
    )
}
year_delta_plots <- lapply(seq_along(plot_vars_sensitivity), function(i) {
  make_year_delta_plot(
    dat = loo_year_plot_data,
    variable_name = plot_vars_sensitivity[i],
    show_season_top_axis = i == 1,
    show_year_labels = i == length(plot_vars_sensitivity)
  )
})
loo_year_plot <- cowplot::plot_grid(
  plotlist = year_delta_plots,
  ncol = 1,
  align = "hv",
  axis = "tblr"
)
loo_year_plot
ggsave(
  filename = file.path(results_dir, "leave_one_year_all_seasons.png"),
  plot = loo_year_plot,
  width = 14,
  height = 8,
  dpi = 300
)
## ------------------------------------------------
## Taxon sensitivity
## ------------------------------------------------

n_top_taxa <- 10
loo_thresholds <- loo_taxon_all %>%
  filter(taxon_removed != "Baseline") %>%
  group_by(season, varname) %>%
  summarise(
    lwr = quantile(delta, 0.05, na.rm = TRUE),
    upr = quantile(delta, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

loo_taxon_plot_data <- loo_taxon_all %>%
  filter(
    taxon_removed != "Baseline",
    varname %in% plot_vars_sensitivity
  ) %>%
  left_join(
    loo_thresholds,
    by = c("season", "varname")
  ) %>%
  mutate(
    season = factor(season, levels = seasons),
    varname = factor(varname, levels = plot_vars_sensitivity),
    
    significant = delta < lwr | delta > upr,
    
    effect_class = case_when(
      !significant ~ "Not significant",
      delta > 0 ~ "Positive",
      delta < 0 ~ "Negative",
      TRUE ~ "Not significant"
    ),
    
    effect_class = factor(
      effect_class,
      levels = c("Positive", "Negative", "Not significant")
    )
  )

if (plot_top_n_taxa) {
  loo_taxon_plot_data <- loo_taxon_plot_data %>%
    group_by(.data$season, .data$varname) %>%
    slice_max(
      order_by = .data$abs_delta,
      n = n_top_taxa,
      with_ties = FALSE
    ) %>%
    ungroup()
}

make_taxon_value_plot <- function(dat,
                                  variable_name,
                                  show_taxon_labels = FALSE,
                                  show_season_top_axis = FALSE,
                                  show_effect_legend = FALSE) {
  
  dat_var <- dat %>%
    filter(.data$varname == variable_name) %>%
    group_by(.data$season) %>%
    arrange(.data$taxon_removed, .by_group = TRUE) %>%
    mutate(
      taxon_rank = row_number()
    ) %>%
    ungroup()
  
  taxon_slots_var <- max(dat_var$taxon_rank)
  
  taxon_bar_width <- 0.8 / taxon_slots_var
  
  dat_var <- dat_var %>%
    mutate(
      season_num = as.numeric(.data$season),
      
      taxon_offset =
        (.data$taxon_rank -
           (taxon_slots_var + 1) / 2) *
        taxon_bar_width,
      
      xmin =
        .data$season_num +
        .data$taxon_offset -
        taxon_bar_width * 0.45,
      
      xmax =
        .data$season_num +
        .data$taxon_offset +
        taxon_bar_width * 0.45,
      
      xmid = (.data$xmin + .data$xmax) / 2
    )
  
  taxon_breaks_var <- dat_var %>%
    distinct(
      .data$xmid,
      .data$taxon_removed
    ) %>%
    arrange(.data$xmid)
  
  baseline_dat <- dat_var %>%
    distinct(
      .data$season,
      .data$baseline
    )
  
  y_min_var <- min(
    c(dat_var$estimate,
      dat_var$baseline),
    na.rm = TRUE
  )
  
  y_max_var <- max(
    c(dat_var$estimate,
      dat_var$baseline),
    na.rm = TRUE
  )
  
  y_pad_var <- (y_max_var - y_min_var) * 0.18
  
  if (y_pad_var == 0)
    y_pad_var <- abs(y_max_var) * 0.18
  
  if (y_pad_var == 0)
    y_pad_var <- 1
  
  ggplot(dat_var) +
    
    geom_vline(
      xintercept = seq(
        1.5,
        length(seasons) - 0.5,
        by = 1
      ),
      linewidth = 0.35,
      color = "gray70"
    ) +
    
    geom_segment(
      data = baseline_dat,
      aes(
        x = as.numeric(.data$season) - 0.42,
        xend = as.numeric(.data$season) + 0.42,
        y = .data$baseline,
        yend = .data$baseline
      ),
      linewidth = 0.5,
      color = "black"
    ) +
    
    geom_point(
      aes(
        x = .data$xmid,
        y = .data$estimate,
        color = .data$effect_class
      ),
      size = 2.5,
      alpha = 0.95
    ) +
    
    scale_color_manual(
      values = taxon_effect_colors,
      limits = c(
        "Positive",
        "Negative",
        "Not significant"
      ),
      drop = FALSE
    ) +
    
    scale_x_continuous(
      breaks = taxon_breaks_var$xmid,
      
      labels = if (show_taxon_labels) {
        taxon_breaks_var$taxon_removed
      } else {
        rep("", nrow(taxon_breaks_var))
      },
      
      sec.axis = dup_axis(
        breaks = seq_along(seasons),
        
        labels = if (show_season_top_axis) {
          unname(season_labels[seasons])
        } else {
          rep("", length(seasons))
        },
        
        name = NULL
      ),
      
      limits = c(
        0.5,
        length(seasons) + 0.5
      ),
      
      expand = c(0, 0)
    ) +
    
    coord_cartesian(
      ylim = c(
        y_min_var - y_pad_var * 0.15,
        y_max_var + y_pad_var
      ),
      clip = "off"
    ) +
    
    labs(
      x = NULL,
      y = metric_plot_labels[[variable_name]]
    ) +
    
    box_theme(base_size = 12) +
    
    theme(
      legend.position =
        if (show_effect_legend)
          "bottom"
      else
        "none",
      
      legend.title = element_blank(),
      
      panel.grid.major.x =
        element_blank(),
      
      panel.grid.minor.x =
        element_blank(),
      
      axis.text.x.bottom =
        element_text(
          angle = 45,
          hjust = 1,
          vjust = 0.5,
          size = 6
        ),
      
      axis.ticks.x.bottom =
        element_line(),
      
      axis.text.x.top =
        element_text(
          angle = 0,
          margin = margin(b = 4)
        ),
      
      axis.ticks.x.top =
        element_line(),
      
      axis.title.y =
        element_text(
          angle = 90,
          vjust = 0.5
        ),
      
      axis.ticks.y.left =
        element_line(),
      
      axis.ticks.length.y =
        unit(3, "pt"),
      
      plot.margin =
        margin(
          t = 4,
          r = 4,
          b = 4,
          l = 4
        )
    )
}

taxon_value_plots <- lapply(
  seq_along(plot_vars_sensitivity),
  function(i) {
    make_taxon_value_plot(
      dat = loo_taxon_plot_data,
      variable_name = plot_vars_sensitivity[i],
      show_season_top_axis = i == 1,
      show_taxon_labels = TRUE,
      show_effect_legend = FALSE
    )
  }
)

taxon_effect_legend <- cowplot::get_legend(
  make_taxon_value_plot(
    dat = loo_taxon_plot_data,
    variable_name = plot_vars_sensitivity[1],
    show_effect_legend = TRUE
  )
)

loo_taxon_plot <- cowplot::plot_grid(
  cowplot::plot_grid(
    plotlist = taxon_value_plots,
    ncol = 1,
    align = "hv",
    axis = "tblr"
  ),
  taxon_effect_legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

loo_taxon_plot

ggsave(
  filename = file.path(
    results_dir,
    "leave_one_taxon_top_effects_all_seasons.png"
  ),
  plot = loo_taxon_plot,
  width = 16,
  height = 9,
  dpi = 300
)
## ------------------------------------------------
## 4. Build metric comparison plots
## ------------------------------------------------
## Build cross-season metric-space plots from bootstrap replicates.
##
## For each pair of metrics, a 95% covariance ellipse summarizes the
## bootstrap cloud:
##   (x - mu)' S^-1 (x - mu) = chi^2_0.95, df = 2
##
## Points show observed seasonal estimates; ellipses show uncertainty
## from timestep resampling.
theta <- seq(0, 2 * pi, length.out = 200)
make_ellipse <- function(data, season_levels, xvar, yvar, theta) {
  purrr::map_dfr(season_levels, function(s) {
    d <- data %>%
      filter(.data$season == s) %>%
      select(all_of(c(xvar, yvar))) %>%
      tidyr::drop_na()
    if (nrow(d) < 3) {
      return(tibble())
    }
    S <- stats::cov(d)
    if (any(!is.finite(S)) || det(S) <= 0) {
      return(tibble())
    }
    mu <- colMeans(d)
    r2 <- stats::qchisq(0.95, df = 2)
    E <- eigen(S)
    A <- E$vectors %*% diag(sqrt(pmax(E$values, 0) * r2))
    pts <- t(mu + A %*% rbind(cos(theta), sin(theta)))
    tibble(season = s, x = pts[, 1], y = pts[, 2])
  }) %>%
    mutate(season = factor(.data$season, levels = seasons))
}
make_ellipse_lims <- function(ellipse_df, point_df, xvar, yvar, pad = 0.08) {
  ## Use a common square plotting window for each metric pair so the
  ## x and y axes have comparable visual scale under coord_equal().
  xr <- range(c(ellipse_df$x, point_df[[xvar]]), na.rm = TRUE)
  yr <- range(c(ellipse_df$y, point_df[[yvar]]), na.rm = TRUE)
  span <- max(diff(xr), diff(yr))
  if (!is.finite(span) || span == 0) span <- max(abs(c(xr, yr)), na.rm = TRUE)
  if (!is.finite(span) || span == 0) span <- 1
  list(
    xlim = mean(xr) + c(-0.5, 0.5) * span * (1 + pad),
    ylim = mean(yr) + c(-0.5, 0.5) * span * (1 + pad)
  )
}
make_center <- function(data, xvar, yvar) {
  data %>%
    summarise(
      x = mean(.data[[xvar]], na.rm = TRUE),
      y = mean(.data[[yvar]], na.rm = TRUE)
    )
}
bootstrap_layers <- function(ellipse_data, xvar, yvar) {
  ## Draw bootstrap ellipses and optionally the replicate points used
  ## to estimate each ellipse.
  if (!show_bootstrap) return(list())
  list(
    geom_path(
      data = ellipse_data,
      aes(x = .data$x, y = .data$y, color = .data$season, group = .data$season),
      linewidth = 0.7,
      alpha = 0.5
    ),
    geom_point(
      data = boot_all,
      aes(x = .data[[xvar]], y = .data[[yvar]], color = .data$season),
      alpha = point_alpha,
      size = 1
    )
  )
}
make_metric_plot <- function(x_boot, y_boot, x_summary, y_summary,
                             ellipse_df, center_df, lims,
                             x_label, y_label) {
  ## Combine uncertainty, observed seasonal estimates, and grand
  ## bootstrap centers for one two-metric comparison.
  ggplot() +
    bootstrap_layers(ellipse_df, x_boot, y_boot) +
    geom_point(
      data = summary_wide,
      aes(x = .data[[x_summary]], y = .data[[y_summary]], fill = .data$season),
      color = "black",
      shape = 21,
      size = 3
    ) +
    geom_text(
      data = summary_wide,
      aes(x = .data[[x_summary]], y = .data[[y_summary]], label = .data$season, color = .data$season),
      fontface = "bold",
      vjust = -0.8,
      show.legend = FALSE
    ) +
    geom_vline(xintercept = center_df$x, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = center_df$y, linetype = "dashed", color = "grey50") +
    labs(x = x_label, y = y_label) +
    coord_equal(xlim = lims$xlim, ylim = lims$ylim) +
    theme_classic(base_size = 15) +
    theme(legend.position = "none") +
    scale_color_manual(values = season_colors) +
    scale_fill_manual(values = season_colors)
}
ell1 <- make_ellipse(boot_all, seasons, "BD_alpha", "CV_alpha", theta)
ell2 <- make_ellipse(boot_all, seasons, "BD_phi", "CV_phi", theta)
ell3 <- make_ellipse(boot_all, seasons, "BD_gamma", "CV_gamma", theta)
ell4 <- make_ellipse(boot_all, seasons, "BD_beta", "BD_gamma", theta)
p_local <- make_metric_plot(
  "BD_alpha", "CV_alpha", "BD_alpha_estimate", "CV_alpha_estimate",
  ell1, make_center(boot_all, "BD_alpha", "CV_alpha"),
  make_ellipse_lims(ell1, summary_wide, "BD_alpha_estimate", "CV_alpha_estimate"),
  bquote("Local compositional variability - " * BD[alpha]^h),
  bquote("Local aggregate variability - " * CV[alpha]^2)
)
p_sync <- make_metric_plot(
  "BD_phi", "CV_phi", "BD_phi_estimate", "CV_phi_estimate",
  ell2, make_center(boot_all, "BD_phi", "CV_phi"),
  make_ellipse_lims(ell2, summary_wide, "BD_phi_estimate", "CV_phi_estimate"),
  bquote("Compositional synchrony - " * BD[phi[1]]^h),
  bquote("Aggregate synchrony - " * phi[1])
)
p_gamma <- make_metric_plot(
  "BD_gamma", "CV_gamma", "BD_gamma_estimate", "CV_gamma_estimate",
  ell3, make_center(boot_all, "BD_gamma", "CV_gamma"),
  make_ellipse_lims(ell3, summary_wide, "BD_gamma_estimate", "CV_gamma_estimate"),
  bquote("Metacommunity compositional variability - " * BD[gamma]^h),
  bquote("Metacommunity aggregate variability - " * CV[gamma]^2)
)
p_spatial <- make_metric_plot(
  "BD_beta", "BD_gamma", "BD_beta_estimate", "BD_gamma_estimate",
  ell4, make_center(boot_all, "BD_beta", "BD_gamma"),
  make_ellipse_lims(ell4, summary_wide, "BD_beta_estimate", "BD_gamma_estimate"),
  bquote(atop("Compositional variability among local communities", "through time (spatial) - " * BD[beta]^h)),
  bquote(atop("Compositional metacommunity", "variability (temporal) - " * BD[gamma]^h))
)
metric_comparison_plot <- plot_grid(
  p_local, p_gamma,
  p_sync, p_spatial,
  ncol = 2,
  labels = c("A", "C", "B", "D"),
  label_size = 14,
  align = "hv"
)
metric_comparison_plot
ggsave(
  filename = file.path(results_dir, "bootstrap_metric_comparison.png"),
  plot = metric_comparison_plot,
  width = 11,
  height = 11,
  dpi = 1200
)
## ------------------------------------------------
## 5. Plot selected variables together
## ------------------------------------------------
## Plot selected BD relationships directly. These panels compare local,
## spatial, and regional compositional variability:
##   BD_alpha^h = local compositional variability
##   BD_beta^h  = spatial compositional variability through time
##   BD_gamma^h = regional compositional variability
make_variable_pair_plot <- function(x_boot, y_boot, x_summary, y_summary, x_label, y_label) {
  ggplot() +
    {
      if (show_bootstrap) {
        geom_point(
          data = boot_all,
          aes(x = .data[[x_boot]], y = .data[[y_boot]], color = .data$season),
          alpha = 0.12,
          size = 1
        )
      }
    } +
    geom_point(
      data = summary_wide,
      aes(x = .data[[x_summary]], y = .data[[y_summary]], fill = .data$season),
      color = "black",
      shape = 21,
      size = 3.5
    ) +
    geom_text(
      data = summary_wide,
      aes(x = .data[[x_summary]], y = .data[[y_summary]], label = .data$season, color = .data$season),
      fontface = "bold",
      vjust = -0.8,
      show.legend = FALSE
    ) +
    scale_color_manual(values = season_colors) +
    scale_fill_manual(values = season_colors) +
    coord_equal() +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") +
    labs(x = x_label, y = y_label)
}
pair_bd_gamma_alpha <- make_variable_pair_plot(
  x_boot = "BD_alpha",
  y_boot = "BD_gamma",
  x_summary = "BD_alpha_estimate",
  y_summary = "BD_gamma_estimate",
  x_label = expression(BD[alpha]^h),
  y_label = expression(BD[gamma]^h)
)
pair_bd_gamma_beta <- make_variable_pair_plot(
  x_boot = "BD_beta",
  y_boot = "BD_gamma",
  x_summary = "BD_beta_estimate",
  y_summary = "BD_gamma_estimate",
  x_label = expression(BD[beta]^h),
  y_label = expression(BD[gamma]^h)
)
variable_pair_plot <- plot_grid(
  pair_bd_gamma_alpha,
  pair_bd_gamma_beta,
  nrow = 1,
  labels = c("A", "B"),
  label_size = 14,
  align = "hv"
)
variable_pair_plot
ggsave(
  filename = file.path(results_dir, "variable_pairs_bd_gamma.png"),
  plot = variable_pair_plot,
  width = 10,
  height = 5,
  dpi = 600
)
## ------------------------------------------------
## 6. Plot spatial compositional variability through time
## ------------------------------------------------
## Read per-timestep spatial compositional variability:
##   BD_t^h = sum_j Var_i(z_itj)
##
## This quantity measures among-station compositional differences at
## each timestep. It is distinct from BD_phi, which is a synchrony ratio
## in the temporal alpha-gamma partition.
spavar_all <- purrr::map_dfr(seasons, function(season) {
  readr::read_csv(
    file.path(results_dir, paste0("spatial_variance_", season, ".csv")),
    show_col_types = FALSE
  ) %>%
    mutate(season = season, year = as.numeric(.data$timestep))
})
spavar_all <- spavar_all %>%
  mutate(
    season = factor(.data$season, levels = seasons, ordered = TRUE),
    ## Place seasons within each year on a quarter-year axis so the
    ## time series preserves seasonal order inside each annual cycle.
    time_index = .data$year + (as.numeric(.data$season) - 1) / 4,
    biomass_scaled = .data$total_metacommunity_biomass / 1e8
  ) %>%
  arrange(.data$time_index)
year_breaks <- seq(
  floor(min(spavar_all$year, na.rm = TRUE)),
  ceiling(max(spavar_all$year, na.rm = TRUE)),
  by = 1
)
legend_plot <- ggplot(spavar_all, aes(x = .data$time_index, y = .data$biomass_scaled, color = .data$season)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = season_colors) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_blank())
p_biomass <- ggplot(spavar_all, aes(x = .data$time_index, y = .data$biomass_scaled)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(color = .data$season), size = 3.5) +
  scale_color_manual(values = season_colors) +
  scale_x_continuous(breaks = year_breaks, labels = year_breaks) +
  labs(
    x = NULL,
    y = expression(paste("Metacommunity biomass (10"^8, " ", mu, "g C L"^{-1}, ")"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_line(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )
p_spatial_time <- ggplot(spavar_all, aes(x = .data$time_index, y = .data$BD)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(color = .data$season), size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = season_colors) +
  scale_x_continuous(breaks = year_breaks, labels = year_breaks) +
  labs(x = "Year", y = expression("Variability among sites (" * BD[beta]^h * ")")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")
spatial_time_plot <- plot_grid(
  plot_grid(
    p_biomass,
    p_spatial_time,
    ncol = 1,
    align = "v",
    labels = c("A", "B"),
    label_size = 14,
    rel_heights = c(1, 1)
  ),
  cowplot::get_legend(legend_plot),
  ncol = 1,
  rel_heights = c(1, 0.08)
)
spatial_time_plot
ggsave(
  filename = file.path(results_dir, "bootstrap_spatial_variance_timeseries.png"),
  plot = spatial_time_plot,
  width = 8,
  height = 7,
  dpi = 1200
)
logger$info("R seasonal-comparison workflow completed")
