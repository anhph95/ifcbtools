## ================================================================
## Generic community plotting helpers
## ================================================================
metric_plot_labels <- c(
  CV_gamma = expression(CV[gamma]^2),
  CV_alpha = expression(CV[alpha]^2),
  CV_phi = expression(phi1),
  BD_gamma = expression(BD[gamma]^h),
  BD_alpha = expression(BD[alpha]^h),
  BD_phi = expression(BD[phi1]^h),
  BD_beta = expression(BD[beta]^h)
)
box_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.4
      ),
      panel.background = ggplot2::element_blank()
    )
}
extract_legend <- function(p, aesthetic = "fill", nrow = 3) {
  guide <- if (aesthetic == "color") {
    ggplot2::guides(color = ggplot2::guide_legend(nrow = nrow))
  } else {
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = nrow, byrow = TRUE))
  }
  ggplot2::ggplotGrob(
    p +
      guide +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_text(face = "bold", size = 12),
        legend.text = ggplot2::element_text(size = 12),
        legend.key.size = grid::unit(0.4, "cm"),
        legend.spacing.y = grid::unit(0.1, "cm"),
        legend.box.spacing = grid::unit(0.1, "cm")
      )
  ) |>
    gtable::gtable_filter("guide-box")
}
plot_partition_bar <- function(metric_table,
                               labels = metric_plot_labels,
                               y_max = NULL) {
  plot_data <- metric_table %>%
    dplyr::mutate(
      varname = factor(.data$varname, levels = community_metric_order),
      group = dplyr::if_else(
        grepl("^CV_", as.character(.data$varname)),
        "orange",
        "skyblue"
      )
    )
  if (is.null(y_max)) {
    y_max <- max(1, max(plot_data$estimate, na.rm = TRUE) * 1.05)
  }
  ggplot2::ggplot(plot_data, ggplot2::aes(.data$varname, .data$estimate, fill = .data$group)) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$estimate)),
      vjust = -0.5,
      size = 3
    ) +
    ggplot2::scale_fill_manual(values = c(orange = "orange", skyblue = "skyblue")) +
    ggplot2::scale_x_discrete(labels = labels) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_max),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    box_theme(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.ticks = ggplot2::element_line(),
      axis.ticks.length = grid::unit(0.1, "cm")
    ) +
    ggplot2::annotation_custom(
      grid::textGrob(
        "Partitioning",
        x = grid::unit(0.02, "npc"),
        y = grid::unit(0.98, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 12)
      )
    )
}
prepare_composition_long <- function(data_with_meta,
                                     taxa_cols,
                                     top_species,
                                     site_col = "nearest_station",
                                     year_col = "year") {
  data_with_meta %>%
    dplyr::select(
      site = dplyr::all_of(site_col),
      year = dplyr::all_of(year_col),
      dplyr::all_of(taxa_cols)
    ) %>%
    tidyr::pivot_longer(
      dplyr::all_of(taxa_cols),
      names_to = "species",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      species_grp = dplyr::if_else(.data$species %in% top_species, .data$species, "Others")
    ) %>%
    dplyr::group_by(.data$site, .data$year, .data$species_grp) %>%
    dplyr::summarise(count = sum(.data$count, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(species_grp = factor(.data$species_grp, levels = c(top_species, "Others")))
}
plot_leave_one_delta <- function(results,
                                 removed_col,
                                 plot_vars = c("BD_gamma", "CV_gamma", "BD_phi", "CV_phi"),
                                 group_col = "season",
                                 group_levels = NULL,
                                 labels = metric_plot_labels) {
  if (is.null(group_levels)) {
    group_levels <- unique(results[[group_col]])
  }
  plot_data <- results %>%
    dplyr::filter(.data$varname %in% plot_vars) %>%
    dplyr::mutate(
      "{group_col}" := factor(.data[[group_col]], levels = group_levels),
      varname = factor(.data$varname, levels = plot_vars)
    ) %>%
    add_baseline_delta(removed_col = removed_col, group_cols = group_col) %>%
    dplyr::filter(.data[[removed_col]] != "Baseline")
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[removed_col]], y = .data$delta, fill = .data[[group_col]])
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::geom_col(position = "dodge", color = "black", linewidth = 0.25) +
    ggplot2::facet_wrap(~ varname, scales = "free_y", labeller = ggplot2::as_labeller(labels)) +
    box_theme(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    ggplot2::labs(x = NULL, y = "Change from baseline", fill = NULL)
}
plot_nmds_paths <- function(ord_sites,
                            path_data,
                            color_values,
                            model_stress,
                            site_col = "nearest_station",
                            metacommunity_label = "Metacommunity") {
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 1) +
    ggplot2::geom_point(
      data = path_data$start_points,
      ggplot2::aes(x = .data$Axis1, y = .data$Axis2, col = .data[[site_col]]),
      size = 3
    ) +
    ggplot2::geom_segment(
      data = path_data$all_segments %>% dplyr::filter(.data[[site_col]] != metacommunity_label),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, col = .data[[site_col]]),
      linewidth = 0.75,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = path_data$end_segments %>% dplyr::filter(.data[[site_col]] != metacommunity_label),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, col = .data[[site_col]]),
      arrow = grid::arrow(length = grid::unit(0.25, "cm"), type = "closed"),
      linewidth = 0.75,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = path_data$all_segments %>% dplyr::filter(.data[[site_col]] == metacommunity_label),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "black",
      linewidth = 1.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = path_data$end_segments %>% dplyr::filter(.data[[site_col]] == metacommunity_label),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      arrow = grid::arrow(length = grid::unit(0.25, "cm"), type = "closed"),
      color = "black",
      linewidth = 1.25,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = color_values) +
    ggplot2::labs(x = "NMDS1", y = "NMDS2") +
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0("stress = ", round(model_stress, 3)),
      hjust = 1.1,
      vjust = 1.2,
      size = 3.5,
      fontface = "bold"
    )
}
