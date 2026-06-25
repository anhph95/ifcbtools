## ================================================================
## IFCB single-season workflow
## ================================================================
##
## This project script prepares one seasonal IFCB metacommunity and
## follows it through three linked analyses:
##   1. observed variability metrics from X[time, site, taxon],
##   2. biomass composition and alpha-gamma-phi partition plots,
##   3. NMDS ordination of community composition.
##
## Ecological meaning:
##   time  = repeated annual surveys t within one season
##   site  = NES-LTER stations i, treated as local communities
##   taxon = IFCB taxa j
rm(list = ls())
library(vegan)
library(dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)
library(cowplot)
library(colorspace)
library(ggrepel)
library(grid)
library(readr)
library(community.variability)
source("ifcb_common.R")
logger <- setup_workflow_logging("ifcb_single_season_R")
logger$info("Starting R single-season workflow")
## ------------------------------------------------
## 1. Project settings
## ------------------------------------------------
## Define the data source, focal season, spatial domain, cruise
## priorities, and plotting variables shared by the workflow.
##
## The station list fixes the local-community set i = 1, ..., n_site
## used to construct a balanced metacommunity for each year.
data_dir <- default_data_dir
input_file <- command_arg_value("--input-file", default_carbon_input_file)
season_filter <- "JAS"
station_list <- c(
  "L1", "L2", "L3", "L4", "L5", "L6",
  "L7", "L8", "L9", "L10", "L11"
)
main_cruise <- c(
  "EN608", "AR28B", "EN617", "AR32",
  "EN627", "AR34B", "EN644", "AR39B",
  "EN649", "AR44", "EN655", "EN657",
  "EN661", "AR52B", "EN668", "AR61B",
  "AT46", "AR66B", "EN687", "AR70B",
  "EN695", "HRS2303", "EN706", "AR77",
  "EN712", "EN715", "EN720", "AE2426",
  "EN727", "AR88", "AR92", "AR95",
  "AR99"
)
ctd_cols <- c(
  "temperature", "salinity", "fleco_afl", "nitrate_nitrite",
  "ammonium", "phosphate", "silicate", "latitude", "year"
)
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#332288",
  "#117733", "#AA4499", "#44AA99", "#999933"
)
site_colors <- c(
  setNames(okabe_ito[seq_along(station_list)], station_list),
  Metacommunity = "black"
)
logger$config(list(
  command = logger$command(commandArgs()),
  working_directory = getwd(),
  data_dir = normalizePath(data_dir, winslash = "/", mustWork = FALSE),
  input_file = normalizePath(input_file, winslash = "/", mustWork = FALSE),
  results_dir = normalizePath(results_dir, winslash = "/", mustWork = FALSE),
  season = season_filter,
  station_list = station_list,
  top_taxa_per_station = 3,
  ctd_columns = ctd_cols
))
## ------------------------------------------------
## 2. Load IFCB carbon data and taxonomy
## ------------------------------------------------
## Read carbon biomass observations and the taxonomy table that
## identifies which columns represent taxon biomass.
##
## The taxon columns form the j dimension of:
##   X[t, i, j] = biomass of taxon j at station i in year t
df <- read.csv(
  input_file,
  stringsAsFactors = FALSE
)
colnames(df) <- gsub("\\.", " ", colnames(df))
df$sample_time <- lubridate::ymd_hms(df$sample_time, tz = "UTC")
taxonomy <- read.csv(
  file.path(data_dir, "ifcb_taxonomy.csv"),
  stringsAsFactors = FALSE
)
taxa_cols <- intersect(names(df), taxonomy$Label)
df[taxa_cols] <- lapply(df[taxa_cols], as.numeric)
## Remove the incomplete current/future year and count observations
## within each cruise x cast.
##
## n_obs is used below as a sampling-coverage criterion: among
## candidate casts for the same station-year, prefer the cast with
## the most depth-resolved IFCB observations.
df <- df %>%
  filter(.data$year != 2026) %>%
  group_by(.data$cruise, .data$cast) %>%
  mutate(n_obs = n()) %>%
  ungroup()
## 2020 has no standard L8 cast for JAS. Use the northernmost
## L9 underway sample as an L8 stand-in for the selected season.
##
## This preserves the balanced site set for t = 2020:
##   n_sites(t) = length(station_list)
df <- df %>%
  bind_rows(
    df %>%
      filter(
        .data$year == 2020,
        .data$season == "JAS",
        .data$nearest_station == "L9",
        .data$sample_type == "underway"
      ) %>%
      slice_max(.data$latitude, n = 1, with_ties = FALSE) %>%
      mutate(
        nearest_station = "L8",
        cast = "L8",
        sample_type = "cast_from_udw"
      )
  )
## ------------------------------------------------
## 3. Apply IFCB-specific seasonal filtering
## ------------------------------------------------
## Define one comparable seasonal metacommunity snapshot per year.
##
## For each year x station x cast, keep the shallowest sample so that
## each local community represents near-surface biomass. If multiple
## casts remain for a station-year, prefer the cast with the greatest
## sample coverage, then prefer known main cruises.
##
## The final completeness filter keeps only years with all stations:
##   n_sites(t) = length(station_list)
ds <- df %>%
  filter(.data$nearest_station %in% station_list) %>%
  filter(.data$sample_type %in% c("cast", "cast_from_udw")) %>%
  filter(.data$season == season_filter) %>%
  group_by(.data$year, .data$nearest_station, .data$cast) %>%
  slice_min(.data$depth, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(.data$year, .data$nearest_station) %>%
  filter(.data$n_obs == max(.data$n_obs, na.rm = TRUE)) %>%
  arrange(desc(.data$cruise %in% main_cruise), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(nearest_station = factor(.data$nearest_station, levels = station_list, ordered = TRUE)) %>%
  group_by(.data$year) %>%
  filter(n_distinct(.data$nearest_station) == length(station_list)) %>%
  ungroup()
logger$output(table(ds$year, ds$nearest_station))
## ------------------------------------------------
## 4. Build community tables and calculate observed metrics
## ------------------------------------------------
## Convert the filtered observation table to the metacommunity array:
##   X[t, i, j] = biomass of taxon j at station i in year t
##
## This is the shared input for aggregate variability
## (CV_gamma, CV_alpha), compositional variability
## (BD_gamma, BD_alpha), and their synchrony partitions.
community_wide <- ds %>%
  transmute(
    site = .data$nearest_station,
    timestep = .data$year,
    across(all_of(taxa_cols))
  )
community_array <- make_community_array(community_wide, taxa_cols)
metrics <- calc_metacommunity_metrics(community_array)
logger$output(metrics)
readr::write_csv(
  wide_metric_table(metrics),
  file.path(results_dir, paste0("estimate_", season_filter, ".csv"))
)
logger$info("Saved metric estimates for ", season_filter)
## Plot the observed alpha-gamma-phi partition, excluding BD_beta
## because BD_beta is a time-weighted spatial summary rather than a
## Lamy et al. synchrony partition.
p_partition <- plot_partition_bar(metrics[-nrow(metrics), ])
## ------------------------------------------------
## 5. Add metacommunity rows for composition plots and NMDS
## ------------------------------------------------
## Add regional metacommunity samples by summing taxon biomass across
## stations within each year:
##   X_.tj = sum_i X_itj
##
## These rows visualize regional composition beside local-community
## station composition. Environmental variables are averaged across
## stations for the corresponding regional row.
totals <- ds %>%
  group_by(.data$year) %>%
  summarise(
    across(all_of(taxa_cols), ~ sum(.x, na.rm = TRUE)),
    sample_time = min(.data$sample_time, na.rm = TRUE),
    across(any_of(ctd_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(nearest_station = "Metacommunity")
ds_with_meta <- bind_rows(ds, totals) %>%
  mutate(
    nearest_station = factor(
      .data$nearest_station,
      levels = c(station_list, "Metacommunity"),
      ordered = TRUE
    )
  )
## ------------------------------------------------
## 6. Plot community composition and metric partitioning
## ------------------------------------------------
## Select dominant taxa for plotting by taking the largest local
## biomass contributors at each station, then pooling those names
## across stations.
##
## Remaining taxa are grouped as "Others" so stacked bars preserve
## total biomass while keeping the legend interpretable.
top_species <- ds %>%
  group_by(.data$nearest_station) %>%
  summarise(across(all_of(taxa_cols), sum, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(-nearest_station, names_to = "species", values_to = "total") %>%
  group_by(.data$nearest_station) %>%
  slice_max(.data$total, n = 3, with_ties = TRUE) %>%
  ungroup() %>%
  group_by(.data$species) %>%
  summarise(total = sum(.data$total), .groups = "drop") %>%
  arrange(desc(.data$total)) %>%
  pull(.data$species)
cb_colors <- c(
  "#E69F00", "#D55E00", "#CC79A7", "#F0E442", "#0072B2",
  "#332288", "#56B4E9", "#AA4499", "#44AA99", "#999933",
  "#88CCEE", "#117733", "#DDCC77", "#661100", "#882255",
  "#6699CC", "#AA4466", "#4477AA", "#228833", "#CC6677"
)
if (length(top_species) > length(cb_colors)) {
  cb_colors <- colorspace::qualitative_hcl(length(top_species), palette = "Set 3")
}
fill_colors <- c(setNames(cb_colors[seq_along(top_species)], top_species), Others = "#BDBDBD")
## Convert wide taxon biomass to long composition data for stacked
## station and metacommunity bars.
ds_long <- prepare_composition_long(ds_with_meta, taxa_cols, top_species)
all_years <- sort(unique(ds$year))
site_labels <- ds_long %>%
  filter(.data$site != "Metacommunity") %>%
  distinct(.data$site)
ymax_shared <- ds_long %>%
  filter(.data$site != "Metacommunity") %>%
  summarise(ymax = max(.data$count, na.rm = TRUE)) %>%
  pull(.data$ymax)
y_force <- ds_long %>%
  filter(.data$site != "Metacommunity") %>%
  distinct(.data$site) %>%
  mutate(ymax = ymax_shared)
p_sites <- ggplot(
  ds_long %>% filter(.data$site != "Metacommunity"),
  aes(factor(.data$year, levels = all_years), .data$count, fill = .data$species_grp)
) +
  geom_col(width = 0.75) +
  geom_blank(data = y_force, aes(y = .data$ymax), inherit.aes = FALSE) +
  facet_wrap(~ site, nrow = 4, ncol = 3, dir = "h") +
  scale_fill_manual(values = fill_colors) +
  scale_x_discrete(limits = as.character(seq(min(all_years), max(all_years), by = 1)), drop = FALSE) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 4),
    labels = function(x) round(x / 1e7, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = NULL, y = expression("Carbon (" * "\u00D7" * 10^7 * " ug C L"^{-1} * ")")) +
  box_theme(base_size = 12) +
  theme(
    panel.spacing = unit(0.75, "lines"),
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(),
    axis.ticks.length = unit(0.1, "cm"),
    legend.position = "none"
  ) +
  geom_text(
    data = site_labels,
    aes(x = -Inf, y = Inf, label = .data$site),
    inherit.aes = FALSE,
    hjust = -0.25,
    vjust = 1.25,
    size = 4
  ) +
  coord_cartesian(clip = "off")
p_meta <- ggplot(
  ds_long %>% filter(.data$site == "Metacommunity"),
  aes(factor(.data$year, levels = all_years), .data$count, fill = .data$species_grp)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = fill_colors) +
  scale_x_discrete(limits = as.character(seq(min(all_years), max(all_years), by = 1)), drop = FALSE) +
  scale_y_continuous(
    labels = function(x) round(x / 1e7, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = expression("Carbon (" * "\u00D7" * 10^7 * " ug C L"^{-1} * ")"),
    fill = "Species"
  ) +
  box_theme(base_size = 12) +
  theme(
    legend.position = "none",
    axis.ticks = element_line(),
    axis.ticks.length = unit(0.1, "cm")
  ) +
  annotation_custom(
    grid::textGrob(
      "Metacommunity",
      x = unit(0.02, "npc"),
      y = unit(0.98, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 12)
    )
  ) +
  coord_cartesian(clip = "off")
composition_plot <- plot_grid(
  plot_grid(p_sites, ncol = 1),
  plot_grid(p_meta, p_partition, ncol = 2, align = "hv"),
  extract_legend(p_meta, aesthetic = "fill", nrow = 3),
  ncol = 1,
  rel_heights = c(2, 1, 0.4)
)
composition_plot
ggsave(
  filename = file.path(results_dir, paste0("community_composition_", season_filter, ".png")),
  plot = composition_plot,
  width = 12,
  height = 12,
  dpi = 300
)
## ------------------------------------------------
## 7. Run NMDS ordination and environmental fit
## ------------------------------------------------
## Ordinate community composition in reduced two-dimensional NMDS
## space.
##
## Rows are station-year samples plus annual metacommunity rows; columns
## are taxon biomass variables. The ordination visualizes compositional
## trajectories through time rather than total biomass variability.
community_matrix <- ds_with_meta[, taxa_cols, drop = FALSE]
model <- run_nmds_ordination(community_matrix, seed = 123, trace = TRUE)
ord_sites <- extract_nmds_sites(
  model,
  metadata = ds_with_meta %>%
    transmute(nearest_station = .data$nearest_station, sample_time = .data$sample_time),
  site_col = "nearest_station",
  time_col = "sample_time"
)
path_data <- extract_nmds_path_segments(
  ord_sites,
  site_col = "nearest_station",
  time_col = "sample_time"
)
nmds_path_plot <- plot_nmds_paths(
  ord_sites = ord_sites,
  path_data = path_data,
  color_values = site_colors,
  model_stress = model$stress,
  site_col = "nearest_station"
)
env_idx <- ds_with_meta$nearest_station != "Metacommunity"
## Fit environmental vectors using only local-community samples.
## Metacommunity rows are regional summaries and should not contribute
## as independent station observations in the environmental fit.
ctd <- ds_with_meta[, ctd_cols]
ctd <- as.data.frame(scale(ctd))
sp_scr <- as.data.frame(vegan::scores(model, display = "species"))
sp_scr$species <- rownames(sp_scr)
sp_scr$len <- sqrt(sp_scr[, 1]^2 + sp_scr[, 2]^2)
sp_top <- sp_scr %>%
  arrange(desc(.data$len)) %>%
  slice_head(n = 10) %>%
  rename(NMDS1 = 1, NMDS2 = 2)
env_fit <- fit_env_vectors(model, ctd, subset = env_idx)
env_df <- env_fit$vectors
lim <- max(abs(c(sp_scr[, 1], sp_scr[, 2])), na.rm = TRUE) * 1.1
p_env <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 1) +
  geom_segment(
    data = env_df,
    aes(x = 0, y = 0, xend = .data$NMDS1, yend = .data$NMDS2),
    arrow = arrow(length = unit(0.08, "in"), type = "closed"),
    colour = "red",
    linewidth = 0.75
  ) +
  geom_text_repel(
    data = env_df,
    aes(x = .data$NMDS1, y = .data$NMDS2, label = .data$var),
    colour = "red",
    fontface = "bold",
    size = 3,
    max.overlaps = Inf
  ) +
  geom_point(data = sp_top, aes(x = .data$NMDS1, y = .data$NMDS2), colour = "blue", size = 2) +
  geom_text_repel(
    data = sp_top,
    aes(x = .data$NMDS1, y = .data$NMDS2, label = .data$species),
    colour = "blue",
    fontface = "italic",
    size = 3,
    max.overlaps = Inf
  ) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_bw() +
  theme(panel.grid = element_blank(), aspect.ratio = 1)
nmds_plot <- plot_grid(
  plot_grid(nmds_path_plot, p_env, ncol = 2, align = "hv"),
  extract_legend(nmds_path_plot, aesthetic = "color", nrow = 2),
  ncol = 1,
  rel_heights = c(1, 0.12)
)
nmds_plot
ggsave(
  filename = file.path(results_dir, paste0("nmds_", season_filter, ".png")),
  plot = nmds_plot,
  width = 10,
  height = 5.5,
  dpi = 300
)
logger$info("R single-season workflow completed")
