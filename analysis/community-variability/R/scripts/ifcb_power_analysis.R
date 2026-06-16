## ================================================================
## IFCB power analysis workflow
## ================================================================
##
## This script runs seasonal bootstrap analyses for the metacommunity
## variability metrics calculated from:
##   X[time, site, taxon]
##
## Each bootstrap replicate resamples timesteps t, then recalculates
## aggregate variability, compositional variability, and synchrony
## partitions. Outputs feed the cross-season comparison plots.
rm(list = ls())
library(dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)
library(readr)
library(community.variability)
source("R/scripts/ifcb_common.R")
## ------------------------------------------------
## 1. Project settings
## ------------------------------------------------
## Define the shared spatial domain, seasons, and bootstrap design.
##
## Resampling occurs along the time dimension:
##   X*[t, i, j] = X[sampled t, i, j]
##
## Sites and taxa remain paired within each resampled year so that the
## spatial and compositional structure of each annual metacommunity is
## preserved.
data_dir <- default_data_dir
seasons <- c("JFM", "AMJ", "JAS", "OND")
n_boot <- 1000
bootstrap_seed <- 123
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
## ------------------------------------------------
## 2. Load IFCB carbon data and taxonomy
## ------------------------------------------------
## Read carbon biomass observations and identify the IFCB taxon
## columns that define the taxon dimension j in X[t, i, j].
df <- read.csv(
  file.path(data_dir, "ifcb_carbon_mix.csv"),
  stringsAsFactors = FALSE
)
colnames(df) <- gsub("\\.", " ", colnames(df))
df$sample_time <- lubridate::ymd_hms(df$sample_time, tz = "UTC")
taxonomy <- read.csv(
  file.path(data_dir, "ifcb_taxonomy_mix.csv"),
  stringsAsFactors = FALSE
)
taxa_cols <- intersect(names(df), taxonomy$Label)
df[taxa_cols] <- lapply(df[taxa_cols], as.numeric)
## Remove the incomplete current/future year and count observations
## within each cruise x cast.
##
## n_obs is used as a sampling-coverage criterion when selecting one
## representative surface cast per station-year.
df <- df %>%
  filter(.data$year != 2026) %>%
  group_by(.data$cruise, .data$cast) %>%
  mutate(n_obs = n()) %>%
  ungroup()
## 2020 has no standard L8 cast for JAS. Use the northernmost
## L9 underway sample as an L8 stand-in for the selected season.
##
## This preserves a balanced station set for JAS 2020:
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
## Pre-filter rows that are shared by all seasonal bootstrap analyses.
## The season loop then only applies the seasonal and completeness
## filters, avoiding repeated station/sample-type filtering.
df_analysis <- df %>%
  filter(.data$nearest_station %in% station_list) %>%
  filter(.data$sample_type %in% c("cast", "cast_from_udw"))
## ------------------------------------------------
## 3. Bootstrap each season and write outputs
## ------------------------------------------------
## For each season, build a balanced annual metacommunity array and
## bootstrap timesteps t. The same workflow is repeated by season so
## seasonal differences reflect data and metric uncertainty rather than
## different preprocessing rules.
for (season_filter in seasons) {
  message("Processing season: ", season_filter)
  ## Define one comparable seasonal metacommunity snapshot per year.
  ##
  ## For each year x station x cast, keep the shallowest sample so the
  ## local community represents near-surface biomass. Then choose one
  ## station-year cast by maximizing n_obs and preferring known main
  ## cruises when ties remain.
  ##
  ## Keep only complete years:
  ##   n_sites(t) = length(station_list)
  ds <- df_analysis %>%
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
  message("Complete years retained: ", n_distinct(ds$year))
  community_wide <- ds %>%
    transmute(
      site = .data$nearest_station,
      timestep = .data$year,
      across(all_of(taxa_cols))
    )
  ## Convert the filtered observations into:
  ##   X[t, i, j] = biomass of taxon j at station i in year t
  ##
  ## Bootstrap replicates resample t, preserving the full site x taxon
  ## community matrix within each selected year.
  community_array <- make_community_array(community_wide, taxa_cols)
  boot_res <- bootstrap_by_dimension(
    X = community_array,
    margin = "timestep",
    n_boot = n_boot,
    seed = bootstrap_seed,
    baseline_in_boot = TRUE,
    show_progress = TRUE
  )
  readr::write_csv(
    boot_res$boot,
    file.path(results_dir, paste0("boot_", season_filter, ".csv"))
  )
  readr::write_csv(
    boot_res$summary,
    file.path(results_dir, paste0("summary_", season_filter, ".csv"))
  )
  readr::write_csv(
    calc_spatial_bd_by_time(community_array),
    file.path(results_dir, paste0("spatial_variance_", season_filter, ".csv"))
  )
}
