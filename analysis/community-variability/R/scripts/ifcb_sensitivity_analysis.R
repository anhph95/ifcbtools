## ================================================================
## IFCB sensitivity analysis workflow
## ================================================================
##
## This script evaluates how strongly each year t and each taxon j
## influences metacommunity variability metrics calculated from:
##   X[time, site, taxon]
##
## Leave-one-year removes one timestep from X. Leave-one-taxon removes
## one taxon slice from X. Metrics are recalculated after each removal
## and compared with the seasonal baseline.
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
## Define the shared seasonal and spatial domain for all sensitivity
## analyses. The station list fixes the local-community dimension i,
## and the season loop repeats the same analysis for each quarter.
data_dir <- default_data_dir
seasons <- c("JFM", "AMJ", "JAS", "OND")
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
## Read carbon biomass observations and identify taxon biomass columns.
##
## These taxon columns form dimension j in:
##   X[t, i, j] = biomass of taxon j at station i in year t
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
## n_obs is used below to prefer the most sampled surface cast when
## multiple casts are available for the same station-year.
df <- df %>%
  filter(.data$year != 2026) %>%
  group_by(.data$cruise, .data$cast) %>%
  mutate(n_obs = n()) %>%
  ungroup()
## 2020 has no standard L8 cast for JAS. Use the northernmost
## L9 underway sample as an L8 stand-in for the selected season.
##
## This keeps the annual JAS metacommunity balanced:
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
## Pre-filter rows that are shared by all seasonal sensitivity analyses.
## The season loop then only applies the seasonal and completeness
## filters, avoiding repeated station/sample-type filtering.
df_analysis <- df %>%
  filter(.data$nearest_station %in% station_list) %>%
  filter(.data$sample_type %in% c("cast", "cast_from_udw"))
## ------------------------------------------------
## 3. Run sensitivity analyses for each season
## ------------------------------------------------
## For each season, build one balanced metacommunity array and then
## recalculate metrics after removing one timestep or one taxon.
##
## The baseline row is retained so each perturbation can be expressed
## as:
##   delta = metric_without_component - metric_baseline
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
  ## Leave-one-year slices remove X[t, , ]. Leave-one-taxon slices
  ## remove X[, , j]. Each result is compared with the full-array
  ## baseline for the same season.
  community_array <- make_community_array(community_wide, taxa_cols)
  message("Leave-one-year")
  readr::write_csv(
    leave_one_out(community_array, margin = "timestep", show_progress = TRUE) %>%
      rename(year_removed = timestep_removed) %>%
      mutate(season = season_filter, .before = 1),
    file.path(results_dir, paste0("leave_one_year_out_", season_filter, ".csv"))
  )
  message("Leave-one-taxon")
  readr::write_csv(
    leave_one_out(community_array, margin = "taxon", show_progress = TRUE) %>%
      mutate(season = season_filter, .before = 1),
    file.path(results_dir, paste0("leave_one_taxon_out_", season_filter, ".csv"))
  )
}
