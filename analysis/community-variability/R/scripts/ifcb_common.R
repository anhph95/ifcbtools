## ================================================================
## Shared IFCB community-variability workflow settings for R scripts
## ================================================================
##
## Source this file from analysis/community-variability. It defines the
## common project paths and constants used by the R analysis workflows.

analysis_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
repo_dir <- normalizePath(file.path(analysis_dir, "..", ".."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(analysis_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

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

default_data_dir <- file.path(repo_dir, "data", "NESLTER_transect")
community_variability_r_dir <- file.path(repo_dir, "community.variability", "R", "community.variability")
