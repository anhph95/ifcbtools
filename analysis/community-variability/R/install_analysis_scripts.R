install_ifcb_r_analysis_scripts <- function(destdir = getwd(), overwrite = FALSE) {
  base_url <- "https://raw.githubusercontent.com/anhph95/ifcbtools/main/analysis/community-variability/R/scripts"
  script_files <- c(
    "ifcb_common.R",
    "ifcb_single_season.R",
    "ifcb_power_analysis.R",
    "ifcb_sensitivity_analysis.R",
    "ifcb_seasonal_comparison.R"
  )

  destdir <- normalizePath(destdir, winslash = "/", mustWork = TRUE)
  message("Installing R analysis scripts into:")
  message("  ", destdir)

  for (script_file in script_files) {
    destfile <- file.path(destdir, script_file)
    if (file.exists(destfile) && !overwrite) {
      message("Skipping existing file: ", script_file)
      next
    }

    url <- paste(base_url, script_file, sep = "/")
    tempfile_path <- tempfile(fileext = ".R")
    on.exit(unlink(tempfile_path), add = TRUE)
    utils::download.file(url, tempfile_path, mode = "wb", quiet = TRUE)
    file.copy(tempfile_path, destfile, overwrite = TRUE)
    message("Downloaded: ", script_file)
  }

  message("Done. Open the copied scripts in this working directory and run them interactively.")
  invisible(script_files)
}

install_ifcb_r_analysis_scripts(
  overwrite = isTRUE(getOption("ifcb.analysis.overwrite", FALSE))
)
