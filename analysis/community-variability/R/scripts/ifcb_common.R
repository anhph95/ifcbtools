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

## ------------------------------------------------
## Shared dependency-free workflow logging
## ------------------------------------------------
## Each workflow writes timestamped informational and error logs under
## results/logs while retaining the same messages in the terminal.
setup_workflow_logging <- function(workflow_name, log_dir = file.path(results_dir, "logs")) {
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_path <- file.path(log_dir, paste0(workflow_name, "_", timestamp, ".out.log"))
  err_path <- file.path(log_dir, paste0(workflow_name, "_", timestamp, ".err.log"))
  file.create(out_path)
  file.create(err_path)

  format_log_line <- function(level, text) {
    paste(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      level,
      workflow_name,
      text,
      sep = " | "
    )
  }

  log_message <- function(..., level = "INFO") {
    text <- paste0(..., collapse = "")
    line <- format_log_line(level, text)
    cat(line, "\n", file = stdout())
    cat(line, "\n", file = out_path, append = TRUE)
    invisible(line)
  }

  redact_command_args <- function(args) {
    secret_terms <- c("password", "token", "secret", "api_key", "apikey")
    redacted <- character()
    redact_next <- FALSE
    for (argument in args) {
      if (redact_next) {
        redacted <- c(redacted, "<redacted>")
        redact_next <- FALSE
        next
      }
      option <- sub("=.*$", "", argument)
      is_secret <- any(vapply(secret_terms, grepl, logical(1), x = tolower(option), fixed = TRUE))
      if (is_secret && grepl("=", argument, fixed = TRUE)) {
        redacted <- c(redacted, paste0(option, "=<redacted>"))
      } else {
        redacted <- c(redacted, argument)
        redact_next <- is_secret
      }
    }
    paste(redacted, collapse = " ")
  }

  log_output <- function(object, level = "INFO") {
    text <- paste(capture.output(print(object)), collapse = "\n")
    log_message(text, level = level)
  }

  log_config <- function(settings) {
    secret_terms <- c("password", "token", "secret", "api_key", "apikey")
    log_message("Run configuration:")
    for (name in names(settings)) {
      value <- settings[[name]]
      if (any(vapply(secret_terms, grepl, logical(1), x = tolower(name), fixed = TRUE))) {
        value <- "<redacted>"
      } else if (length(value) > 1) {
        value <- paste(value, collapse = ", ")
      }
      log_message("  ", name, ": ", value)
    }
    invisible(settings)
  }

  log_error <- function(error) {
    trace <- paste(capture.output(traceback(2)), collapse = "\n")
    detail <- conditionMessage(error)
    if (nzchar(trace)) {
      detail <- paste(detail, trace, sep = "\n")
    }
    line <- format_log_line("ERROR", detail)
    cat(line, "\n", file = stderr())
    cat(line, "\n", file = out_path, append = TRUE)
    cat(line, "\n", file = err_path, append = TRUE)
    invisible(line)
  }

  log_message("Logging to: ", out_path)
  log_message("Errors to: ", err_path)
  list(
    info = log_message,
    output = log_output,
    config = log_config,
    command = redact_command_args,
    error = log_error,
    out_path = out_path,
    err_path = err_path
  )
}
