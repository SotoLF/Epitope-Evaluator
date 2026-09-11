# deploy.R -- publish Epitope-Evaluator 2 to shinyapps.io
#
#   Rscript deploy.R --staging     # publish to Epitope-Evaluator-v2 (safe)
#   Rscript deploy.R --production  # overwrite the live Epitope-Evaluator
#   Rscript deploy.R --check       # bundle contents and size, publish nothing
#
# Credentials are read from the environment (see ensure_account() below), so
# put them in a gitignored .Renviron beside this script:
#
#   SHINYAPPS_ACCOUNT=fuxmanlab
#   SHINYAPPS_TOKEN=...
#   SHINYAPPS_SECRET=...
#
# Get them from https://www.shinyapps.io/admin/#/tokens (Show -> Show secret).

STAGING    <- "Epitope-Evaluator-v2"
PRODUCTION <- "Epitope-Evaluator"

# ---------------------------------------------------------------------------
# What goes in the bundle
# ---------------------------------------------------------------------------
# Everything the running app touches, and nothing else. Tests, the benchmark,
# the screenshot tool and the changelog are development files: shipping them
# only slows the upload and lengthens the cold start.

app_files <- function() {
  keep <- c(
    "app.R", "global.R",
    list.files("utils",   pattern = "[.]R$",   full.names = TRUE),
    list.files("modules", pattern = "[.]R$",   full.names = TRUE),
    list.files("ui",      pattern = "[.]R$",   full.names = TRUE),
    # www/ top level carries the stylesheet and the favicon set.
    list.files("www", pattern = "[.](css|svg|png|ico|webmanifest)$", full.names = TRUE),
    list.files("www/images", pattern = "[.]png$", full.names = TRUE),
    list.files("data", pattern = "[.](xls|txt|fasta)$", full.names = TRUE),
    # The paper's datasets keep the top-level path they were published under, so
    # links into the repository from the 2022 paper still resolve.
    list.files("Biological_Applications_Data", pattern = "[.](xls|txt|fasta)$",
               full.names = TRUE)
  )
  keep <- keep[file.exists(keep)]

  # Fail loudly rather than shipping a bundle whose examples or assets 404 at
  # runtime -- a missing data file breaks a tool, a missing favicon breaks the
  # browser tab, and neither shows up until someone opens the deployed app.
  need <- character(0)

  src <- readLines("modules/data_input.R", warn = FALSE)
  pat <- '"(data|Biological_Applications_Data)/[^"]+"'
  need <- c(need, gsub('"', "", unlist(regmatches(src, gregexpr(pat, src)))))

  # href=/src= assets in app.R resolve relative to www/.
  app <- readLines("app.R", warn = FALSE)
  assets <- unlist(regmatches(app, gregexpr('(href|src) = "[^"/][^"]*"', app)))
  assets <- gsub('^(href|src) = "|"$', "", assets)
  assets <- assets[grepl("[.](css|svg|png|ico|webmanifest)$", assets)]
  need <- c(need, file.path("www", assets))

  missing <- setdiff(unique(need), keep)
  if (length(missing)) {
    stop("The app references files that the bundle would not include:\n  ",
         paste(missing, collapse = "\n  "), call. = FALSE)
  }
  keep
}

report <- function(files) {
  sz <- file.size(files)
  cat(sprintf("\n%d files, %.1f MB total\n", length(files), sum(sz) / 1024^2))
  by_dir <- tapply(sz, dirname(files), sum)
  for (d in names(sort(by_dir, decreasing = TRUE))) {
    cat(sprintf("  %-14s %6.1f MB\n", d, by_dir[[d]] / 1024^2))
  }
  big <- files[sz > 3 * 1024^2]
  if (length(big)) {
    cat("\nLargest files:\n")
    for (f in big[order(-file.size(big))]) {
      cat(sprintf("  %-34s %5.1f MB\n", f, file.size(f) / 1024^2))
    }
  }
  invisible(files)
}

# ---------------------------------------------------------------------------
# Pre-flight
# ---------------------------------------------------------------------------

preflight <- function() {
  stopifnot(file.exists("app.R"), file.exists("global.R"))

  pkgs <- c("shiny", "bslib", "plotly", "DT", "data.table", "matrixStats", "stringi")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing package(s): ", paste(miss, collapse = ", "), call. = FALSE)

  cat("R ", as.character(getRversion()), "\n", sep = "")
  for (p in pkgs) cat(sprintf("  %-12s %s\n", p, packageVersion(p)))

  # A broken app must fail here, not after a five-minute upload.
  cat("\nSourcing the app... ")
  e <- new.env()
  app <- suppressPackageStartupMessages(source("app.R", local = e)$value)
  if (!inherits(app, "shiny.appobj")) stop("app.R did not return a Shiny app object.", call. = FALSE)
  cat("OK\n")

  if (!file.exists("www/images/tool_distribution.png")) {
    warning("Documentation screenshots are missing; run tools/screenshots.py first.",
            call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mode <- if ("--production" %in% args) {
  "production"
} else if ("--staging" %in% args) {
  "staging"
} else {
  "check"
}

files <- app_files()
preflight()
report(files)

if (mode == "check") {
  cat("\nCheck only. Re-run with --staging or --production to publish.\n")
  quit(save = "no")
}

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("install.packages(\"rsconnect\") first.", call. = FALSE)
}

# Account credentials come from the environment, normally set in a gitignored
# .Renviron beside this script:
#
#   SHINYAPPS_ACCOUNT=fuxmanlab
#   SHINYAPPS_TOKEN=...
#   SHINYAPPS_SECRET=...
#
# They are deliberately NOT written into this file, which is meant to be
# committed. rsconnect caches the account in ~/.config/R/rsconnect after the
# first successful call, so the variables are only needed once per machine.
ensure_account <- function() {
  if (nrow(rsconnect::accounts())) return(invisible(TRUE))

  acct   <- Sys.getenv("SHINYAPPS_ACCOUNT")
  token  <- Sys.getenv("SHINYAPPS_TOKEN")
  secret <- Sys.getenv("SHINYAPPS_SECRET")

  if (!nzchar(acct) || !nzchar(token) || !nzchar(secret)) {
    stop(
      "No shinyapps.io account configured.\n",
      "  Either create a .Renviron next to deploy.R containing\n",
      "      SHINYAPPS_ACCOUNT=fuxmanlab\n",
      "      SHINYAPPS_TOKEN=...\n",
      "      SHINYAPPS_SECRET=...\n",
      "  or run rsconnect::setAccountInfo() once in an interactive R session.\n",
      "  Tokens: https://www.shinyapps.io/admin/#/tokens",
      call. = FALSE)
  }
  cat(sprintf("Registering shinyapps.io account '%s' from the environment...\n", acct))
  rsconnect::setAccountInfo(name = acct, token = token, secret = secret)
  invisible(TRUE)
}
ensure_account()

app_name <- if (mode == "production") PRODUCTION else STAGING
cat(sprintf("\nDeploying to '%s'...\n", app_name))

rsconnect::deployApp(
  appDir      = ".",
  appFiles    = files,
  appName     = app_name,
  appTitle    = "Epitope-Evaluator",
  appPrimaryDoc = "app.R",
  forceUpdate = TRUE,
  launch.browser = FALSE
)
