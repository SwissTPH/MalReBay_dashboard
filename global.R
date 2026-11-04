if (file.exists(".Renviron")) readRenviron(".Renviron")

# MalReBay ships a precompiled Stan model and never compiles at runtime, but
# {instantiate}/{cmdstanr} still hard-error if no CmdStan installation can be
# resolved at all (see cmdstanr::set_cmdstan_path / stan_package_model()), and
# {cmdstanr}'s $sample() shells out to a real `bin/stanc` under that path for
# model introspection even with compile = FALSE. Deploy targets like
# shinyapps.io have no real CmdStan, so point at a tiny stub directory (just a
# makefile with a version line, plus a real portable `stanc` binary) that
# satisfies both needs without an actual toolchain.
if (!nzchar(Sys.getenv("CMDSTAN"))) {
  Sys.setenv(CMDSTAN = normalizePath(file.path(getwd(), "stub_cmdstan")))
}

# rsconnect's own bundling does not preserve the executable bit on plain
# files uploaded from the app directory (unlike MalReBay's own binaries,
# which travel via git -> GitHub -> R's package install and keep it fine) -
# confirmed via shinyapps.io logs showing stanc at mode 644 post-deploy.
# Restore it defensively on every startup rather than trust the bundle.
stanc_path <- file.path(Sys.getenv("CMDSTAN"), "bin", "stanc")
if (file.exists(stanc_path)) Sys.chmod(stanc_path, mode = "0755")

library(shiny)
library(bslib)
library(DT)
library(shinyjs)
library(MalReBay)
library(ggplot2)
library(mirai)
EXTDATA_PATH <- system.file("extdata", package = "MalReBay")

# Site drop-down entries: the pooled view, and (probability histogram only)
# the package's faceted one-panel-per-site view
ALL_SITES   <- "All sites"
EACH_SITE   <- "One panel per site"

# Pick the figure for a site drop-down value from a MalReBay plot result of
# the form list(all_sites = <plot>, by_site = list(<site> = <plot>, ...)).
# Falls back to the pooled figure for a value that isn't a site of the
# current data (e.g. briefly, right after another dataset is loaded).
select_site_plot <- function(plots, site) {
  if (!is.null(site) && site %in% names(plots$by_site)) plots$by_site[[site]] else plots$all_sites
}

# Several MalReBay plot functions also print() their plots when given no
# output_folder. Evaluate `expr` with those prints sent to a null device, so
# they don't land on (and overwrite) the renderPlot() device, then restore it.
without_printing <- function(expr) {
  old_dev <- grDevices::dev.cur()
  grDevices::pdf(NULL)
  on.exit({
    grDevices::dev.off()
    if (old_dev > 1) grDevices::dev.set(old_dev)
  })
  expr
}

# Default MCMC configuration for MalReBay's Stan-based sampler
# (matches the parameter names read by MalReBay::classify_infections() /
# run_stan_sites(): n_chains, iter, burn_in_frac, adapt_delta, random_seed)
# Most MCMC chains a dashboard run may use: chains run in parallel, and the
# hosted app has limited CPU and memory
MAX_MCMC_CHAINS <- 3

# Largest analysis the dashboard accepts, so a single run stays within about
# 15-20 minutes on a shared instance (measured locally: ~1.5-2 s per patient
# at 2000 iterations and 3 chains, plus ~20-25 s per site). Larger studies
# should use the MalReBay R package directly.
MAX_RECURRENT_PATIENTS <- 200     # Day 0 + recurrence pairs, all sites together
MAX_ADDITIONAL_SAMPLES <- 1000    # extra Day 0 samples (allele frequencies only)
MAX_ITERATIONS         <- 5000    # per chain
MAX_PATIENT_ITERATIONS <- 400000  # patients x iterations, e.g. 200 x 2000

# Why a run of this size can't be started, or NULL if it can
run_size_problem <- function(n_patients, n_additional, iter) {
  if (n_patients > MAX_RECURRENT_PATIENTS)
    return(sprintf("The dataset has %d patients with a recurrence; the dashboard accepts at most %d.",
                   n_patients, MAX_RECURRENT_PATIENTS))
  if (n_additional > MAX_ADDITIONAL_SAMPLES)
    return(sprintf("The dataset has %d additional Day 0 samples; the dashboard accepts at most %d.",
                   n_additional, MAX_ADDITIONAL_SAMPLES))
  if (iter > MAX_ITERATIONS)
    return(sprintf("%d iterations per chain requested; the dashboard accepts at most %d.",
                   iter, MAX_ITERATIONS))
  if (n_patients * iter > MAX_PATIENT_ITERATIONS)
    return(sprintf("%d patients x %d iterations is too large a run for the dashboard: lower the iterations to %d or fewer.",
                   n_patients, iter, floor(MAX_PATIENT_ITERATIONS / n_patients / 500) * 500))
  NULL
}

# MalReBay runs in background R processes (mirai daemons), so a run doesn't
# freeze the app for everyone else connected to it. N_RUN_WORKERS runs can
# execute at the same time on an instance; further runs wait for a free
# worker. Each run uses up to MAX_MCMC_CHAINS cores and ~300-400 MB.
N_RUN_WORKERS <- 2
mirai::daemons(N_RUN_WORKERS)
# The workers need the app's packages and CmdStan location
mirai::everywhere({
  .libPaths(lib_paths)
  Sys.setenv(CMDSTAN = cmdstan)
}, lib_paths = .libPaths(), cmdstan = Sys.getenv("CMDSTAN"))
shiny::onStop(function() mirai::daemons(0))

# Start MalReBay(run_args) on a background worker; resolves to the results
# and the run time in seconds (only the run itself, not the wait for a
# worker). The worker creates `started_flag` (a file path) when it picks the
# run up, so the app can tell a queued run from a running one.
run_malrebay_background <- function(run_args, started_flag) {
  mirai::mirai({
    file.create(started_flag)
    suppressPackageStartupMessages(library(MalReBay))
    start   <- Sys.time()
    results <- do.call(MalReBay::MalReBay, run_args)
    list(results = results, seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
  }, run_args = run_args, started_flag = started_flag)
}


DEFAULT_MCMC <- list(
  n_chains     = 3,
  iter         = 2000,
  burn_in_frac = 0.5,
  adapt_delta  = 0.85,
  random_seed  = 42
)

# Input templates + example marker metadata, copied from the MalReBay
# package (inst/templates/, inst/extdata/makers_details.xlsx) into
# www/templates/ so they're served as plain static files: no dependence on
# which MalReBay install is on the library path, and no Shiny session needed
# (downloadHandler links break inside hidden dropdowns and some iframes).
# Re-copy them when the package's templates change.
static_download_link <- function(file, label, class = NULL) {
  tags$a(href = paste0("templates/", file), download = file, class = class, label)
}

# Placeholder shown inside a chart/table card instead of leaving it blank
empty_panel_msg <- function(icon_name, text) {
  div(class = "panel-empty-msg", icon(icon_name), " ", text)
}

# Length-polymorphic examples, both cut from the package's 3-site TES dataset
# so the bundled demos' MCMC runtime stays reasonable:
# - 1 site:  a single site ("Lunda Sul", 13 patients)
# - 2 sites: two sites ("Benguela" + "Lunda Sul"), dropping "Zaire", the
#          largest site at 76/140 late-failure rows
LP_EXAMPLE_SOURCE <- "Angola_2021_TES_7NMS.xlsx"
LP_EXAMPLE_SHORT_SITES <- "Lunda Sul"
LP_EXAMPLE_LONG_SITES  <- c("Benguela", "Lunda Sul")

build_lp_example <- function(sites, file_name) {
  raw_path <- file.path(EXTDATA_PATH, LP_EXAMPLE_SOURCE)
  sheets   <- readxl::excel_sheets(raw_path)

  sheet_data <- lapply(sheets, function(s) {
    df <- as.data.frame(readxl::read_excel(raw_path, sheet = s))
    df[df$Site %in% sites, ]
  })
  names(sheet_data) <- sheets

  out_path <- file.path(tempdir(), file_name)
  writexl::write_xlsx(sheet_data, out_path)
  out_path
}
