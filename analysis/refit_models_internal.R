# Internal implementation used by the top-level refit_models.R entry point.
if (!exists("project_path")) source(file.path("R", "project.R"))
source_project(file.path("R", "metadata.R"))
source_project(file.path("R", "data_io.R"))
source_project(file.path("R", "models.R"))

if (!requireNamespace("nimble", quietly = TRUE)) stop("Package 'nimble' is required for model fitting.", call. = FALSE)
suppressPackageStartupMessages(library(nimble))

read_integer_setting <- function(name, default) {
  value <- Sys.getenv(name, unset = as.character(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) stop(name, " must be a positive integer.", call. = FALSE)
  parsed
}

requested <- strsplit(Sys.getenv("LBN_MODELS", unset = "static,dynamic,hidden_markov"), ",", fixed = TRUE)[[1L]]
requested <- trimws(requested)
allowed <- c("static", "dynamic", "hidden_markov")
if (length(setdiff(requested, allowed)) > 0L) stop("Unknown model(s): ", paste(setdiff(requested, allowed), collapse = ", "), call. = FALSE)

settings <- list(
  niter = read_integer_setting("LBN_NITER", paper_mcmc_settings$niter),
  nburnin = read_integer_setting("LBN_NBURNIN", paper_mcmc_settings$nburnin),
  thin = read_integer_setting("LBN_THIN", paper_mcmc_settings$thin),
  nchains = read_integer_setting("LBN_NCHAINS", paper_mcmc_settings$nchains)
)
if (settings$nburnin >= settings$niter) stop("LBN_NBURNIN must be smaller than LBN_NITER.", call. = FALSE)

inputs <- make_model_inputs()
ensure_output_directories()
for (model in requested) {
  message("Fitting ", model, " LBN with ", settings$nchains, " chain(s)...")
  fit <- nimble::nimbleMCMC(
    code = model_code(model), constants = inputs$constants, data = inputs$data,
    inits = model_inits(model), monitors = model_monitors(model),
    niter = settings$niter, nburnin = settings$nburnin, thin = settings$thin,
    nchains = settings$nchains, setSeed = paper_seed + seq_len(settings$nchains),
    summary = TRUE, WAIC = TRUE
  )
  object_name <- paste0("refit_", model)
  assign(object_name, fit)
  save(list = object_name, file = project_path("results", "mcmc", paste0(object_name, ".RData")))
  message(model, " WAIC: ", fit$WAIC$WAIC)
}
