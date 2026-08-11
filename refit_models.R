# OPTIONAL ENTRY POINT
#
# Run this script only to estimate the three NIMBLE models again. The default
# manuscript settings are computationally expensive. Most users should run
# reproduce.R instead, which uses the archived posterior samples.

if (!file.exists("longitudinal-bns-paper.Rproj")) {
  stop(
    "Open longitudinal-bns-paper.Rproj and run refit_models.R from the project root.",
    call. = FALSE
  )
}

source(file.path("analysis", "refit_models_internal.R"))
