# MAIN ENTRY POINT
#
# Run this script to reproduce all tables and figures from the archived
# posterior samples. A full model refit is not required.

if (!file.exists("longitudinal-bns-paper.Rproj")) {
  stop(
    "Open longitudinal-bns-paper.Rproj and run reproduce.R from the project root.",
    call. = FALSE
  )
}

source(file.path("R", "project.R"))

steps <- c(
  "analysis/validate_data.R",
  "analysis/make_tables.R",
  "analysis/make_dag_figures.R",
  "analysis/make_posterior_figures.R",
  "analysis/validate_results.R"
)

for (step in steps) {
  message("\n--- Running ", step, " ---")
  sys.source(project_path(step), envir = globalenv())
}

writeLines(capture.output(sessionInfo()), project_path("results", "session_info.txt"))
message("\nReproduction completed successfully. See results/ for all outputs.")
