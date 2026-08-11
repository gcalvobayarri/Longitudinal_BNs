.lbn_project_cache <- new.env(parent = emptyenv())

find_project_root <- function(start = getwd()) {
  if (!is.null(.lbn_project_cache$root)) {
    return(.lbn_project_cache$root)
  }

  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    project_files <- list.files(current, pattern = "[.]Rproj$", full.names = TRUE)
    if (length(project_files) > 0L) {
      .lbn_project_cache$root <- current
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the RStudio project root.", call. = FALSE)
    }
    current <- parent
  }
}

project_path <- function(...) {
  file.path(find_project_root(), ...)
}

ensure_output_directories <- function() {
  dirs <- c(
    project_path("results"),
    project_path("results", "tables"),
    project_path("results", "figures"),
    project_path("results", "mcmc")
  )
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

source_project <- function(file) {
  sys.source(project_path(file), envir = globalenv())
}
