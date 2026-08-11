load_rdata <- function(path) {
  if (!file.exists(path)) {
    stop("Required file does not exist: ", path, call. = FALSE)
  }
  env <- new.env(parent = baseenv())
  object_names <- load(path, envir = env)
  stats::setNames(lapply(object_names, function(name) env[[name]]), object_names)
}

load_raw_data <- function() {
  required <- c("C1", "C2", "C3", "T1", "T2", "T3", "MP", "foulsR", "home", "C", "PG", "SG", "PF", "SF")
  values <- load_rdata(project_path("data", "raw", "matrix_data_PHI_13players.RData"))
  missing <- setdiff(required, names(values))
  if (length(missing) > 0L) {
    stop("Raw data are missing objects: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  values[required]
}

load_fit <- function(model = c("static", "dynamic", "hidden_markov"), refit = FALSE) {
  model <- match.arg(model)
  if (refit) {
    path <- project_path("results", "mcmc", paste0("refit_", model, ".RData"))
    objects <- load_rdata(path)
    if (length(objects) != 1L) {
      stop("Expected exactly one fit object in ", path, call. = FALSE)
    }
    fit <- objects[[1L]]
  } else {
    filenames <- c(
      static = "static_posterior.rds",
      dynamic = "dynamic_posterior.rds",
      hidden_markov = "hidden_markov_posterior.rds"
    )
    path <- project_path("results", "mcmc", filenames[[model]])
    if (!file.exists(path)) {
      stop("Required posterior file does not exist: ", path, call. = FALSE)
    }
    fit <- readRDS(path)
  }

  if (!is.list(fit$samples) || length(fit$samples) == 0L || is.null(fit$WAIC$WAIC)) {
    stop("Invalid posterior fit structure in ", path, call. = FALSE)
  }
  fit
}

validate_raw_data <- function(data) {
  expected_dim <- c(13L, 82L)
  matrix_names <- names(data)
  bad_dims <- matrix_names[!vapply(data, function(x) identical(dim(x), expected_dim), logical(1))]
  if (length(bad_dims) > 0L) {
    stop("Unexpected dimensions in: ", paste(bad_dims, collapse = ", "), call. = FALSE)
  }

  for (k in seq_len(3L)) {
    made <- data[[paste0("C", k)]]
    attempted <- data[[paste0("T", k)]]
    if (any(made < 0 | attempted < 0 | made > attempted, na.rm = TRUE)) {
      stop("Invalid made/attempted counts for shot type ", k, call. = FALSE)
    }
  }

  if (any(!data$home %in% c(0, 1))) stop("The home indicator is not binary.", call. = FALSE)
  if (any(data$MP < 0, na.rm = TRUE)) stop("Negative minutes detected.", call. = FALSE)

  absent <- data$MP == 0
  count_names <- c("C1", "C2", "C3", "T1", "T2", "T3", "foulsR")
  inconsistent <- count_names[vapply(count_names, function(name) any(data[[name]][absent] != 0), logical(1))]
  if (length(inconsistent) > 0L) {
    stop("Non-zero outcomes when minutes are zero in: ", paste(inconsistent, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}
