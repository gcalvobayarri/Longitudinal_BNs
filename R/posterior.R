posterior_matrix <- function(fit) {
  if (!is.list(fit$samples) || length(fit$samples) == 0L) {
    stop("The fit does not contain a non-empty list of MCMC chains.", call. = FALSE)
  }
  columns <- lapply(fit$samples, colnames)
  if (!all(vapply(columns[-1L], identical, logical(1), columns[[1L]]))) {
    stop("MCMC chains have different parameter columns.", call. = FALSE)
  }
  do.call(rbind, lapply(fit$samples, as.matrix))
}

parameter_draws <- function(draws, parameter) {
  if (!parameter %in% colnames(draws)) {
    stop("Parameter not present in posterior draws: ", parameter, call. = FALSE)
  }
  draws[, parameter]
}

summarise_vector <- function(x) {
  c(
    mean = mean(x),
    sd = stats::sd(x),
    q2.5 = unname(stats::quantile(x, 0.025)),
    q97.5 = unname(stats::quantile(x, 0.975))
  )
}

summarise_parameters <- function(draws, parameters, labels = parameters) {
  if (length(parameters) != length(labels)) stop("parameters and labels must have equal length.", call. = FALSE)
  values <- t(vapply(parameters, function(p) summarise_vector(parameter_draws(draws, p)), numeric(4)))
  data.frame(parameter = labels, values, row.names = NULL, check.names = FALSE)
}

wide_shot_summary <- function(draws, parameter_roots, labels) {
  rows <- vector("list", length(parameter_roots))
  for (row in seq_along(parameter_roots)) {
    output <- data.frame(parameter = labels[[row]], stringsAsFactors = FALSE)
    for (k in seq_len(3L)) {
      stats <- summarise_vector(parameter_draws(draws, sprintf("%s[%d]", parameter_roots[[row]], k)))
      suffix <- shot_labels[[k]]
      output[[paste0("mean_", suffix)]] <- stats[["mean"]]
      output[[paste0("sd_", suffix)]] <- stats[["sd"]]
      output[[paste0("q2.5_", suffix)]] <- stats[["q2.5"]]
      output[[paste0("q97.5_", suffix)]] <- stats[["q97.5"]]
    }
    rows[[row]] <- output
  }
  do.call(rbind, rows)
}

extract_waic <- function(fit) {
  value <- fit$WAIC$WAIC
  if (length(value) != 1L || !is.finite(value)) stop("Could not extract a scalar WAIC.", call. = FALSE)
  as.numeric(value)
}
