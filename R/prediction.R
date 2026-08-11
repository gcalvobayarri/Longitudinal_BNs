inv_logit <- function(x) stats::plogis(x)

position_value <- function(player_id, data, prefix, k, draws) {
  value <- rep(0, nrow(draws))
  for (position in c("SG", "PG", "SF", "PF")) {
    indicator <- data[[position]][player_id, 1]
    if (indicator == 1) value <- value + draws[, sprintf("%s%s[%d]", prefix, position, k)]
  }
  value
}

success_probability_draws <- function(draws, player_id, home = 1, data = load_raw_data()) {
  result <- matrix(NA_real_, nrow = nrow(draws), ncol = 3L, dimnames = list(NULL, shot_labels))
  for (k in seq_len(3L)) {
    eta <- draws[, sprintf("beta0[%d]", k)] +
      position_value(player_id, data, "beta", k, draws) +
      draws[, sprintf("b0[%d, %d]", player_id, k)] +
      draws[, sprintf("betaH[%d]", k)] * home
    result[, k] <- inv_logit(eta)
  }
  result
}

simulate_new_match <- function(draws, player_ids = seq_len(13L), home = 1, previous_minutes = NULL, seed = paper_seed, data = load_raw_data()) {
  set.seed(seed)
  output <- vector("list", length(player_ids))

  for (index in seq_along(player_ids)) {
    i <- player_ids[[index]]
    n <- nrow(draws)
    plays <- stats::rbinom(n, 1, 1 - draws[, sprintf("p[%d]", i)])
    minutes_eta <- draws[, "delta0"] + draws[, sprintf("d0[%d]", i)]
    if (!is.null(previous_minutes)) minutes_eta <- minutes_eta + draws[, "deltaw"] * log1p(previous_minutes[[index]])
    minutes <- stats::rpois(n, exp(minutes_eta) * plays + 1.0E-7)
    fouls <- stats::rpois(n, exp(draws[, "gamma0"] + draws[, sprintf("c0[%d]", i)] + draws[, "gamma"] * minutes) * plays + 1.0E-7)

    attempts <- made <- matrix(0L, nrow = n, ncol = 3L)
    probabilities <- success_probability_draws(draws, i, home = home, data = data)
    for (k in seq_len(3L)) {
      exposure <- if (k == 1L) fouls else minutes
      eta <- draws[, sprintf("alpha0[%d]", k)] +
        position_value(i, data, "alpha", k, draws) +
        draws[, sprintf("a0[%d, %d]", i, k)] +
        draws[, sprintf("alpha[%d]", k)] * exposure
      attempts[, k] <- stats::rpois(n, exp(eta) * plays + 1.0E-7)
      made[, k] <- stats::rbinom(n, attempts[, k], probabilities[, k])
    }

    output[[index]] <- data.frame(
      draw = seq_len(n), player_id = i, player = player_metadata$player[[i]], home = home,
      plays = plays, minutes = minutes, fouls_drawn = fouls,
      attempts_1pt = attempts[, 1], made_1pt = made[, 1],
      attempts_2pt = attempts[, 2], made_2pt = made[, 2],
      attempts_3pt = attempts[, 3], made_3pt = made[, 3],
      points = made[, 1] + 2L * made[, 2] + 3L * made[, 3],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}
