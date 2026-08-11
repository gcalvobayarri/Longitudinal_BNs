# Internal step used by reproduce.R: create Figures 4-6 using base R only.
if (!exists("project_path")) source(file.path("R", "project.R"))
source_project(file.path("R", "metadata.R"))
source_project(file.path("R", "data_io.R"))
source_project(file.path("R", "posterior.R"))
source_project(file.path("R", "prediction.R"))

ensure_output_directories()
draws <- posterior_matrix(load_fit("dynamic"))
selected_players <- c(1L, 6L)

save_figure_pair <- function(draw, stem, width, height = 4.2) {
  pdf_path <- project_path("results", "figures", paste0(stem, ".pdf"))
  png_path <- project_path("results", "figures", paste0(stem, ".png"))

  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  draw()
  grDevices::dev.off()

  grDevices::png(
    png_path, width = width, height = height, units = "in", res = 180,
    bg = "white"
  )
  draw()
  grDevices::dev.off()
}

success_probabilities <- lapply(selected_players, function(player_id) {
  success_probability_draws(draws, player_id, home = 1)
})
names(success_probabilities) <- player_metadata$player[selected_players]

shot_colours <- c("1PT" = "#E69F00", "2PT" = "#0072B2", "3PT" = "#009E73")

draw_figure_4 <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1.0), las = 1)

  for (player_index in seq_along(success_probabilities)) {
    player_draws <- success_probabilities[[player_index]]
    density_curves <- lapply(seq_len(3L), function(k) {
      stats::density(player_draws[, k], from = 0, to = 1, n = 512)
    })

    graphics::plot(
      NA, xlim = c(0, 1), ylim = c(0, 30), xlab = "Probability",
      ylab = if (player_index == 1L) "Density" else "", axes = FALSE,
      main = names(success_probabilities)[[player_index]]
    )
    graphics::grid(col = "#E5E5E5", lty = 1)
    graphics::axis(1, at = seq(0, 1, 0.25))
    graphics::axis(2, at = seq(0, 30, 10))
    graphics::box(col = "#BFBFBF")

    for (k in seq_len(3L)) {
      curve <- density_curves[[k]]
      colour <- shot_colours[[k]]
      graphics::polygon(
        c(curve$x, rev(curve$x)), c(curve$y, rep(0, length(curve$y))),
        col = grDevices::adjustcolor(colour, alpha.f = 0.65), border = NA
      )
    }

    if (player_index == 2L) {
      graphics::legend(
        "topright", legend = names(shot_colours), fill = shot_colours,
        title = "Shots scored", border = NA, bty = "n"
      )
    }
  }
}

prediction <- simulate_new_match(
  draws, player_ids = selected_players, home = 1, seed = paper_seed
)
saveRDS(prediction, project_path("results", "posterior_predictive.rds"))

relative_frequency <- function(data, variable) {
  counts <- as.data.frame(
    table(player = data$player, value = data[[variable]]),
    stringsAsFactors = FALSE
  )
  counts <- counts[counts$Freq > 0, ]
  counts$value <- as.numeric(counts$value)
  counts$proportion <- ave(
    counts$Freq, counts$player, FUN = function(x) x / sum(x)
  )
  counts
}

player_colours <- c("Allen Iverson" = "#6A3D9A", "Kyle Korver" = "#FFDF4D")

draw_frequency_plot <- function(frequencies, x_label, x_limits) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4.2, 4.2, 1.0, 1.0), las = 1)

  y_max <- max(frequencies$proportion) * 1.08
  graphics::plot(
    NA, xlim = x_limits, ylim = c(0, y_max), xlab = x_label, ylab = "",
    axes = FALSE, xaxs = "i", yaxs = "i"
  )
  graphics::grid(col = "#E5E5E5", lty = 1)
  graphics::axis(1)
  graphics::axis(2)
  graphics::box(col = "#BFBFBF")

  # Draw the yellow distribution first so both overlapping distributions remain visible.
  draw_order <- c("Kyle Korver", "Allen Iverson")
  for (player in draw_order) {
    values <- frequencies[frequencies$player == player, ]
    graphics::rect(
      values$value - 0.5, 0, values$value + 0.5, values$proportion,
      col = grDevices::adjustcolor(player_colours[[player]], alpha.f = 0.70),
      border = NA
    )
  }

  graphics::legend(
    "topright", legend = names(player_colours), fill = player_colours,
    title = "Player", border = NA, bty = "n"
  )
}

points_over_30 <- prediction[prediction$minutes > 30 & prediction$minutes <= 48, ]
points_frequency <- relative_frequency(points_over_30, "points")

minutes_at_most_10 <- prediction[prediction$points <= 10, ]
minutes_frequency <- relative_frequency(minutes_at_most_10, "minutes")

save_figure_pair(draw_figure_4, "figure_4_success_probabilities", width = 8)
save_figure_pair(
  function() draw_frequency_plot(points_frequency, "Points", c(0, 60)),
  "figure_5_points_given_minutes", width = 7
)
save_figure_pair(
  function() draw_frequency_plot(minutes_frequency, "Minutes", c(0, 40)),
  "figure_6_minutes_given_points", width = 7
)

message(
  "Posterior and posterior-predictive figures 4-6 created with base R and seed ",
  paper_seed, "."
)
