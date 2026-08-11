draw_node <- function(x, y, label, radius = 0.055) {
  symbols(x, y, circles = radius, inches = FALSE, add = TRUE, bg = "white", fg = "black")
  text(x, y, labels = label, cex = 0.8)
}

draw_edge <- function(from, to, radius = 0.055, curvature = 0) {
  delta <- to - from
  distance <- sqrt(sum(delta^2))
  unit <- delta / distance
  start <- from + unit * radius
  end <- to - unit * radius
  if (curvature == 0) {
    arrows(start[1], start[2], end[1], end[2], length = 0.08, angle = 20)
  } else {
    control <- (start + end) / 2 + c(0, curvature)
    t <- seq(0, 1, length.out = 80)
    x <- (1 - t)^2 * start[1] + 2 * (1 - t) * t * control[1] + t^2 * end[1]
    y <- (1 - t)^2 * start[2] + 2 * (1 - t) * t * control[2] + t^2 * end[2]
    lines(x, y)
    arrows(x[length(x) - 2L], y[length(y) - 2L], end[1], end[2], length = 0.08, angle = 20)
  }
}

time_slice_nodes <- function(offset = 0, suffix = "j") {
  list(
    M = c(0.50 + offset, 0.88), F = c(0.24 + offset, 0.66),
    T1 = c(0.18 + offset, 0.42), T2 = c(0.50 + offset, 0.42), T3 = c(0.82 + offset, 0.42),
    C1 = c(0.18 + offset, 0.16), C2 = c(0.50 + offset, 0.16), C3 = c(0.82 + offset, 0.16),
    suffix = suffix
  )
}

draw_slice_edges <- function(nodes) {
  coordinates <- nodes[setdiff(names(nodes), "suffix")]
  edges <- list(c("M", "F"), c("F", "T1"), c("M", "T2"), c("M", "T3"), c("T1", "C1"), c("T2", "C2"), c("T3", "C3"))
  for (edge in edges) draw_edge(coordinates[[edge[1]]], coordinates[[edge[2]]])
  invisible(coordinates)
}

draw_slice_nodes <- function(nodes) {
  coordinates <- nodes[setdiff(names(nodes), "suffix")]
  labels <- c(M = paste0("M[", nodes$suffix, "]"), F = paste0("F[", nodes$suffix, "]"), T1 = paste0("T1[", nodes$suffix, "]"), T2 = paste0("T2[", nodes$suffix, "]"), T3 = paste0("T3[", nodes$suffix, "]"), C1 = paste0("C1[", nodes$suffix, "]"), C2 = paste0("C2[", nodes$suffix, "]"), C3 = paste0("C3[", nodes$suffix, "]"))
  for (name in names(coordinates)) draw_node(coordinates[[name]][1], coordinates[[name]][2], labels[[name]])
  invisible(coordinates)
}

draw_slice <- function(nodes) {
  draw_slice_edges(nodes)
  draw_slice_nodes(nodes)
}

plot_static_dag <- function() {
  par(mar = rep(0.3, 4)); plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
  draw_slice(time_slice_nodes())
}

plot_dynamic_dag <- function(hidden = FALSE) {
  par(mar = rep(0.3, 4)); plot.new(); plot.window(xlim = c(0, 2.5), ylim = c(0, 1.05), asp = 1)
  previous_nodes <- time_slice_nodes(0.25, "j-1")
  current_nodes <- time_slice_nodes(1.45, "j")
  previous <- draw_slice_edges(previous_nodes)
  current <- draw_slice_edges(current_nodes)
  if (!hidden) draw_edge(previous$M, current$M, curvature = 0.25)
  if (hidden) {
    z_previous <- c(0.12, 0.60); z_current <- c(1.32, 0.60)
    draw_edge(z_previous, z_current, curvature = 0.40)
    for (target in c("C1", "C2", "C3")) { draw_edge(z_previous, previous[[target]]); draw_edge(z_current, current[[target]]) }
  }
  draw_slice_nodes(previous_nodes)
  draw_slice_nodes(current_nodes)
  if (hidden) { draw_node(z_previous[1], z_previous[2], "Z[j-1]"); draw_node(z_current[1], z_current[2], "Z[j]") }
}

save_base_figure <- function(plot_function, stem, width = 7, height = 4.2) {
  grDevices::pdf(project_path("results", "figures", paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  plot_function(); grDevices::dev.off()
  grDevices::png(project_path("results", "figures", paste0(stem, ".png")), width = width, height = height, units = "in", res = 180)
  plot_function(); grDevices::dev.off()
}
