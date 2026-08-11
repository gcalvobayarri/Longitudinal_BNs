# Internal step used by reproduce.R: create Tables 2-7.
if (!exists("project_path")) source(file.path("R", "project.R"))
source_project(file.path("R", "metadata.R"))
source_project(file.path("R", "data_io.R"))
source_project(file.path("R", "posterior.R"))

ensure_output_directories()
fits <- lapply(c("static", "dynamic", "hidden_markov"), load_fit)
names(fits) <- c("Static LBN", "Dynamic LBN", "Hidden Markov LBN")
table_2 <- data.frame(model = names(fits), WAIC = vapply(fits, extract_waic, numeric(1)), row.names = NULL)
utils::write.csv(table_2, project_path("results", "tables", "table_2_waic.csv"), row.names = FALSE)

draws <- posterior_matrix(fits[["Dynamic LBN"]])
made_roots <- c("beta0", "betaH", "betaPF", "betaPG", "betaSF", "betaSG", "sigmab")
made_labels <- c("Intercept (center)", "Home", "Power forward", "Point guard", "Small forward", "Shooting guard", "Player SD")
table_3 <- wide_shot_summary(draws, made_roots, made_labels)
utils::write.csv(table_3, project_path("results", "tables", "table_3_shots_made.csv"), row.names = FALSE)

positive_roots <- made_roots[1:6]
table_4 <- data.frame(parameter = made_labels[1:6], stringsAsFactors = FALSE)
for (k in seq_len(3L)) {
  table_4[[paste0("probability_positive_", shot_labels[[k]])]] <- vapply(positive_roots, function(root) mean(draws[, sprintf("%s[%d]", root, k)] > 0), numeric(1))
}
utils::write.csv(table_4, project_path("results", "tables", "table_4_probability_positive.csv"), row.names = FALSE)

attempt_roots <- c("alpha0", "alpha", "alphaPF", "alphaSF", "alphaPG", "alphaSG", "sigmaa")
attempt_labels <- c("Intercept (center)", "Fouls/minutes effect", "Power forward", "Small forward", "Point guard", "Shooting guard", "Player SD")
table_5 <- wide_shot_summary(draws, attempt_roots, attempt_labels)
utils::write.csv(table_5, project_path("results", "tables", "table_5_shots_attempted.csv"), row.names = FALSE)

table_6 <- summarise_parameters(draws, c("delta0", "deltaw", "sigmad"), c("Overall mean", "Autoregressive coefficient", "Player SD"))
utils::write.csv(table_6, project_path("results", "tables", "table_6_minutes.csv"), row.names = FALSE)

participation <- lapply(seq_len(13L), function(i) summarise_vector(1 - draws[, sprintf("p[%d]", i)]))
table_7 <- cbind(player_metadata, as.data.frame(do.call(rbind, participation), row.names = NULL))
utils::write.csv(table_7, project_path("results", "tables", "table_7_participation.csv"), row.names = FALSE)
message("Tables 2-7 created from 3,000 retained dynamic-model posterior draws.")
