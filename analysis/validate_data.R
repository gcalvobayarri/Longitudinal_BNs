# Internal step used by reproduce.R: validate inputs and create Table 1.
if (!exists("project_path")) source(file.path("R", "project.R"))
source_project(file.path("R", "metadata.R"))
source_project(file.path("R", "data_io.R"))

data <- load_raw_data()
validate_raw_data(data)

table_1 <- data.frame(
  player_metadata,
  minutes_played = rowSums(data$MP),
  fouls_received = rowSums(data$foulsR),
  ft_attempted = rowSums(data$T1), ft_made = rowSums(data$C1),
  two_pt_attempted = rowSums(data$T2), two_pt_made = rowSums(data$C2),
  three_pt_attempted = rowSums(data$T3), three_pt_made = rowSums(data$C3),
  check.names = FALSE
)

ensure_output_directories()
utils::write.csv(table_1, project_path("results", "tables", "table_1_player_summary.csv"), row.names = FALSE)
message("Raw data validation passed: 13 players x 82 games.")
