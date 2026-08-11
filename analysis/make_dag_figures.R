# Internal step used by reproduce.R: create Figures 1-3.
if (!exists("project_path")) source(file.path("R", "project.R"))
source_project(file.path("R", "dag.R"))
ensure_output_directories()

save_base_figure(plot_static_dag, "figure_1_static_dag", width = 4.5, height = 4.5)
save_base_figure(function() plot_dynamic_dag(FALSE), "figure_2_dynamic_dag", width = 8, height = 4.5)
save_base_figure(function() plot_dynamic_dag(TRUE), "figure_3_hidden_markov_dag", width = 8, height = 4.5)
message("Conceptual DAG figures 1-3 created.")
