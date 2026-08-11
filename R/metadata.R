player_metadata <- data.frame(
  player_id = seq_len(13L),
  player = c(
    "Allen Iverson", "Andre Iguodala", "Chris Webber", "John Salmons",
    "Kevin Ollie", "Kyle Korver", "Lee Nailon", "Lou Williams",
    "Matt Barnes", "Michael Bradley", "Samuel Dalembert",
    "Shavlik Randolph", "Steven Hunter"
  ),
  position = c("PG", "PG", "PF", "SG", "SG", "SF", "PF", "PG", "PF", "PF", "C", "PF", "C"),
  stringsAsFactors = FALSE
)

shot_labels <- c("1PT", "2PT", "3PT")
paper_seed <- 260809824L
paper_mcmc_settings <- list(niter = 1000000L, nburnin = 500000L, thin = 500L, nchains = 3L)
