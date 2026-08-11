# Longitudinal Bayesian Networks: reproducibility repository

This repository contains the data, archived posterior samples, and R code used
to reproduce the empirical results in:

> Calvo et al. *Longitudinal Bayesian networks for assessing team performance
> in the National Basketball Association* (arXiv:2608.09824v1).

## Start here

For the standard reproduction, only run **`reproduce.R`**.

1. Open `longitudinal-bns-paper.Rproj` in RStudio.
2. Reproduce every table and figure:

   ```r
   source("reproduce.R")
   ```

The generated files are written to `results/tables/` and `results/figures/`.
The archived posterior samples are stored in package-independent RDS files, so
this workflow does not require NIMBLE, ggplot2, or any other contributed R
package. It does not refit the models and normally finishes quickly.

From a terminal, run:

```sh
Rscript reproduce.R
```

## Which script should I use?

| Goal | Script | Typical use |
|---|---|---|
| Reproduce the paper results | `reproduce.R` | **Use this in almost all cases.** |
| Re-estimate all NIMBLE models | `refit_models.R` | Optional; computationally expensive. |

Files inside `analysis/` and `R/` are internal building blocks called by these
two entry points. They are provided for transparency and development, but users
do not need to run them individually.

## What `reproduce.R` does

The script executes the complete workflow in the correct order:

1. validates the 13-player by 82-game input data and creates Table 1;
2. creates Tables 2-7 from the archived posterior samples;
3. creates the three network diagrams (Figures 1-3);
4. runs the posterior-predictive simulation and creates Figures 4-6;
5. checks the numerical outputs against the values reported in the paper;
6. records the R session used for the reproduction.

Posterior-predictive simulation uses the fixed seed `260809824`.

## Optional full model refit

Run the following only if you want to estimate the models again:

```r
renv::restore()
source("refit_models.R")
```

The default settings match the manuscript:

- 3 chains per model;
- 1,000,000 iterations per chain;
- 500,000 burn-in iterations;
- thinning interval of 500.

This can take a long time and requires a working C++ toolchain because NIMBLE
compiles each model. On Windows, install the Rtools version matching your R
version before starting a refit.

For a smaller diagnostic run, set the environment variables before sourcing
the script:

```r
Sys.setenv(
  LBN_MODELS = "dynamic",
  LBN_NITER = 2000,
  LBN_NBURNIN = 1000,
  LBN_THIN = 10,
  LBN_NCHAINS = 1
)
source("refit_models.R")
```

`LBN_MODELS` accepts `static`, `dynamic`, `hidden_markov`, or a comma-separated
combination. New fits are saved as `results/mcmc/refit_<model>.RData`; the
archived paper fits are never overwritten.

## Repository structure

```text
.
|-- reproduce.R                 # Main script: reproduce the paper
|-- refit_models.R              # Optional script: re-estimate the models
|-- longitudinal-bns-paper.Rproj
|-- data/raw/                   # Player-game input data
|-- results/mcmc/               # Package-independent posterior samples
|-- results/tables/             # Generated Tables 1-7
|-- results/figures/            # Generated Figures 1-6
|-- analysis/                   # Internal workflow steps
|-- R/                          # Internal functions and model definitions
|-- renv.lock                   # Packages required only for model refitting
`-- README.md
```

## Output files

| Paper item | Generated file |
|---|---|
| Table 1 | `results/tables/table_1_player_summary.csv` |
| Table 2 | `results/tables/table_2_waic.csv` |
| Table 3 | `results/tables/table_3_shots_made.csv` |
| Table 4 | `results/tables/table_4_probability_positive.csv` |
| Table 5 | `results/tables/table_5_shots_attempted.csv` |
| Table 6 | `results/tables/table_6_minutes.csv` |
| Table 7 | `results/tables/table_7_participation.csv` |
| Figures 1-6 | `results/figures/figure_*.pdf` and `.png` |
| Numerical check | `results/validation_report.txt` |
| R session | `results/session_info.txt` |

## Software

The reproduction workflow was verified with R 4.5.1 and R 4.6.0 on Windows and
uses base R only. Full model refitting was verified with NIMBLE 1.3.0. Exact
refitting package versions are recorded in `renv.lock`.
