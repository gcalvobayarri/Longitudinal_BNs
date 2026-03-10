# Code Repository for: Longitudinal Bayesian Networks

This repository contains the code and resources used to implement the analyses presented in the associated research paper on **Longitudinal Bayesian Networks**. The repository provides the scripts required to reproduce the posterior analyses and predictive simulations derived from the fitted Bayesian model.

## Repository Structure

The repository is organised to separate data, model outputs, and analysis scripts:

```
paper/
├── data/ # Data used in the analysis
├── models/ # Model specifications and definitions
├── results/ # Stored outputs from model estimation (e.g., MCMC samples)
└── scripts/ # Scripts used for posterior analysis and predictive simulations
```




## Description of the Code

The scripts in the `paper/scripts/` directory perform the post–estimation analyses based on posterior samples generated from the Bayesian model.

### Posterior Analysis

**`posterior_results.R`**

This script loads the stored posterior samples and computes summary statistics of the posterior distributions. In particular, it derives quantities of interest such as posterior means, credible intervals, and posterior probabilities associated with model parameters.

### Posterior Predictive Simulation

**`predictive_post_dynamic.R`**

This script performs posterior predictive simulations based on the fitted model. The predictive procedure uses the stored posterior samples to generate simulated outcomes for the variables of interest and assess the predictive behaviour of the model.

## Model Output

The scripts rely on posterior samples stored in:

```
paper/results/nimbleMCMCDynamic_samples.RData
```

These samples correspond to the Markov Chain Monte Carlo (MCMC) output obtained during the Bayesian model estimation.

## Reproducibility

To reproduce the analyses reported in the paper:

1. Clone this repository.
2. Install the required R packages.
3. Run the scripts located in `paper/scripts/`.

The scripts will load the stored MCMC samples and reproduce the posterior summaries and predictive simulations used in the study.

## Software

All analyses are implemented in **R**, using Bayesian modelling tools including:

- `nimble`
- standard R packages for data manipulation and analysis.

## Notes

This repository is intended to facilitate **transparency and reproducibility** of the results presented in the associated research article.
