# Posterior results
library(nimble, warn.conflicts = FALSE)
library(mcmcplots)
library(ggplot2)
library(MCMCvis)
library(coda)

load('paper/results/nimbleMCMCDynamic_samples.RData')

# 0. Posterior summary-----------

nimbleMCMC_samples$summary$all.chains

# 1. Parameters and random effects-----------

library(coda)
samples_mcmc_list <- mcmc.list(lapply(nimbleMCMC_samples$samples, mcmc))


all_params <- varnames(samples_mcmc_list)

## Unir todas las cadenas en una sola matriz
posterior_mat <- do.call(rbind, samples_mcmc_list)

posterior <- lapply(all_params, function(p) posterior_mat[, p])
names(posterior) <- all_params

#2. Probabilidad posterior P(parametro > 0) para un subconjunto


# 1) Parámetros de interés 
params_interes <- c(
  paste0("beta0[", 1:3, "]"),
  paste0("betaH[", 1:3, "]"),
  paste0("betaPF[", 1:3, "]"),
  paste0("betaPG[", 1:3, "]"),
  paste0("betaSF[", 1:3, "]"),
  paste0("betaSG[", 1:3, "]"),
  paste0("sigmab[", 1:3, "]")
)


prob_mayor_0 <- sapply(posterior[params_interes], function(x) mean(x > 0))

# Resultado en data.frame ordenado
resultado <- data.frame(
  parametro = names(prob_mayor_0),
  P_mayor_0 = as.numeric(prob_mayor_0),
  row.names = NULL
)

resultado <- resultado[order(resultado$parametro), ]

print(resultado)
