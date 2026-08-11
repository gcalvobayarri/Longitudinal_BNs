make_model_inputs <- function(data = load_raw_data()) {
  participation <- (data$MP > 0) * 1
  model_data <- list(
    C1 = data$C1, C2 = data$C2, C3 = data$C3,
    T1 = data$T1, T2 = data$T2, T3 = data$T3,
    foulsR = data$foulsR, MP = round(data$MP), z = participation
  )
  constants <- list(
    N = nrow(data$MP), T = ncol(data$MP),
    home = data$home, PF = data$PF, SF = data$SF, PG = data$PG, SG = data$SG
  )
  list(data = model_data, constants = constants)
}

static_lbn_code <- function() {
  nimble::nimbleCode({
    for (i in 1:N) {
      for (j in 1:T) {
        C1[i, j] ~ dbin(p1[i, j], T1[i, j])
        C2[i, j] ~ dbin(p2[i, j], T2[i, j])
        C3[i, j] ~ dbin(p3[i, j], T3[i, j])
        logit(p1[i, j]) <- beta0[1] + betaSG[1] * SG[i, j] + betaPG[1] * PG[i, j] + betaSF[1] * SF[i, j] + betaPF[1] * PF[i, j] + betaH[1] * home[i, j] + b0[i, 1]
        logit(p2[i, j]) <- beta0[2] + betaSG[2] * SG[i, j] + betaPG[2] * PG[i, j] + betaSF[2] * SF[i, j] + betaPF[2] * PF[i, j] + betaH[2] * home[i, j] + b0[i, 2]
        logit(p3[i, j]) <- beta0[3] + betaSG[3] * SG[i, j] + betaPG[3] * PG[i, j] + betaSF[3] * SF[i, j] + betaPF[3] * PF[i, j] + betaH[3] * home[i, j] + b0[i, 3]

        z[i, j] ~ dbern(psi[i])
        T1[i, j] ~ dpois(lambda1[i, j])
        T2[i, j] ~ dpois(lambda2[i, j])
        T3[i, j] ~ dpois(lambda3[i, j])
        lambda1[i, j] <- exp(alpha0[1] + alphaSG[1] * SG[i, j] + alphaPG[1] * PG[i, j] + alphaSF[1] * SF[i, j] + alphaPF[1] * PF[i, j] + a0[i, 1] + alpha[1] * foulsR[i, j]) * z[i, j] + 1.0E-7
        lambda2[i, j] <- exp(alpha0[2] + alphaSG[2] * SG[i, j] + alphaPG[2] * PG[i, j] + alphaSF[2] * SF[i, j] + alphaPF[2] * PF[i, j] + a0[i, 2] + alpha[2] * MP[i, j]) * z[i, j] + 1.0E-7
        lambda3[i, j] <- exp(alpha0[3] + alphaSG[3] * SG[i, j] + alphaPG[3] * PG[i, j] + alphaSF[3] * SF[i, j] + alphaPF[3] * PF[i, j] + a0[i, 3] + alpha[3] * MP[i, j]) * z[i, j] + 1.0E-7

        foulsR[i, j] ~ dpois(lambdaF[i, j])
        lambdaF[i, j] <- exp(gamma0 + c0[i] + gamma * MP[i, j]) * z[i, j] + 1.0E-7
        MP[i, j] ~ dpois(lambdaM[i, j])
        lambdaM[i, j] <- exp(delta0 + d0[i]) * z[i, j] + 1.0E-7
      }
    }

    for (k in 1:3) {
      beta0[k] ~ dnorm(0, 0.01); betaSG[k] ~ dnorm(0, 0.01); betaPG[k] ~ dnorm(0, 0.01); betaSF[k] ~ dnorm(0, 0.01); betaPF[k] ~ dnorm(0, 0.01); betaH[k] ~ dnorm(0, 0.01)
      alpha0[k] ~ dnorm(0, 0.01); alpha[k] ~ dnorm(0, 0.01); alphaSG[k] ~ dnorm(0, 0.01); alphaPG[k] ~ dnorm(0, 0.01); alphaSF[k] ~ dnorm(0, 0.01); alphaPF[k] ~ dnorm(0, 0.01)
      for (i in 1:N) { b0[i, k] ~ dnorm(0, lambdab[k]); a0[i, k] ~ dnorm(0, lambdaa[k]) }
      lambdab[k] <- pow(sigmab[k], -2); sigmab[k] ~ dunif(0, 100)
      lambdaa[k] <- pow(sigmaa[k], -2); sigmaa[k] ~ dunif(0, 100)
    }
    gamma0 ~ dnorm(0, 0.01); gamma ~ dnorm(0, 0.01); delta0 ~ dnorm(0, 0.01)
    for (i in 1:N) { c0[i] ~ dnorm(0, lambdac); d0[i] ~ dnorm(0, lambdad); psi[i] ~ dunif(0, 1); p[i] <- 1 - psi[i] }
    lambdac <- pow(sigmac, -2); sigmac ~ dunif(0, 100)
    lambdad <- pow(sigmad, -2); sigmad ~ dunif(0, 100)
  })
}

dynamic_lbn_code <- function() {
  # The explicit model below differs from the static LBN only in the minutes node.
  nimble::nimbleCode({
    for (i in 1:N) {
      for (j in 1:T) {
        C1[i, j] ~ dbin(p1[i, j], T1[i, j]); C2[i, j] ~ dbin(p2[i, j], T2[i, j]); C3[i, j] ~ dbin(p3[i, j], T3[i, j])
        logit(p1[i, j]) <- beta0[1] + betaSG[1] * SG[i, j] + betaPG[1] * PG[i, j] + betaSF[1] * SF[i, j] + betaPF[1] * PF[i, j] + betaH[1] * home[i, j] + b0[i, 1]
        logit(p2[i, j]) <- beta0[2] + betaSG[2] * SG[i, j] + betaPG[2] * PG[i, j] + betaSF[2] * SF[i, j] + betaPF[2] * PF[i, j] + betaH[2] * home[i, j] + b0[i, 2]
        logit(p3[i, j]) <- beta0[3] + betaSG[3] * SG[i, j] + betaPG[3] * PG[i, j] + betaSF[3] * SF[i, j] + betaPF[3] * PF[i, j] + betaH[3] * home[i, j] + b0[i, 3]
        z[i, j] ~ dbern(psi[i])
        T1[i, j] ~ dpois(lambda1[i, j]); T2[i, j] ~ dpois(lambda2[i, j]); T3[i, j] ~ dpois(lambda3[i, j])
        lambda1[i, j] <- exp(alpha0[1] + alphaSG[1] * SG[i, j] + alphaPG[1] * PG[i, j] + alphaSF[1] * SF[i, j] + alphaPF[1] * PF[i, j] + a0[i, 1] + alpha[1] * foulsR[i, j]) * z[i, j] + 1.0E-7
        lambda2[i, j] <- exp(alpha0[2] + alphaSG[2] * SG[i, j] + alphaPG[2] * PG[i, j] + alphaSF[2] * SF[i, j] + alphaPF[2] * PF[i, j] + a0[i, 2] + alpha[2] * MP[i, j]) * z[i, j] + 1.0E-7
        lambda3[i, j] <- exp(alpha0[3] + alphaSG[3] * SG[i, j] + alphaPG[3] * PG[i, j] + alphaSF[3] * SF[i, j] + alphaPF[3] * PF[i, j] + a0[i, 3] + alpha[3] * MP[i, j]) * z[i, j] + 1.0E-7
        foulsR[i, j] ~ dpois(lambdaF[i, j]); lambdaF[i, j] <- exp(gamma0 + c0[i] + gamma * MP[i, j]) * z[i, j] + 1.0E-7
      }
      MP[i, 1] ~ dpois(lambdaM[i, 1]); lambdaM[i, 1] <- exp(delta0 + d0[i]) * z[i, 1] + 1.0E-7
      for (j in 2:T) { MP[i, j] ~ dpois(lambdaM[i, j]); lambdaM[i, j] <- exp(delta0 + d0[i] + deltaw * log(1 + MP[i, j - 1])) * z[i, j] + 1.0E-7 }
    }
    for (k in 1:3) {
      beta0[k] ~ dnorm(0, 0.01); betaSG[k] ~ dnorm(0, 0.01); betaPG[k] ~ dnorm(0, 0.01); betaSF[k] ~ dnorm(0, 0.01); betaPF[k] ~ dnorm(0, 0.01); betaH[k] ~ dnorm(0, 0.01)
      alpha0[k] ~ dnorm(0, 0.01); alpha[k] ~ dnorm(0, 0.01); alphaSG[k] ~ dnorm(0, 0.01); alphaPG[k] ~ dnorm(0, 0.01); alphaSF[k] ~ dnorm(0, 0.01); alphaPF[k] ~ dnorm(0, 0.01)
      for (i in 1:N) { b0[i, k] ~ dnorm(0, lambdab[k]); a0[i, k] ~ dnorm(0, lambdaa[k]) }
      lambdab[k] <- pow(sigmab[k], -2); sigmab[k] ~ dunif(0, 100); lambdaa[k] <- pow(sigmaa[k], -2); sigmaa[k] ~ dunif(0, 100)
    }
    gamma0 ~ dnorm(0, 0.01); gamma ~ dnorm(0, 0.01); delta0 ~ dnorm(0, 0.01); deltaw ~ dnorm(0, 0.01)
    for (i in 1:N) { c0[i] ~ dnorm(0, lambdac); d0[i] ~ dnorm(0, lambdad); psi[i] ~ dunif(0, 1); p[i] <- 1 - psi[i] }
    lambdac <- pow(sigmac, -2); sigmac ~ dunif(0, 100); lambdad <- pow(sigmad, -2); sigmad ~ dunif(0, 100)
  })
}

hidden_markov_lbn_code <- function() {
  nimble::nimbleCode({
    for (i in 1:N) {
      for (j in 1:T) {
        C1[i, j] ~ dbin(p1[i, j], T1[i, j]); C2[i, j] ~ dbin(p2[i, j], T2[i, j]); C3[i, j] ~ dbin(p3[i, j], T3[i, j])
        logit(p1[i, j]) <- beta0[1, Z[j]] + betaSG[1] * SG[i, j] + betaPG[1] * PG[i, j] + betaSF[1] * SF[i, j] + betaPF[1] * PF[i, j] + betaH[1] * home[i, j] + b0[i, 1]
        logit(p2[i, j]) <- beta0[2, Z[j]] + betaSG[2] * SG[i, j] + betaPG[2] * PG[i, j] + betaSF[2] * SF[i, j] + betaPF[2] * PF[i, j] + betaH[2] * home[i, j] + b0[i, 2]
        logit(p3[i, j]) <- beta0[3, Z[j]] + betaSG[3] * SG[i, j] + betaPG[3] * PG[i, j] + betaSF[3] * SF[i, j] + betaPF[3] * PF[i, j] + betaH[3] * home[i, j] + b0[i, 3]
        z[i, j] ~ dbern(psi[i])
        T1[i, j] ~ dpois(lambda1[i, j]); T2[i, j] ~ dpois(lambda2[i, j]); T3[i, j] ~ dpois(lambda3[i, j])
        lambda1[i, j] <- exp(alpha0[1] + alphaSG[1] * SG[i, j] + alphaPG[1] * PG[i, j] + alphaSF[1] * SF[i, j] + alphaPF[1] * PF[i, j] + a0[i, 1] + alpha[1] * foulsR[i, j]) * z[i, j] + 1.0E-7
        lambda2[i, j] <- exp(alpha0[2] + alphaSG[2] * SG[i, j] + alphaPG[2] * PG[i, j] + alphaSF[2] * SF[i, j] + alphaPF[2] * PF[i, j] + a0[i, 2] + alpha[2] * MP[i, j]) * z[i, j] + 1.0E-7
        lambda3[i, j] <- exp(alpha0[3] + alphaSG[3] * SG[i, j] + alphaPG[3] * PG[i, j] + alphaSF[3] * SF[i, j] + alphaPF[3] * PF[i, j] + a0[i, 3] + alpha[3] * MP[i, j]) * z[i, j] + 1.0E-7
        foulsR[i, j] ~ dpois(lambdaF[i, j]); lambdaF[i, j] <- exp(gamma0 + c0[i] + gamma * MP[i, j]) * z[i, j] + 1.0E-7
        MP[i, j] ~ dpois(lambdaM[i, j]); lambdaM[i, j] <- exp(delta0 + d0[i]) * z[i, j] + 1.0E-7
      }
    }
    for (k in 1:3) {
      beta0.p[k, 1] ~ dnorm(0, 0.1); beta0.p[k, 2] ~ dnorm(0, 0.1)
      beta0[k, 1] <- min(beta0.p[k, 1:2]); beta0[k, 2] <- max(beta0.p[k, 1:2])
      betaSG[k] ~ dnorm(0, 0.01); betaPG[k] ~ dnorm(0, 0.01); betaSF[k] ~ dnorm(0, 0.01); betaPF[k] ~ dnorm(0, 0.01); betaH[k] ~ dnorm(0, 0.01)
      alpha0[k] ~ dnorm(0, 0.01); alpha[k] ~ dnorm(0, 0.01); alphaSG[k] ~ dnorm(0, 0.01); alphaPG[k] ~ dnorm(0, 0.01); alphaSF[k] ~ dnorm(0, 0.01); alphaPF[k] ~ dnorm(0, 0.01)
      for (i in 1:N) { b0[i, k] ~ dnorm(0, lambdab[k]); a0[i, k] ~ dnorm(0, lambdaa[k]) }
      lambdab[k] <- pow(sigmab[k], -2); sigmab[k] ~ dunif(0, 100); lambdaa[k] <- pow(sigmaa[k], -2); sigmaa[k] ~ dunif(0, 100)
    }
    gamma0 ~ dnorm(0, 0.01); gamma ~ dnorm(0, 0.01); delta0 ~ dnorm(0, 0.01)
    for (i in 1:N) { c0[i] ~ dnorm(0, lambdac); d0[i] ~ dnorm(0, lambdad); psi[i] ~ dunif(0, 1); p[i] <- 1 - psi[i] }
    lambdac <- pow(sigmac, -2); sigmac ~ dunif(0, 100); lambdad <- pow(sigmad, -2); sigmad ~ dunif(0, 100)
    Z[1] ~ dcat(p0[1:2]); p0[1] ~ dbeta(1, 1); p0[2] <- 1 - p0[1]
    for (j in 2:T) { Z[j] ~ dcat(omega[Z[j - 1], 1:2]) }
    omega[1, 1] <- 1 - pCH; omega[1, 2] <- pCH; omega[2, 1] <- pHC; omega[2, 2] <- 1 - pHC
    pCH ~ dbeta(1, 1); pHC ~ dbeta(1, 1)
  })
}

model_monitors <- function(model = c("static", "dynamic", "hidden_markov")) {
  model <- match.arg(model)
  common <- c("beta0", "betaH", "b0", "alpha0", "alpha", "a0", "gamma0", "gamma", "c0", "p", "delta0", "d0", "betaSG", "betaPG", "betaSF", "betaPF", "alphaSG", "alphaPG", "alphaSF", "alphaPF", "sigmab", "sigmaa", "sigmac", "sigmad")
  if (model == "dynamic") common <- c(common, "deltaw")
  if (model == "hidden_markov") common <- c(common, "beta0.p", "pCH", "pHC", "p0")
  common
}

model_inits <- function(model = c("static", "dynamic", "hidden_markov"), n_games = 82L, seed = paper_seed) {
  model <- match.arg(model)
  set.seed(seed)
  common <- list(betaH = c(0, 0, 0), alpha0 = c(0.5, 0.5, -5), alpha = c(0, 0, 0), gamma0 = -1, gamma = 0, delta0 = 3, sigmab = rep(1, 3), sigmaa = rep(1, 3), sigmac = 1, sigmad = 1)
  if (model != "hidden_markov") common$beta0 <- c(0.5, 0.5, -3)
  if (model == "dynamic") common$deltaw <- 0
  if (model == "hidden_markov") {
    common$beta0.p <- matrix(c(0.45, 0.45, -3.5, 0.5, 0.5, -3), nrow = 3, ncol = 2)
    common$pCH <- 0.5; common$pHC <- 0.5; common$p0 <- c(0.5, 0.5); common$Z <- sample(1:2, n_games, replace = TRUE)
  }
  common
}

model_code <- function(model = c("static", "dynamic", "hidden_markov")) {
  model <- match.arg(model)
  switch(model, static = static_lbn_code(), dynamic = dynamic_lbn_code(), hidden_markov = hidden_markov_lbn_code())
}
