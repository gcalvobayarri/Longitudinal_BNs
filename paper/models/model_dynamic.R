library(nimble, warn.conflicts = FALSE)
load('paper/data/matrix_data_PHI_13players.RData')

z <- matrix(data=NA, nrow=13, ncol= 82)
z[MP==0]<-0
z[MP!=0]<-1


code_logit <- nimbleCode({
  for (i in 1:N) {
    for (j in 1:T) {
      C1[i, j] ~ dbin(p1[i, j], T1[i, j])
      C2[i, j] ~ dbin(p2[i, j], T2[i, j])
      C3[i, j] ~ dbin(p3[i, j], T3[i, j])
      
      logit(p1[i, j]) <- beta0[1] + betaSG[1]*SG[i, j] + betaPG[1]*PG[i, j] +
        betaSF[1]*SF[i, j] + betaPF[1]*PF[i, j] +
        betaH[1]*home[i, j] + b0[i, 1]
      
      logit(p2[i, j]) <- beta0[2] + betaSG[2]*SG[i, j] + betaPG[2]*PG[i, j] +
        betaSF[2]*SF[i, j] + betaPF[2]*PF[i, j] +
        betaH[2]*home[i, j] + b0[i, 2]
      
      logit(p3[i, j]) <- beta0[3] + betaSG[3]*SG[i, j] + betaPG[3]*PG[i, j] +
        betaSF[3]*SF[i, j] + betaPF[3]*PF[i, j] +
        betaH[3]*home[i, j] + b0[i, 3]
    }
  }
  
  for (i in 1:N) {
    for (j in 1:T) {
      z[i, j] ~ dbern(psi[i])
      
      T1[i, j] ~ dpois(lambda1[i, j])
      log(lambda1[i, j]) <- alpha0[1] + alphaSG[1]*SG[i, j] + alphaPG[1]*PG[i, j] +
        alphaSF[1]*SF[i, j] + alphaPF[1]*PF[i, j] +
        a0[i, 1] + alpha[1]*foulsR[i, j] + log(z[i, j] + 1.0E-7) # revisar 
      
      T2[i, j] ~ dpois(lambda2[i, j])
      log(lambda2[i, j]) <- alpha0[2] + alphaSG[2]*SG[i, j] + alphaPG[2]*PG[i, j] +
        alphaSF[2]*SF[i, j] + alphaPF[2]*PF[i, j] +
        a0[i, 2] + alpha[2]*MP[i, j] + log(z[i, j] + 1.0E-7)
      
      T3[i, j] ~ dpois(lambda3[i, j])
      log(lambda3[i, j]) <- alpha0[3] + alphaSG[3]*SG[i, j] + alphaPG[3]*PG[i, j] +
        alphaSF[3]*SF[i, j] + alphaPF[3]*PF[i, j] +
        a0[i, 3] + alpha[3]*MP[i, j] + log(z[i, j] + 1.0E-7)
    }
  }
  
  for (i in 1:N) {
    for (j in 1:T) {
      foulsR[i, j] ~ dpois(lambdaF[i, j])
      log(lambdaF[i, j]) <- gamma0 + c0[i] + gamma*MP[i, j] + log(z[i, j] + 1.0E-7)
    }
  }
  
  for (i in 1:N) {
    MP[i, 1] ~ dpois(lambdaM[i, 1])
    log(lambdaM[i, 1]) <- delta0 + d0[i] + log(z[i, 1] + 1.0E-7)
    for (j in 2:T) {
      MP[i, j] ~ dpois(lambdaM[i, j])
      log(lambdaM[i, j]) <- delta0 + d0[i] + deltaw*log(1 + MP[i, j-1]) + log(z[i, j] + 1.0E-7)
    }
  }
  
  ## Priors 
  for (k in 1:3) {
    beta0[k]  ~ dnorm(0, 0.01)
    betaSG[k] ~ dnorm(0, 0.01)
    betaPG[k] ~ dnorm(0, 0.01)
    betaSF[k] ~ dnorm(0, 0.01)
    betaPF[k] ~ dnorm(0, 0.01)
    betaH[k]  ~ dnorm(0, 0.01)
    
    alpha0[k] ~ dnorm(0, 0.01)
    alpha[k]  ~ dnorm(0, 0.01)
    alphaSG[k] ~ dnorm(0, 0.01)
    alphaPG[k] ~ dnorm(0, 0.01)
    alphaSF[k] ~ dnorm(0, 0.01)
    alphaPF[k] ~ dnorm(0, 0.01)
  }
  gamma0 ~ dnorm(0, 0.01)
  gamma  ~ dnorm(0, 0.01)
  delta0 ~ dnorm(0, 0.01)
  deltaw ~ dnorm(0, 0.01)
  
  for (k in 1:3) {
    for (i in 1:N) {
      b0[i, k] ~ dnorm(0, lambdab[k])
    }
    lambdab[k] <- 1 / (sigmab[k] * sigmab[k])
    sigmab[k] ~ dunif(0, 100)
  }
  for (k in 1:3) {
    for (i in 1:N) {
      a0[i, k] ~ dnorm(0, lambdaa[k])
    }
    lambdaa[k] <- 1 / (sigmaa[k] * sigmaa[k])
    sigmaa[k] ~ dunif(0, 100)
  }
  for (i in 1:N) {
    c0[i] ~ dnorm(0, lambdac)
    d0[i] ~ dnorm(0, lambdad)
    psi[i] ~ dunif(0, 1)
    p[i] <- 1 - psi[i]
  }
  lambdac <- 1 / (sigmac * sigmac)
  sigmac ~ dunif(0, 100)
  lambdad <- 1 / (sigmad * sigmad)
  sigmad ~ dunif(0, 100)
})



# Computing model-------------------

data_model <- list(N = dim(MP)[1], T = dim(MP)[2],
                   C1 = C1,
                   C2 = C2,
                   C3 = C3,
                   T1 = T1,
                   T2 = T2,
                   T3 = T3,
                   foulsR = foulsR,
                   MP = round(MP),
                   home = home,
                   PF = PF,
                   SF = SF,
                   PG = PG,
                   SG = SG, z=z)

parameters <- c( "beta0", "betaH",  "b0", #w1", "w2", "w3",
                 "alpha0", "alpha", "a0",
                 "gamma0", "gamma", "c0", "p", 
                 "delta0", "d0", "deltaw", 'betaSG', 'betaPG', 'betaSF', 'betaPF',
                 'alphaSG', 'alphaPG', 'alphaSF', 'alphaPF',
                 "sigmab", "sigmaa", "sigmac", "sigmad")

inits <- list( beta0=c(0.5, 0.5, -3), betaH=c(0, 0, 0), alpha0=c(0.5, 0.5, -5),
        alpha=c(0, 0, 0), gamma0=-1, gamma=0, delta0 = 3, deltaw = 0, sigmab=c(1, 1, 1),
        sigmaa=c(1, 1, 1), sigmac=1, sigmad=1)


nimbleMCMC_samples <- nimbleMCMC(code = code_logit, 
                                 constants = data_model, 
                                 
                                 inits = inits,
                                 nburnin = 500000, niter = 1000000,
                                 summary = T, thin = 500,
                                 monitors = parameters,
                                 nchains = 3,
                                 WAIC = T)

save(nimbleMCMC_samples, file = 'paper/results/nimbleMCMCDynamic_samples.RData')

load('paper/results/nimbleMCMCDynamic_samples.RData')

nimbleMCMC_samples$summary$all.chains 

nimbleMCMC_samples$WAIC$WAIC
