# Posterior predictive distribution
library(nimble, warn.conflicts = FALSE)
library(mcmcplots)
library(ggplot2)
library(MCMCvis)
library(coda)

load('paper/results/nimbleMCMCDynamic_samples.RData')

# 1. Parameters and random effects-----------

library(coda)
samples_mcmc_list <- mcmc.list(lapply(nimbleMCMC_samples$samples, mcmc))


all_params <- varnames(samples_mcmc_list)

## Unir todas las cadenas en una sola matriz
posterior_mat <- do.call(rbind, samples_mcmc_list)

posterior <- lapply(all_params, function(p) posterior_mat[, p])
names(posterior) <- all_params


## Pasar a data.frame (opcional, pero suele ser cómodo)
#posterior <- as.data.frame(samples_mat)

## Comprobaciones
str(posterior)
head(posterior)
colnames(posterior)   # deberían coincidir con all_params


# 2. Posterior predictive distributions-------------------

n_players <- 13
n_iter <- 3000
J <- 2

a  <- array(NA, dim = c(n_players, n_iter, J))
m  <- array(NA, dim = c(n_players, n_iter, J))
f  <- array(NA, dim = c(n_players, n_iter, J))
t1 <- array(NA, dim = c(n_players, n_iter, J))
t2 <- array(NA, dim = c(n_players, n_iter, J))
t3 <- array(NA, dim = c(n_players, n_iter, J))
c1 <- array(NA, dim = c(n_players, n_iter, J))
c2 <- array(NA, dim = c(n_players, n_iter, J))
c3 <- array(NA, dim = c(n_players, n_iter, J))

p1 <- array(NA, dim = c(n_players, n_iter, J))
p2 <- array(NA, dim = c(n_players, n_iter, J))
p3 <- array(NA, dim = c(n_players, n_iter, J))

#covariates
PG <- c(1,1,0,0,0,0,0,1,0,0,0,0,0)
SG <- c(0,0,0,1,1,0,0,0,0,0,0,0,0)
SF <- c(0,0,0,0,0,1,0,0,0,0,0,0,0)
PF <- c(0,0,1,0,0,0,1,0,1,1,0,1,0)
home <- c(1,1
)




for(i in 1:n_players){
  for(s in 1:n_iter){
    for(j in 1){
      # indicador de participación
      a[i, s, j] <- rbinom(1, 1, prob = (1 - posterior[[sprintf("p[%d]", i)]][s])) #fallo
      
      # minutos jugados
      m[i, s, j] <- rpois(1, lambda = exp(posterior$delta0[s] + posterior[[sprintf("d0[%d]", i)]][s]) * a[i, s, j] + 1e-7)
      
      # faltas cometidas
      f[i, s, j] <- rpois(1, lambda = exp(posterior$gamma0[s] +
                                            posterior[[sprintf("c0[%d]", i)]][s] +
                                            posterior$gamma[s] * m[i, s, j]) * a[i, s, j] + 1e-7)
      
      # tiros
      # de 1
      t1[i, s, j] <- rpois(1, lambda = exp(posterior$`alpha0[1]`[s] + 
                                             posterior$`alphaSG[1]`[s] * SG[i] +
                                             posterior$`alphaPG[1]`[s] * PG[i] +
                                             posterior$`alphaSF[1]`[s] * SF[i] +
                                             posterior$`alphaPF[1]`[s] * PF[i] +
                                             posterior[[sprintf("a0[%d, 1]", i)]][s] +
                                             posterior$`alpha[1]`[s] * f[i, s, j]) * a[i, s, j] + 1e-7)
      
      # de 2
      t2[i, s, j] <- rpois(1, lambda = exp(posterior$`alpha0[2]`[s] + 
                                             posterior$`alphaSG[2]`[s] * SG[i] +
                                             posterior$`alphaPG[2]`[s] * PG[i] +
                                             posterior$`alphaSF[2]`[s] * SF[i] +
                                             posterior$`alphaPF[2]`[s] * PF[i] +
                                             posterior[[sprintf("a0[%d, 2]", i)]][s] +
                                             posterior$`alpha[2]`[s] * m[i, s, j]) * a[i, s, j] + 1e-7)
      
      # de 3
      t3[i, s, j] <- rpois(1, lambda = exp(posterior$`alpha0[3]`[s] + 
                                             posterior$`alphaSG[3]`[s] * SG[i] +
                                             posterior$`alphaPG[3]`[s] * PG[i] +
                                             posterior$`alphaSF[3]`[s] * SF[i] +
                                             posterior$`alphaPF[3]`[s] * PF[i] +
                                             posterior[[sprintf("a0[%d, 3]", i)]][s] +
                                             posterior$`alpha[3]`[s] * m[i, s, j]) * a[i, s, j] + 1e-7)
      
      #canastas
      #de 1
      p1[i, s, j] <- 1 / (1 + exp(-(posterior$`beta0[1]`[s] + 
                                      posterior$`betaSG[1]`[s] * SG[i] +
                                      posterior$`betaPG[1]`[s] * PG[i] +
                                      posterior$`betaSF[1]`[s] * SF[i] +
                                      posterior$`betaPF[1]`[s] * PF[i] +
                                      posterior[[sprintf("b0[%d, 1]", i)]][s] +
                                      posterior$`betaH[3]`[s] * home[j])))
      
      c1[i, s, j] <- rbinom(1, size = t1[i,s, j], prob = p1[i, s, j])
      
      #de 2
      p2[i, s, j] <- 1 / (1 + exp(-(posterior$`beta0[2]`[s] + 
                                      posterior$`betaSG[2]`[s] * SG[i] +
                                      posterior$`betaPG[2]`[s] * PG[i] +
                                      posterior$`betaSF[2]`[s] * SF[i] +
                                      posterior$`betaPF[2]`[s] * PF[i] +
                                      posterior[[sprintf("b0[%d, 2]", i)]][s] +
                                      posterior$`betaH[2]`[s] * home[j])))
      
      c2[i, s, j] <- rbinom(1, size = t2[i,s,j], prob = p2[i, s, j])
      
      #de 3
      p3[i, s, j] <- 1 / (1 + exp(-(posterior$`beta0[3]`[s] + 
                                      posterior$`betaSG[3]`[s] * SG[i] +
                                      posterior$`betaPG[3]`[s] * PG[i] +
                                      posterior$`betaSF[3]`[s] * SF[i] +
                                      posterior$`betaPF[3]`[s] * PF[i] +
                                      posterior[[sprintf("b0[%d, 3]", i)]][s] +
                                      posterior$`betaH[3]`[s] * home[j])))
      
      c3[i, s, j] <- rbinom(1, size = t3[i,s,j], prob = p3[i, s, j])
    }
  }
}


# j = 2.....
for(i in 1:n_players){
  for(s in 1:n_iter){
    for(j in J){
    # indicador de participación
    a[i, s, j] <- rbinom(1, 1, prob = (1 - posterior[[sprintf("p[%d]", i)]][s]))
    
    # minutos jugados
    m[i, s, j] <- rpois(1, lambda = exp(posterior$delta0[s] + posterior$deltaw[s] * m[i, s, j-1] + posterior[[sprintf("d0[%d]", i)]][s]) * a[i, s, j] + 1e-7)
    
    # faltas cometidas
    f[i, s, j] <- rpois(1, lambda = exp(posterior$gamma0[s] +
                                       posterior[[sprintf("c0[%d]", i)]][s] +
                                       posterior$gamma[s] * m[i, s, j]) * a[i, s, j] + 1e-7)
    
    # tiros
    # de 1
    t1[i, s, j] <- rpois(1, lambda = exp(posterior$`alpha0[1]`[s] + 
                                        posterior$`alphaSG[1]`[s] * SG[i] +
                                        posterior$`alphaPG[1]`[s] * PG[i] +
                                        posterior$`alphaSF[1]`[s] * SF[i] +
                                        posterior$`alphaPF[1]`[s] * PF[i] +
                                        posterior[[sprintf("a0[%d, 1]", i)]][s] +
                                        posterior$`alpha[1]`[s] * f[i, s, j]) * a[i, s, j] + 1e-7)
    
    # de 2
    t2[i, s, j] <- rpois(1, lambda = exp(posterior$`alpha0[2]`[s] + 
                                        posterior$`alphaSG[2]`[s] * SG[i] +
                                        posterior$`alphaPG[2]`[s] * PG[i] +
                                        posterior$`alphaSF[2]`[s] * SF[i] +
                                        posterior$`alphaPF[2]`[s] * PF[i] +
                                        posterior[[sprintf("a0[%d, 2]", i)]][s] +
                                        posterior$`alpha[2]`[s] * m[i, s, j]) * a[i, s, j] + 1e-7)
    
    # de 3
    t3[i, s, j] <- rpois(1, lambda = exp(posterior$`alpha0[3]`[s] + 
                                        posterior$`alphaSG[3]`[s] * SG[i] +
                                        posterior$`alphaPG[3]`[s] * PG[i] +
                                        posterior$`alphaSF[3]`[s] * SF[i] +
                                        posterior$`alphaPF[3]`[s] * PF[i] +
                                        posterior[[sprintf("a0[%d, 3]", i)]][s] +
                                        posterior$`alpha[3]`[s] * m[i, s, j]) * a[i, s, j] + 1e-7)
    
    #canastas
    #de 1
    p1[i, s, j] <- 1 / (1 + exp(-(posterior$`beta0[1]`[s] + 
                                 posterior$`betaSG[1]`[s] * SG[i] +
                                 posterior$`betaPG[1]`[s] * PG[i] +
                                 posterior$`betaSF[1]`[s] * SF[i] +
                                 posterior$`betaPF[1]`[s] * PF[i] +
                                 posterior[[sprintf("b0[%d, 1]", i)]][s] +
                                 posterior$`betaH[3]`[s] * home[j])))
    
    c1[i, s, j] <- rbinom(1, size = t1[i,s, j], prob = p1[i, s, j])
    
    #de 2
    p2[i, s, j] <- 1 / (1 + exp(-(posterior$`beta0[2]`[s] + 
                                 posterior$`betaSG[2]`[s] * SG[i] +
                                 posterior$`betaPG[2]`[s] * PG[i] +
                                 posterior$`betaSF[2]`[s] * SF[i] +
                                 posterior$`betaPF[2]`[s] * PF[i] +
                                 posterior[[sprintf("b0[%d, 2]", i)]][s] +
                                 posterior$`betaH[2]`[s] * home[j])))
    
    c2[i, s, j] <- rbinom(1, size = t2[i,s,j], prob = p2[i, s, j])
    
    #de 3
    p3[i, s, j] <- 1 / (1 + exp(-(posterior$`beta0[3]`[s] + 
                                 posterior$`betaSG[3]`[s] * SG[i] +
                                 posterior$`betaPG[3]`[s] * PG[i] +
                                 posterior$`betaSF[3]`[s] * SF[i] +
                                 posterior$`betaPF[3]`[s] * PF[i] +
                                 posterior[[sprintf("b0[%d, 3]", i)]][s] +
                                 posterior$`betaH[3]`[s] * home[j])))
    
    c3[i, s, j] <- rbinom(1, size = t3[i,s,j], prob = p3[i, s, j])
    }
  }
}


library(dplyr)

# Base de combinaciones: jugador más lento, iteración más rápido
base_df <- expand.grid(
  iteracion = 1:n_iter,
  jugador = 1:n_players,
  j =1:2
) %>% arrange(j, jugador, iteracion)


posterior_predictive <- base_df %>%
  mutate(
    minutos    = c(as.vector(t(m[,,1])), as.vector(t(m[,,2]))),  # t() transpone para que se alineen, lo añada por filas
    faltas     = c(as.vector(t(f[,,1])),  as.vector(t(f[,,2]))),
    tiros1     = c(as.vector(t(t1[,,1])), as.vector(t(t1[,,2]))),
    canastas1  = c(as.vector(t(c1[,,1])), as.vector(t(c1[,,2]))),
    tiros2     = c(as.vector(t(t2[,,1])), as.vector(t(t2[,,2]))),
    canastas2  = c(as.vector(t(c2[,,1])), as.vector(t(c2[,,2]))),
    tiros3     = c(as.vector(t(t3[,,1])), as.vector(t(t3[,,2]))),
    canastas3  = c(as.vector(t(c3[,,1])), as.vector(t(c3[,,2]))),
    puntos     = canastas1 + 2 * canastas2 + 3 * canastas3
  )


# el tiempo 2 se descalabra es normal 

# los puntos del jugador 1 y 6 dado que juegan más de 30 minutos 

library(dplyr)
sim_j1_30m <- posterior_predictive %>%
  filter(jugador == 1, minutos <= 48, minutos >=30, j==1)

ggplot(data=sim_j1_30m, aes(x = puntos)) +
  geom_bar(aes(y = ..prop.., group = 1)) +
  labs(
    title = "Distribución predictiva de puntos",
    subtitle = "Condicionada a que el jugador 1 jugara más de 30 minutos",
    x = "Puntos", y = "Frecuencia relativa"
  ) + xlim(-1, 60) +
  theme_minimal()


sim_j6_30m <- posterior_predictive %>%
  filter(jugador == 6, minutos <= 48, minutos >=30, j==1)



ggplot(data=sim_j6_30m, aes(x = puntos)) +
  geom_bar(aes(y = ..prop.., group = 1)) +
  labs(
    title = "Distribución predictiva de puntos",
    subtitle = "Condicionada a que el jugador 6 jugara más de 30 minutos",
    x = "Puntos", y = "Frecuencia relativa"
  ) + xlim(-1, 60) +
  theme_minimal()


# ejemplo 1: 2 jugadores --------------

library(dplyr)
library(ggplot2)
library(viridis)


# Filter and rename players
sim_j1j6_30m <- posterior_predictive %>%
  filter(jugador %in% c(1, 6), minutos >= 30, minutos <= 48, j==1) %>%
  mutate(player = case_when(
    jugador == 1 ~ "Allen Iverson",
    jugador == 6 ~ "Kyle Korver"
  ))

# Plot
ggplot(sim_j1j6_30m, aes(x = puntos, fill = player)) +
  geom_bar(aes(y = ..prop.., group = player), position = "dodge") +
  labs(
    x = "Points", y = "Relative Frequency",
    fill = "Player"
  ) +
  xlim(-1, 60) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

ggplot(sim_j1j6_30m, aes(x = puntos, fill = player)) +
  geom_bar(aes(y = ..prop.., group = player), position = "identity", alpha = 0.6) +
  labs(
    x = "Points",
    fill = "Player", y = NULL
  ) +
  xlim(-1, 60) +
  scale_fill_viridis(discrete = TRUE, option = "D") +  # Paleta estética
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16)
  )





# al revés

sim_j1_m10pts <- posterior_predictive %>%
  filter(jugador == 1, puntos <= 10, j==1)

sim_j6_m10pts <- posterior_predictive %>%
  filter(jugador == 6, puntos <= 10, j == 1)


ggplot(data=sim_j6_m10pts, aes(x = minutos)) +
  geom_bar(aes(y = ..prop.., group = 1)) +
  labs(
    title = "Predictive distribution of minutes played",
    subtitle = "Conditioned on the Kyle Korver scoring fewer than 10 points",
    x = "Minutes", y = "Relative frequency"
  ) +
  theme_minimal()

library(ggplot2)
library(viridis)

ggplot(data = sim_j6_m10pts, aes(x = minutos)) +
  geom_bar(aes(y = ..prop.., group = 1), fill = viridis(1, option = "C")) +
  labs(
    x = "Minutes",
    y = "Relative frequency"
  ) +
  theme_minimal(base_size = 16) +   # Ajusta tamaño base
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16, face = "italic"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14)
  )



ggplot(data=sim_j1_m10pts, aes(x = minutos)) +
  geom_bar(aes(y = ..prop.., group = 1)) +
  labs(
    title = "Predictive distribution of minutes played",
    subtitle = "Conditioned on the Alen Iverson scoring fewer than 10 points",
    x = "Minutes", y = "Relative frequency"
  ) +
  theme_minimal()

ggplot(data = sim_j1_m10pts, aes(x = minutos)) +
  geom_bar(aes(y = ..prop.., group = 1), fill = viridis(1, option = "C")) +
  labs(
    x = "Minutes",
    y = "Relative frequency"
  ) +
  theme_minimal(base_size = 16) +   # tamaño base mayor
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16, face = "italic"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14)
  )


# Filter and label players
sim_combined <- posterior_predictive %>%
  filter(jugador %in% c(1, 6), puntos <= 10) %>%
  mutate(player_name = case_when(
    jugador == 1 ~ "Allen Iverson",
    jugador == 6 ~ "Kyle Korver"
  ))

# Combined plot with relative frequencies
ggplot(sim_combined, aes(x = minutos, fill = player_name)) +
  geom_bar(aes(y = ..prop.., group = player_name), position = "identity",
           alpha = 0.6) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Minutes", y = NULL,
    fill = "Player"
  ) +scale_fill_viridis(discrete = TRUE, option = "D") +  # Paleta estética
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16)
  )


# 3. Success probabilities estimation plots-----------
#Allen
library(ggplot2)

df1 <- data.frame(
  p = c(p1[1,,1], p2[1,,1], p3[1,,1]),
  shot = factor(
    rep(c("1PT", "2PT", "3PT"),
        each = length(p1[1,,1])),
    levels = c("1PT", "2PT", "3PT")
  )
)

ggplot(df1, aes(x = p, fill = shot)) +
  geom_density(alpha = 0.6, colour=NA) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 30)) +
  labs(
    title = NULL,
    x = "Probability",
    y = "Density",
    fill = "Shots scored"   # título de la leyenda
  ) +
  scale_fill_manual(
    name = "Shots scored",
    values = c(
      "1PT" = "#E69F00",  # naranja
      "2PT" = "#0072B2",  # azul
      "3PT" = "#009E73"   # verde
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16)
  )

#Kyle
library(ggplot2)

df6 <- data.frame(
  p = c(p1[6,,1], p2[6,,1], p3[6,,1]),
  shot = factor(
    rep(c("1PT", "2PT", "3PT"),
        each = length(p1[6,,1])),
    levels = c("1PT", "2PT", "3PT")
  )
)

ggplot(df6, aes(x = p, fill = shot)) +
  geom_density(alpha = 0.6, colour=NA) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 30)) +
  labs(
    title = NULL,
    x = "Probability",
    y = "Density",
    fill = "Shots scored"   # título de la leyenda
  ) +
  scale_fill_manual(
    name = "Shots scored",
    values = c(
      "1PT" = "#E69F00",  # naranja
      "2PT" = "#0072B2",  # azul
      "3PT" = "#009E73"   # verde
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16)
  )

