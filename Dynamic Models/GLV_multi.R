##########################################################################################
# Generalized Lotka-Volterra model ---------------------------------------
##########################################################################################

# Fit to multiple time series of equally configured pops.

# Notes -------------------------------------------------------------------
# The competition variant has logistic growth! (exponential in predator-prey model)
# There are solutions for one equilibrium, multiple equilibria, no equilibria

## The generalised Lotka-Volterra equations can represent both competition (a_ij and aji negative) and predation (one of a_ij and a_ji negative), depending on the values of the parameters, as described below.

## r: intrinsic (at absence of others) birth or death rates of the species.
## A positive value for r means that a species is able to reproduce in the absence of any other species, whereas a negative value means that its population will decline unless the appropriate other species are present (e.g. a herbivore that cannot survive without plants to eat, or a predator that cannot persist without its prey).

## A: The values of matrix A represent the relationships between the species.
## The value of <math>a_{ij}</math> represents the effect that species j has upon species i.  The effect is proportional to the populations of both species, as well as to the value of <math>a_{ij}</math>.  Thus, if both <math>a_{ij}</math> and <math>a_{ji}</math> are negative then the two species are said to be in direct competition with one another, since they each have a direct negative effect on the other's population.  If <math>a_{ij}</math> is positive but <math>a_{ji}</math> is negative then species i is considered to be a predator (or parasite) on species j, since i's population grows at j's expense.
## Positive values for both <math>a_{ij}</math> and <math>a_{ji}</math> would be considered mutualism.  However, this is not often used in practice, because it can make it possible for both species' populations to grow indefinitely.
## Indirect negative and positive effects are also possible.  For example, if two predators eat the same prey then they compete indirectly, even though they might not have a direct competition term in the community matrix.
## The DIAGONAL terms <math>a_{ii}</math> are usually taken to be negative (i.e. species i's population has a negative effect on itself).  This self-limitation prevents populations from growing indefinitely.

# Library -----------------------------------------------------------------
library(deSolve)
library(tidyverse)

# Multi-species matrix model -------------------------------------------------------------------

## In the COMPETITION MODEL, the system is expressed by equations, here for one species in a two species system
## dm <-  m * r (1 - (A %*% m)) # matrix version

## dn/dt ==  N r_n * (1 - (a_nn N + a_nm M) # two species multiplication

## - There, the growth rate r_n N is the observed growth rate.
## - And it is logistically limited by the added effects of inter- and intraspecific competition/mutualism.
## - a act multiplicatively on r, are all competition terms (intra/inter!): negative mutualism, positive competition
## - Species compete for the same "resource"; model is asymptotically stable when all but one species are extinct

## In G-LV, whose parameters have nothing to do with the parameters above (here for one species in a two species system)
## dm <- m * (r + (A %*% m))
## which translates for one species in a two species system to
## dN/dt == N * (rn_intr + (a_nn * N + a_nm * M))
## - Here, the growth rate r_n N is the intrinsic (at absence of others) growth (or decline) rate. "Fundamental".
## - The sum of interactions are are added. Thus, at constant competitors, this would just be exponential growth.
## - a act additively on r: both negative competition, both positive mutualism, negative and positive: predator-prey
## - The limitation of r is not logistic, AFAIS.
## - There are coexistence equilibria!

calcGLV <- function(t,
                    m, # A vector of species states.
                    par){
  r <- par[[1]] # Length n vector of growth rates.
  A <- par[[2]] # Full n x n matrix of competition factors.
  ## This would be the competition model
  # dm <- r * m * (1 - (A %*% m)) # Here, a logistic limitation is built in: the bigger m gets
  dm <- m * (r + (A %*% m)) # Here, limitation emerges only from the diagonal.
  return(list(c(dm)))
}



# Simulation -------------------------------------------------------------------
set.seed(1)

#### Set parameters.
## Generate some random parameters.
n_species <- 3
r <- rlnorm(n_species, log(2), 0.1) # Vector of intrinsic growth rates.

A <- -matrix(rlnorm(n_species^2, log(0.2), 0.1), nrow = n_species) # A full matrix of mutual competition factors.
diag(A) <- -rlnorm(n_species, log(0.3), 0.01)  # Intraspecific interaction.

par <- list(r, A) # Parameters list, including a matrix of alpha values.

#### Integrate model for one sim
simulateSeries <- function(times = seq(4, 40, by = 2), sigma = 0) {
  m0 <- runif(n_species, 0.5, 4) # Initial state matrix.
  Sim <- ode(m0, times, calcGLV, par, method = 'ode45')
  Sim[,2:(n_species+1)] <- matrix(rlnorm(Sim, log(Sim), sigma), nrow = nrow(Sim))[, 2:(n_species+1)]

  ## only necessary for rnorm
  # Sim[, 2:(n_species+1)][Sim[, 2:(n_species+1)] < 0] <- 0

  return(Sim)
}

Sim <- simulateSeries(sigma = 0.05)
matplot(Sim[, 1], Sim[, -1], type="l", ylab="N") # log='y'



#### Multiple simulations
simulateMultipleSeries <- function(n_series = 15, n_times = 20, sigma = 0.05, format = c("long", "wide", "list")) {
  Sims <- replicate(n_series,
                    simulateSeries(1:(n_times), sigma = sigma),
                    simplify = F)

  if (match.arg(format) %in% c("wide", "long")) {
    Sims <- cbind(do.call(rbind, Sims), series = rep(1:n_series, each = n_times))
  }

  if (match.arg(format) == "long") {
    Sims <- tidyr::pivot_longer(as.data.frame(Sims),
                                cols = all_of(paste(1:n_species)),
                                names_to = "pop",
                                values_to = "abundance") %>%
      mutate(species = rep(1:n_species, length.out = nrow(.)))
  }
  return(Sims)
}


#### Sim of multiple series, just for visual exploration
S <- simulateMultipleSeries(sigma = 0.05)

ggplot(S, mapping = aes(x = time,
                        y = abundance,
                        color = species,
                        group = interaction(species, series))) +
  geom_line() +
  facet_wrap(facets = c("series"))



# Stan model --------------------------------------------------------------
library(cmdstanr)
# install_cmdstan(cores = 3)

## Get a list of lists of data, where each second-level list belongs to one time series
## The lists will be unpacked as arrays of doubles, matrices, vectors etc. in stan
getStanData <- function(sims){
  stanlist <- list(
    N_series = length(sims),

    ## For now, n_obs, and n_species have to be the same over all series. Thus, only first element extraction
    N_obs = sapply(sims, function(s) nrow(s)-1)[[1]], # observations, other than at t_0
    N_species = sapply(sims, function(s) ncol(s) - 1)[[1]],

    time_init = sapply(sims, function(s) s[1, "time"]), # stan expects a one-dimensional array
    times = t(sapply(sims, function(s) s[2:nrow(s), "time"])), # stan expects a two-dimensional array real times[N_series,N_obs]; // arrays are row major!
    y_init = t(sapply(sims, function(s) s[1, -1])), # stan expects two-dimensional array: real y_init[N_series,N_species];
    Y = aperm(sapply(sims, function(s) s[-1, -1], simplify = "array"), c(3, 1, 2)) # stan expects two-dimensional array of vectors: vector[N_species] Y[N_series,N_obs] ## provided indexing order is opposite in R and stan!
  )
  return(stanlist)
}

## Simulate multiple time series in list format for stan fit
sims <- simulateMultipleSeries(format = "list", sigma = 0.05)
data <- getStanData(sims)

getInits <- function() {
  list(
    nA = matrix(rlnorm(data$N_species^2, log(0.2), 0.01), nrow = data$N_species, ncol = data$N_species), # very similar values

    ## different intrinsic growth rates
    r = rnorm(n_species, 2, 0.01),

    state_init_hat = data$y_init, # stan expects vector<lower=0>[N_species] state_init[N_series]
    sigma = 1 # rep(0.1, data$N_species) # Big error at the start.
  )
}


#### Compile model
model <- cmdstan_model("Dynamic/GLV_multi.stan")
# model$exe_file() ## path to compiled model


#### Optimize model.
fit_optim <- model$optimize(data = data,
                            output_dir = "Dynamic/Fits",
                            init = getInits,
                            iter = 1000)

fit_optim$summary() %>% View
fit_optim$summary(variables = 'A')
A

# fit_var <- model$variational(data = data,
#                              output_dir = "Dynamic/Fits",
#                              init = getInits)

#### Sample.
n_chains <- 3
fit <- model$sample(data = data,
                    output_dir = "Dynamic/Fits",
                    init = getInits,
                    iter_warmup = 200, iter_sampling = 500,
                    chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))

fit$summary(variables = 'state_init')
data$y_init

# fit$output()
# fit <- rstan::read_stan_csv(glue::glue("Dynamic/Fits/GLV_multi-202010151422-{1:3}-efa638.csv"))

fit_rstan <- rstan::read_stan_csv(fit$output_files()) ## This is the easiest workaround to get an rstan object
fit_rstan
shinystan::launch_shinystan(fit_rstan)
# saveRDS(fit_rstan, "Dynamic/fit_GLV.rds")

# gauseR --------------------------------------------------------------------
## Some strategies with autocorrelation are explained with library(gauseR)


