##########################################################################################
# Generalized Lotka-Volterra model ---------------------------------------
##########################################################################################

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
set.seed(5)
n_species <- 3
time <- 0:100

## Generate some random parameters.
r <- rlnorm(n_species, log(1), 0.01) # Vector of intrinsic growth rates.

A <- -matrix(rlnorm(n_species^2, log(0.2), 0.01), nrow = n_species) # A full matrix of mutual competition factors.
# diag(A) < runif(length(diag(A)), 0, 0.2)  # Intraspecific interaction.

par <- list(r, A) # Parameters list, including a matrix of alpha values.

m_0 <- runif(n_species, 1, 2) # Initial state matrix.

Sim_m <- ode(m_0, time, calcGLV, par)
matplot(Sim_m[, 1], Sim_m[, -1], type="l", ylab="N") # log='y'

# add error
sigma <- 0.05
Sim_m[,(1:n_species)+1] <- matrix(rlnorm(Sim_m, log(Sim_m), sigma), nrow = nrow(Sim_m))[,(1:n_species)+1]
matplot(Sim_m[,-1], type = "l")

# Stan model --------------------------------------------------------------

library(cmdstanr)
# install_cmdstan(cores = 3, overwrite = T)

data <- list(
  N = nrow(Sim_m) - 1,
  N_species = ncol(Sim_m) - 1,

  time_init = Sim_m[1,1],
  times = Sim_m[-1,1],
  y_init = Sim_m[1,-1],
  Y = as.data.frame(Sim_m[-1,-1]) # Stan expects an array of vectors, i.e. list of vectors.
)

any(data$Y < 0)


getInits <- function() {
  list(
    ##
    nA = matrix(rlnorm(data$N_species, log(0.2), 0.001), nrow = data$N_species, ncol = data$N_species),

    ## different intrinsic growth rates
    r = rlnorm(n_species, log(1), 0.001),
    state_init = data$y_init - 0.5*(0.1)^2,
    sigma = 1 # rep(0.1, data$N_species)
    )
}


model <- cmdstan_model("Dynamic/GLV.stan")
# model$exe_file() ## path to compiled model

#### Optimize.
fit_optim <- model$optimize(data = data,
                            output_dir = "Dynamic/Fits",
                            init = getInits)
fit_optim$summary() %>% View()

#### Sample.
n_chains <- 3
fit <- model$sample(data = data,
                    init = getInits,
                    output_dir = "Dynamic/Fits",
                    iter_sampling = 1000, iter_warmup = 500,
                    chains = n_chains, parallel_chains = getOption("mc.cores", n_chains))

# fit_rstan <- rstan::read_stan_csv(fit$output_files())
# fit_rstan <- rstan::read_stan_csv(glue::glue("Dynamic/Fits/GLV-202010142133-{1:3}-918672.csv"))

fit_rstan <- rstan::read_stan_csv(fit$output_files()) ## This is the easiest workaround to get an rstan object
fit_rstan
saveRDS(fit_rstan, "Dynamic/fit_GLV.rds")
shinystan::launch_shinystan(fit_rstan)


# gauseR --------------------------------------------------------------------
## Some strategies with autocorrelation are explained with library(gauseR)


