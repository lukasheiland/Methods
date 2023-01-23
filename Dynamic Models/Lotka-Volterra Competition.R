##########################################################################################
# Lotka-Volterra model for competition system  ---------------------------------------
##########################################################################################

# Notes -------------------------------------------------------------------
## The competition variant has logistic growth! (exponential in predator-prey model)
## There are solutions for one equilibrium, multiple equilibria, no equilibria
## - Species compete for the same "resource"; model is asymptotically stable when all but one species are extinct


# Library -----------------------------------------------------------------
library(deSolve)


# Two-species model -------------------------------------------------------------------
## states N, M of two species
## dn/dt ==  r_n N * (1 - a_nn N - a_nm M)
## increment of rate r_n*N is at max for both very small intraspecific and interspecific competition
## a12 is the competition of M on N
## as states N or M grow bigger, r_n*N approaches 0

calcComp2 <- function (time, state, par)
{
  with(as.list(par), {
    dn <- r_n * state[1] * (1 - a11 * state[1] - a12 * state[2])
    dm <- r_m * state[2] * (1 - a22 * state[2] - a21 * state[1])
    list(c(dn, dm))
  })
}


# Two-species simulation --------------------------------------------------------------
par2 <- c(r_n = 0.05,
         r_m = 0.05,
         a11 = 0.2,
         a21 = 0.2,
         a22 = 0.2,
         a12 = 0.1)
state_init <- c(0.1, 1)
time = 1:1000

Sim_2 <- ode(y = state_init, times = time, func = calcComp2, parms = par2)
matplot(Sim_2[,-1], type = "l", xlab = "time", ylab = "population")
legend("topright", c("Phragmites australis", "Valeriana dioica"), lty = c(1,2), col = c(1,2), box.lwd = 0)


# Multi-species matrix model -------------------------------------------------------------------

## For one species the model equation is:
## dn/dt ==  r_n N * (1 - (a_nn N + a_nm M))

## Same for multiple species:
## dm/dt ==  r m * (1 - (A %*% m))
## first term (r m): here, r, m are just vectorized
## second term (1 - (A %*% m)): for a two species model equivalent to
## c(A_11 m1 + A_12 m2, A_21 m1 + A_22 m2)

calcCompM <- function(t,
                      m, # A vector of species states.
                      par){
  r <- par[[1]] # Length n vector of growth rates.
  A <- par[[2]] # Full n x n matrix of competition factors.
  dm <- r * m * (1 - (A %*% m)) # This is logistic as well.
  return(list(c(dm)))
}

# Multi-species simulation -------------------------------------------------------------------
n_species <- 6
time <- 1:2000

## Generate some random parameters.
r <- runif(n_species)*2 # Vector of growth rates.

alpha <- .01 # The mean competition.
A <- matrix(rnorm(n_species^2, alpha, alpha/10), nrow = n_species) # A full matrix of mutual competition factors.

par <- list(r, A) # Parameters list, including a matrix of alpha values.

m0 <- runif(n_species)/(n_species*alpha) # Initial state matrix.

Sim_m <- ode(m0, time, calcCompM, par)
matplot(time, Sim_m[,-1], type="l", ylab="N") # log='y'



# Real data Bayesian tools fit -----------------------------------------------------------

## Questionable if it works, see if rlnorm parameterization is right.
# library(BayesianTools)
#
# ## Data
# library(dave) # Package from WILDI text book.
# data(lveg)
# ?lveg
# ?ltim
#
# Veg <- cbind(t = ltim, lveg)
# v1963 <- unlist(Veg[1, -1])
# V <- as.matrix(Veg[-1,-1])
#
# n_spec <- length(v1963)
# time <- Veg$Year - 1963
#
# n_par <- n_spec + n_spec^2 + 1
#
# ## Likelihood
# likelihood <- function(p){
#
#   r <- p[1:n_spec]
#   A <- matrix(p[(n_spec+1):(n_par-1)], nrow = n_spec)
#   A <- matrix(rnorm(A, A, 3), nrow = n_spec) # Shrinkage
#   sigma <- p[n_par]
#
#   par <- list(r, A)
#
#   Sim_hat <- ode(v1963, time, calcCompM, par)[,-1] # Just dropping the first column of times.
#   Sim <- rnorm(Sim_hat, Sim_hat, sigma) # Stochasticity.
#
#   L <- dlnorm(V, mean = Sim, sd = 1, log = T)
#   return(sum(L, na.rm = T))
# }
#
#
# ## BayesianTools fit
# setup <- createBayesianSetup(likelihood,
#                              lower = c(rep(-1, n_par-1), 0), # The last par is the error sd.
#                              upper = c(rep(2, n_par-1), 5),
#                              parallel = T)
# runsettings <- list(iterations = 10^4)
# draws <- runMCMC(setup, settings = runsettings)
# saveRDS(draws, "Bayesian/Lotka-Volterra-MCMC.rds")
# # draws <- readRDS("Bayesian/Lotka-Volterra-MCMC.rds")
#
# summary(draws)
# plot(draws, start = 0)
#
# mle <- DEoptim::DEoptim(function(x) -likelihood(x),
#                         lower = c(rep(-1, n_par-1), 0),
#                         upper = c(rep(2, n_par-1), 5))




# gauseR --------------------------------------------------------------------
library(gauseR)

