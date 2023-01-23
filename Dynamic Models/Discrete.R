############################################################################
# Discrete models ---------------------------------------------------------
############################################################################

## The problem with discrete models is often, that with time steps that aren't infinitesimally small, states can become negative. (C.f. Euler discretization example below)
## Solution: Model on the states, e.g. Ricker and Gompertz model.

## The process error only acts on the difference, is thus a function ot delta time.

## Ricker model -----------------------------------------------------------

## Note: ``log(a/b) == log(a) - log(b)``
## => The difference between states at two time steps: log(N_t1) - log(N_t1) = log(N_t1/N_t0)
##
## Expoential growth:
##    ```
##    log(N_t1) = log(N_t) + r + normal(0, s)
##    log(N_t1/N_t) = r + normal(0, s)
##    ```
##
## Logistic growth (Ricker):
##    ```
##    log(N_t1/N_t) = r_max - b*N_t + normal(0, s)
##    log(N_t1) = log(N_t) + r_max - b*N_t + normal(0, s)
##    ```
## - there is a constant linear decrease in the growth rate (r) as population size increases
## - Here, r_max is the growth rate when the population size is at its smallest possible value (i.e, 0).
## - note the second formulation can easily be estimated through linear regression with log.
##    - log state ~ time will yield a line


simRicker <- function(times,
                      initialstate, # A vector of species states.
                      pars) {


  ## Just unpacking for readability.
  r <- pars$r # Length n_species vector of input rates
  b <- pars$b # Length n_species vector of shading rates
  s <- pars$sigma_process

  ## Set the count variables
  times_intern <- 1:max(times)
  n_times <- length(times_intern)

  state <- rep(initialstate, times = n_times)
  state_log <- log(state)

  ## Here comes the model.
  for (t in 2:n_times) {
    state_log[t] <- state_log[t-1] + r - b*state[t-1] + rnorm(1, 0, s)
    ## or:  # state_log[t] <- state_log[t-1] + r - b*exp(state_log[t-1]) + rnorm(1, 0, s)
    state[t] <- exp(state_log[t]) ## necessary within loop for state[t-1]
  }

  return(state[which(times_intern %in% times)])
}


times <- 1:100
initialstate <- 0.5
pars <- list(r = 0.4, b = 0.002, sigma_process = 0.02)

y <- simRicker(times, initialstate, pars)
plot(y ~ times, type = "b")

increment <- y[2:length(y)] - y[1:(length(y)-1)]
plot(y[2:length(y)] ~ y[1:(length(y)-1)])
abline(0, 1)
plot(increment ~ y[1:(length(y)-1)])


## GOMPERTZ model -----------------------------------------------------------

## - Similar to the Ricker model, except that the constant linear decrease in growth rate (r) does not depend on the state but on log-state.
##
## - for states < 1, there is actually facilitation
## - with increasing state, density-effects grows slower than in the Ricker model.
## - This means that the density-limitation growth through the logstate. at small population sizes but as population size increases, the effect becomes less and less pronounced.
##
##    ```
##    log(N_t1) =  log(N_t) + r_max + b*log(N_t) + normal(0, s)
##    ```
## - Here, r_max is the growth rate when the population size equals 1
## - Interestingly, under the Gompertz model, the growth rate goes to infinity as the population size goes to 0
##    - log state ~ time will not yield a line

simGompertz <- function(times,
                      initialstate, # A vector of species states.
                      pars) {


  ## Just unpacking for readability.
  r <- pars$r # Length n_species vector of input rates
  b <- pars$b # Length n_species vector of shading rates
  s <- pars$sigma_process

  ## Set the count variables
  times_intern <- 1:max(times)
  n_times <- length(times_intern)

  ## there is only state_log here
  state_log <- log(rep(initialstate, times = n_times))

  ## Here comes the model.
  for (t in 2:n_times) {
    state_log[t] <- state_log[t-1] + r - b*state_log[t-1] + rnorm(1, 0, s)
  }

  state <- exp(state_log)
  return(state[which(times_intern %in% times)])
}


times <- 1:300
initialstate <- 0.5

## NOTE: Both r and b increased compared to Rickert
pars <- list(r = 0.2, b = 0.02, sigma_process = 0.02)

y <- simGompertz(times, initialstate, pars)

increment <- y[2:length(y)] - y[1:(length(y)-1)]
plot(y[2:length(y)] ~ y[1:(length(y)-1)])
abline(0, 1)
plot(increment ~ y[1:(length(y)-1)])



## NOTE: Euler-discretization of an ODE system -----------------------------

## from:
##    dN <- r - (m + c*N)*N

## to (with additional process error):
##    e_j <- rnorm(n, 0, s)
##    N_t1 <- N_t + (r - (m + c*N_t)*N_t + e_j) * dt[t-1]



# New ---------------------------------------------------------------------


simRicker <- function(times,
                      initialstate, # A vector of species states.
                      pars) {


  ## Just unpacking for readability.
  r <- pars$r # Length n_species vector of input rates
  b <- pars$b # Length n_species vector of shading rates
  s <- pars$sigma_process

  ## Set the count variables
  times_intern <- 1:max(times)
  n_times <- length(times_intern)

  state <- rep(initialstate, times = n_times)
  state_log <- log(state)

  ## Here comes the model.
  for (t in 2:n_times) {
    state_log[t] <- state_log[t-1] + r - b*state[t-1] + rnorm(1, 0, s)
    ## or:  # state_log[t] <- state_log[t-1] + r - b*exp(state_log[t-1]) + rnorm(1, 0, s)
    state[t] <- exp(state_log[t]) ## necessary within loop for state[t-1]
  }

  return(state[which(times_intern %in% times)])
}


times <- 1:100
initialstate <- 0.5
pars <- list(r = 0.4, b = 0.002, sigma_process = 0.02)

y <- simRicker(times, initialstate, pars)
plot(y ~ times, type = "b")

increment <- y[2:length(y)] - y[1:(length(y)-1)]
plot(y[2:length(y)] ~ y[1:(length(y)-1)])
abline(0, 1)
plot(increment ~ y[1:(length(y)-1)])


