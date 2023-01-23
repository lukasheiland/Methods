##########################################################################################
# Multinomial model ------------------------------------------------------------
##########################################################################################

# Library -----------------------------------------------------------------
library(rstan)

# Simulation --------------------------------------------------------------
set.seed(1)
n <- 100 #sample size
n_species <- 10

Data_x <- data.frame(ph = runif(n, -1, 1), moist = runif(n, -1, 1))
X <- model.matrix(~ ph + moist, data = Data_x) # adds an intercept
Beta <- cbind(intercept = rnorm(n_species, 0.7),
              ph = -1:8,
              moist = -7:2) # rows are species!

Y_linear <- apply(Beta, MARGIN = 1, function(beta) X %*% beta)
Y_hat <- plogis(Y_linear) # columns are species
Y <- t(apply(Y_hat, MARGIN = 1, function(prob) rmultinom(1, size = 1000, prob = prob)))



# Bind data and plot ----------------------------------------------------------------
Data <- data.frame(Data_x, y = Y)

library(ggplot2)
Data_long <- tidyr::gather(Data, key = 'species', value = 'y', -ph, -moist)
ggplot(Data_long, aes(x = ph, y = y, fill = species)) +
  geom_area()

# Fit ---------------------------------------------------------------------
library(nnet)
fit_nnet <- multinom(Y ~ ph + moist, data = Data) # Y is a matrix with K columns
summary(fit_nnet)
summary(fit_nnet)


fit_stan <- stan('Hierarchical.stan',
                 data = data,
                 init = 'random',
                 warmup = 500, # if not specified: automatically half of the iterations
                 iter = 1000,
                 cores = parallel::detectCores())

bayesplot::mcmc_areas(as.matrix(fit_stan)[, 3:6])
launch_shinystan(fit_stan)
