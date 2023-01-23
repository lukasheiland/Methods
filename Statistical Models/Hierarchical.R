# A hierarchical model is a particular multilevel model where parameters are nested within one another.
# Apply if observations of one individual have the same value for one measure (e.g. treatment) but differen values for individual measures.
# This results in pseudor- eplication if we do not account for the hierarchical data structure in the model.
library(rstan)
library(shinystan)

# Simulation --------------------------------------------------------------
set.seed(1)
n <- 1000 #sample size
eggsize <- runif(n, -1, 1)
nesttemp <- runif(n, -1, 1)

n_groups <- 10 #number of groups
group <- rep(1:n_groups, each = n/n_groups)

#population-level regression coefficient
beta_pop <- c(intercept = 2, eggsize = -1, nesttemp = 0.5, nesttemp_sq = 0.2)
sigma_pop <- 0.9

#standard deviations of effects, used for random group-level coefficients
sd_beta <- c(3, 15, 0.7, 1.2)

# Group-level regression coefficients
# Draw groupwise random betas around the population wide parameters,
# so that nrow(Beta_groups) == n_groups.
# Random intercept and slopes!
Beta_groups <- mapply(function(c, sd) rnorm(n_groups, c, sd), c = beta_pop, sd = sd_beta)

X <- model.matrix( ~ eggsize + nesttemp + I(nesttemp^2),
                   data = data.frame(eggsize = eggsize, nesttemp = nesttemp))
y_hat <- sapply(1:n, function(n) X[n, ] %*% Beta_groups[group[n], ])
y <- rnorm(n, y_hat, sigma_pop)

plot(y ~ nesttemp, col = group)

# Bind data ----------------------------------------------------------------
data <- list(X = X,
             group = group,
             y = y,
             N_beta = length(beta_pop),
             N = n,
             N_groups = n_groups)

# Fit ---------------------------------------------------------------------
fit_lmer <- lme4::lmer(y ~ X[,2] + X[,3] + X[,4] + ( 1 + X[,2] + X[,3] + X[,4] | group),
           data = data)
summary(fit_lmer)

fit_stan <- stan('Hierarchical.stan',
                 data = data,
                 init = 'random',
                 warmup = 500, # if not specified: automatically half of the iterations
                 iter = 1000,
                 cores = parallel::detectCores())

bayesplot::mcmc_areas(as.matrix(fit_stan)[, 3:6])

launch_shinystan(fit_stan)
