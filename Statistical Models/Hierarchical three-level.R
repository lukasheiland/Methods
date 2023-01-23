library(tensorflow)
# install_tensorflow()
library(tfprobability)
# install_tfprobability()
library(cmdstanr)
# install_cmdstan()

# Simulate data from three-level hierarchical process ---------------------------------
## model structure with two nested groups: groups1/groups2/observations
set.seed(100)

n_groups1 <- 15
n_groupswithin1 <- 10 # groups2
n_groups2 <- n_groups1 * n_groupswithin1
n_obswithin2 <- 5
n <- n_groups1 * n_groupswithin1 * n_obswithin2

groups1 <- rep(1:n_groups1, each = n_groupswithin1*n_obswithin2)
groups2 <- rep(1:n_groups2, each = n_obswithin2)
id <- 1:n

G <- cbind(groups1, groups2, id)
  
x1 <- runif(n_groups1, -1, 1)
x2 <- runif(n_groups2, -1, 1)

## level 0
m0 <- 3 # overall mean

#                      m0
# (sd = 8)
# m1                              m1              
# (sd = 4)
# m2               m2             m2             m2
# (sd = 2)
# obs obs obs     obs obs obs     obs obs obs    obs obs obs

## level 1, dependent on x1
m1 <- rnorm(n_groups1, mean = m0, sd = 8)
b1 <- 1
m1 <- m1 + b1*x1
m1 <- rep(m1, each = n_groupswithin1) # length == n_groups2

## level 2, independent across groups1
m2 <- rnorm(n_groups2, mean = m1, sd = 4)
b2 <- -1
m2 <- m2 + b2 * x2
m2 <- rep(m2, each = n_obswithin2) # length == n

## level 3
y <- rnorm(n, mean = m2, sd = 2) ## within group2 error

D <- data.frame(G, y, x1 = x1[groups1], x2 = x2[groups2])
boxplot(y ~ groups1:groups2, data = D[1:175,], col = groups1)

# Stan fit ----------------------------------------------------------------
standata <- list(y = D$y,
                 x2 = x2,
                 x1 = x1,
                 
                 N_groups1 = length(unique(D$groups1)),
                 N_groups2 = length(unique(D$groups2)),
                 N = nrow(D),
                 
                 lookup2in3 = D$groups2,
                 lookup1in2 = unique(D[c('groups2', 'groups1')])$groups1,
                 
                 hypersigma = 1
                 )

## Compilation
stanmodel <- cmdstan_model('Hierarchical/Hierarchical.stan')

#### Sampling/Optimization
## 0. Optimization
stanfit_optim <- stanmodel$optimize(data = standata)
stanfit_optim$summary() # %>% View()

## 1. NUTS sampling
stanfit_nuts <- stanmodel$sample(data = standata,
                                chains = 3,
                                parallel_chains = getOption("mc.cores", 3))



stanfit_nuts$summary(variables = c("m0", "b1", "b2", "sigma1", "sigma2", "sigma3"))
draw <- stanfit_nuts$draws()
bayesplot::mcmc_trace(draw, pars = c("m0", "b1", "b2", "sigma1", "sigma2", "sigma3"))
bayesplot::mcmc_areas(draw, pars = c("m0", "b1", "b2", "sigma1", "sigma2", "sigma3"))
bayesplot::ppc_dens_overlay(y = unique(m2),
                 yrep = posterior::as_draws_matrix(stanfit_nuts$draws(variables = "m3_hat"))[1:500,],
                 alpha = 0.2)

## 2. Variational Bayes
stanfit_variational <- stanmodel$variational(data = standata)
bayesplot::mcmc_areas(stanfit_variational$draws(), pars = c("m0", "b1", "b2", "sigma1", "sigma2", "sigma3"))


## 3. TMB with Laplace sampling
# library(tmbstan)


# stanfit <- rstan::read_stan_csv(stanfit_nuts$output_files())
# shinystan::launch_shinystan(stanfit)


# mle fit -------------------------------------------------------------
fit <- lme4::lmer(y ~ x2 + x1 + (1 | groups1/groups2), data = D)
summary(fit)

# brms

