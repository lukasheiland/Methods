library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE) # To avoid recompilation of unchanged Stan programs

# library(bayesplot)

# Simulation --------------------------------------------------------------
# y_i = a + (β_g)*x_i + ε_i
# with groups g and observations i

set.seed(1)
n <- 2000

## The data generating process
x <- runif(n, -3, 3)
X <- matrix(c(rep(1, n), x, x^2), nrow = n) # 
beta <- c(66, 33, -6.6)

y_hat <- X %*% beta
plot(y_hat ~ x)

## Random effects
n_groups <- 3
g <- sort(sample.int(n_groups, size = n, replace = T))
group <-  as.factor(g)

sigma_r <- 66 ### beautiful thing: only this one parameter is fit for all of the random effects!
Sigma_r <- sigma_r * diag(nrow = length(beta)) # covariance matrix as sigma_r * identity matrix
                                               # this way, with sigma on the diagonal, the structure generates homogeneous variance
                                               # 1 on diag for correlation,
                                               # varying sigmas for heterogeneity
# strucchange::root.matrix(Sigma_r)

# mehrdimensionaler zentraler Grenzwertsatz:
# Grenzwerte bestimmter Summen unabhängiger mehrdimensionaler Zufallsvariablen -> mvnorm
# Parameter: 1) vector of expected values, 2) how do the three distributions covary
Beta_r <- MASS::mvrnorm(n_groups,
                        mu = beta, # vector of expected values
                        Sigma_r) # analog to sd

mu <- sapply(1:n, FUN = function(i) X[i,] %*% Beta_r[g[i],]) 
plot(mu ~ x, col = g)

## The normal model
sigma <- 33
y <- rnorm(n, mu, sd = sigma) 
plot(y ~ x, col = g)


# MLE fit -----------------------------------------------------------------
library(lme4) # vs. nlme: lme4 more stable

# Each random effect could also be plotted as fixed, should even, if random effects would not have been invented ;)
# only problem: degrees of freedom are not really known, AIC selection problematic.

# Main difference: variable gets a bias towards 0 -> random effects are distributed tighter than fitted fixed effects.

fit_lm <- lm(y ~ x + I(x^2))
summary(fit_lm)
points(predict(fit_lm) ~ x, col = 5)

fit_lmm <- lmer(y ~ x + I(x^2) + (1 + x + I(x^2) | group)) # only beta[c(1,2)] are also random
points(predict(fit_lmm) ~ x, col = g)
points(mu ~ x, col = g)
res <-  y - predict(fit_lmm)
qqnorm(res)

summary(fit_lmm)
# what the summary returns as fixed effects …
fixef(fit_lmm) # are the parameter estimates for an average individual, similar (how different?) to calling lm
# the random effects section returns estimates for Std. deviation sigma (and Variance = sigma^2)
# for (1) between individuals (x), and (2) residual sd …
ranef(fit_lmm) # returns the contrasts to fixef for the groups
# absolute estimates:
matrix(rep(fixef(fit_lmm), 3), nrow = 3, byrow = T) + as.matrix(ranef(fit_lmm)$group)
Beta_r # pretty good


# brms ----------------------------------------------------------------
# library(brms)
# fit_rstanarm <- stan_lmer(y ~ x + I(x^2) + (1 + x + I(x^2)| group))
# rstan::get_stanmodel(fit_rstanarm$stanfit)

# data --------------------------------------------------------------------
# data for Stan needs to be named list.
data <- list(N = n,
             N_groups = n_groups,
             N_beta = length(beta),
             x = x,
             y = y,
             group = as.numeric(group))

# Stan model --------------------------------------------------------------------
## actually, in a Bayesian framework, one would rather do it hierarchically,
## with nested stuff, but here it is: the usual method recreated in Stan

fit_stan <- stan('Bayesian/LMM.stan',
                 data = data,
                 init = 'random',
                 warmup = 1000, # if not specified: automatically half of the iterations
                 iter = 1500,
                 cores = parallel::detectCores())
fit_stan # Rhat near 1: the chains have converged.
fixef(fit_lmm)

posterior <- extract(fit_stan)
plotpars <- c("beta0", "beta1", "beta2", "sigma", 'sigma_beta0', 'sigma_beta1', 'sigma_beta2')
traceplot(fit_stan, pars = plotpars)
pairs(fit_stan, pars = plotpars, las = 1) # Red points indicate divergent transitions.

# Plots -------------------------------------------------------------------
# parameter plots
plot(fit_stan,
     show_density = T,
     ci_level = 0.5,
     outer_level = 0.95,
     fill_color = 4,
     pars = plotpars[1:3])


plot(y ~ x, col = g)
for (j in 1:n_groups) {
        for (i in 1:500) { # draw a thousand random lines with parameters from the posterior
                curve((posterior$beta0[i] + posterior$R_beta[i, 1, j]) +
                              (posterior$beta1[i] + posterior$R_beta[i, 2, j]) * x +
                              (posterior$beta2[i] + posterior$R_beta[i, 3, j]) * x^2,
                      col = scales::alpha(j, 0.1),
                      lty = 1,
                      lwd = 0.3,
                      add = T)
                }
}
library(shinystan)
sso <- launch_shinystan(fit_stan)
