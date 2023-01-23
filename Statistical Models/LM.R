library(rstan)
library(bayesplot)

# Simulation --------------------------------------------------------------
n <- 100
x <- runif(n, -30, 30)
X <- matrix(c(rep(1, n), x, x^2), nrow = n) #
sigma <- 22
beta <-  c(333, 3.3, -0.33)
y_hat <- X %*% beta
y <- rnorm(n, y_hat , sd = sigma)

plot(y ~ x)
points(y_hat ~ x, col = 3)


# MLE fit -----------------------------------------------------------------
fit_lm <- lm(y ~ x + I(x^2))
summary(fit_lm)
beta_lm <- summary(fit_lm)$coeff[,1]  # intercept, slope, quadratic slope
sigma_lm <- sigma(fit_lm)

fitted_lm <- predict(fit_lm)
points(fitted_lm ~ x, col = 2, cex = 0.5)

res <-  y - y_hat
plot(res ~ fitted_lm)


# Data --------------------------------------------------------------------
# Data for Stan needs to be named list.
data <- list(N = n,
             x = x,
             y = y)

# Stan model --------------------------------------------------------------------
fit_stan <- stan('Bayesian/LM.stan',
                 data = data,
                 warmup = 200, # if not specified: automatically half of the iterations
                 iter = 1000,
                 cores = parallel::detectCores())
fit_stan # Rhat near 1: the chains have converged.
str(fit_stan)
summary(fit_stan)
posterior <- extract(fit_stan)


# Plots -------------------------------------------------------------------
# parameter plots
plot(fit_stan,
     show_density = T,
     ci_level = 0.5,
     outer_level = 0.95,
     fill_color = 'salmon',
     pars = c('beta', 'sigma'))

plot(density(posterior$beta[,1]), main = 'Beta0')
abline(v = beta_lm[1], col = 4, lty = 2)
abline(v = beta[1], col = 2, lty = 2)


plot(density(posterior$beta[,2]), main = 'Beta1')
abline(v = beta_lm[2], col = 4, lty = 2)
abline(v = beta[2], col = 2, lty = 2)

plot(density(posterior$sigma), main = 'Sigma')
abline(v = sigma_lm, col = 4, lty = 2)
abline(v = sigma, col = 2, lty = 2)


# plot lines
plot(y ~ x, pch = 20)
for (i in 1:1000) { # draw a thousand random lines with parameters from the posterior
  curve(posterior$beta[i,1] + posterior$beta[i,2] * x + posterior$beta[i,3] * x^2,
        col = scales::alpha('grey', 0.5),
        lty = 1,
        lwd = 0.2,
        add = T)
}

# traceplots

plot(posterior$sigma, type = 'l')
traceplot(fit_stan, pars = c('beta', 'sigma'))
stan_dens(fit_stan, pars = c('beta', 'sigma'))

# hypothesis tests
sum(posterior$beta[,1] > 330) / length(posterior$beta[,1]) # prob that beta1 > 0.2

# get the randomly generated predictions, random estimate posterior?
y_rep <- as.matrix(fit_stan, pars = 'y_rep') # each row is an iteration (single posterior estimate) from the model.
                                             # the columns correspond to the y[i]

error <- apply(y_rep[seq(from = 1, to = 1000, by = 10), ], MARGIN = 1,
               function(iteration) iteration - y) # Bayesian framework: simulated error distribution. (residuals would be difference with the estimate)
# wahrscheinlich macht man das mit dem Fehler anders

# bayesplot: posterior predictive checks
ppc_dens_overlay(y, y_rep[1:200, ])
#available_ppc()

# residuals
library(shinystan)
sso <- launch_shinystan(fit_stan)
