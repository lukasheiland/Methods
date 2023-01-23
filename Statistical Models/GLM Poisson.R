library(rstan)

# Simulate ----------------------------------------------------------------
n <- 1000
x <- runif(n, -1,1)
X <- cbind(1, x)
beta <- c(-1, 1.2)
y_hat <- X %*% beta
mu <- exp(y_hat) # poisson is only defined >= 0 => transform 'y_hat' to that range
y <-  rpois(n, mu)
plot(y ~ x, col = 'green')
points(mu ~ x, col = 'red')

# Bind data ---------------------------------------------------------------
data <- list(N = n,
             y = y,
             x = x)

# MLE fit -----------------------------------------------------------------
fit_mle <- glm(y ~ x, family = poisson)
model.matrix(fit_mle) == X

# Stan fit ----------------------------------------------------------------
fit_stan <- stan('Bayesian/GLM Poisson.stan',
                 data = data,
                 warmup = 200, # if not specified: automatically half of the iterations
                 iter = 500,
                 cores = parallel::detectCores())
posterior <- extract(fit_stan)

# scale the predictions back to linear
for (i in 1:500) { # draw the first thousand random lines with parameters from the posterior
  curve(exp(posterior$beta0[i] + posterior$beta1[i] * x),
        col = scales::alpha('red', 0.5),
        lty = 1,
        lwd = 0.2,
        add = T)
}

# scale parameters back to original linear
for (i in 1:1000) {
  curve(posterior$beta0_exp[i] + posterior$beta1_exp[i] * x,
        col = scales::alpha('green', 0.5),
        lty = 1,
        lwd = 0.2,
        add = T)
}

library(shinystan)
sso <- launch_shinystan(fit_stan)