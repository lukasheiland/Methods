##########################################################################################
# Regularization ---------------------------------------------------------------
##########################################################################################

# L1 (lasso), i.e. exp prior, leads to either 0 or something different (massive weight around 0, but tapers off less strictly in the tails)
# L2 (ridge), i.e. normal prior, leads something (like sd) closeer to 0 (mass broader around 0, but strictöy tapering tails)

# after shrinkage, there will be no p-value!

# Library -----------------------------------------------------------------

library(glmnet)
library(rstan)
library(shinystan)
library(rjags)
library(bayesplot)


# Simulate data with xes with and without effect --------------------------
n <- 100
n_effects <-  30
X <- matrix(rnorm(n * n_effects), n, n_effects)
effects <- c(runif(n_effects - 20, -10, 10), rep(0, 20)) # only ten real effects, rest 0
y_hat <- X %*% effects
y <- rnorm(n, y_hat)


# glmnet ---------------------------------------------------------
# vignette: https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html

fit_lm <- lm(y ~ X)
summary(fit_lm)

fit_lasso <- glmnet(X, y, alpha = 0.8) # alpha = 0: lasso, alpha = 1: ridge
  # penalty: lambda*abs(beta + beta), intercept wont get ridge regression penalty
  # It is known that the ridge penalty shrinks the coefficients of correlated predictors towards each other while the lasso tends to pick one of them and discard the others.

  # The glmnet algorithms use cyclical coordinate descent, which successively optimizes the objective function over each parameter with others fixed, and cycles repeatedly until convergence.

plot(fit_lasso, label = T) # A coefficient profile plot against the L1-norm
  # Each curve corresponds to a variable.
  # It shows the path of its coefficient against the l_1-norm of the whole coefficient vector at as λ varies.
  # The axis above indicates the number of nonzero coefficients at the current λ, which is the effective degrees of freedom (df) for the lasso.
summary(fit_lasso)

coef(fit_lasso)
coef_lasso <- coef(fit_lasso, s = 1) # returns coefs of the sequence at a certain penalty lambda
names(coef_lasso[which(abs(coef_lasso) > 0),])

glmnet_cv <- cv.glmnet(X, y) # cross validation
plot(glmnet_cv) # x-axis: 'rubber band strength'
  # marked in the plot: two lambdas
  # lambda.min is the value of λ that gives minimum mean cross-validated error.
  # lambda.1se, is the most regularized model such that error is within one standard error of the minimum.

## Btw., there is no good unbiased estimator which can cope with colineratiy
# but biased ones, such as in lasso regression.


# Data --------------------------------------------------------------------
data <- list(N = n,
             N_effects = n_effects,
             X = X,
             y = y)

# Stan model --------------------------------------------------------------------
fit_stan <- stan('Bayesian/Regularization.stan',
                 data = data,
                 init = 'random',
                 warmup = 500, # if not specified: automatically half of the iterations
                 iter = 900,
                 cores = parallel::detectCores())
fit_stan

posterior <- extract(fit_stan)
plot(fit_stan,
     show_density = T,
     ci_level = 0.5,
     outer_level = 0.95,
     fill_color = 4)


# Fit JAGS model ----------------------------------------------------------
jagscode <- "model{
# Likelihood
for(i in 1:N){
  mu[i] <- b + inprod(X[i,], a) # this
  y[i] ~ dnorm(mu[i], tau)
}

# Prior distributions

for(i in 1:N_effects){
  a[i] ~ dnorm(0, lambda)
}
b ~ dnorm(0, 0.0001)

tau ~ dgamma(0.001, 0.001)
lambda ~ dgamma(0.001, 0.001)
sigma <- 1/sqrt(tau) # this line is optional, just in case you want to observe sigma or set sigma (e.g. for inits)

}
"

jagsmodel <- jags.model(file = textConnection(jagscode), data = data, n.chains = 3)

parnames <- c("a","b","sigma")
samples_jags <- coda.samples(jagsmodel, variable.names = parnames, n.iter = 10^4)
launch_shinystan(as.shinystan(samples_jags))
summary(samples_jags)
