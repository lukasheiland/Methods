##########################################################################################
# Error in variables models  -------------------------------------------------------------------
##########################################################################################

# Simulate regression dilution -------------------------------------------------
n <- 500
x <- runif(n, -1, 1)
x_obs <- x + rnorm(n, sd = 2)
b <- 2
y <- b * x + rnorm(n) # The process on y happens before measurement error on x

data <- data.frame(y, x_obs)
m_true <- lm(y ~ x)
m2 <- lm(y ~ x_obs)

plot(x, y, col = "red")
abline(m_true, col = "red")
points(x_obs, y, col = "blue")
abline(m2, col = "blue")

##########################################################################################
# Methods  -------------------------------------------------------------------
##########################################################################################

# brms -------------------------------------------------
library(brms)

## include a fixed measurement error with me
fit <- brm(y ~ me(x_obs, 2), data = data)
summary(fit)

# stan -------------------------------------------------
library(rstan)

stancode <- "
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x_obs;
}

parameters {
  real beta0;
  real beta;
  real mu_x;
  real<lower=0> sigma_y;
  real<lower=0> sigma_x;
  real<lower=0> sigma_x_true;
  vector[N] x_true;
}

model {
  // priors
  sigma_y ~ normal(0, 5);
  sigma_x ~ normal(0, 5);
  global_sigma ~ normal(0, 5);

  // sampling x_true from a regularizing distribution
  x_true ~ normal(mu_x, sigma_x_true);

  // sampling x_true from a regularizing distribution
  x_obs ~ normal(x_true, sigma_x);

  // likelihood
  y ~ normal( x_true*beta + beta0, sigma_y);
}"

