data {
  int <lower=0> N;
  int <lower=0> N_effects;
  vector[N] y;
  matrix[N, N_effects] X;
}

parameters {
  real <lower=0> tau;
  real <lower=0> inverselambda;
  real intercept;
  vector[N_effects] effects;
}

transformed parameters {
  vector[N] y_hat;

  y_hat = intercept + X * effects; // there is no penalty on the intercept
}

// Tibshirani (1996) suggested that Lasso estimates can be interpreted as posterior mode estimates
// when the regression parameters have independent and identical Laplace (i.e., double-exponential) priors.
// In fact, when you place a Laplace prior over the parameters, the MAP solution should be identical (not merely similar) to regularization with the L1 penalty
// and the Laplace prior will produce an identical shrinkage effect to the L1 penalty.
// Due to numerical issues, solutions may not actually be identical.


model {
  // priors
  // intercept ~ normal(0, 0.0001); // i guess this is just for faster convergence
  tau ~ double_exponential(0, 0.001);
  inverselambda ~ double_exponential(0, 0.001); // double_exponential() == Laplace(mu, sd) this is extremely tight around 0
  effects ~ normal(0, inverselambda); // 1/sqrt(lambda) = sd!; effects are priored to zero, and inversedlambda will be sampled!

  // likelihood
  y ~ normal(y_hat, tau);
}

generated quantities{
  real<lower=0> sigma;

  sigma = 1/sqrt(tau); // why?
}

