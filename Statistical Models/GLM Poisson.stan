data {
  int<lower=0> N; 
  int<lower=0> y[N];
  vector[N] x;
}

parameters {
  real beta0; 
  real<lower=0> beta1;
}

model {
  // some prior
  beta1 ~ uniform(0, 10); 
  y ~ poisson_log(beta0 + beta1 * x);
    //  parameterization of the Poisson using the log rate α = log λ as a parameter.
    // equivalent to poisson(log(y_hat))
    // also see poisson_log_glm() which additionally accepts a model matrix
}

generated quantities {
  vector[2] beta_exp;
  // re-calculate intercept on original scale:
  beta_exp = [exp(beta0), exp(beta1)]'; //' transformws the row_vector to a vector
}
