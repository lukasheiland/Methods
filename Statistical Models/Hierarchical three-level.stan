
data {
  int<lower=0> N_groups1; // level 1
  int<lower=0> N_groups2; // level 2
  int<lower=0> N; // level 3

  // int<lower=0> id[N];
  // int<lower=1> groups1[N];
  // int<lower=1> groups2[N];
  
  int<lower=1> lookup1in2[N_groups2];
  int<lower=1> lookup2in3[N];

  vector[N_groups1] x1;
  vector[N_groups2] x2;
  vector[N] y;
  
  real hypersigma;
}

parameters {
  real b1;
  real b2;
  
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0> sigma3;
  
  real m0;
  vector[N_groups1] m1;
  vector[N_groups2] m2;
}

transformed parameters {
  vector[N_groups1] m2_hat;
  vector[N_groups2] m3_hat;
  
  m2_hat = m1 + b1*x1;
  m3_hat = m2 + b2*x2;
}

model {
  // Priors
  sigma3 ~ cauchy(0, hypersigma); // cauchy(0, 1)
  
  // 'Random effects'
  m1 ~ normal(m0, sigma1);
  m2 ~ normal(m2_hat[lookup1in2], sigma2);
  
  // Likelihood
  y ~ normal(m3_hat[lookup2in3], sigma3);
}

