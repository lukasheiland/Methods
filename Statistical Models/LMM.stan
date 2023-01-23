data {
  int <lower=0> N;
  int <lower=0> N_groups;
  int <lower=1> N_beta;
  int <lower=1,upper=N_groups> group[N]; // this is an int array (vectors are always real)
  vector[N] y;
  vector[N] x;
}

parameters { // block defines sampling space
  real beta0;
  real beta1;
  real beta2;
  
  real<lower=0> sigma;
  
  vector[N_groups] r; // actually only this will get sanpled
  vector[N_beta] sigma_beta; // sigmas (or rather root of the covariance matrix?) regarding the random effect for level 'group'
    // all just a multiplied with r,
    // in order to receive completely correlated r_betas
}

transformed parameters {
  vector[N] y_hat;
  matrix[N_beta, N_groups] R_beta;
  
  for (g in 1:N_groups)
    R_beta[,g] = square(sigma_beta) * r[g];

  for (n in 1:N)
    y_hat[n] = beta0 + R_beta[1, group[n]] +
              (beta1 + R_beta[2, group[n]]) * x[n] +
              (beta2 + R_beta[3, group[n]]) * square(x[n]);
}

model {
  // random effect
  r ~ normal(0, 1);

  // Likelihood
  y ~ normal(y_hat, sigma);
}

generated quantities{
  matrix[N_groups, N_beta] Beta_r;
  for (g in 1:N_groups)
    Beta_r[,g] = [beta0, beta1, beta2]' + R_beta[,g]; // matrix of group x beta0-2
}

