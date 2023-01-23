// for now, this is just a hierarchical linear model

data {
  int N;
  int N_groups; // the number of groups
  int N_beta; // number of columns in the model matrix
  int group[N]; // array of group ids
  matrix[N, N_beta] X; //the model matrix

  vector[N] y; //the response
}

parameters {
  vector[N_beta] beta; //population-level regression coefficients
  vector<lower=0>[N_beta] sigma_groups; //the standard deviation for the random regression coefficients
  matrix[N_beta, N_groups] Beta; //matrix of group-level regression coefficients
  real<lower=0> sigma; //standard deviation of the individual observations
}

transformed parameters{

}

model {
  vector[N] mu; //linear predictor
  //priors
  // beta ~ normal(0,5); //weakly informative priors on the regression coefficients
  //sigma_groups ~ cauchy(0,2.5); //weakly informative priors, see section 6.9 in STAN user guide
  //sigma ~ gamma(2,0.1); //weakly informative priors, see section 6.9 in STAN user guide

  for(g in 1:N_groups){
   Beta[, g] ~ normal(beta, sigma_groups); //fill the matrix of group-level regression coefficients
                                           // this happens groupwise, but always within a specific sigma around beta.
  }

  for(n in 1:N){
    mu[n] =  X[n,] * Beta[, group[n]]; //compute the linear predictor using relevant group-level regression coefficients
  }

  //likelihood
  y ~ normal(mu, sigma);
}
