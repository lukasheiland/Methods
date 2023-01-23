//  ~ indicates sampling
//
// optional blocks:
// “functions”
// “transformed data”
// “transformed parameters”
// “generated quantities”


// The input data is a vector 'y' of length 'N'.
// Declare variables, their data types, their dimensions, any restrictions (i.e. upper = or lower =)

data {
 int <lower = 1> N; // Sample size
 vector[N] x; // Predictor
 vector[N] y; // Outcome
 // real y[N]
}


// The parameters.
// datatype, dimensions, restrictions, name
// implicitly using uniform(-infinity, +infinity)

parameters {
 vector[3] beta;

 real <lower = 0> sigma; // Error SD
}

// The model.
//
// If no prior is defined, Stan uses default priors with the specifications uniform(-infinity, +infinity).

model {
  // prior
  beta[2] ~ normal(10, 20);
  
  y ~ normal(beta[1] + x * beta[2] + square(x) * beta[3], sigma);
}

// The posterior predictive distribution.
// block can be used to get info about posterior.
// randomly generate predicted values for each data point
// Typically, the data generating functions will be the distributions you used in the model block but with an _rng suffix (random number generator).

generated quantities {
 real y_rep[N];
 // real res[n]; // but don't because it can easily be subtracted in R post hoc
  y_rep = normal_rng(beta[1] + x * beta[2] + square(x) * beta[3], sigma); // generate randomly estimated posterior
  // res[n] = y[n] - y_rep[n];
}
