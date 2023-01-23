// A generic implementation of the generalized Lotka-Volterra system.
functions {

  vector ds_dt(real time, vector state, vector r, matrix A, int n_spec) {

    vector[n_spec] ds = state .* (r + (A * state));

    return ds;
  }
}

data {
  int<lower=0> N_series;
  int<lower=0> N_obs; // not including t0, for now N_obs have to be the same in every series :(
  int<lower=0> N_species; // for now N_species have to be the same in every series :(

  real time_init[N_series];
  real times[N_series,N_obs]; // arrays are row major!
  real y_init[N_series,N_species];
  vector[N_species] Y[N_series,N_obs]; // two-dimensional array of vectors[N_species]
}

parameters {
  matrix<lower=0>[N_species, N_species] nA; // here A is assumed to be all competition, thus A negative. nA = -A (positive) exists only for being able to directly sample from an exponential prior
  vector<lower=0>[N_species] r;

  vector<lower=0>[N_species] state_init[N_series];
  real<lower=0> sigma; // vector<lower=0>[N_species] sigma;
}

transformed parameters {
  matrix[N_species, N_species] A = -nA;

  vector<lower=0>[N_species] State[N_series,N_obs];

  // integrate ode returns an array of column vectors, 'columns' accessed vie State[i]
  for(z in 1:N_series) {
      State[z,] = ode_rk45(ds_dt,
                           state_init[z], time_init[z], times[z,],
                           r, A, // model parameters
                           N_species);
  }

}


model {
  // Priors

  // for (s in 1:N_species) nA[, s] ~ cauchy(0, 1);
  r ~ lognormal(log(2), 0.1);

  // Regularizing prior for sigma
  // Already small over-estimation of sigma will lead to severely small parameter estimates.
  // The bias gets bigger with greater sigma relative to mu, the size of the median/values fitted to. But too steep prior will bias towards
  // sigma ~ cauchy(0, 0.01); // Is half-cauchy, not strong enough for lognormal??? Here, sigma is constrained to be > 0, thus half-cauchy
  sigma ~ exponential(1000);

  // Model
  for(s in 1:N_series) {
    y_init[s, ] ~ lognormal(log(state_init[s]), sigma); // Separately fitting initial state to feed to integrator.

    for (t in 1:N_obs) {
      Y[s,t] ~ lognormal(log(State[s,t]), sigma); // Lognormal model of states. // Note.: mean parameterization through mu = log(mean) - 0.5*sigma^2
    }
  }
}
