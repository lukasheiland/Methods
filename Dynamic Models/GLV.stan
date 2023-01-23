// A generic implementation of the generalized Lotka-Volterra system.
functions {

  vector ds_dt(real time, vector state, vector r, matrix A, int n_spec) {

    vector[n_spec] ds = state .* (r + (A * state));

    return ds;
  }
}

data {
  int<lower=0> N;
  int<lower=0> N_species;

  real time_init;
  real<lower=0> times[N];
  real y_init[N_species];
  vector<lower=0>[N_species] Y[N];
}

parameters {
  matrix<lower=0>[N_species, N_species] nA; // here A is assumed to be all competition, thus A negative and nA positive
  vector<lower=0>[N_species] r;

  vector<lower=0>[N_species] state_init;
  real<lower=0> sigma; // vector<lower=0>[N_species] sigma;
}

transformed parameters {
  matrix<upper=0>[N_species, N_species] A = -nA;
  vector<lower=0>[N_species] State[N];

  // integrate ode returns an array of column vectors, 'columns' accessed vie State[i]
  State = ode_rk45_tol(ds_dt,
                       state_init, time_init, times,
                       1e-5, 1e-3, 10000, // old integrator tolerances settings: real rel_tol, real abs_tol, int max_num_steps
                       r, A, // model parameters
                       N_species);

}


model {
  // Priors

  // slightly informative
  // for (s in 1:N_species) nA[, s] ~ exponential(10); // negative decaying prior, because all-competitive is assumed
  // r ~ lognormal(0, 0.1); // positive prior around 1; 0 = log(1)

  // Regularizing: half-cauchy isn't enough:
  // sigma ~ cauchy(0, 0.1); // sigma is constrained to be > 0, thus half-cauchy
  sigma ~ exponential(1000); //


  // Model
  y_init ~ lognormal(log(state_init), sigma); // Separately fitting initial state to feed to integrator.

  for (t in 1:N) {
    Y[t,] ~ lognormal(log(State[t,]), sigma); // Lognormal model of states.
    // Y[t,] ~ normal(State[t,], sigma); // Lognormal model of states.
  }
}
