//Growth function
functions{
  //Growth function for use with Runge-Kutta method
  real growth(real y, real g_max, real S_max, real k){
    return g_max * exp(-0.5 * pow(log(y / S_max) / k, 2));
  }
}

// Data structure
data {
  int N_obs;
  real S_obs[N_obs];
  real delta_obs[N_obs];
  real census_interval[N_obs];
}

// The parameters accepted by the model.
parameters {
  //Species level
  real<lower=0> species_max_growth;
  real<lower=0> species_diameter_at_max_growth;
  real<lower=0> species_K;

  //Global level
  real<lower=0> global_error_sigma;
}

// The model to be estimated.
model {
  real delta_hat[N_obs];
  
  for(i in 1:N_obs){
    delta_hat[i] = growth(S_obs[i], species_max_growth, 
            species_diameter_at_max_growth,
            species_K) * census_interval[i];
  }

  //Likelihood
  delta_obs ~ normal(delta_hat, global_error_sigma);

  //Priors
  //Species level
  species_max_growth ~lognormal(0, 5);
  species_diameter_at_max_growth ~lognormal(1,5);
  species_K ~lognormal(0, 5);

  //Global level
  global_error_sigma ~cauchy(0, 5);
}


generated quantities{
  real delta_hat[N_obs];
  
  for(i in 1:N_obs){
    delta_hat[i] = growth(S_obs[i], species_max_growth, 
            species_diameter_at_max_growth,
            species_K) * census_interval[i];
  }
}
