//Growth function
functions{
  //Growth function for use with Runge-Kutta method
  real growth(real y, real beta){
    return beta;
  }
}

// Data structure
data {
  int int_method;
  int N_obs;
  int N_ind;
  real S_obs[N_obs];
  int treeid_factor[N_obs];
  int census[N_obs];
  real census_interval[N_obs];
  real S_0_obs[N_ind];
  int tree_id_vec[N_ind];
}

// The parameters accepted by the model.
parameters {
  //Individual level
  real<lower=0> ind_S_0[N_ind];
  real<lower=0> ind_beta[N_ind];
  
  real species_beta_mu;
  real<lower=0> species_beta_sigma;
  
  //Global level
  real<lower=0> global_error_sigma;
}

// The model to be estimated. 
model {
  real S_hat[N_obs];
  real G_hat[N_obs];
  real ind_pars; // (beta)

  for(i in 1:N_obs){
    //get parameters
    ind_pars = ind_beta[treeid_factor[i]];
    
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0[treeid_factor[i]];
    }
    
    //Assign next size
    if(i < N_obs){ //Avoid writing outside the bounds of the data
      if(treeid_factor[i+1]==treeid_factor[i]){ //Don't overwrite next individual
        //Estimate growth
        G_hat[i] = growth(S_hat[i], ind_pars) * census_interval[i+1];
    
        S_hat[i+1] = S_hat[i] + G_hat[i];
      } 
    }
  }
  
  //Likelihood
  S_obs ~ normal(S_hat, global_error_sigma);
  
  //Priors
  //Individual level
  ind_S_0 ~ normal(S_0_obs, global_error_sigma);
  ind_beta ~ lognormal(species_beta_mu, 
                    species_beta_sigma);
  
  //Species level
  species_beta_mu ~ normal(0.1, 1);
  species_beta_sigma ~cauchy(0.1, 1);
  
  //Global level
  global_error_sigma ~cauchy(0.1, 1);
}

generated quantities{
  real S_hat[N_obs];
  real G_hat[N_obs];
  real ind_pars; // (beta)

  for(i in 1:N_obs){
    //get parameters
    ind_pars = ind_beta[treeid_factor[i]];
    
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0[treeid_factor[i]];
    }
    
    //Assign next size
    if(i < N_obs){ //Avoid writing outside the bounds of the data
      if(treeid_factor[i+1]==treeid_factor[i]){ //Don't overwrite next individual
        //Estimate growth
        G_hat[i] = growth(S_hat[i], ind_pars) * census_interval[i+1];
    
        S_hat[i+1] = S_hat[i] + G_hat[i];
      } else { #Uses previous interval to estimate the growth over an equivalent time forwards
        G_hat[i] = growth(S_hat[i], ind_pars) * census_interval[i];
      }
    } else {
      G_hat[i] = growth(S_hat[i], ind_pars) * census_interval[i];
    }
  }
}
