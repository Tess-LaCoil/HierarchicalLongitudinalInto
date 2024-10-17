//Growth function
functions{
  //Growth function for use with Runge-Kutta method
  real growth(real y, vector ind_pars[1]){
    return ind_pars[1];
  }
}

// Data structure
data {
  int int_method;
  int N_obs;
  real S_obs[N_obs];
  int census[N_obs];
  real census_interval[N_obs];
  real S_0_obs;
}

// The parameters accepted by the model.
parameters {
  //Individual level
  real<lower=0> ind_S_0;
  real<lower=0> ind_beta;
  
  //Global level
  real<lower=0> global_error_sigma;
}

// The model to be estimated. 
model {
  real S_hat[N_obs];
  real G_hat[N_obs];
  vector ind_pars[1];
  
  for(i in 1:N_obs){
    ind_pars[1] = ind_beta;
    
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0;
    }
      
    if(i < N_obs){
      G_hat[i] = growth(S_hat[i], ind_pars) * census_interval[i+1];
      S_hat[i+1] = S_hat[i] + G_hat[i];
    } 
  }
  
  //Likelihood
  S_obs ~ normal(S_hat, global_error_sigma);
  
  //Priors
  //Individual level
  ind_S_0 ~ normal(S_0_obs, global_error_sigma);
  ind_beta ~lognormal(0, 4);
  
  //Global level
  global_error_sigma ~cauchy(1, 1);
}

generated quantities{
  real S_hat[N_obs];
  real G_hat[N_obs];
  vector ind_pars[1];
  
  for(i in 1:N_obs){
    ind_pars[1] = ind_beta;
    
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0;
    }
      
    if(i < N_obs){
      G_hat[i] = growth(S_hat[i], ind_pars) * census_interval[i+1];
      S_hat[i+1] = S_hat[i] + G_hat[i];
    } else {
      G_hat[i] = growth(S_hat[i], ind_pars) * census_interval[i];
    }
  }
}
