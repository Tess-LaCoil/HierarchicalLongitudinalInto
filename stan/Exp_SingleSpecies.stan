//Growth function
functions{
  //Growth function for use with Runge-Kutta method
  real growth(real y, real beta){
    return beta * y;
  }
  
  real euler(real y, real beta, real interval){
    real y_hat;
    
    y_hat = y + growth(y, beta)*interval;
    
    return y_hat;
  }
  
  real midpoint(real y, real beta, real interval){
    real mid;
    
    mid = y + 0.5 * interval * growth(y, beta);
    return growth(mid, beta);
  }
  
  real rk4(real y, real beta, real interval){
    real k1;
    real k2;
    real k3;
    real k4;
    real g_est;
    
    k1 = growth(y, beta);
    k2 = growth(y+interval*k1/2.0, beta);
    k3 = growth(y+interval*k2/2.0, beta);
    k4 = growth(y+interval*k3, beta);
    g_est = (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
    
    return g_est;
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
  real beta_temp; // (beta)

  for(i in 1:N_obs){
    //get parameters
    beta_temp = ind_beta[treeid_factor[i]];
    
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0[treeid_factor[i]];
    }
    
    //Assign next size
    if(i < N_obs){ //Avoid writing outside the bounds of the data
      if(treeid_factor[i+1]==treeid_factor[i]){ //Don't overwrite next individual
    
    //Estimate growth
      if(int_method == 1){ //Euler method
        G_hat[i] = euler(S_hat[i], beta_temp,
                        census_interval[i+1]);

      } else if(int_method == 2){ //Midpoint method
        G_hat[i] = midpoint(S_hat[i], beta_temp,
                        census_interval[i+1]);
        
      } else if(int_method == 3){ //RK4 method
        G_hat[i] = rk4(S_hat[i], beta_temp,
                        census_interval[i+1]);
      }
    
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
  real beta_temp; // (beta)
  
  
  for(i in 1:N_obs){
    //get parameters
    beta_temp = ind_beta[treeid_factor[i]];
    //Assign next size
    if(i < N_obs){ //Avoid writing outside the bounds of the data
      if(census[i]==1){//Fits the first size
        S_hat[i] = ind_S_0[treeid_factor[i]];
      }
      
      if(treeid_factor[i+1]==treeid_factor[i]){ //Don't overwrite next individual
      
      //Estimate growth
        if(int_method == 1){ //Euler method
          G_hat[i] = euler(S_hat[i], beta_temp,
                          census_interval[i+1]);
  
        } else if(int_method == 2){ //Midpoint method
          G_hat[i] = midpoint(S_hat[i], beta_temp,
                          census_interval[i+1]);
          
        } else if(int_method == 3){ //RK4 method
          G_hat[i] = rk4(S_hat[i], beta_temp,
                          census_interval[i+1]);
        }
    
        S_hat[i+1] = S_hat[i] + G_hat[i];
      } else { #If final value, use previous census interval instead.
        if(int_method == 1){ //Euler method
          G_hat[i] = euler(S_hat[i], beta_temp,
                          census_interval[i]);
  
        } else if(int_method == 2){ //Midpoint method
          G_hat[i] = midpoint(S_hat[i], beta_temp,
                          census_interval[i]);
          
        } else if(int_method == 3){ //RK4 method
          G_hat[i] = rk4(S_hat[i], beta_temp,
                          census_interval[i]);
        }
      }
    } else {
      if(int_method == 1){ //Euler method
        G_hat[i] = euler(S_hat[i], beta_temp,
                        census_interval[i]);
  
      } else if(int_method == 2){ //Midpoint method
        G_hat[i] = midpoint(S_hat[i], beta_temp,
                        census_interval[i]);
        
      } else if(int_method == 3){ //RK4 method
        G_hat[i] = rk4(S_hat[i], beta_temp,
                        census_interval[i]);
      }
    }
  }
}
