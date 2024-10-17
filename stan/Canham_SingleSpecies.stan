//Growth function
functions{
  //Growth function for use with Runge-Kutta method
  real growth(real y, real g_max, real S_max, real k){
    return g_max * exp(-0.5 * pow(log(y / S_max) / k, 2));
  }
  
  real euler(real y, real g_max, real S_max, real k, real interval){
    real y_hat;
    
    y_hat = y + growth(y, g_max, S_max, k)*interval;
    
    return y_hat;
  }

  real midpoint(real y, real g_max, real S_max, real k, real interval){
    real mid;
    real y_hat;

    mid = y + 0.5 * interval * growth(y, g_max, S_max, k);
    
    y_hat = y + growth(mid, g_max, S_max, k) * interval;
    
    return y_hat;
  }

  real rk4_step(real y, real g_max, real S_max, real k, real interval){
    real k1;
    real k2;
    real k3;
    real k4;
    real y_hat;

    k1 = growth(y, g_max, S_max, k);
    k2 = growth(y+interval*k1/2.0, g_max, S_max, k);
    k3 = growth(y+interval*k2/2.0, g_max, S_max, k);
    k4 = growth(y+interval*k3, g_max, S_max, k);
    
    y_hat = y + (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4) * interval;

    return y_hat;
  }
  
  real rk4(real y, real g_max, real S_max, real k, real interval, real step_size){
    int steps;
    real duration;
    real y_hat;
    real step_size_temp;
    
    duration = 0;
    y_hat = y;
    
    while(duration < interval){
      //Determine the relevant step size
      step_size_temp = min([step_size, interval-duration]);
      
      //Get next size estimate
      y_hat = rk4_step(y_hat, g_max, S_max, k, step_size_temp);
      
      //Increment observed duration
      duration = duration + step_size_temp;
    }
    
    return y_hat;
  }
}

// Data structure
data {
  int int_method;
  real step_size;
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
  real<lower=0> ind_max_growth[N_ind];
  real<lower=0> ind_diameter_at_max_growth[N_ind];
  real<lower=0> ind_K[N_ind];

  //Species level
  real species_max_growth_mean;
  real<lower=0> species_max_growth_sd;
  real species_diameter_at_max_growth_mean;
  real<lower=0> species_diameter_at_max_growth_sd;
  real species_K_mean;
  real<lower=0> species_K_sd;

  //Global level
  real<lower=0> global_error_sigma;
}

// The model to be estimated.
model {
  real S_hat[N_obs];

  for(i in 1:N_obs){
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0[treeid_factor[i]];
    }
    
    //Estimate growth
    //Assign next size
    if(i < N_obs){ //Avoid writing outside the bounds of the data
      if(treeid_factor[i+1]==treeid_factor[i]){
        if(int_method == 1){ //Euler method
          S_hat[i+1] = euler(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i+1]);
    
        } else if(int_method == 2){ //Midpoint method
          S_hat[i+1] = midpoint(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i+1]);
                          
        } else if(int_method == 3){ //RK4 method
          S_hat[i+1] = rk4(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i+1], step_size);
        }
      }
    }
  }

  //Likelihood
  S_obs ~ normal(S_hat, global_error_sigma);

  //Priors
  //Individual level
  ind_S_0 ~ normal(S_0_obs, global_error_sigma);
  ind_max_growth ~lognormal(species_max_growth_mean, species_max_growth_sd);
  ind_diameter_at_max_growth ~lognormal(species_diameter_at_max_growth_mean,
                                        species_diameter_at_max_growth_sd);
  ind_K ~lognormal(species_K_mean, species_K_sd);

  //Species level
  species_max_growth_mean ~normal(0, 1);
  species_max_growth_sd ~cauchy(0, 1);
  species_diameter_at_max_growth_mean ~normal(0, 1);
  species_diameter_at_max_growth_sd ~cauchy(0, 1);
  species_K_mean ~normal(0, 1);
  species_K_sd ~cauchy(0, 1);

  //Global level
  global_error_sigma ~cauchy(0, 2);
}

generated quantities{
  real S_hat[N_obs];
  real G_hat[N_obs];
  
  for(i in 1:N_obs){
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0[treeid_factor[i]];
    }
    
    //Estimate growth
    //Assign next size
    if(i < N_obs){ //Avoid writing outside the bounds of the data
      if(treeid_factor[i+1]==treeid_factor[i]){
        if(int_method == 1){ //Euler method
          S_hat[i+1] = euler(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i+1]);
    
        } else if(int_method == 2){ //Midpoint method
          S_hat[i+1] = midpoint(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i+1]);
                          
        } else if(int_method == 3){ //RK4 method
          S_hat[i+1] = rk4(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i+1], step_size);
        }
        
        G_hat[i] = S_hat[i+1] - S_hat[i];
        
      } else { #Uses previous census interval to predict G_hat for final size
        if(int_method == 1){ //Euler method
          S_hat_temp = euler(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i]);
    
        } else if(int_method == 2){ //Midpoint method
          S_hat_temp = midpoint(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i]);
                          
        } else if(int_method == 3){ //RK4 method
          S_hat_temp = rk4(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i], step_size);
        }
        G_hat[i] = S_hat_temp - S_hat[i];
      }
    } else {
      if(int_method == 1){ //Euler method
          S_hat_temp = euler(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i]);
    
        } else if(int_method == 2){ //Midpoint method
          S_hat_temp = midpoint(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i]);
                          
        } else if(int_method == 3){ //RK4 method
          S_hat_temp = rk4(S_hat[i], ind_max_growth[treeid_factor[i]], 
            ind_diameter_at_max_growth[treeid_factor[i]], 
            ind_K[treeid_factor[i]], census_interval[i], step_size);
        }
        G_hat[i] = S_hat_temp - S_hat[i];
    }
  }
}
