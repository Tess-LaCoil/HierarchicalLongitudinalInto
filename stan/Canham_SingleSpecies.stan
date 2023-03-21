//Growth function
functions{
  //Growth function for use with Runge-Kutta method
  real growth(real y, real max_growth, real size_max_growth, real K){
    return max_growth * 
    exp(-0.5 * pow(log(y / size_max_growth) / K, 2));
  }
  
  real midpoint(real y, real max_growth, real size_max_growth, real K, real interval){
    real mid;
    
    mid = y + 0.5 * interval * growth(y, max_growth, size_max_growth, K);
    return growth(mid, max_growth, size_max_growth, K);
  }
  
  real rk4(real y, real max_growth, real size_max_growth, real K, real interval){
    real k1;
    real k2;
    real k3;
    real k4;
    real g_est;
    
    k1 = growth(y, max_growth, size_max_growth, K);
    k2 = growth(y+interval*k1/2.0, max_growth, size_max_growth, K);
    k3 = growth(y+interval*k2/2.0, max_growth, size_max_growth, K);
    k4 = growth(y+interval*k3, max_growth, size_max_growth, K);
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
  real G_hat[N_obs];
  
  for(i in 1:N_obs){
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0[treeid_factor[i]];
    }
    
    if(i < N_obs){
      if(treeid_factor[i+1]==treeid_factor[i]){
        //Estimate growth rate
        if(int_method == 1){ //Euler method
          G_hat[i] = growth(S_hat[i], 
                          ind_max_growth[treeid_factor[i]],
                          ind_diameter_at_max_growth[treeid_factor[i]], 
                          ind_K[treeid_factor[i]]);
                            
        } else if(int_method == 2){ //Midpoint method
          G_hat[i] = midpoint(S_hat[i], 
                              ind_max_growth[treeid_factor[i]],
                              ind_diameter_at_max_growth[treeid_factor[i]], 
                              ind_K[treeid_factor[i]], 
                              census_interval[i+1]);
        } else if(int_method == 3){ //Midpoint method
          G_hat[i] = rk4(S_hat[i], 
                          ind_max_growth[treeid_factor[i]],
                          ind_diameter_at_max_growth[treeid_factor[i]], 
                          ind_K[treeid_factor[i]], 
                          census_interval[i+1]);
        }
                          
        S_hat[i+1] = S_hat[i] + G_hat[i]*census_interval[i+1];
      } else {
        G_hat[i] = 0; //Gives 0 as the growth estimate for the last obs of particular tree.
      }
    } else {
      G_hat[i] = 0; //Gives 0 as the growth estimate for the last data point.
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
  species_max_growth_mean ~normal(0, 2);
  species_max_growth_sd ~cauchy(4, 2);
  species_diameter_at_max_growth_mean ~normal(5, 1);
  species_diameter_at_max_growth_sd ~cauchy(4, 2);
  species_K_mean ~normal(-0.6694307, 1);
  species_K_sd ~cauchy(0.1, 1);
  
  //Global level
  global_error_sigma ~cauchy(1, 5);
}

generated quantities{
  real S_hat[N_obs];
  real G_hat[N_obs];
  
  for(i in 1:N_obs){
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0[treeid_factor[i]];
    }
    
    if(i < N_obs){
      if(treeid_factor[i+1]==treeid_factor[i]){
        //Estimate growth rate
        if(int_method == 1){ //Euler method
          G_hat[i] = growth(S_hat[i], 
                          ind_max_growth[treeid_factor[i]],
                          ind_diameter_at_max_growth[treeid_factor[i]], 
                          ind_K[treeid_factor[i]]);
                            
        } else if(int_method == 2){ //Midpoint method
          G_hat[i] = midpoint(S_hat[i], 
                              ind_max_growth[treeid_factor[i]],
                              ind_diameter_at_max_growth[treeid_factor[i]], 
                              ind_K[treeid_factor[i]], 
                              census_interval[i+1]);
        } else if(int_method == 3){ //Midpoint method
          G_hat[i] = rk4(S_hat[i], 
                          ind_max_growth[treeid_factor[i]],
                          ind_diameter_at_max_growth[treeid_factor[i]], 
                          ind_K[treeid_factor[i]], 
                          census_interval[i+1]);
        }
                          
        S_hat[i+1] = S_hat[i] + G_hat[i]*census_interval[i+1];
      } else {
        G_hat[i] = 0; //Gives 0 as the growth estimate for the last obs of particular tree.
      }
    } else {
      G_hat[i] = 0; //Gives 0 as the growth estimate for the last data point.
    }
  }
}
