//Growth function
functions{
  //Growth function for use with Runge-Kutta method
  real growth(real y, vector ind_pars[3]){
    return ind_pars[1] *
    exp(-0.5 * pow(log(y / ind_pars[2]) / ind_pars[3], 2));
  }

  real midpoint(real y, vector ind_pars[3], real interval){
    real mid;

    mid = y + 0.5 * interval * growth(y, ind_pars);
    return growth(mid, ind_pars);
  }

  real rk4(real y, vector ind_pars[3], real interval){
    real k1;
    real k2;
    real k3;
    real k4;
    real g_est;

    k1 = growth(y, ind_pars);
    k2 = growth(y+interval*k1/2.0, ind_pars);
    k3 = growth(y+interval*k2/2.0, ind_pars);
    k4 = growth(y+interval*k3, ind_pars);
    g_est = (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);

    return g_est;
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
  real<lower=0> ind_max_growth;
  real<lower=0> ind_diameter_at_max_growth;
  real<lower=0> ind_K;
  
  //Global level
  real<lower=0> global_error_sigma;
}

// The model to be estimated. 
model {
  real S_hat[N_obs];
  real G_hat[N_obs];
  vector ind_pars[3];
  
  for(i in 1:N_obs){
    //get parameters
    ind_pars[1] = ind_max_growth;
    ind_pars[2] = ind_diameter_at_max_growth;
    ind_pars[3] = ind_K;
    
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0;
    }
    
    //Estimate growth
    if(int_method == 1){ //Euler method
        G_hat[i] = growth(S_hat[i], ind_pars);

      } else if(int_method == 2){ //Midpoint method
        G_hat[i] = midpoint(S_hat[i], ind_pars,
                            census_interval[i+1]);
      } else if(int_method == 3){ //RK4 method
        G_hat[i] = rk4(S_hat[i], ind_pars,
                        census_interval[i+1]);
      }
    
    //Assign next size
    if(i < N_obs){ //Avoid writing outside the bounds of the data\
      S_hat[i+1] = S_hat[i] + G_hat[i]*census_interval[i+1];
    }
  }
  
  //Likelihood
  S_obs ~ normal(S_hat, global_error_sigma);
  
  //Priors
  //Individual level
  ind_S_0 ~ normal(S_0_obs, global_error_sigma);
  ind_max_growth ~lognormal(0, 1);
  ind_diameter_at_max_growth ~lognormal(3, 1);
  ind_K ~lognormal(0, 0.5);
  
  //Global level
  global_error_sigma ~cauchy(1,5);
}

generated quantities{
  real S_hat[N_obs];
  real G_hat[N_obs];
  vector ind_pars[3];
  
  for(i in 1:N_obs){
    //get parameters
    ind_pars[1] = ind_max_growth;
    ind_pars[2] = ind_diameter_at_max_growth;
    ind_pars[3] = ind_K;
    
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0;
    }
    
    
    if(i < N_obs){ //Avoid writing outside the bounds of the data\
      //Estimate growth
     if(int_method == 1){ //Euler method
        G_hat[i] = growth(S_hat[i], ind_pars)*census_interval[i+1];

      } else if(int_method == 2){ //Midpoint method
        G_hat[i] = midpoint(S_hat[i], ind_pars,
                            census_interval[i+1]);
      } else if(int_method == 3){ //RK4 method
        G_hat[i] = rk4(S_hat[i], ind_pars,
                        census_interval[i+1]);
      }
    
    //Assign next size
      S_hat[i+1] = S_hat[i] + G_hat[i];
    } else {
      if(int_method == 1){ //Euler method
        G_hat[i] = growth(S_hat[i], ind_pars)*census_interval[i];

      } else if(int_method == 2){ //Midpoint method
        G_hat[i] = midpoint(S_hat[i], ind_pars,
                            census_interval[i]);
      } else if(int_method == 3){ //RK4 method
        G_hat[i] = rk4(S_hat[i], ind_pars,
                        census_interval[i]);
      }
    
    }
  }
}
