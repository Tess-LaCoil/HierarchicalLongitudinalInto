//Growth function
functions{
  //Growth function for use with Runge-Kutta method
  real growth(real y, real beta){
    return beta * y;
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
  
  for(i in 1:N_obs){
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0;
    }
    
    if(i < N_obs){
      if(int_method == 1){ //Euler method
        G_hat[i] = growth(S_hat[i], ind_beta);
                          
      } else if(int_method == 2){ //Midpoint method
        G_hat[i] = midpoint(S_hat[i], 
                            ind_beta,
                            census_interval[i+1]);
                          
      } else if(int_method == 3){ //RK4 method
        G_hat[i] = rk4(S_hat[i], 
                        ind_beta,
                        census_interval[i+1]);
      }
  
      S_hat[i+1] = S_hat[i] + G_hat[i]*census_interval[i+1];
    } else {
      G_hat[i] = 0; //Gives 0 as the growth estimate for the last data point.
    }
  }
  
  //Likelihood
  S_obs ~ normal(S_hat, global_error_sigma);
  
  //Priors
  //Individual level
  ind_S_0 ~ normal(S_0_obs, global_error_sigma);
  ind_beta ~lognormal(0, 0.5);
  
  //Global level
  global_error_sigma ~cauchy(1, 1);
}

generated quantities{
  real S_hat[N_obs];
  real G_hat[N_obs];
  
  for(i in 1:N_obs){
    if(census[i]==1){//Fits the first size
      S_hat[i] = ind_S_0;
    }
    
    if(i < N_obs){
      if(int_method == 1){ //Euler method
        G_hat[i] = growth(S_hat[i], ind_beta);
                          
      } else if(int_method == 2){ //Midpoint method
        G_hat[i] = midpoint(S_hat[i], 
                            ind_beta,
                            census_interval[i+1]);
                          
      } else if(int_method == 3){ //Midpoint method
        G_hat[i] = rk4(S_hat[i], 
                        ind_beta,
                        census_interval[i+1]);
      }
  
      S_hat[i+1] = S_hat[i] + G_hat[i]*census_interval[i+1];
    } else {
      G_hat[i] = 0; //Gives 0 as the growth estimate for the last data point.
    }
  }
}
