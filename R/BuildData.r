#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                           Data Construction Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Build datasets
build_data <- function(data_construction_controls){
  for(i in data_construction_controls$data_build_models){
    print(paste("Building data for model:", i$model_name, sep=" "))
    
    source(i$model_spec, local=TRUE) #Load in spec
    
    filename <- paste(data_construction_controls$rstan_data_locations$single_species_data,
                      i$model_name, sep="_")
    filename <- paste(filename, ".rds", sep="")
    build_sim_data(growth_function, growth_function_pars, 
                   data_construction_controls$sim_pars, 
                   sim_spec, filename) #Produce data
    
    if(data_construction_controls$build_data_single_individual){ #Sample single individual
      data <- readRDS(filename)$sim_data
      
      samp_id <- sample(data$treeid, size=2)[1] #Trick to get around sample function issues
      ind <- data %>% filter(treeid == samp_id)
      
      
      
      rstan_data_ind <- build_rstan_data(N_obs = data_construction_controls$sim_pars$N_obs, 
                                         S_obs = ind$S_obs, 
                                         census = ind$census, 
                                         census_interval = ind$census_interval)
      save_data <- list(sim_data = ind, rstan_data = rstan_data_ind)
      
      filename <- paste(data_construction_controls$rstan_data_locations$single_individual_data,
                        i$model_name, sep="_")
      filename <- paste(filename, ".rds", sep="")
      saveRDS(save_data, file=filename)
    }
  }
}

#Add error based on Ruger et al 2011 error model and parameter estimates
add_error <- function(size){
  SD1 <- 0.927 + 0.0038*size
  SD2 <- 2.56
  size_with_error <- size + rnorm(1, mean=0, sd=SD1) + rbinom(1, size=1, prob=0.0276)*2.56
  return(size_with_error)
}

#Produce list of RStan data lists for each integration method
build_rstan_data <- function(N_obs, S_obs, census, census_interval){
  rstan_data <- list()
  for(i in 1:3){
    rstan_data[[i]] <- list(
      int_method = i,
      N_obs = N_obs,
      S_obs = S_obs,
      census = census,
      census_interval = census_interval,
      S_0_obs = S_obs[1]
    )
  }
  return(rstan_data)
}

#Produce list of RStan data lists for each integration method with multiple individuals
build_rstan_data_multi_ind <- function(N_obs, N_ind, sim_data_list){
  rstan_data <- list()
  sim_data <- sim_data_list$sim_data
  for(i in 1:3){
    rstan_data[[i]] <- list(
      int_method = i,
      N_obs = N_obs,
      N_ind = N_ind,
      S_obs = sim_data$S_obs,
      treeid_factor = sim_data$treeid, #Vector indexed by N_obs
      census = sim_data$census, #Vector indexed by N_obs
      census_interval = sim_data$census_interval, #Vector indexed by N_obs
      S_0_obs = sim_data_list$S_0_obs, #Vector indexed by N_ind
      tree_id_vec = sim_data_list$tree_id_vec #Vector indexed by N_ind
    )
  }
  return(rstan_data)
}

#Build simulated data based on specified growth function and parameter distributions
build_sim_data <- function(growth_function, growth_function_pars, sim_pars, sim_spec, filename){
  #Extract values to easier variable names
  time <- sim_pars$time
  N_step <- sim_pars$N_step
  N_ind <- sim_pars$N_ind
  
  #Construct initial tibble
  sim_data <- tibble(
    treeid = seq_len(N_ind),
    S_0 = exp(rnorm(N_ind, #Get initial size.
                    sim_pars$S0_mean, 
                    sim_pars$S0_sd)) +1, 
    sim_spec(sim_pars$N_ind)
  ) %>%
    expand_grid(time) %>%
    arrange(treeid, time) %>% 
    group_by(treeid) %>%
    mutate(census=rank(time)) %>%
    group_by(time) %>%
    mutate(S_true = ifelse(time==min(time), S_0, NA),
           S_obs = NA,
           census_interval = sim_pars$census_interval) %>%
    ungroup()
  S_0_obs <- c()
  tree_id_vec <- c(1:sim_pars$N_ind)
  
  #Uses the RK4 algorithm to give a simulated growth based on the true values.
  #Adds measurement error based on Ruger et al 2011.
  for(i in 1:sim_pars$N_ind){
    #The number of steps is determined by the level of precision required in step_size
    N_step <- sim_pars$N_obs * sim_pars$census_interval / sim_pars$step_size
    
    #We start at a random size +1cm in order to ensure a minimum size of 1
    runge_kutta_int <- rk4_est(S_0 = (sim_data$S_0[((i-1)*sim_pars$N_obs+1)]),
                               growth = growth_function,
                               pars = sim_data[i,growth_function_pars],
                               sim_pars$step_size, N_step)
    
    #Take a subset of the estimates which are the observations
    runge_kutta_obs <- runge_kutta_int[seq(from=1, 
                           to = length(runge_kutta_int), 
                           by = (length(runge_kutta_int)/sim_pars$N_obs))]
    
    #Check that the treeid lines up and add error
    if(sim_data$treeid[((i-1)*sim_pars$N_obs+1)] == sim_data$treeid[(i*sim_pars$N_obs)]){
      sim_data$S_true[((i-1)*sim_pars$N_obs+1): (i*sim_pars$N_obs)] <- runge_kutta_obs
      
      for(j in 1:sim_pars$N_obs){ #Add error, abs() to prevent negative size
        sim_data$S_obs[((i-1)*sim_pars$N_obs+j)] <- abs(add_error(runge_kutta_obs[j]))
      }
      
    } else {
      print("Indexing error: check treeid")
      return(-1)
    }
    S_0_obs[i] <- sim_data$S_obs[((i-1)*(sim_pars$N_obs)+1)]
  }
  
  sim_data <- sim_data %>% arrange(treeid, census) %>%
    group_by(treeid) %>%
    mutate(
      S_true_next = lead(S_true),
      S_obs_next = lead(S_obs),
      delta_S = (S_true_next - S_true)/census_interval,
      delta_S_obs = (S_obs_next - S_obs)/census_interval
    ) %>%
    ungroup()
  
  sim_data_list <- list(
    sim_data = sim_data,
    S_0_obs = S_0_obs,
    tree_id_vec = c(1:N_ind)
  )
  
  sim_data_list$rstan_data <- build_rstan_data_multi_ind(N_obs = (sim_pars$N_obs*sim_pars$N_ind), 
                                                         N_ind = sim_pars$N_ind, 
                                                         sim_data_list)
  
  saveRDS(sim_data_list, file=filename)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                     Numerical Integration Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Euler method
euler_est <- function(S_0, growth, pars, step_size, N_step){
  euler_int <- c(S_0)
  for(i in 2:N_step){
    euler_int[i] <- euler_int[i-1] + growth(euler_int[i-1], pars) * step_size
  }
  return(euler_int)
}


#Midpoint method
midpoint_est <- function(S_0, growth, pars, step_size, N_step){
  midpoint_int <- c(S_0)
  for(i in 2:N_step){
    midpoint <- midpoint_int[i-1] + 0.5 * step_size * growth(midpoint_int[i-1], pars)
    midpoint_int[i] <- midpoint_int[i-1] + growth(midpoint, pars) * step_size
  }
  return(midpoint_int)
}


#Runge-Kutta 4th order
rk4_est <- function(S_0, growth, pars, step_size, N_step){
  runge_kutta_int <- c(S_0)
  for(i in 2:N_step){
    k1 <- growth(runge_kutta_int[i-1], pars)
    k2 <- growth((runge_kutta_int[i-1] + step_size*k1/2), pars)
    k3 <- growth((runge_kutta_int[i-1] + step_size*k2/2), pars)
    k4 <- growth((runge_kutta_int[i-1] + step_size*k3), pars)
    
    runge_kutta_int[i] <- runge_kutta_int[i-1] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*step_size
  }
  return(runge_kutta_int)
}