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
    build_sim_data(data_construction_controls$error_type,
                   growth_function, growth_function_pars,
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

#Add error based on chosen model
add_error <- function(size, error_type="Norm", error_size = 0.001){
  if(error_type == "Ruger"){ #Model from Ruger 2011
    size_with_error <- size + ruger_error_mm(size*10)/10

  } else if(error_type == "Norm") { #Normally distributed error
    size_with_error <- size + rnorm(1, mean=0, sd=error_size)

  } else if(error_type == "none") { #No error
    size_with_error <- size

  } else { #No error type or wrong error type
    print("Please input valid error type.")
    size_with_error <- -1
  }

  return(size_with_error)
}

#Model and values from Ruger et al 2011 Growth Strategies of Tropica Tree Species
ruger_error_mm <- function(size_mm){#Model from Ruger 2011 p. 2 based on mm measurement
  SD1 <- 0.927 #0.927mm
  SD2 <- 0.0038*size_mm #Size-dependent bit
  SD3 <- 25.6 #25.6mm
  error_mm <- rnorm(1, mean=0, sd=(SD1 + SD2))
  + rbinom(1, size=1, prob=0.0276) * rnorm(1, 0, SD3)

  return(error_mm)
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
build_sim_data <- function(error_type, growth_function, growth_function_pars,
                           sim_pars, sim_spec, filename){
  #Extract values to easier variable names
  time <- sim_pars$time
  N_step <- sim_pars$N_step
  N_ind <- sim_pars$N_ind

  #Construct initial tibble
  sim_ind_info <- tibble(
    treeid = seq_len(N_ind),
    S_0 = exp(rnorm(N_ind, #Get initial size.
                    sim_pars$S0_mean,
                    sim_pars$S0_sd)) +1,
    sim_spec(sim_pars$N_ind) #Returns parameter values for each individual
  )

  if(sim_pars$model == "canham"){ #Control for identifiability problems in Canham
    for(i in 1:N_ind){
      while((sim_ind_info$S_0[i] >= 15)||(sim_ind_info$S_0[i] <= 1)){ #Constrain initial sizes for identifiability
        sim_ind_info$S_0[i] <- exp(rnorm(1, #Get initial size.
                                         sim_pars$S0_mean,
                                         sim_pars$S0_sd)) +1
      }
    }
  }

  sim_data <- sim_ind_info %>%
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
    runge_kutta_int <- rk4_est(S_0 = sim_ind_info$S_0[i],
                               growth = growth_function,
                               pars = sim_ind_info[i, growth_function_pars],
                               sim_pars$step_size, N_step)

    #Take a subset of the estimates which are the observations
    runge_kutta_obs <- runge_kutta_int[seq(from=1,
                                           to = length(runge_kutta_int),
                                           by = (length(runge_kutta_int)/sim_pars$N_obs))]

    #Check that the treeid lines up and add error
    if(sim_data$treeid[((i-1)*sim_pars$N_obs+1)] == sim_data$treeid[(i*sim_pars$N_obs)]){
      sim_data$S_true[((i-1)*sim_pars$N_obs+1): (i*sim_pars$N_obs)] <- runge_kutta_obs

      for(k in 1:sim_pars$N_obs){ #Add error, abs() to prevent negative size
        sim_data$S_obs[((i-1)*sim_pars$N_obs+k)] <- abs(add_error(runge_kutta_obs[k], error_type))
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
      delta_S = (lead(S_true) - S_true),
      delta_S_obs = (lead(S_obs) - S_obs),
      average_growth_obs = sum(delta_S_obs, na.rm=TRUE)/max(time),
      S_first = first(S_true),
      S_final = last(S_true),
      total_growth = S_final - S_first
    ) %>%
    ungroup()

  sim_ind_info$average_growth_obs <- sim_data$average_growth_obs[seq(from=1,
                                                                     by=sim_pars$N_obs,
                                                                     length.out=sim_pars$N_ind)]

  sim_data_list <- list(
    sim_data = sim_data,
    sim_ind_info = sim_ind_info,
    S_0_obs = S_0_obs,
    tree_id_vec = c(1:N_ind)
  )

  sim_data_list$rstan_data <- build_rstan_data_multi_ind(N_obs = (sim_pars$N_obs*sim_pars$N_ind),
                                                         N_ind = sim_pars$N_ind,
                                                         sim_data_list)

  if(sim_pars$model == "canham"){ #Add step size
    rstan_data <- c(sim_data_list$rstan_data[1],
                    step_size = 1,
                    sim_data_list$rstan_data[-1])
    sim_data_list$rstan_data <- rstan_data
  }

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
