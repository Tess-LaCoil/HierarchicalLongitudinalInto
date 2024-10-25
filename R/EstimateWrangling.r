#Currently only operates for multi-species models as the single species
# and single individual models are used for testing.

#---------------------------------------------------------------------------#
#                       Model-level wrangling
#---------------------------------------------------------------------------#
build_extracted_samples <- function(data_extract_controls){
  for(i in data_extract_controls$model_list){ #For each model, do the data wrangling
    #Multiple individuals in the same species
    if(data_extract_controls$level_of_models$single_species){
      sim_data_loc <- paste(data_extract_controls$rstan_data_locations$single_species_data, 
                         i$model_name, sep="_")
      sim_data_loc <- paste(sim_data_loc, ".rds", sep="")
      
      sim_data <- readRDS(sim_data_loc)
      
      data_extraction(data_extract_controls, i, sim_data, level="multi_ind")
      
      
    } 
    
    #Single individual
    if(data_extract_controls$level_of_models$single_individual){
      sim_data_loc <- paste(data_extract_controls$rstan_data_locations$single_individual_data, 
                         i$model_name, sep="_")
      sim_data_loc <- paste(sim_data_loc, ".rds", sep="")
      
      sim_data <- readRDS(sim_data_loc)
      
      data_extraction(data_extract_controls, i, sim_data, level="single_ind")
      
    } 
  }
}

#Collate data from multiple individuals
data_extraction <- function(data_extract_controls, model_list_item, 
                            sim_data, level, fit=NULL){
  source(model_list_item$model_spec, local=TRUE)
  
  #Initialise data structures
  #Extract true values for the measurement level
  measurement_data <- sim_data$sim_data %>%
    select(treeid, time, census, S_true, S_obs, delta_S, delta_S_obs)
  
  #Extract true individual parameters
  ind_sim_data <- sim_data$sim_data %>%
    select(treeid, growth_function_pars) %>%
    unique()
  
  individual_data <- tibble(treeid = ind_sim_data$treeid)
  for(j in 1:length(growth_pars)){
    name <- paste("true", growth_par_names[j], sep="_")
    individual_data[[name]] <- ind_sim_data[[growth_function_pars[j]]]
  }
  
  #Get rstan data
  rstan_data <- sim_data$rstan_data[[data_extract_controls$int_method]]
  
  #Load fit if not provided
  if(is.null(fit)){
    if(level == "multi_ind"){
      filename <- "output/data/Single_Species_Fit"
    } else if(level == "single_ind"){
      filename <- "output/data/Single_Ind_Fit"
    }
    #Add the model details
    filename <- paste(filename,
                      model_list_item$model_file,
                      data_extract_controls$est_method,
                      data_extract_controls$int_method,
                      sep="_")
    filename <- paste(filename, ".rds", sep="")
    
    if(file.exists(filename)){ #Check if file exists
      fit <- readRDS(filename)
    } else {
      print(paste("Model fit not found:", filename))
      return()
    }
  }
  
  #Extract estimates at different levels
  #Extraction depends on which method of estimation was used.
  if(data_extract_controls$est_method =="samp"){
    est_data <- rstan::extract(fit[[data_extract_controls$int_method]], permuted=TRUE)
    size_estimates <- apply(est_data$S_hat, 2, mean)
    growth_estimates <- apply(est_data$G_hat, 2, mean)
    
  } else if(data_extract_controls$est_method == "opt"){
    est_data <- fit[[data_extract_controls$int_method]]$par
    size_estimates <- est_data[grepl("^S_hat", names(est_data))]
    names(size_estimates) <- NULL #Remove names
    growth_estimates <- est_data[grepl("^G_hat", names(est_data))]
    names(growth_estimates) <- NULL #Remove names
    
  } else {
    print("Please identify which estimation method was used: samp or opt.")
    return()
  }

  #Measurement level
  measurement_data$S_hat <- size_estimates #Adds size estimates for the particular model
  measurement_data$G_hat <- growth_estimates #Adds growth estimates
  
  #Calculate MSE for growth estimates, uses pairwise index for estimates to remove imputed 0
  MSE_growth <- build_MSE_vec(measurement_data$delta_S_obs[!is.na(measurement_data$delta_S_obs)],
                              growth_estimates[!is.na(measurement_data$delta_S_obs)],
                              measurement_data$delta_S[!is.na(measurement_data$delta_S)],
                              c("Pairwise Difference", "Model"))
  
  #Calculate MSE for sizes
  MSE_size <- build_MSE_vec(measurement_data$S_obs,
                            size_estimates,
                            measurement_data$S_true,
                            c("Observed Size", "Model"))
  
  #Individual level
  individual_data <- extract_ind(est_data, 
                                 data_extract_controls$est_method,
                                 ind_pars,
                                 ind_pars_names,
                                 individual_data,
                                 obs_tree_id_vec = rstan_data$treeid)
  
  #MSE for ind-level parameters
  MSE_ind <- ind_par_MSE(individual_data, growth_pars, growth_par_names,
                         obs_tree_id_vec = rstan_data$treeid)
  
  #Export data structure
  save_data <- list(individual_data = individual_data,
                    measurement_data = measurement_data,
                    MSE_growth = MSE_growth,
                    MSE_size = MSE_size,
                    MSE_ind = MSE_ind)
  
  #Construct data to output to file and filename
  if(level=="multi_ind"){
    #Species level data
    species_data <- extract_species(est_data,
                                    data_extract_controls$est_method,
                                    sp_pars, 
                                    sp_pars_names)
    
    save_data$species_data <- species_data
    
    #95% CIs for species hyperparameters
    hyperpar_CI <- build_par_CI(est_data, sp_pars, sp_pars_names)
    save_data$hyperpar_CI <- hyperpar_CI
    
    save_filename <- "output/data/Single_species_compiled_data"
  } else if(level =="single_ind"){
    save_filename <- "output/data/Single_individual_compiled_data"
  }
  
  save_filename <- paste(save_filename,
                         model_list_item$model_name,
                         data_extract_controls$int_method, sep="_")
  
  if(data_extract_controls$est_method == "samp"){
    save_filename <- paste(save_filename, "_sampling.rds", sep="")
  } else if(data_extract_controls$est_method == "opt"){
    save_filename <- paste(save_filename, "_optimizing.rds", sep="")
  }
  
  saveRDS(save_data, file=save_filename)
}

#Produce a vector which calculates MSE for model compared to pairwise difference
build_MSE_vec <- function(pairwise, model, true_val, names){
  pairwise_diff_MSE <- est_MSE(pairwise, true_val)
  growth_model_MSE <- est_MSE(model, true_val)
  
  MSE_vec <- c(pairwise_diff_MSE, growth_model_MSE)
  MSE_vec <- setNames(MSE_vec, names)
  
  return(MSE_vec)
}

#Calculate mean squared error
est_MSE <- function(estimates, true_val){
  MSE <- sum((estimates-true_val)^2)/length(estimates)
  return(MSE)
}

#---------------------------------------------------------------------------#
#                       Individual level wrangling
#---------------------------------------------------------------------------#
#Extract individual parameters for each model
#Requires separate vector of tree ids per observation, obs_tree_id_vec
extract_ind <- function(est_data, 
                        est_method,
                        pars, 
                        pars_names,
                        individual_data,
                        obs_tree_id_vec){
  #Extract mean of parameter posterior distributions
  for(i in pars){
    if(est_method == "samp"){
      individual_data[[pars_names[which(pars == i)]]] <- apply(est_data[[i]], 2, mean)
      
    } else if (est_method == "opt"){
      individual_data[[pars_names[which(pars == i)]]] <- est_data[grepl(pars[which(pars == i)], 
                                                                        names(est_data))]
    }
    
    #Ratio and difference of true value to estimated for growth pars
    if(i != "ind_S_0"){ #Don't build ratio for the initial size value
      name <- paste("true", pars_names[which(pars == i)], sep="_")
      ratio_name <- paste("ratio", pars_names[which(pars == i)], sep="_")
      diff_name <- paste("difference", pars_names[which(pars == i)], sep="_")
      individual_data[[ratio_name]] <- individual_data[[name]]/individual_data[[pars_names[which(pars == i)]]]
      individual_data[[diff_name]] <- individual_data[[name]] - individual_data[[pars_names[which(pars == i)]]] 
    }
  }
  
  return(individual_data)
}

#Calculates the MSE for the ind-level growth parameters
ind_par_MSE<- function (individual_data, growth_pars, growth_par_names, obs_tree_id_vec){
  ind_MSE <- list()
  
  for(j in 1:length(growth_pars)){ #Iterate through list of parameters
    name <- paste("true", growth_par_names[j], sep="_")
    true <- individual_data[[name]]
    est <- individual_data[[growth_par_names[j]]]
    
    ind_MSE[[growth_pars[j]]] <- est_MSE(est, true)
  }
  
  return(ind_MSE)
}

#---------------------------------------------------------------------------#
#                         Species level wrangling
#---------------------------------------------------------------------------#
#Extract species parameter estimates
extract_species <- function(est_data,
                            est_method,
                            pars, 
                            pars_names){
  species_data <- c()
  
  for(i in 1:length(pars)){ #Extract estimates
    if(est_method == "samp"){
      species_data[i] <- mean(est_data[[pars[i]]])
    } else if (est_method == "opt"){
      species_data[i] <- est_data[grepl(pars[which(pars[i])], names(est_data))]
    }
  }
  
  #Name the values in vector
  species_data <- setNames(species_data, pars_names)
  
  return(species_data)
}

#Build a 95% CI for each parameter in a list
build_par_CI <- function(est_data, pars, pars_names){
  ci_list <- list()
  
  for(i in 1:length(pars)){ #Extract estimates, build CI
    data <- est_data[[pars[i]]]
    ci_list[[pars_names[i]]] <- bayes_ci(data)
  }
  
  return(ci_list)
}

#Build a 95% CI from the posterior samples
bayes_ci <-function(vec){
  #Get empirical quantiles from vec
  ci <- quantile(vec, probs=c(0.025,0.975))
  
  #Name elements of the vector
  ci <- setNames(ci, c("Lower", "Upper"))
  
  return(ci)
}

#-----------------------------------------------------------------------------#
#                 Predictive power estimation
#-----------------------------------------------------------------------------#
S_6_obs <- tibble_6th_obs$S_obs
census_interval <- tibble_6th_obs$census_interval
model_spec_name <- "R/ConstSpec.r"

extract_est_tibble <- function(fit, model_spec_name, dataset, S_6_obs, census_interval){
  source(model_spec_name, local=TRUE)
  #Extract samples
  est_data <- rstan::extract(fit, permuted=TRUE)
  #Building data frames
  measurement_data <- tibble(treeid_factor = dataset$treeid_factor,
                             census = dataset$census)
  measurement_data$S_hat <- apply(est_data$S_hat, 2, mean)
  temp_starting_size <- measurement_data %>%
    filter(census == 5) %>%
    mutate(S_5 = S_hat) %>%
    select(treeid_factor, S_5)
  temp_starting_size$S_6_obs <- S_6_obs
  temp_starting_size$census_interval <- census_interval
  
  temp_par_list <- list()
  par_list_vec <- c()
  for(j in 1:nrow(temp_starting_size)){
    temp_list <- list()
    for(k in ind_pars){
      temp_list[[ind_pars_names[which(ind_pars == k)]]] <- apply(est_data[[k]], 2, mean)
    }
    par_list_vec[j] <- temp_list
  }
}