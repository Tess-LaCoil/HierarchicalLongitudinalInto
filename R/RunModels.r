build_all_models <- function(model_build_controls){
  #Iterate through the models in the provided list
  for(i in 1:length(model_build_controls$model_list)){
    #Load in model specs
    source(model_build_controls$model_list[[i]]$model_spec, local=TRUE)
    
    #Build models at each level required
      if(model_build_controls$level_of_models$single_species){ #Single species model
      temp_filename <- model_build_controls$rstan_data_locations$single_species_data
      temp_filename <- paste(temp_filename,
                             model_build_controls$model_list[[i]]$model_name, 
                             sep="_")
      single_species_controls$rstan_data_loc <- paste(temp_filename, ".rds", sep="")
      
      build_fit(single_species_controls, 
                model_build_controls$plot_diagnostics, 
                model_build_controls$est_method,
                model_build_controls$int_method,
                model_build_controls$inc_warmup)
    }
    
    if(model_build_controls$level_of_models$single_individual){ #Single individual model
      temp_filename <- model_build_controls$rstan_data_locations$single_individual_data
      temp_filename <- paste(temp_filename,
                             model_build_controls$model_list[[i]]$model_name, 
                             sep="_")
      single_individual_controls$rstan_data_loc <- paste(temp_filename, ".rds", sep="")
      
      build_fit(single_individual_controls, 
                model_build_controls$plot_diagnostics, 
                model_build_controls$est_method,
                model_build_controls$int_method,
                model_build_controls$inc_warmup)
    }
  }
}


#-----------------------------------------------------------------------#
#                        Generic model functions
#-----------------------------------------------------------------------#
build_fit <- function(model_controls, plot_diagnostics, 
                      est_method, int_method, inc_warmup){
  if(file.exists(model_controls$rstan_data_loc)){ #Check that the data file exists
    rstan_data <- readRDS(model_controls$rstan_data_loc)$rstan_data
  } else {
    print(paste("Data file not found:", model_controls$rstan_data_loc, ". Please build data."))
    return(NA)
  }
  
  fit <- list()
  
  for(i in int_method){
    fit[[i]] <- build_rstan_model(stan_controls = model_controls$stan_controls, 
                             model_name = model_controls$model_name,
                             stan_file_name = model_controls$stan_file_path,
                             est_method,
                             rstan_data[[i]])
    
    #Plot diagnostics
    if(plot_diagnostics && (est_method == "samp")){
      plot_all_diagnostics(fit[[i]], int_method = i, 
                           model_controls$diagnostic_list,
                           inc_warmup = inc_warmup, 
                           RHat_name = model_controls$model_name)
    }
  }
  
  filename <- paste(model_controls$fit_save_path, est_method, 
                    paste(int_method, collapse="_"), sep="_")
  filename <- paste(filename, ".rds", sep="")
  saveRDS(fit, file=filename) #Output to file
  
  rm(fit)#Delete object from memory after it has been saved.
}

#Function to load stan file and build model
build_rstan_model <- function(stan_controls, 
                              model_name, 
                              stan_file_name,
                              est_method, 
                              rstan_data){
  print(paste("Fitting", model_name))
  
  model <- stan_model(file=stan_file_name)
  
  start_time <- Sys.time()
  print(paste("Estimation method:", est_method))
  print(paste("Integration method:", rstan_data$int_method))
  print(paste("Model fit start at:", start_time))
  
  if(est_method == "samp"){
    fit <- sampling(model, data=rstan_data, 
                    iter=stan_controls$iter,
                    chains=stan_controls$chains,
                    cores=stan_controls$cores,
                    control=stan_controls$control)
  } else if(est_method == "opt"){
    fit <- optimizing(model, data=rstan_data)
  }
  
  
  end_time <- Sys.time()
  #duration <- end_time - start_time
  print(paste("Model fit end at:", end_time))
  return(fit)
}

#-----------------------------------------------------------------------#
#                        Diagnostic plotting
#-----------------------------------------------------------------------#
#Diagnostic plots requested for particular fit
# list_of_data is a list of lists, each entry is a list consisting 
# of a vector of parameters to be plotted, and the name of the plot.
plot_all_diagnostics <- function(fit, int_method, list_of_data, inc_warmup, RHat_name){
  for(i in 1:length(list_of_data)){
    plot_diagnostic_trace(fit, int_method, 
                          pars=list_of_data[[i]]$pars, 
                          name=list_of_data[[i]]$name,
                          inc_warmup)
  } 
  
  #Histogram of R_hat
  #Expect to see NaNs for each final size as growth is not estimated.
  fit_summary <- summary(fit)
  Rhat <- fit_summary$summary[,10]
  n_NaN <- length(which(is.na(Rhat)))
  n <- length(Rhat)
  
  filename <- paste0("output/figures/diagnostic/", RHat_name, int_method, "_RHat_Hist.png")
  png(filename, width=500, height=300)
  hist(Rhat, main=paste(n_NaN,"NaN values from", n, "Rhats"))
  dev.off()
}

#Save diagnostic plot to file
plot_diagnostic_trace <- function(fit, int_method, pars, name, inc_warmup, size=800){
  filename <- paste0("output/figures/diagnostic/", name, int_method, ".png")
  png(filename, width=size, height=size)
  print(traceplot(fit, pars=pars, inc_warmup=inc_warmup))
  dev.off()
}
