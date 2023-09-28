#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                     Numerical Integration Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Euler method
euler_est <- function(S_0, growth, pars, step_size, nstep){
  euler_int <- c(S_0)
  for(i in 2:n_step){
    euler_int[i] <- euler_int[i-1] + growth(euler_int[i-1], pars) * step_size
  }
  return(euler_int)
}


#Midpoint method
midpoint_est <- function(S_0, growth, pars, step_size, nstep){
  midpoint_int <- c(S_0)
  for(i in 2:n_step){
    midpoint <- midpoint_int[i-1] + 0.5 * step_size * growth(midpoint_int[i-1], pars)
    midpoint_int[i] <- midpoint_int[i-1] + growth(midpoint, pars) * step_size
  }
  return(midpoint_int)
}


#Runge-Kutta 4th order
rk4_est <- function(S_0, growth, pars, step_size, nstep){
  runge_kutta_int <- c(S_0)
  for(i in 2:n_step){
    k1 <- growth(runge_kutta_int[i-1], pars)
    k2 <- growth((runge_kutta_int[i-1] + step_size*k1/2), pars)
    k3 <- growth((runge_kutta_int[i-1] + step_size*k2/2), pars)
    k4 <- growth((runge_kutta_int[i-1] + step_size*k3), pars)
    
    runge_kutta_int[i] <- runge_kutta_int[i-1] + (step_size/6)*(k1 + 2*k2 + 2*k3 + k4)
  }
  return(runge_kutta_int)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             Plotting Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Output plot to file
save_plot_of_estimates <- function(lines, points, names, colours, pch_vals, 
                                   lty_vals, title, file_name, lwd){
  png(file_name, width=500, height=500)
  plot_estimates(lines, points, names, colours, pch_vals, lty_vals, title, lwd)
  dev.off()
}

#Produce plot for save_plot_of_estimates
plot_estimates <- function(lines, points, names, colours, pch_vals, lty_vals, title, lwd){
  par(mai =  c(bottom=1, left=1, top=1, right=0.5))
  output_plot <- plot(x=time, y=lines[[1]], type="l", col=colours[1], 
       ylim=c(0,exp(4)),
       xlab="Time",
       ylab="Size",
       main=title,
       lwd=lwd,
       cex.lab = 2,
       cex.axis = 2,
       cex.main = 2)
  
  start <- 1
  if(length(lines) > 1){ start <- 2 }
  for(j in start:length(lines)){
    lines(time, lines[[j]], col=colours[j], lwd=lwd, cex=2)
  }
  
  for(j in 1:length(points)){
    points(time, points[[j]], pch=pch_vals[j], col=colours[j], lwd=lwd, cex=2)
  }
  curve(exp(x), from=0, to=4, col="black", add=TRUE, cex=2)
  legend(x=0.5, y=50, 
         legend=names,
         col=colours,
         lty=lty_vals,
         lwd=rep(2, times=length(lty_vals)),
         pch=pch_vals,
         bty="n", 
         cex=2)
  
  return(output_plot)
}

#Plot results from estimates
plot_results <- function(estimates, 
                         rstan_data,
                         est_names, 
                         title, 
                         file_name,
                         names,
                         colours,
                         lty_list,
                         pch_list){
  lines_list <- list()
  for(i in 1:length(estimates$S_hat)){
    lines_list[[i]] <- estimates$S_hat[[i]]
  }
  
  points_list <- lines_list
  points_list[[length(points_list)+1]] <- rstan_data[[1]]$S_obs
  save_plot_of_estimates(lines_list, points_list, names, colours, pch_list, 
                         lty_list, title, file_name, lwd=1.5)
}

#Produces side-by-side violin plot of beta samples for each model
produce_estplot <- function(model_list, model_names){
  #Initialise data frame
  plot_data <- data.frame(val = c(), model = c())
  plots <- list()
  means <- data.frame(val = c(1,1,1), model = model_names, 
                      lower= c(1,1,1),
                      upper = c(1,1,1))
  
  #Build plot data
  for(i in 1:length(model_list)){
    samp <- rstan::extract(model_list[[i]])
    temp <- data.frame(val = samp$ind_beta, 
                       model = rep(model_names[i], times=length(samp$ind_beta)))
    means$val[i] <- round(mean(temp$val), digits=4)
    means$model[i] <- model_names[i]
    means$lower[i] <- quantile(temp$val, 0.05)
    means$upper[i] <- quantile(temp$val, 0.95)
    
    plot_data <- rbind(plot_data, temp)
  }

  #Produce a plot to compare the different estimates of beta
  png("output/figures/RungeKuttaDemo/BetaPlot.png", width=500, height=500)
  build_estplot(means)
  dev.off()
}

#Plots the beta estimates
build_estplot <- function(means){
  pch_list <- c(3,4,5)
  col_list <- c("blue", "red", "green4")
  par(mai =  c(bottom=1, left=1, top=1, right=0.5))
  output_plot <- plot(x=means$model, y=means$val, xlab="Integration Method",
       ylab="Beta", lwd=0.2, main="Posterior Beta Estimates",
       cex.lab = 2,
       cex.axis = 2,
       cex.main = 2)
  for(i in 1:3){
    points(x=means$model[i], y=means$val[i], pch=pch_list[i], col=col_list[i], 
           cex=2, lwd=1.1)
  }
  abline(a=1, b=0, lty="dashed")
  
  return(output_plot)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                           Data Construction Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Add error based on Ruger et al 2011 error model and parameter estimates
add_error <- function(size, error_type="Norm"){
  if(error_type == "Ruger"){ #Model from Ruger 2011
    SD1 <- 0.927 + 0.0038*size
    SD2 <- 2.56
    size_with_error <- size + rnorm(1, mean=0, sd=SD1) + rbinom(1, size=1, prob=0.0276)*2.56
    
  } else if(error_type == "Norm") { #Normally distributed error
    size_with_error <- size + rnorm(1, mean=0, sd=1)
    
  } else { #No error type or wrong error type
    print("Please input valid error type.")
    size_with_error <- -1
  }
  
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
build_rstan_data_multi_ind <- function(N_obs, N_ind, sim_data, int_method, 
                                       S_0_obs, tree_id_vec){
  rstan_data <- list()
  for(i in 1:3){
    rstan_data[[i]] <- list(
      int_method = i,
      N_obs = N_obs,
      N_ind = N_ind,
      S_obs = sim_data$S_obs,
      treeid_factor = sim_data$treeid, #Vector indexed by N_obs
      census = sim_data$census, #Vector indexed by N_obs
      census_interval = sim_data$census_interval, #Vector indexed by N_obs
      S_0_obs = S_0_obs, #Vector indexed by N_ind
      tree_id_vec = tree_id_vec #Vector indexed by N_ind
    )
  }
  return(rstan_data)
}

#Build simulated data based on specified growth function and parameter distributions
build_sim_data <- function(growth_function, sim_pars){
  sim_data <- tibble(
    treeid = seq_len(sim_pars$n_ind),
    S_0 = exp(rnorm(sim_pars$n_ind, #Get initial size.
                    sim_pars$s0_mean, 
                    sim_pars$s0_sd)) 
    ) %>%
    mutate(#Get values for parameters
      for(i in sim_pars$spec){
        i$name = i$funct(i$pars) #CORRECT THIS
      }
    ) %>%
    expand_grid(time) %>%
    arrange(treeid, time) %>% 
    group_by(treeid, time) %>%
    mutate(S_true = ifelse(time==min(time), S_0, NA),
           S_obs = NA,
           census = time+1,
           census_interval = step_size)
  S_0_obs <- c()
  tree_id_vec <- c(1:n_ind)
  
  #Uses the RK4 algorithm to give a simulated growth based on the true values.
  #Adds measurement error based on Ruger et al 2011.
  for(i in 1:n_ind){
    runge_kutta_int <- rk4_est(S_0=c(sim_data$S_0[((i-1)*n_step+1)]),
                               growth = growth_function,
                               pars = c(sim_data$g_max[i], sim_data$s_max[i], sim_data$k[i]),
                               step_size,
                               nstep)
    if(sim_data$treeid((i-1)*n_step+1) == sim_data$treeid(i*n_step)){
      sim_data$S_true[((i-1)*n_step+1): (i*n_step)] <- runge_kutta_int
      for(j in 1:length(runge_kutta_int)){ #Add error to each term
        sim_data$S_obs[((i-1)*n_step+j)] <- add_error(runge_kutta_int[j])
      }
    } else {
      print("Indexing error: check treeid")
      return(-1)
    }
    S_0_obs[i] <- runge_kutta_int[1]
  }
  
  return_data <- list(
    sim_data = sim_data,
    S_0_obs = S_0_obs,
    tree_id_vec = c(1:n_ind)
  )
  
  return(return_data)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                           Estimate Wrangling Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Extract samples from models
est_wrangling <- function(models, names, est_method){
  S_hat <- list()
  beta_hat <- list()
  
  for(i in 1:length(models)){
    if(est_method == "samp"){
      samples <- rstan::extract(models[[i]])
      S_hat[[names[i]]] <- apply(samples$S_hat, 2, mean)
      beta_hat[[names[i]]] <- mean(samples$ind_beta)
      
    } else if(est_method == "opt")
    {
      S_hat[[names[i]]] <- models[[i]]$par[grepl("^S_hat", names(models[[i]]$par))]
      names(S_hat) <- NULL #Remove names
      
      beta_hat[[names[i]]] <- models[[i]]$par[grepl("^ind_beta", names(models[[i]]$par))]
      names(beta_hat) <- NULL #Remove names
    }
  }
  
  estimates <- list(S_hat=S_hat, beta_hat=beta_hat)
  return(estimates)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             Running Models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Run models, plot output for chosen growth function, estimation method
#Default values are for data with error
build_output <- function(rstan_file, 
                         rstan_data, 
                         est_names, 
                         est_method,
                         int_methods = c(3),
                         title, 
                         file_name,
                         names,
                         colours = c("green4", "black", "black"),
                         lty_list = c(1,0,1),
                         pch_list = c(5,1,NA),
                         alg="LBFGS",
                         estplot = FALSE){
  #Run model
  model_output <- perform_estimation(rstan_file, rstan_data, 
                                     est_names, est_method, alg,
                                     int_methods, estplot)
  
  #Plot results
  plot_results(model_output$estimates, 
               rstan_data,
               est_names, 
               title, 
               file_name,
               names,
               colours,
               lty_list,
               pch_list)
  
  #Return values
  return(model_output)
}


#Build models with chosen estimation method
perform_estimation <- function(rstan_file, rstan_data, est_names, est_method, alg,
                               int_methods, estplot){
  model <- stan_model(file=rstan_file)
  model_list <- list()
  
  if(est_method == "samp"){
    if(1 %in% int_methods){
      model_euler <- sampling(model, data=rstan_data[[1]], chains=3)
      model_list$model_euler <- model_euler
    }
    if(2 %in% int_methods){
      model_midpoint <- sampling(model, data=rstan_data[[2]], chains=3)
      model_list$model_midpoint <- model_midpoint
      }
    if(3 %in% int_methods){
      model_rk4 <- sampling(model, data=rstan_data[[3]], chains=3)
      model_list$model_rk4 <- model_rk4
      }
    
  } else if(est_method == "opt"){
    if(1 %in% int_methods){
      model_euler <- optimizing(model, data=rstan_data[[1]], algorithm=alg)
      model_list$model_euler <- model_euler
    }
    if(2 %in% int_methods){
      model_midpoint <- optimizing(model, data=rstan_data[[2]], algorithm=alg)
      model_list$model_midpoint <- model_midpoint
    }
    if(3 %in% int_methods){
      model_rk4 <- optimizing(model, data=rstan_data[[3]], algorithm=alg)
      model_list$model_midpoint <- model_midpoint
    }
  }
  
  #Produces side-by-side plot of beta samples for each model
  if(estplot){ produce_estplot(model_list, est_names) }
  
  #Sample wrangling
  estimates <- est_wrangling(model_list, est_names, est_method)
  
  return(list(model_list = model_list, 
              estimates = estimates, 
              est_method = est_method))
}
