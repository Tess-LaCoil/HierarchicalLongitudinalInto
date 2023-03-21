#Top-level management function
model_analysis <- function(analysis_controls){
  for(i in analysis_controls$model_list){ #For each model, do the data wrangling
    filename <- "output/data/Single_species_compiled_data"
    
    filename <- paste(filename,
                      i$model_name,
                      analysis_controls$int_method, sep="_")
    
    if(data_extract_controls$est_method == "samp"){
      filename <- paste(filename, "_sampling.rds", sep="")
    } else if(data_extract_controls$est_method == "opt"){
      filename <- paste(filename, "_optimizing.rds", sep="")
    }
    
    #Load in required data
    if(file.exists(filename)){
      data <- readRDS(filename)
    } else {
      print("Compiled data not found, please wrangle samples.")
      return()
    }
    
    run_analysis(data, i, analysis_controls)
  }
}

#Produce plots for each model
run_analysis <- function(data, model_list_item, analysis_controls){
  source(model_list_item$model_spec, local=TRUE)
  
  #Scatter plots
  for(j in scatterplot_list){
    build_scatterplot(data, j, model_list_item$model_name)
  }
  
  #Parameter plots
  for(j in parplot_list){
    build_beeswarm_boxplot(data$individual_data, j, model_list_item$model_name)
  }
}

#~~~~~~   Scatter plot   ~~~~~~#
#Build data for scatterplots
build_scatterplot <- function(data, scatterplot_list_item, model_name){
  #Size plot
  if(scatterplot_list_item$par == "S_hat"){
    #Get model values
    plot_data <- data.frame(x = data$measurement_data$S_true, 
                            y = data$measurement_data$S_hat,
                            Model = rep("Bayesian Model", 
                                       times=length(data$measurement_data$S_true)))
    
    #Attach values without model
    plot_data <- rbind(data.frame(x = data$measurement_data$S_true, 
                                  y = data$measurement_data$S_obs,
                                  Model = rep("Observed Size", 
                                              times=length(data$measurement_data$S_true))
                                  ),
                       plot_data
                       )
    
    plot_data <- plot_data %>% na.omit() #Remove NAs
    
    #Build spec object
    filename <- paste("output/figures/Scatter_Size_", model_name, ".png", sep="")
    plot_spec <- list(xlab="True size cm",
                      ylab="Estimated size cm",
                      filename = filename,
                      log_log = scatterplot_list_item$log_log)
    
    #Call plotting function
    plot_scatterplot(plot_data, scatterplot_list_item, plot_spec)
    
  #Growth plot
  } else if(scatterplot_list_item$par == "G_hat"){
    #Build dataset
    plot_data <- data.frame(x = data$measurement_data$delta_S, 
                            y = data$measurement_data$G_hat,
                            Model = rep("Bayesian Model", 
                                        times=length(data$measurement_data$delta_S)))
    plot_data <- rbind(data.frame(x = data$measurement_data$delta_S, 
                                  y = data$measurement_data$delta_S_obs,
                                  Model = rep("Pairwise Difference", 
                                              times=length(data$measurement_data$S_true))
                       ),
                       plot_data
    )
    
    plot_data <- plot_data %>% na.omit() #Remove NAs
    
    #Build spec object
    filename <- paste("output/figures/Scatter_Growth", model_name, ".png", sep="")
    plot_spec <- list(xlab="True growth cm/yr",
                      ylab="Estimated growth cm/yr",
                      filename = filename,
                      log_log = scatterplot_list_item$log_log)
    
    #Call plotting function
    plot_scatterplot(plot_data, scatterplot_list_item, plot_spec)
  
  #Other parameters
  } else {
    plot_data <- data.frame(x = data$individual_data[[paste("true", scatterplot_list_item$par, sep="_")]], 
                            y = data$individual_data[[scatterplot_list_item$par]],
                            Model = rep("Bayesian", 
                                        times=length(data$measurement_data$delta_S)))
    plot_data <- plot_data %>% na.omit() #Remove NAs
    
    #Build spec object
    filename <- paste("output/figures/Scatter_", scatterplot_list_item$par_name, "_", model_name, ".png", sep="")
    plot_spec <- list(xlab=paste("True", scatterplot_list_item$par_name, sep=" "),
                      ylab=paste("Est.", scatterplot_list_item$par_name, sep=" "),
                      filename = filename,
                      log_log = scatterplot_list_item$log_log)
    
    plot_scatterplot(plot_data, scatterplot_list_item, plot_spec)
  }
}

#Save scatterplot to file
plot_scatterplot <- function(plot_data, scatterplot_list_item, plot_spec){
  plot <- ggplot_scatterplot(plot_data, plot_spec)
  ggsave(plot_spec$filename, plot=plot, scale=1, width=100, height=100, units="mm")
}

#Function for producing plot
ggplot_scatterplot <- function(plot_data, plot_spec){
  plot <- ggplot(data=plot_data, aes(x,y, col=Model, alpha=Model))+
    geom_point(aes(x=x, y=y, shape=Model), size=2) +
    scale_shape_manual(values = c(16, 3)) +
    scale_color_manual(values = c("#4A4A4A", "green4")) +
    scale_alpha_discrete(range = c(0.2, 1)) +
    labs(x = plot_spec$xlab, y = plot_spec$ylab) +
    geom_abline(intercept=0, slope=1, linetype="dashed") +
    theme_classic() +
    theme(legend.position = c(0.3,0.9))
  
  if(min(plot_data$y) < 0){
    plot <- plot +
      geom_hline(yintercept=0)
  }
  
  if(plot_spec$log_log){
    plot <- plot +
      scale_x_continuous(trans = "log10") +
      scale_y_continuous(trans = "log10")
  }
  
  return(plot)
}

#~~~~~~   Parameter plots   ~~~~~~#
build_beeswarm_boxplot <- function(individual_data, parplot_list_element, model_name){
  #Build data frame
  plot_data <- data.frame()
  
  for(j in parplot_list_element$pars){
    plot_data_temp <- data.frame(par = rep(j$par_name, 
                                           times=length(individual_data[[j$par]])),
                                 val = individual_data[[j$par]])
    plot_data <- rbind(plot_data, plot_data_temp)
  }
  
  #Beeswarm plot
  filename <- paste("output/figures/",
                    model_name,
                    parplot_list_element$plot_name,
                    "_Beeswarm.png", 
                    sep="")
  beeswarm_plot <- ggplot_parameter_beeswarm(plot_data, 
                                             parplot_list_element$xlab, 
                                             parplot_list_element$ylab,
                                             parplot_list_element$hline_height)
  ggsave(filename, plot=beeswarm_plot, scale=1, width=(60*length(parplot_list_element$pars)), 
         height=100, units="mm")
  
  #Box plot
  filename <- paste("output/figures/",
                    model_name,
                    parplot_list_element$plot_name,
                    "_Boxplot.png", 
                    sep="")
  box_plot <- ggplot_parameter_boxplot(plot_data,
                                       parplot_list_element$xlab, 
                                       parplot_list_element$ylab,
                                       parplot_list_element$hline_height)
  ggsave(filename, plot=box_plot, scale=1, width=(60*length(parplot_list_element$pars)), 
         height=100, units="mm")
}

#~~~~~~   Beeswarm plot   ~~~~~~#
#Produces beeswarm plot object for build_beeswarmplot
ggplot_parameter_beeswarm <- function(plot_data, xlab, ylab, hline_height){
  plot <- ggplot(data=plot_data, aes(x=par, y=val)) +
    geom_quasirandom(alpha=0.5, groupOnX=TRUE) +
    xlab(xlab) +
    ylab(ylab) +
    geom_hline(yintercept=hline_height, linetype="dashed") + 
    theme_classic()
  
  return(plot)
}

#~~~~~~   Box plot   ~~~~~~#
ggplot_parameter_boxplot <- function(plot_data,
                                     xlab, 
                                     ylab,
                                     hline_height){
  plot <- ggplot(plot_data, aes(x=par, y=val))+
    geom_boxplot() +
    xlab(xlab) +
    ylab(ylab) +
    geom_hline(yintercept=hline_height, linetype="dashed") +
    theme_classic()
  
  return(plot)
}

#~~~~~~   Time series plots   ~~~~~~#
#Plots a sample of individuals to show how the model estimates behave against the observed size
plot_obs_and_est_life_history <-function(plotting_data, n_ind, model_name){
  for(j in 1:n_ind){
    sample_id <- unique(plotting_data$treeid)[j] #sample(plotting_data$treeid, size=2)[1]
    
    sample <- plotting_data %>%
      filter(treeid == sample_id) %>%
      select(S_obs, time, treeid, S_hat)
    
    
    #Build data frame
    data <- data.frame(size=sample$S_obs, 
                       time=sample$time, 
                       cond=rep("Observed", times=length(sample$time)))
    data <- rbind (data, data.frame(size=sample$S_hat,
                                    time=sample$time, 
                                    cond=rep(model_name, times=length(sample$time))))
    
    #Produce plot
    file_name <- paste("output/figures/sampled/", model_name,
                       "_Actual_And_Est_Sample_dot_", 
                       j,".png", sep="")
    
    plot <- ggplot_obs_and_est_life_history(data, title = "")
    
    ggsave(file_name, plot=plot, width=130, height=100, units="mm")
  }
}

#Produces plot for plot_obs_and_est_life_history()
ggplot_obs_and_est_life_history<- function(data, title){
  plot <- ggplot(data=data, aes(x, y)) +
    geom_line(aes(x=(time+1990), y=size, color=as.factor(cond),
                  group=cond, linetype=as.factor(cond)), size=1.1) +
    geom_point(aes(x=(time+1990), y=size, color=as.factor(cond), 
                   shape=as.factor(cond),
                   group=cond), size=2.5) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_color_manual(values = c("black", "green4")) +
    xlab("Years") +
    ylab("Size (cm)") +
    ggtitle(title) +
    labs(color = NULL) +
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  
  return(plot)
}

#Plots individuals to show how the model estimates behave
plot_est_life_history <-function(measurement_data, model_name){
  #Build data frames
  size_data <- data.frame(y=measurement_data$S_hat,
                          x=measurement_data$time,
                          treeid=measurement_data$treeid)
  
  growth_data <- data.frame(y=measurement_data$G_hat,
                          x=measurement_data$S_hat,
                          treeid=measurement_data$treeid)
  
  #Produce plots
  file_name_size <- paste("output/figures/", model_name,
                     "_EstimatedSizes.png", sep="")
  size_plot <- ggplot_est_life_history(size_data, x_lab = "Years", y_lab = "Size",
                                       logAxis="y")
  ggsave(file_name_size, plot=size_plot, width=180, height=100, units="mm")
  
  file_name_growth <- paste("output/figures/", model_name,
                          "_EstimatedGrowth.png", sep="")
  growth_plot <- ggplot_est_life_history(growth_data, x_lab = "Size", y_lab = "Growth",
                                         logAxis="x")
  ggsave(file_name_growth, plot=growth_plot, width=180, height=100, units="mm")
}

#Produces plot for plot_obs_and_est_life_history()
ggplot_est_life_history <- function(data, x_lab, y_lab, logAxis = NA){
  plot <- ggplot(data=data, aes(x=x, y=y)) +
    geom_line(aes(group=treeid), color="green4", alpha=0.1) +
    geom_point(color="green4", size=1.5, alpha=0.1) +
    xlab(x_lab) +
    ylab(y_lab) +
    theme_classic() +
    theme(legend.position = "none", axis.text=element_text(size=12))
  
  if("x" %in% logAxis){
    plot <- plot + scale_x_continuous(trans='log10')
  }
  if("y" %in% logAxis){
    plot <- plot + scale_y_continuous(trans='log10')
  }
  
  return(plot)
}

#~~~~~~   Lifetime size   ~~~~~~#
#Plots the predicted life of an individual based on mean and median growth parameters
plot_pair_life_history <-function(pop_pars_est, ind_pars_est, growth_function,
                             name, S_0, step_size, lifespan){
  #Build predicted data
  time <- seq(from=1, to=lifespan, length.out=(lifespan/step_size))
  
  #Build individuals for 5%, mean, 95% of species parameter dist.
  individual_pop <- data.frame(S_true = rk4_est(S_0,
                                                 growth = growth_function,
                                                 pars = pop_pars_est,
                                                 step_size,
                                                 N_step=lifespan/step_size),
                                time = time,
                                cond = rep("Population", times = (lifespan/step_size))
  )
  individual_ind <- data.frame(S_true = rk4_est(S_0,
                                                growth = growth_function,
                                                pars = ind_pars_est,
                                                step_size,
                                                N_step=lifespan/step_size),
                               time = time,
                               cond = rep("Individual", times = (lifespan/step_size))
  )
  avg_plot_data <- rbind(individual_pop, individual_ind)

  #Produce plot
  file_name <- paste("output/figures/", name,
                     "_Avg_Life_History.png", sep="")
  
  plot <- ggplot_pair_life_history(avg_plot_data)
  
  ggsave(file_name, plot=plot, width=150, height=100, units="mm")
}

#Produces plot for plot_pair_life_history()
ggplot_pair_life_history<- function(plot_data){
  plot <- ggplot(data=plot_data, aes(x=time, y=S_true, group=cond)) +
    geom_line(aes(linetype=cond), size=0.8) +
    scale_linetype_manual(values=c("longdash", "solid")) +
    xlab("Years") +
    ylab("Size (cm)") +
    labs(color = NULL) +
    theme_classic() +
    theme(legend.position=c(0.2,0.8), legend.title=element_blank())
  
  return(plot)
}

#Plot multiple life histories based on data frame of parameters
plot_sample_life_history <- function(post_pars,
                                     growth_function,
                                     name = "BCIAvg", 
                                     n_ind, S_0, 
                                     S_final = NA,
                                     step_size, lifespan,
                                     min_growth_size = 0,
                                     max_growth_size = 50){
  time <- seq(from=0, to=lifespan, length.out=(lifespan/step_size))
  plot_data <- data.frame(S_true=c(), time=c(), cond=c())
  
  for(i in 1:n_ind){
    temp <- data.frame(S_true = rk4_est(S_0 = 3,
                                        growth = growth_function,
                                        pars = post_pars[i,],
                                        step_size,
                                        N_step=lifespan/step_size),
                       time = time,
                       individual = rep(i, times = (lifespan/step_size))
    )
    plot_data <- rbind(temp, plot_data)
  }
  
  
  #Produce plot of sizes over time
  file_name <- paste("output/figures/", name,
                     "_Sampled_Life_History.png", sep="")
  plot <- ggplot_sample_life_history(plot_data)
  ggsave(file_name, plot=plot, width=150, height=150, units="mm")
  
  
  #Produce plot of growth rates by size
  file_name <- paste("output/figures/", name,
                     "_Sampled_Growth_History.png", sep="")
  if(is.na(S_final[1])){S_final <- rep(max_growth_size, times=length(S_0))}
  plot <- ggplot_sample_growth_trajectories(post_pars, growth_function,  
                                            max_growth_size, min_growth_size, 
                                            S_0, S_final)
  ggsave(file_name, plot=plot, width=150, height=150, units="mm")
  
}

ggplot_sample_life_history <- function(plot_data){
  plot <- ggplot(data=plot_data, aes(x=time, y=S_true, group=as.factor(individual))) +
    geom_line(alpha=0.2, color="#003300", size=1) +
    xlab("Years") +
    ylab("Size (cm)") +
    labs(color = NULL) +
    theme_classic() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18,face="bold"))
  
  return(plot)
}

ggplot_sample_growth_trajectories <- function(post_pars, growth_function, 
                                              max_growth_size, min_growth_size,
                                              S_0, S_final){
  plot <- ggplot() +
    xlim(min_growth_size, min(c(max_growth_size, max(S_final)))) +
    xlab("Size") +
    ylab("Growth rate") +
    theme_classic() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18,face="bold"))
  
  for(i in 1:nrow(post_pars)){
    args_list <- list(pars=post_pars[i,])
    plot <- plot +
      geom_function(fun=growth_function, args=args_list, alpha=0.2, 
                    color="#003300", size=1, xlim=c(S_0[i], S_final[i]))
  }
  
  return(plot)
}
