# Code to extract results from simulations and corresponding model fits. 
# Assumes that the saved fits are stored in a directory called input
# with specific naming structure.

#Load in model fits, extract info, export to data structures
{
  library(tidyverse)
  library(lubridate)
  library(rstan)
  library(grid)
  library(ggridges)
  library(cowplot)
}

structure <- "25Yr_ObsPeriod_6Obs"
error_type <- "Ruger" #"Norm"
batches <- 1100

#simulated data available
sim_data_list <- readRDS("data/50ind_25Yr_ObsPeriod_6Obs_CanhamSim.rds")

sp_pars <- c("species_max_growth_mean",
             "species_max_growth_sd",
             "species_diameter_at_max_growth_mean",
             "species_diameter_at_max_growth_sd",
             "species_K_mean",
             "species_K_sd")
sp_pars_true_vals <- c(mean(log(sim_data_list$sim_ind_info$g_max)),
                       sd(log(sim_data_list$sim_ind_info$g_max)),
                       mean(log(sim_data_list$sim_ind_info$s_max)),
                       sd(log(sim_data_list$sim_ind_info$s_max)),
                       mean(log(sim_data_list$sim_ind_info$k)),
                       sd(log(sim_data_list$sim_ind_info$k)))
ind_par_names <- c("ind_max_growth", "ind_diameter_at_max_growth", "ind_K")
growth_function_pars <- c("g_max", "s_max", "k")

#Diagnostics
diag_data <- read.csv(paste0("input/",error_type,"Diag.csv"))
diag_data_reduced <- diag_data %>%
  filter(Rhat > 1.05, DivTrans >50)

#Exclusion lists come from the diag_data_reduced frame with double checking
#exclude <- c(192, 227, 85, 888, 910) #Norm
#exclude <- c(1013, 523, 655, 657) #Ruger

load_data <- FALSE
if(load_data){
  #Initialise data structures
  error_reduction_data <- tibble()
  ind_data <- tibble()
  species_data <- tibble()
  error_par_data <- tibble()
  sim_data <- tibble()
  
  for(i in 1:batches){ #Iterate through simulations
    if((i-1) %% 25 == 0){
      print(paste0("Model: ", i))
    }
    
    filename <- paste0("input/fits/Fit_", error_type, "_", i, "_1_25Yr_ObsPeriod_6Obs_FullModel.rds")
    #Load model if it exists, otherwise create dummy row
    if(!file.exists(filename) || (i %in%exclude)){
      temp_sim_data <- tibble(
        batch = i,
        start_time = 0,
        end_time = 0,
        runtime = -1
      )
      
      sim_data <- rbind(sim_data, temp_sim_data)
    } else {
      model_save <- readRDS(filename)
      sim_data_list <- readRDS("data/50ind_25Yr_ObsPeriod_6Obs_CanhamSim.rds")
      
      if(length(model_save$error_save_data)>0){
        sim_data_list$sim_data$S_obs <- model_save$error_save_data$S_obs
        sim_data_list$sim_data <- sim_data_list$sim_data %>%
          group_by(treeid) %>%
          arrange(time) %>%
          mutate(delta_S_obs = lead(S_obs) - S_obs) %>%
          ungroup() %>%
          arrange(treeid, time)
        
        #Extract samples
        est_data <- rstan::extract(model_save$fit, permuted=TRUE)
        
        #Extract error parameter data
        error_par_temp <- tibble(
          batch = i,
          sigma_e_mean = mean(est_data$global_error_sigma),
          sigma_e_median = median(est_data$global_error_sigma),
          sigma_e_CI_lower = quantile(est_data$global_error_sigma, 0.025),
          sigma_e_CI_upper = quantile(est_data$global_error_sigma, 0.975)
        )
        error_par_data <- rbind(error_par_data, error_par_temp)
        
        #Extract size and growth measurement error reduction
        size_data <- tibble(
          true_size = sim_data_list$sim_data$S_true,
          obs_size = sim_data_list$sim_data$S_obs,
          ind_fit = apply(est_data$S_hat, 2, mean)
        )
        
        size_data <- size_data %>%
          mutate(
            obs_diff = obs_size - true_size,
            ind_diff = ind_fit - true_size,
          )
        
        growth_data <- tibble(
          true_growth = sim_data_list$sim_data$delta_S,
          obs_growth = sim_data_list$sim_data$delta_S_obs,
          ind_fit = apply(est_data$G_hat, 2, mean)
        )
        
        growth_data <- growth_data %>%
          mutate(
            obs_diff = obs_growth - true_growth,
            ind_diff = ind_fit - true_growth
          ) %>%
          drop_na()
        
        temp_error_reduction_data <- tibble(
          batch = i,
          error_type = error_type,
          obs_size_RMSE = sqrt(sum(size_data$obs_diff^2)/nrow(size_data)),
          ind_size_RMSE = sqrt(sum(size_data$ind_diff^2)/nrow(size_data)),
          obs_growth_RMSE = sqrt(sum(growth_data$obs_diff^2)/nrow(growth_data)),
          ind_growth_RMSE = sqrt(sum(growth_data$ind_diff^2)/nrow(growth_data)),
          size_RMSE_percent_reduction = (1-(ind_size_RMSE/obs_size_RMSE))*100,
          growth_RMSE_percent_reduction = (1-(ind_growth_RMSE/obs_growth_RMSE))*100
        )
        
        error_reduction_data <- rbind(error_reduction_data, temp_error_reduction_data)
        
        #Extract individual parameter estimates
        for(j in 1:length(ind_par_names)){
          #Get percentage of ind_par CIs that contains true value, par estimate
          ind_CI_contain_true <- c()
          for(k in 1:50){
            CI_lower <- as.numeric(quantile(est_data[[ind_par_names[j]]][,k], probs=0.025))
            CI_upper <- as.numeric(quantile(est_data[[ind_par_names[j]]][,k], probs=0.975))
            true_val <- sim_data_list$sim_ind_info[[growth_function_pars[j]]][k]
            
            ind_CI_contain_true[k] <- (true_val >= CI_lower) && (true_val <= CI_upper)
            
            ind_par_temp <- tibble(batch = i,
                                   treeid = k,
                                   par_name = ind_par_names[j],
                                   true_val = true_val,
                                   est = mean(est_data[[ind_par_names[j]]][,k]),
                                   ci_lower = CI_lower,
                                   ci_upper = CI_upper)
            ind_par_temp$ind_ci_contain_true <- ind_CI_contain_true[k]
            ind_par_temp$error <- ind_par_temp$est - ind_par_temp$true_val
            ind_data <- rbind(ind_data, ind_par_temp)
          }
        }
        
        #Extract species-level estimates
        for(j in 1:length(sp_pars)){
          temp_species_data <- tibble(
            batch = i,
            error_type = error_type,
            par_name = sp_pars[j],
            par_true_val = sp_pars_true_vals[j],
            par_est = mean(est_data[[sp_pars[j]]]),
            par_ci_lower = as.numeric(quantile(est_data[[sp_pars[j]]], probs=0.025)),
            par_ci_upper = as.numeric(quantile(est_data[[sp_pars[j]]], probs=0.975))
          )
          temp_species_data <- temp_species_data %>%
            mutate(ci_contain_true = 
                     ((par_true_val > par_ci_lower) && 
                        (par_true_val < par_ci_upper)))
          species_data <- rbind(species_data, temp_species_data)
        }
        
        #Extract sim-level info
        sim_data_temp <- tibble(batch = i,
                                start_time = model_save$model_start_time,
                                end_time = model_save$model_end_time) %>%
          mutate(runtime = difftime(end_time, start_time, units = "mins"))
        sim_data <- rbind(sim_data, sim_data_temp)
      }
    }
  }
  
  save_data <- list(
    error_reduction_data = error_reduction_data,
    ind_data = ind_data,
    species_data = species_data,
    sim_data = sim_data,
    error_par_data = error_par_data
  )
  
  saveRDS(save_data, file= paste0("output/data/", error_type, "_summary_data.rds"))
}

#-----------------------------------------------------------------------------#
# Start here if data already processed
#-----------------------------------------------------------------------------#
#load in saved data
save_data <- readRDS(file= paste0("output/data/", error_type, "_summary_data.rds"))
#Initialise data structures
error_reduction_data <- save_data$error_reduction_data
ind_data <- save_data$ind_data
species_data <- save_data$species_data
sim_data <- save_data$sim_data
sim_data_filtered <- sim_data %>%
  filter(runtime > 0)

mean(sim_data_filtered$runtime)
hist(error_reduction_data$size_RMSE_percent_reduction, main = "", col = "lightblue",
     xlab = "RMSE % reduction for size", breaks = 20)
hist(error_reduction_data$growth_RMSE_percent_reduction, main = "", col = "lightgreen",
     xlab = "RMSE % reduction for growth", breaks = 20)

mean(error_par_data$sigma_e_mean) #Mean estimate
error_par_data <- error_par_data %>% #Construct CIs for N(0,0.1)
  mutate(true_geq_lower = (0.1 >= sigma_e_CI_lower),
         true_leq_upper = (0.1 <= sigma_e_CI_upper))
error_par_data$true_in_CI <- error_par_data$true_geq_lower & error_par_data$true_leq_upper
sum(error_par_data$true_in_CI)/ nrow(error_par_data)

#Extract information on how much error in size and growth was reduced
error_reduction_summary <- error_reduction_data %>%
  summarise(meanSizeRed = mean(size_RMSE_percent_reduction),
            meanGrowthRed = mean(growth_RMSE_percent_reduction),
            meanSizeRMSEobs = mean(obs_size_RMSE),
            meanSizeRMSEfit =mean(ind_size_RMSE),
            meanGrowthRMSEobs = mean(obs_growth_RMSE),
            meanGrowthRMSEfit = mean(ind_growth_RMSE))

ind_data_summary <- ind_data %>%
  group_by(par_name) %>%
  summarise(CIs_contain_true = sum(ind_ci_contain_true),
            mean_CI_width = mean(ci_upper - ci_lower),
            RMSE = sqrt(sum((true_val - est)^2)/nrow(sim_data))) %>%
  mutate(CI_percent = (CIs_contain_true/(50*nrow(sim_data_filtered))) * 100)

#plots for individual data
ind_par_names <- c("ind_max_growth", "ind_diameter_at_max_growth", "ind_K")
ind_par_names_fancy <- c("g_max", "S_max", "k")
plots <- list()
for(i in 1:3){
  temp <- ind_data %>%
    filter(par_name == ind_par_names[i]) %>%
    mutate(CI_width = ci_upper - ci_lower,
           error = est - true_val)
  
  plots[[i]] <- ggplot(data=temp, aes(x = true_val,y = error))+
    geom_point(size=2, colour = "green4", fill="#1EB300", alpha = 0.1) +
    labs(x = paste("True Value: ", ind_par_names_fancy[i], sep=""), 
         y = "Error (est - true)") +
    geom_abline(intercept=0, slope=0, linetype="dashed", linewidth = 1, alpha=0.5) +
    scale_x_log10() +
    #scale_y_log10() +
    theme_classic()
}
ind_par_plot_grid <- plot_grid(plots[[1]], plots[[2]], plots[[3]], 
                               align = "h", nrow = 1, labels = c("(d)", "(e)", "(f)"))

#Extract info on species-level parameter performance
species_data_summary <- species_data %>%
  group_by(par_name) %>%
  summarise(
    true_val = mean(par_true_val),
    mean_est = mean(par_est),
    per_contain_true = (sum(ci_contain_true)/nrow(sim_data))*100,
    mean_ci_width = mean(par_ci_upper - par_ci_lower),
    RMSE = sqrt(sum((par_true_val - par_est)^2)/nrow(sim_data))
  )

#plots for species-level
species_data <- species_data %>%
  group_by(par_name) %>%
  mutate(
    Fancy_name = case_when(
      par_name == "species_max_growth_mean" ~ "ln(g_max) Mean",
      par_name == "species_max_growth_sd" ~ "ln(g_max) SD",
      par_name == "species_diameter_at_max_growth_mean" ~ "ln(S_max) Mean",
      par_name == "species_diameter_at_max_growth_sd" ~ "ln(S_max) SD",
      par_name == "species_K_mean" ~ "ln(k) Mean",
      par_name == "species_K_sd" ~ "ln(k) SD",
      .default = "Whoops"
  )) %>%
  mutate(error = par_est - par_true_val) %>%
  ungroup()


#create density ridge plot
ggplot(species_data, aes(x = error, y = Fancy_name)) +
  geom_density_ridges(aes(fill=Fancy_name), alpha = 0.8) +
  geom_vline(xintercept=0, linetype = "dashed", color = "black", linewidth = 1) +
  ylab("Parameter") +
  xlab("Error (est - true)") +
  theme_classic() +
  theme(legend.position = "none")


#------------------------------------------------------------------------------#
#Analysis of least growth individual
growth_function <- function(x, pars){
  growth <- pars[[1]] * 
    exp(-0.5 * (log(x / pars[[2]]) / pars[[3]])^2 )
  return(growth)
}

least_growth_ind <- ind_data %>%
  filter(treeid == 35)

hist_data1 <- least_growth_ind %>%
  filter(par_name == "ind_max_growth") %>%
  mutate(err = est - true_val)
hist_gmax <- ggplot(data = hist_data1, aes(err)) +
  geom_histogram(fill = "#72b000",
                 color = "black") +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed") +
  labs(x = "Max growth rate g_max error",
       y = "Count") +
  theme_classic()

hist_data2 <- least_growth_ind %>%
  filter(par_name == "ind_diameter_at_max_growth") %>%
  mutate(err = est - true_val)
hist_Smax <- ggplot(data = hist_data2, aes(err)) +
  geom_histogram(fill = "#72b000",
                 color = "black") +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed") +
  labs(x = "DBH at max growth rate S_max error",
       y = "Count") +
  theme_classic()

hist_data3 <- least_growth_ind %>%
  filter(par_name == "ind_K") %>%
  mutate(err = est - true_val)
hist_k <- ggplot(data = hist_data3, aes(err)) +
  geom_histogram(fill = "#72b000",
                 color = "black") +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dashed") +
  labs(x = "Spread K error",
       y = "Count") +
  theme_classic()

rstan_data <- readRDS("data/50ind_25Yr_ObsPeriod_6Obs_CanhamSim.rds")
ind_35_data <- rstan_data$sim_data %>%
  filter(treeid == 35)
ind_35_pars <- list(pars = c(ind_35_data$g_max[1], 
                             ind_35_data$s_max[1], 
                             ind_35_data$k[1]))

S_0 <- ind_35_data$S_0[1]
S_final <- ind_35_data$S_final[1]

ind35_growth <- ggplot() +
  geom_function(fun=growth_function, args=ind_35_pars, alpha=0.5, 
                color="#72b000", linewidth=1, xlim=c(1, max(rstan_data$sim_data$S_true)),
                linetype = "solid") +
  geom_function(fun=growth_function, args=ind_35_pars, alpha=1, 
                color="#72b000", linewidth=2.5, xlim=c(S_0, S_final),
                linetype = "solid") +
  labs(x = "Size (cm)", y = "Growth rate (cm/yr)") +
  theme_classic()

hist_grid <- plot_grid(hist_gmax, hist_Smax, hist_k, align = "hv", ncol = 1)
least_growth_plot <- plot_grid(ind35_growth, hist_grid, align = "h", nrow = 1, 
                               rel_widths = c(0.6, 0.4))

#Take sample of posterior parameter estimates to plot growth functions
sim_ids <- sample(sim_data_filtered$batch, size=20)

ind35_growth_fns <- ggplot() 
for(i in 1:20){
  ests <- least_growth_ind %>% 
    filter(batch == sim_ids[i])
  temp_pars <- list(pars = c(ests$est[which(ests$par_name == "ind_max_growth")], 
                             ests$est[which(ests$par_name == "ind_diameter_at_max_growth")], 
                             ests$est[which(ests$par_name == "ind_K")]))
  
  ind35_growth_fns <- ind35_growth_fns +
    geom_function(fun=growth_function, args=temp_pars, alpha=0.2, 
                  color="#000000", linewidth=0.5, xlim=c(1, max(rstan_data$sim_data$S_true)),
                  linetype = "solid")
}

ind35_growth_fns <- ind35_growth_fns +
  geom_function(fun=growth_function, args=ind_35_pars, alpha=0.5, 
                color="#72b000", linewidth=1, xlim=c(1, max(rstan_data$sim_data$S_true)),
                linetype = "solid") +
  geom_function(fun=growth_function, args=ind_35_pars, alpha=1, 
                color="#72b000", linewidth=2.5, xlim=c(S_0, S_final),
                linetype = "solid") +
  labs(x = "Size (cm)", y = "Growth rate (cm/yr)") +
  theme_classic()

ind35_growth_fns_grid <- plot_grid(ind35_growth_fns)

least_growth_plot_fns <- plot_grid(ind35_growth_fns_grid, hist_grid, align = "h", nrow = 1, 
                                   rel_widths = c(0.6, 0.4))
