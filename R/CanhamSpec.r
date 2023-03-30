#Control file for Canham growth model
#Build list for diagnostics
diagnostic_list_single_species <- list(
  #Global
  plot4=list(
    pars=c("global_error_sigma"),
    name="canham_global_error_sigma_single_species"
  ),
  
  #Species
  plot5=list(
    pars=c("species_max_growth_mean",
           "species_max_growth_sd",
           "species_diameter_at_max_growth_mean",
           "species_diameter_at_max_growth_sd",
           "species_K_mean",
           "species_K_sd"),
    name="canham_pars_single_species"
  )
)

diagnostic_list_single_individual <- list(
  #Global
  plot1=list(
    pars=c("global_error_sigma"),
    name="canham_global_error_sigma_single_individual"
  ),
  
  #Species
  plot2=list(
    pars=c("ind_max_growth", "ind_diameter_at_max_growth", "ind_K"),
    name="canham_ind_single_individual"
  )
)

#List of controls for single-species model
single_species_controls <- list(
  model_name = "Single species Canham model",
  stan_file_path = "stan/Canham_SingleSpecies.stan",
  fit_save_path = "output/data/Single_Species_Fit_Canham",
  stan_controls = list(
    iter = 2000,
    chains = 4,
    cores = 4,
    control = list(max_treedepth = 20)
  ),
  diagnostic_list = diagnostic_list_single_species
)

#List of controls for single-individual model
single_individual_controls <- list(
  model_name = "Single individual Canham model",
  stan_file_path = "stan/Canham_SingleIndividual.stan",
  fit_save_path = "output/data/Single_Ind_Fit_Canham",
  stan_controls = list(
    iter = 2000,
    chains = 4,
    cores = 1,
    control = list(max_treedepth = 15)
  ),
  diagnostic_list = diagnostic_list_single_individual
)

#Individual-level parameter controls
ind_pars <- c("ind_S_0", 
              "ind_max_growth", 
              "ind_diameter_at_max_growth",
              "ind_K")
ind_pars_names <- c("canham_S_0_hat", 
                    "canham_max_growth_hat", 
                    "canham_diameter_at_max_growth_hat",
                    "canham_K_hat") #Param names with model
mean_size_name <- "canham_S_hat_mean"

#Species-level parameter controls
sp_pars <- c("species_max_growth_mean", 
             "species_max_growth_sd", 
             "species_diameter_at_max_growth_mean", 
             "species_diameter_at_max_growth_sd", 
             "species_K_mean",
             "species_K_sd")
sp_pars_names <- c("canham_species_max_growth_mean", 
                   "canham_species_max_growth_sd", 
                   "canham_species_diameter_at_max_growth_mean",
                   "canham_species_diameter_at_max_growth_sd",
                   "canham_species_K_mean",
                   "canham_species_K_sd")

#Plot parameter controls
growth_pars <- c("ind_max_growth", 
                 "ind_diameter_at_max_growth",
                 "ind_K")
growth_par_names <- c("canham_max_growth_hat", 
                      "canham_diameter_at_max_growth_hat",
                      "canham_K_hat")
plot_par_names <- c("Max Growth", "Diameter at Max Growth", "K")
size_name <- "canham_S_hat"
growth_name <- "canham_G_hat"
cond <- "Canham Pred"

#Lists for plot functions
scatterplot_list <- list(
  size_plot = list(
    par = "S_hat",
    name = "S_hat",
    par_name = "Estimated Size",
    log_log = FALSE
  ),
  
  growth_plot = list(
    par = "G_hat",
    name = "G_hat",
    par_name = "Estimated Growth",
    log_log = FALSE
  ),
  
  canham_plot1 = list(
    par = "canham_max_growth_hat",
    name= "canham_max_growth",
    par_name= "Max Growth",
    log_log = FALSE
  ),
  canham_plot2 = list(
    par = "canham_diameter_at_max_growth_hat",
    name= "canham_diameter_at_max_growth",
    par_name= "Diameter at Max Growth",
    log_log = FALSE
  ),
  canham_plot3 = list(
    par = "canham_K_hat",
    name= "canham_K",
    par_name= "K",
    log_log = FALSE
  )
)

#List of plots for beeswarm/box
parplot_list <- list(
  canham_plot1 = list(
    plot_name = "_canham_difference_smax_",
    pars = list(
      par_1 = list(
        par = "difference_canham_diameter_at_max_growth_hat",
        par_name= "Size at Max Growth"
      )
    ),
    xlab="Parameter",
    ylab="Difference true - est",
    hline_height=0
  ),
  
  canham_plot1 = list(
    plot_name = "_canham_difference_gmax_k_",
    pars = list(
      par_1 = list(
        par = "difference_canham_max_growth_hat",
        par_name= "Max Growth"
      ),

      par_2 = list(
        par = "difference_canham_K_hat",
        par_name= "k"
      )
    ),
    xlab="Parameter",
    ylab="Difference true - est",
    hline_height=0
  )
)

#Define growth function
growth_function <- function(x, pars){
  growth <- pars[[1]] * 
    exp(-0.5 * (log(x / pars[[2]]) / pars[[3]])^2 )
  return(growth)
}

growth_function_pars <- c("g_max", "s_max", "k")

sim_spec <- function(n){
  tibble(
    g_max = exp(rnorm(n, mean=0.5, sd=0.1)),
    s_max = exp(rnorm(n, mean=2.5, sd=0.1)),
    k = exp(rnorm(n, mean=log(0.512), sd=0.1))
  )
}


#Define the growth function for plotting
growth_function_plot <- function(x, pars=list(canham_max_growth_hat, 
                                              canham_diameter_at_max_growth_hat, 
                                              canham_K_hat)){
  growth <- pars$canham_max_growth_hat * 
    exp(-0.5 * (log(x / pars$canham_diameter_at_max_growth_hat) / pars$canham_K_hat)^2 )
  return(growth)
}

model_median_growth_history_list <- list(
  canham_med_history = list(
    growth_name = "canham_G_hat",
    size_name = "canham_S_hat",
    name = "Canham",
    growth_pars = growth_par_names,
    growth_function = growth_function_plot
  )
)
