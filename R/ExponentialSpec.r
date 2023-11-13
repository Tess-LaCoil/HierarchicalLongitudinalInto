#Control file for exponential growth model
#Build list for diagnostics
diagnostic_list_single_species <- list(
  plot1=list(
    pars=c(
      "species_beta_mu",
      "species_beta_sigma"
    ),
    name="exp_species_hyperparams_single_species"
  ),
  plot2=list(
    pars=c("global_error_sigma"), 
    name="exp_global_error_sigma_single_species"
  )
)

diagnostic_list_single_individual <- list(
  plot1=list(
    pars=c("ind_beta"),
    name="exp_species_hyperparams_single_individual"
  ),
  plot2=list(
    pars=c("global_error_sigma"), 
    name="exp_global_error_sigma_single_individual"
  )
)

#List of controls for single-species model
single_species_controls <- list(
  model_name = "Single species exponential model",
  stan_file_path = "stan/Exp_SingleSpecies.stan",
  fit_save_path = "output/data/Single_Species_Fit_Exp",
  stan_controls = list(
    iter = 2000,
    chains = 4,
    cores = 4,
    control = list(max_treedepth = 15)
  ),
  diagnostic_list = diagnostic_list_single_species
)

single_individual_controls <- list(
  model_name = "Single individual exponential model",
  stan_file_path = "stan/Exp_SingleIndividual.stan",
  fit_save_path = "output/data/Single_Ind_Fit_Exp",
  stan_controls = list(
    iter = 2000,
    chains = 4,
    cores = 4,
    control = list(max_treedepth = 15)
  ),
  diagnostic_list = diagnostic_list_single_individual
)

#Individual-level parameter controls
ind_pars <- c("ind_S_0", "ind_beta") #Internal model parameter names
ind_pars_names <- c("exp_S_0_hat", "exp_beta_hat") #Param names with model
mean_size_name <- "exp_S_hat_mean"

#Species-level parameter controls
sp_pars <- c("species_beta_mu", "species_beta_sigma")
sp_pars_names <- c("exp_sp_beta_mean", "exp_sp_beta_sigma")

#Plot parameter controls
growth_pars <- c("ind_beta")
growth_par_names <- c("exp_beta_hat")
size_name <- "exp_S_hat"
growth_name <- "exp_G_hat"
cond <- "Exp Pred"

#Lists for plot functions
scatterplot_list <- list(
  size_plot = list(
    par = "S_hat",
    name = "S_hat",
    par_name = "Estimated Size",
    log_log = TRUE
  ),
  
  growth_plot = list(
    par = "G_hat",
    name = "G_hat",
    par_name = "Estimated Growth",
    log_log = TRUE
  ),
  
  exp_plot1 = list(
    par = "exp_beta_hat",
    name = "exp_beta",
    par_name = "Beta",
    log_log = FALSE
  )
)

#List of plots for beeswarm/box
parplot_list <- list(
  exp_plot1 = list(
    plot_name = "_beta_ratio_",
    pars = list(
      par_1 = list(
        par = "ratio_exp_beta_hat",
        par_name = "Beta"
      )
    ),
    xlab="Parameter",
    ylab="Ratio true/est",
    hline_height=1
  ),
  
  exp_plot2 = list(
    plot_name="_beta_difference_",
    pars = list(
      par_1 = list(
        par = "difference_exp_beta_hat",
        par_name = "Beta"
      )
    ),
    xlab="Parameter",
    ylab="Difference true - est",
    hline_height=0
  )
)

#Define growth function
growth_function <- function(x, pars){
  return(pars[[1]]*x)
}

growth_function_pars <- c("beta")

sim_spec <- function(n) {
  tibble(
    beta = exp(rnorm(n, mean=-4, sd=0.3))
  )
}

#Define the growth function to be used to plot the median behaviour
growth_function_plot <- function(x, pars=list(exp_beta)){
  growth <- pars$exp_beta * x
  return(growth)
}

model_median_growth_history_list <- list(
  exp_med_history = list(
    growth_name = "exp_G_hat",
    size_name = "exp_S_hat",
    name = "exponential",
    growth_pars = growth_par_names,
    growth_function = growth_function_plot
  )
)