#Control file for constant growth model
#Build list for diagnostics
diagnostic_list_single_species <- list(
  plot1=list(
    pars=c(
      "species_beta_mu",
      "species_beta_sigma"
    ),
    name="const_species_hyperparams_single_species"
  ),
  plot2=list(
    pars=c("global_error_sigma"), 
    name="const_global_error_sigma_single_species"
  )
)

diagnostic_list_single_individual <- list(
  plot1=list(
    pars=c("ind_beta"),
    name="const_species_hyperparams_single_individual"
  ),
  plot2=list(
    pars=c("global_error_sigma"), 
    name="const_global_error_sigma_single_individual"
  )
)

#List of controls for single species model
single_species_controls <- list(
  model_name = "Single species constant model",
  stan_file_path = "stan/ConstGrowth_SingleSpecies.stan",
  fit_save_path = "output/data/Single_Species_Fit_Const",
  stan_controls = list(
    iter = 2000,
    chains = 4,
    cores = 4
  ),
  diagnostic_list = diagnostic_list_single_species
)

#List of controls for single individual model
single_individual_controls <- list(
  model_name = "Single individual constant model",
  stan_file_path = "stan/ConstGrowth_SingleIndividual.stan",
  fit_save_path = "output/data/Single_Ind_Fit_Const",
  stan_controls = list(
    iter = 2000,
    chains = 4,
    cores = 1
  ),
  diagnostic_list = diagnostic_list_single_individual
)

#Individual-level parameter controls
ind_pars <- c("ind_S_0", "ind_beta") #Internal model parameter names
ind_pars_names <- c("const_S_0_hat", "const_beta_hat") #Param names with model
mean_size_name <- "const_S_hat_mean"

#Species-level parameter controls
sp_pars <- c("species_beta_mu", "species_beta_sigma")
sp_pars_names <- c("const_sp_beta_mean", "const_sp_beta_sigma")

#Plot parameter controls
growth_pars <- c("ind_beta")
growth_par_names <- c("const_beta_hat")
size_name <- "const_S_hat"
growth_name <- "const_G_hat"
cond <- "Const Pred"

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
    log_log = FALSE
  ),
  
  const_plot1 = list(
    par = "const_beta_hat",
    name = "const_beta",
    par_name = "Beta",
    log_log = FALSE
  )
)

#List of plots for beeswarm/box
parplot_list <- list(
  const_plot1 = list(
    plot_name = "_beta_ratio_",
    pars = list(
      par_1 = list(
        par = "ratio_const_beta_hat",
        par_name = "Beta"
      )
    ),
    xlab="Parameter",
    ylab="Ratio true/est",
    hline_height=1
  ),
  
  const_plot2 = list(
    plot_name="_beta_difference_",
    pars = list(
      par_1 = list(
        par = "difference_const_beta_hat",
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
  return(pars[[1]])
}

growth_function_pars <- c("beta")

sim_spec <- function(n) {
  tibble(
    beta = exp(rnorm(n, mean=-2, sd=1))
  )
}


#Define the growth function for plotting
growth_function_plot <- function(x, pars=list(const_beta)){
  growth <- pars$const_beta
  return(growth)
}

model_median_growth_history_list <- list(
  const_med_history = list(
    growth_name = "const_G_hat",
    size_name = "const_S_hat",
    name = "Constant",
    growth_pars = growth_par_names,
    growth_function = growth_function_plot
  )
)