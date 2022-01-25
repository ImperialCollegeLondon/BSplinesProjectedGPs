
cat("\n Begin postprocessing_assess_mising.R \n")

library(rstan)
library(data.table)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(extraDistr)
library(bayesplot)
library(scales)

indir = "/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/" # path to the repo
outdir = '/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/results/'
# states = strsplit('CA,FL,NY,TX,PA,IL,OH,GA,NC,MI',',')[[1]]
states = strsplit('CA,FL,NY,TX',',')[[1]]
stan_model = "220120a"
JOBID = 2455

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-states')
  stopifnot(args_line[[7]]=='-stan_model')
  stopifnot(args_line[[9]]=='-JOBID')
  indir <- args_line[[2]]
  outdir <- args_line[[4]]
  states <-  strsplit((args_line[[6]]),',')[[1]]
  stan_model <- args_line[[8]]
  JOBID <- as.numeric(args_line[[10]])
}


# load functions
source(file.path(indir, "functions", "postprocessing-plotting_functions.R"))
source(file.path(indir, "functions", "postprocessing-summary_functions.R"))
source(file.path(indir, "functions", "postprocessing-utils.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir.fit.post = file.path(outdir, run_tag, "fits")
outdir.data = file.path(outdir, run_tag, "data")
outdir.fig.post = file.path(outdir, run_tag, "figure", run_tag)
outdir.table = file.path(outdir, run_tag, "table", run_tag)

if( !dir.exists( dirname(outdir.fig.post) ) ) dir.create( dirname(outdir.fig.post), recursive = T)
if( !dir.exists( dirname(outdir.table) ) ) dir.create( dirname(outdir.table), recursive = T)

# code
locations = readRDS( file.path(outdir.fit.post, paste0("location_", run_tag,".rds")) )
loc_name = locations[code %in% states,]$loc_label
Code = locations[code %in% states, ]$code
df_state = data.table(loc_label = loc_name, code = Code, state_index = 1:length(Code))

# load image 
load(file.path(outdir.data, paste0("stanin_",run_tag,".RData")))
outdir.fig = outdir.fig.post
outdir.fit = outdir.fit.post

# load fit cumulative deaths
cat("Load fits \n")
file = file.path(outdir.fit.post, paste0("fit_cumulative_deaths_", run_tag,".rds"))
fit_cum <- readRDS(file=file)
fit_samples <- rstan::extract(fit_cum)

# Convergence diagnostics
cat("\nMake convergence diagnostics \n")
summary  <- make_convergence_diagnostics_stats(fit_cum, fit_samples, outdir.table)

# Make predictive checks table
cat("\nMake posterior predictive checks table \n")
predictive_checks_table = make_predictive_checks_table(fit_cum, fit_samples, df_week, df_age_reporting, 
                                                       data, 'deaths_predict_state_age_strata', outdir.table)

# plot predictive checks table
cat("\nMake posterior predictive checks plots \n")
plot_posterior_predictive_checks(predictive_checks_table, 
                                 variable = "weekly.deaths", 
                                 outdir = outdir.fig)

# plot predictive check cumulative
cat("\nMake posterior predictive checks cumulative plots \n")
tmp1 = find_cumsum_nonr_deaths_state_age(fit_samples, df_week, df_age_continuous, unique(df_age_reporting$age), stan_data, 'deaths_predict')
plot_sum_missing_deaths(tmp1, outdir.fig)

tmp1 = find_sum_nonr_deaths_state_age(fit_samples, df_age_continuous, unique(df_age_reporting$age), stan_data, 'deaths_predict', outdir.table)
plot_sum_bounded_missing_deaths(tmp1, outdir.fig)

# trace and paris plots
cat("\n Trace plot weekly deaths params \n")
p <- bayesplot::mcmc_trace(fit_cum, regex_pars = c('nu', 'alpha_gp', 'rho_gp'))
ggsave(p, file = paste0(outdir.fig, '-mcmc_trace_parameters.png'), h = 20, w = 20, limitsize = F)

cat("\n Pairs plot weekly deaths params \n")
p <- bayesplot::mcmc_pairs(fit_cum, regex_pars = c('nu', 'alpha_gp', 'rho_gp'))
ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_parameters.png'), h = 20, w = 20, limitsize = F)

# forest plots
samples <- rstan::extract(fit_cum)
names_samples <- names(samples)
names_fit <- names(fit_cum)

## base model
lambda_table <- make_lambda_table(fit_samples, stan_data, df_week, df_state, outdir.table)
plot_lambda_table(lambda_table, outdir.fig)
var_base_model_table <- make_var_base_model_table(fit_samples, stan_data, df_state, outdir.table)
plot_var_base_model_table(var_base_model_table, outdir.fig)

## vaccination model parameters
if(any(c('varphi', 'psi', 'chi', 'kappa') %in% names_samples)){

  cat("\n Trace plot vaccination params \n")
  p <- bayesplot::mcmc_trace(fit_cum, regex_pars = c('varphi', 'psi', 'chi', 'kappa'))
  ggsave(p, file = paste0(outdir.fig, '-mcmc_trace_vaccine_parameters.png'), h = 10, w = 10, limitsize = F)
  
  cat("\n Pairs plot vaccination params \n")
  p <- bayesplot::mcmc_pairs(fit_cum, regex_pars = c('varphi', 'psi', 'chi', 'kappa'))
  ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_vaccine_parameters.png'), h = 20, w = 20, limitsize = F)

}

if(any(c('slope_resurgence0', 'vaccine_effect_slope') %in% names_samples)){
  
  ## trace plots
  cat("\n Trace plot vaccination slope params \n")
  p <- bayesplot::mcmc_trace(fit_cum, regex_pars = c('slope_resurgence', 'vaccine_effect_slope'))
  ggsave(p, file = paste0(outdir.fig, '-mcmc_trace_vaccine_slope_parameters.png'), h = 20, w = 20, limitsize = F)
  
  ## pairs plot 
  cat("\n Pairs plot vaccination slope params \n")
  names_var = c('slope_resurgence0', 'vaccine_effect_slope')
  tmp <- data.table(name= names_fit[ grepl(paste(paste0('^',names_var),collapse = '|'),names_fit) ])
  p <- bayesplot::mcmc_pairs(fit_cum, pars = tmp[grepl('\\[1', name)]$name)
  ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_vaccine_slope_1864_parameters.png'), h = 20, w = 20, limitsize = F)
  p <- bayesplot::mcmc_pairs(fit_cum, pars = tmp[grepl('\\[2', name)]$name)
  ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_vaccine_slope_65p_parameters.png'), h = 20, w = 20, limitsize = F)
  
  ## interval plots
  cat("\n Interval plot vaccination slope params \n")
  names_var = c('slope_resurgence0', 'slope_resurgence_re', 'slope_resurgence', 'vaccine_effect_slope')
  tmp <- data.table(name= names_fit[ grepl(paste(paste0('^',names_var),collapse = '|'),names_fit) ])
  p <- bayesplot::mcmc_intervals(fit_cum, pars=tmp[grepl('\\[1', name)]$name,prob = .95, prob_outer = 0.95)+
    theme_bw() + 
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'grey50')
  ggsave(p, file = paste0(outdir.fig, '-mcmc_interval_vaccine_slope_1864_parameters.png'), h = 6, w = 4, limitsize = F)
  p <- bayesplot::mcmc_intervals(fit_cum, pars=tmp[grepl('\\[2', name)]$name,prob = .95, prob_outer = 0.95)+
    theme_bw() + 
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'grey50')
  ggsave(p, file = paste0(outdir.fig, '-mcmc_interval_vaccine_slope_65p_parameters.png'), h = 6, w = 4, limitsize = F)
  
  names_var = c('sigma_slope_resurgence')
  tmp <- data.table(name= names_fit[ grepl(paste(paste0('^',names_var),collapse = '|'),names_fit) ])
  if(nrow(tmp) > 0){
    p <- bayesplot::mcmc_intervals(fit_cum, pars=tmp$name,prob = .95, prob_outer = 0.95)+
      theme_bw() + 
      geom_vline(xintercept = 0, linetype = 'dashed', col = 'grey50')
    ggsave(p, file = paste0(outdir.fig, '-mcmc_interval_vaccine_parameters_var_slope.png'), h = 6, w = 4, limitsize = F)
    
  }

}


if(any(c('intercept_resurgence0', 'vaccine_effect_intercept') %in% names_samples)){
  
  ## trace plots
  cat("\n Trace plot vaccination intercept params \n")
  p <- bayesplot::mcmc_trace(fit_cum, regex_pars = c('intercept_resurgence', 'vaccine_effect_intercept'))
  ggsave(p, file = paste0(outdir.fig, '-mcmc_trace_vaccine_intercept_parameters.png'), h = 20, w = 20, limitsize = F)
  
  ## pairs plot 
  cat("\n Pairs plot vaccination intercept params \n")
  names_var = c('intercept_resurgence0', 'vaccine_effect_intercept')
  tmp <- data.table(name= names_fit[ grepl(paste(paste0('^',names_var),collapse = '|'),names_fit) ])
  p <- bayesplot::mcmc_pairs(fit_cum, pars = tmp[grepl('\\[1', name)]$name)
  ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_vaccine_intercept_1864_parameters.png'), h = 20, w = 20, limitsize = F)
  p <- bayesplot::mcmc_pairs(fit_cum, pars = tmp[grepl('\\[2', name)]$name)
  ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_vaccine_intercept_65p_parameters.png'), h = 20, w = 20, limitsize = F)

  ## interval plots
  cat("\n Interval plot vaccination intercept params \n")
  names_var = c('intercept_resurgence0', 'intercept_resurgence_re', 'intercept_resurgence', 'vaccine_effect_intercept')
  tmp <- data.table(name= names_fit[ grepl(paste(paste0('^',names_var),collapse = '|'),names_fit) ])
  p <- bayesplot::mcmc_intervals(fit_cum, pars=tmp[grepl('\\[1', name)]$name,prob = .95, prob_outer = 0.95)+
    theme_bw() + 
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'grey50')
  ggsave(p, file = paste0(outdir.fig, '-mcmc_interval_vaccine_intercept_1864_parameters.png'), h = 6, w = 4, limitsize = F)
  p <- bayesplot::mcmc_intervals(fit_cum, pars=tmp[grepl('\\[2', name)]$name,prob = .95, prob_outer = 0.95)+
    theme_bw() + 
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'grey50')
  ggsave(p, file = paste0(outdir.fig, '-mcmc_interval_vaccine_intercept_65p_parameters.png'), h = 6, w = 4, limitsize = F)
  
  names_var = c('sigma_intercept_resurgence')
  tmp <- data.table(name= names_fit[ grepl(paste(paste0('^',names_var),collapse = '|'),names_fit) ])
  if(nrow(tmp) >0){
    p <- bayesplot::mcmc_intervals(fit_cum, pars=tmp$name,prob = .95, prob_outer = 0.95)+
      theme_bw() + 
      geom_vline(xintercept = 0, linetype = 'dashed', col = 'grey50')
    ggsave(p, file = paste0(outdir.fig, '-mcmc_interval_vaccine_parameters_var_intercept.png'), h = 6, w = 4, limitsize = F)
    
  }

}


names = c('slope_resurgence0', 'slope_resurgence_re', 'vaccine_effect_slope', 'intercept_resurgence0', 'intercept_resurgence_re', 'vaccine_effect_intercept')
if(any(names %in% names_samples)){
  
  min_age_index_vac = 3
  df_age_vaccination2 = df_age_vaccination[age_index >= 3]
  df_age_vaccination2[, age_index := age_index - min_age_index_vac + 1]
  
  math_name = c('psi^"base"*""', 'psi^"state"*""', 'psi^"vacc"*""', 'chi^"base"*""', 'chi^"state"*""', 'chi^"vacc"*""')
  groups = c('slope', 'slope', 'indirect\nvaccine\neffects', 'baseline', 'baseline', 'indirect\nvaccine\neffects')
  groups_levels = c('baseline', 'slope', 'indirect\nvaccine\neffects')
  
  if(!('vaccine_effect_intercept_cross' %in% names_samples)){

    tmp <- make_forest_plot_table(summary, df_age_vaccination2, df_state, names, math_name, groups, groups_levels)

    
  } else{
    names = c('slope_resurgence0', 'slope_resurgence_re', 'vaccine_effect_slope_cross', 'intercept_resurgence0', 'intercept_resurgence_re', 'vaccine_effect_intercept_cross')
    math_name = c('psi^"base"*""', 'psi^"state"*""', 'psi^"vacc-cross"*""', 'chi^"base"*""', 'chi^"state"*""', 'chi^"vacc-cross"*""')
    
    tmp <- make_forest_plot_table2(summary, df_age_vaccination2, df_state, names, math_name, groups, groups_levels)

  }
  
  forest_plot <- plot_forest_plot(tmp, outdir.fig)
  
}

names = c('vaccine_effect_intercept_diagonal', 'vaccine_effect_slope_diagonal')
if(any(names %in% names_samples)){

  math_name = c('chi^"vacc"*""', 'psi^"vacc"*""')

  tmp <- as.data.table(summary[ grepl(paste(paste0('^',names),collapse = '|'),rownames(summary)) ,])
  setnames(tmp, c('50%', '2.5%', '97.5%'), c('M', 'CL', "CU"))
  tmp[, variable := math_name]
  
  plot_forest_plot_with_common_effect(forest_plot, tmp, outdir.fig)

}



  
  
cat("\n End postprocessing_assess_mixing.R \n")



