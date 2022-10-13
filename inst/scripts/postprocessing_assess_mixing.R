
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
library(facetscales)
library(truncnorm)
library(invgamma)

if(0){
  indir = "~/git/BSplinesProjectedGPs/inst/" # path to the repo
  outdir = '~/Downloads/results/'
}
indir = "/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/" # path to the repo
outdir = '/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/results/'
states <- 'CA'
stan_model = "220209a"
JOBID = 3541

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
file = file.path(outdir.fit.post, paste0("posterior_samples_", run_tag,".rds"))
if(file.exists(file)){
  fit_samples <- readRDS(file=file)
}else{
  fit_samples <- rstan::extract(fit_cum)
}

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
tryCatch(
  # This is what I want to do...
  {
    p <- bayesplot::mcmc_trace(fit_cum, regex_pars = c('nu', 'zeta_gp', 'gamma_gp'))
    ggsave(p, file = paste0(outdir.fig, '-mcmc_trace_parameters.png'), h = 50, w = 50, limitsize = F)
    
    cat("\n Pairs plot weekly deaths params \n")
    p <- bayesplot::mcmc_pairs(fit_cum, regex_pars = c('nu', 'zeta_gp', 'gamma_gp'))
    ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_parameters.png'), h = 50, w = 50, limitsize = F)
    
  },
  # ... but if an error occurs, tell me what happened: 
  error=function(error_message) {
    message("This is my custom message.")
    message("And below is the error message from R:")
    message(error_message)
    return(NA)
  }
)

# forest plots
names_samples <- names(fit_samples)
names_fit <- names(fit_cum)

## base model
if(!with_cmdstan){
  model = rstan::stan_model(path.to.stan.model)
}
lambda_table <- make_lambda_table(fit_samples, stan_data, df_week, df_state, outdir.table)
plot_lambda_table(lambda_table, outdir.fig)
var_base_model_table <- make_var_base_model_table(fit_samples, stan_data, df_state, outdir.table)
plot_var_base_model_table(var_base_model_table, outdir.fig)

## vaccination model parameters
names = c('slope_resurgence0', 'slope_resurgence_re', 'intercept_resurgence0', 'intercept_resurgence_re', 
          'vaccine_effect_intercept_cross', 'vaccine_effect_intercept_diagonal', 'vaccine_effect_intercept0',
          'vaccine_effect_slope_cross', 'vaccine_effect_slope_diagonal', 'vaccine_effect_slope0'
          )
if(any(names %in% names_samples)){
  
  # summary <- summary(fit_cum)$summary
  vars <- c('vaccine_effect_intercept_cross', 'vaccine_effect_intercept_diagonal', 'vaccine_effect_intercept0',
            'vaccine_effect_slope_cross', 'vaccine_effect_slope_diagonal', 'vaccine_effect_slope0')
  vars <- vars[vars %in% names_samples]
  p <- bayesplot::mcmc_trace(fit_cum, regex_pars = vars)
  ggsave(p, file = paste0(outdir.fig, '-mcmc_trace_parameters_vaccination.png'), h = 20, w = 20, limitsize = F)
  
  p <- bayesplot::mcmc_pairs(fit_cum, regex_pars = vars)
  ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_parameters_vaccination.png'), h = 50, w = 50, limitsize = F)
  
  
  min_age_index_vac = 3
  df_age_vaccination2 = df_age_vaccination[age_index >= 3]
  df_age_vaccination2[, age_index := age_index - min_age_index_vac + 1]
  
  math_name = c('psi^"base"*""', 'psi^"state"*""', 'chi^"base"*""', 'chi^"state"*""', 
                'chi^"vacc-cross"*""', 'chi^"vacc"*""',  'chi^"vacc0"*""',
                'psi^"vacc-cross"*""', 'psi^"vacc"*""', 'psi^"vacc0"*""')
  groups = c('slope', 'slope', 'baseline', 'baseline', 'baseline', 'baseline', 'baseline','slope', 'slope', 'slope')
  groups_levels = c('baseline', 'slope')
  
  tmp <- make_forest_plot_table(summary, df_age_vaccination2, df_state, names, math_name, groups, groups_levels)
  
  forest_plot <- plot_forest_plot(tmp, outdir.fig)
  
}


cat("\n End postprocessing_assess_mixing.R \n")



