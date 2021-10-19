
cat("\n Begin postprocessing_assess_mising.R \n")

library(rstan)
library(data.table)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(extraDistr)
library(bayesplot)

indir = "/rds/general/user/mm3218/home/git/covid19Vaccination/inst/" # path to the repo
outdir = '/rds/general/user/mm3218/home/git/covid19Vaccination/inst/results/'
states = strsplit('CA,TX',',')[[1]]
stan_model = "211014b"
JOBID = 7259


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
source(file.path(indir, "functions-new", "postprocessing-plotting_functions.R"))
source(file.path(indir, "functions-new", "postprocessing-summary_functions.R"))
source(file.path(indir, "functions-new", "postprocessing-utils.R"))

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


# Convergence diagnostics
cat("\nMake convergence diagnostics \n")
make_convergence_diagnostics_stats(fit_cum, outdir.table)

# Make predictive checks table
cat("\nMake posterior predictive checks table \n")
predictive_checks_table = make_predictive_checks_table(fit_cum, df_week, df_age_reporting, 
                                                       data, 'deaths_predict_state_age_strata', outdir.table)

# plot predictive checks table
cat("\nMake posterior predictive checks plots \n")
plot_posterior_predictive_checks(predictive_checks_table, 
                                 variable = "weekly.deaths", 
                                 outdir = outdir.fig)

# plot predictive check cumulative
cat("\nMake posterior predictive checks cumulative plots \n")
tmp1 = find_cumsum_nonr_deaths_state_age(fit_cum, df_week, df_age_continuous, unique(df_age_reporting$age), stan_data, 'deaths_predict')
plot_sum_missing_deaths(tmp1, outdir.fig)

tmp1 = find_sum_nonr_deaths_state_age(fit_cum, df_age_continuous, unique(df_age_reporting$age), stan_data, 'deaths_predict', outdir.table)
plot_sum_bounded_missing_deaths(tmp1, outdir.fig)

# trace and paris plots
p <- bayesplot::mcmc_trace(fit_cum, regex_pars = c('nu', 'alpha_gp', 'rho_gp'))
ggsave(p, file = paste0(outdir.fig, '-mcmc_trace_parameters.png'), h = 20, w = 20, limitsize = F)

p <- bayesplot::mcmc_pairs(fit_cum, regex_pars = c('nu', 'alpha_gp', 'rho_gp'))
ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_parameters.png'), h = 20, w = 20, limitsize = F)

samples <- rstan::extract(fit_cum)
names_samples <- names(samples)
names_fit <- names(fit_cum)

if(any(c('varphi', 'psi', 'chi', 'kappa') %in% names_samples)){

  p <- bayesplot::mcmc_trace(fit_cum, regex_pars = c('varphi', 'psi', 'chi', 'kappa'))
  ggsave(p, file = paste0(outdir.fig, '-mcmc_trace_vaccine_parameters.png'), h = 10, w = 10, limitsize = F)
  
  p <- bayesplot::mcmc_pairs(fit_cum, regex_pars = c('varphi', 'psi', 'chi', 'kappa'))
  ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_vaccine_parameters.png'), h = 20, w = 20, limitsize = F)
  
  
}

if(any(c('intercept_resurgence0', 'slope_resurgence0', 'vaccine_effect_intercept') %in% names_samples)){
  
  p <- bayesplot::mcmc_trace(fit_cum, regex_pars = c('intercept_resurgence', 'slope_resurgence', 'vaccine_effect'))
  ggsave(p, file = paste0(outdir.fig, '-mcmc_trace_vaccine_parameters.png'), h = 20, w = 20, limitsize = F)
  
  # p <- bayesplot::mcmc_pairs(fit_cum, regex_pars = c('intercept_resurgence', 'slope_resurgence', 'vaccine_effect'))
  # ggsave(p, file = paste0(outdir.fig, '-mcmc_pair_vaccine_parameters.png'), h = 20, w = 20, limitsize = F)
  
  names_var = c('intercept_resurgence_re', 'slope_resurgence_re', 'vaccine_effect_intercept', 'vaccine_effect_slope')
  tmp <- data.table(name= names_fit[ grepl(paste(paste0('^',names_var),collapse = '|'),names_fit) ])
  p <- bayesplot::mcmc_intervals(fit_cum, pars=tmp$name,prob = .95, prob_outer = 0.95)+
    theme_bw()
  ggsave(p, file = paste0(outdir.fig, '-mcmc_interval_vaccine_parameters.png'), h = 10, w = 7, limitsize = F)
  
}


cat("\n End postprocessing_assess_mixing.R \n")



