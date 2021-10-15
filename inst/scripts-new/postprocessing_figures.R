
cat("\n Begin postprocessing_figures.R \n")

library(rstan)
library(data.table)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(extraDistr)
library(bayesplot)
library(jcolors)

indir ="~/git/covid19Vaccination/inst" # path to the repo
outdir = file.path('~/Downloads/', "results")
states = strsplit('CA,TX',',')[[1]]
stan_model = "211014b"
JOBID = 3541

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-location.index')
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
source(file.path(indir, "functions-new", "summary_functions.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir.fit.post = file.path(outdir, run_tag, "fits")
outdir.data = file.path(outdir, run_tag, "data")
outdir.fig.post = file.path(outdir, run_tag, "figure", run_tag)
outdir.table = file.path(outdir, run_tag, "table", run_tag)

# code
locations = readRDS( file.path(outdir.fit.post, paste0("location_", run_tag,".rds")) )
Code = locations[location.index,]$code

# load image 
load(file.path(outdir.data, paste0("stanin_",run_tag,".RData")))
outdir.fig = outdir.fig.post
outdir.fit = outdir.fit.post

# load fit cumulative deaths
cat("Load fits \n")
file = file.path(outdir.fit, paste0("fit_cumulative_deaths_", run_tag,".rds"))
fit_cum <- readRDS(file=file)

# find date for the first 10th cumulative deaths
date_10thcum = subset(data, !is.na(weekly.deaths))
date_10thcum = date_10thcum[, list(weekly.deaths = sum(weekly.deaths)), by = c('date', 'code')]
date_10thcum[, cum.deaths := cumsum(weekly.deaths), by = 'code']
date_10thcum = date_10thcum[ cum.deaths >=10, list(date_10thcum = min(date)), by = 'code']

# age groups
fouragegroups = c('0-24', '25-54', '55-79', '80+')
fiveagegroups = c('0-24', '25-54', '55-74', '75-84', '85+')

# week when resurgence started
# start_resurgence <- as.Date('2021-07-03')

# age group vaccination
if(!is.null(stan_data$max_age_not_vaccinated )){
  df_age_vac_effects = copy(df_age_continuous)
  df_age_vac_effects[, age_index := c(rep(0, stan_data$max_age_not_vaccinated ), stan_data$map_A_to_C)]
  df_age_vac_effects = df_age_vac_effects[, list(age_from = min(age), age_to = max(age), 
                                                 age = paste0(min(age), '-', max(age))), by = 'age_index']
  df_age_vac_effects[age_to == max(df_age_continuous$age), age := paste0(age_from, '+')]
  
}


# Plot estimate B-splines parameters plane 
plot_posterior_plane(fit_cum, df_week, df_age_continuous, stan_data, outdir = outdir.fig)


# Plots continuous and aggregated age distribution phi
age_contribution_continuous_table = make_var_by_age_table(fit_cum, df_week, df_age_continuous, 'phi', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_continuous_table, 'phi', outdir = outdir.fig)
age_contribution_discrete_table = make_var_by_age_table(fit_cum, df_week, df_age_reporting, 'phi_reduced', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_discrete_table, 'phi_reduced', outdir = outdir.fig, discrete = T)


# baseline contribution adjusted and non-adjusted for population composition
make_contribution_ref(fit_cum, date_10thcum, fiveagegroups, data, df_week, df_age_continuous, outdir.table)
make_contribution_ref_adj(fit_cum, date_10thcum, fiveagegroups, df_week, pop_data, outdir.table)


# contribution over time per age groups
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '0-64', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '0-74', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '0-54', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '20-64', date_10thcum, pop_data, data, outdir.table, with_empirical = F)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '55-74', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '65-79', date_10thcum, pop_data, data, outdir.table)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '65-74', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '75-84', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '65+', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '75+', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '80+', date_10thcum, pop_data, data, outdir.table)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '85+', date_10thcum, pop_data, data, outdir.table, with_empirical = T)


# mortality rate
mortality_rate_table = make_mortality_rate_table(fit_cum, fiveagegroups, date_10thcum, df_week, pop_data,
                                                 JHUData, df_age_continuous, 'cumulative_deaths' , outdir.table)
plot_mortality_rate(mortality_rate_table, outdir.fig)


# predicted weekly deaths by various age groups 
deatht = make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table)
tmp = make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'alpha_reduced', df_age_reporting, outdir.table, withempirical = T)
make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                    age_groups = c('0-54', '55-74', '75+'), lab = '3agegroups', withempirical = T,
                                    reduction = c(vaccine_data[date %in% df_week$date & prop > 0, min(date)], start_resurgence-7))
make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                    age_groups = c('0-74', '75+'), lab = '2agegroups', withempirical = T,
                                    reduction = c(vaccine_data[date %in% df_week$date & prop > 0, min(date)], start_resurgence-7))
make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                    age_groups = unique(df_age_vaccination$age), lab = 'vacagegroups', withempirical = F,
                                    reduction = c(vaccine_data[date %in% df_week$date & prop > 0, min(date)], start_resurgence-7))


# vaccine effect
if(!is.null(stan_data$prop_vac)){
  min_age_index_vac = 3
  df_age_vaccination2 = df_age_vaccination[age_index >= 3]
  df_age_vaccination2[, age_index := age_index - min_age_index_vac + 1]
  
  df_week2 = df_week[week_index >= stan_data$w_start_resurgence & week_index <= stan_data$w_stop_resurgence]
  df_week2[, week_index := 1:nrow(df_week2)]
  
  chi_table = make_var_by_age_age_table(fit_cum, df_age_vaccination2, 'chi', outdir.table)
  psi_table = make_var_by_age_age_table(fit_cum, df_age_vaccination2, 'psi', outdir.table)
  plot_estimate_vaccine(psi_table, 'pre-resurgence proportion', 'pre_resurgence', outdir.fig)
  plot_estimate_vaccine(chi_table, 'trends in the proportion', 'trends', outdir.fig)
  
  E_pdeaths = make_var_by_age_table(fit_cum, df_week, df_age_vaccination2, 'E_pdeaths', outdir.table)
  E_pdeaths_counterfactual = make_var_by_age_table(fit_cum, df_week2, df_age_vaccination2, 'E_pdeaths_counterfactual', outdir.table)
  diff_E_pdeaths_counterfactual = make_var_by_age_table(fit_cum, df_week2, df_age_vaccination2, 'diff_E_pdeaths_counterfactual', outdir.table)
  plot_vaccine_effects_counterfactual(E_pdeaths_counterfactual, E_pdeaths, diff_E_pdeaths_counterfactual, outdir.fig)
  
  prop_vac = prepare_prop_vac_table(stan_data)
  r_pdeaths = make_var_by_age_table(fit_cum, df_week2, df_age_vaccination2, 'r_pdeaths', outdir.table)
  plot_relative_resurgence_vaccine(r_pdeaths, prop_vac, df_age_vaccination2, df_week2, outdir.fig)
  find_stats_vaccine_effects(start_resurgence, pick_resurgence, diff_E_pdeaths_counterfactual, prop_vac, outdir.table)
  
}


# compare to DoH data
if(nrow(subset(scrapedData, code %in% Code)) > 0 ){

  tmp = find_cumulative_deaths_prop_givensum_state_age(fit_cum, date_10thcum, df_week, df_age_continuous, scrapedData, 'cum.deaths', outdir.table)
  
  compare_CDCestimation_DoH_age_plot(CDC_data = copy(tmp), scraped_data = scrapedData,
                                     var.cum.deaths.CDC = c('M_abs_cum', 'CL_abs_cum', 'CU_abs_cum'), outdir = outdir.fig)
  compare_CDCestimation_DoH_age_prop_plot(copy(tmp), outdir.fig)
  compare_CDCestimation_DoH_age_weekly_plot(copy(tmp), outdir.fig)

}

scrapedData = subset(scrapedData, code == Code)

if(Code == 'GA')
  scrapedData = reduce_agebands_scrapedData_GA(scrapedData)



# make panel figure
age_contribution_continuous_table$method = 'GP-BS-SE'
deatht$method = 'GP-BS-SE'
p2=plot_contribution_continuous_comparison_method(copy(age_contribution_continuous_table), copy(deatht), copy(data),
                                                  'GP-BS-SE', 'GP-BS-SE',
                                                  show.method = F,
                                                  heights = c(1,1), outdir.fig)



cat("\n End postprocessing_figures.R \n")



