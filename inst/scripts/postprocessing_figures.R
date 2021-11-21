
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

indir ="~/git/BSplinesProjectedGPs/inst" # path to the repo
outdir = file.path('/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst', "results")
# states = strsplit('CA,FL,NY,TX,PA,IL,OH,GA,NC,MI',',')[[1]]
states = strsplit('CA,FL,NY,TX',',')[[1]]
stan_model = "211031d2"
JOBID = 20145

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
source(file.path(indir, "functions", "postprocessing-statistics_functions.R"))
source(file.path(indir, "functions", "postprocessing-utils.R"))
source(file.path(indir, "functions", "summary_functions.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir.fit.post = file.path(outdir, run_tag, "fits")
outdir.data = file.path(outdir, run_tag, "data")
outdir.fig.post = file.path(outdir, run_tag, "figure", run_tag)
outdir.table = file.path(outdir, run_tag, "table", run_tag)

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

# selected states
selected_code = c('CA', 'FL', 'NY', 'TX')

####


# Plots continuous and aggregated age distribution phi
age_contribution_continuous_table = make_var_by_age_by_state_table(fit_cum, df_week, df_age_continuous, df_state, 'phi', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_continuous_table, 'phi', outdir = outdir.fig)
age_contribution_discrete_table = make_var_by_age_by_state_table(fit_cum, df_week, df_age_reporting, df_state, 'phi_reduced', outdir.table)
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
                                    reduction = c(vaccine_data[date %in% df_week$date & prop > 0, min(date)], min(resurgence_dates$start_resurgence)-7))
make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                    age_groups = c('0-74', '75+'), lab = '2agegroups', withempirical = T,
                                    reduction = c(vaccine_data[date %in% df_week$date & prop > 0, min(date)], min(resurgence_dates$start_resurgence)-7))
make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                    age_groups = unique(df_age_vaccination$age), lab = 'vacagegroups', withempirical = F,
                                    reduction = c(vaccine_data[date %in% df_week$date & prop > 0, min(date)], min(resurgence_dates$start_resurgence)-7))

# vaccine effect
if(!is.null(stan_data$prop_vac)){
  
  min_age_index_vac = 3
  df_age_vaccination2 = df_age_vaccination[age_index >= 3]
  df_age_vaccination2[, age_index := age_index - min_age_index_vac + 1]
  
  df_week[, dummy := 1]; resurgence_dates[, dummy := 1]
  df_week2 = merge(df_week, resurgence_dates, by = 'dummy', allow.cartesian=TRUE)
  df_week2 = df_week2[date >= start_resurgence & date <= stop_resurgence]
  df_week2[, week_index := 1:length(date), by = 'code']
  df_week2 = select(df_week2, - start_resurgence, - stop_resurgence, -dummy)
  df_week2 = df_week2[order(code)]
  df_week = select(df_week, -dummy)
  # df_state =  select(df_state, -dummy)
  
  stan_data1 = add_vaccine_prop(stan_data, df_week, Code, vaccine_data, resurgence_dates)
  prop_vac = prepare_prop_vac_table(stan_data1, df_week2, df_age_vaccination2)
  
  make_var_by_age_by_state_table(fit_cum, df_week, df_age_vaccination2, df_state, 'phi_reduced_vac', outdir.table)
  
  E_pdeaths_cum_predict = make_var_cum_by_age_table(fit_cum, df_week, df_age_vaccination2, 'E_pdeaths_predict', outdir.table)
  E_pdeaths_cum_counterfactual = make_var_cum_by_age_table_counterfactual(fit_cum, df_week, df_week2, resurgence_dates, df_age_vaccination2, c('E_pdeaths_predict', 'E_pdeaths_counterfactual'), outdir.table)
  plot_vaccine_effects_counterfactual(E_pdeaths_cum_counterfactual, E_pdeaths_cum_predict, resurgence_dates, 'cumulative', 'cumulative', outdir.fig)
  
  E_pdeaths_predict = make_var_by_age_by_state_table(fit_cum, df_week, df_age_vaccination2, df_state, 'E_pdeaths_predict', outdir.table)
  E_pdeaths_counterfactual = make_var_by_age_by_state_table(fit_cum, df_week2, df_age_vaccination2, df_state, 'E_pdeaths_counterfactual', outdir.table)
  plot_vaccine_effects_counterfactual(E_pdeaths_counterfactual, E_pdeaths_predict, resurgence_dates, 'weekly', 'weekly', outdir.fig)
  
  E_pdeaths_predict_resurgence_cumulative = make_var_by_age_by_state_table(fit_cum, df_week2, df_age_vaccination2, df_state, 'E_pdeaths_predict_resurgence_cumulative', outdir.table)
  E_pdeaths_counterfactual_resurgence_cumulative = make_var_by_age_by_state_table(fit_cum, df_week2, df_age_vaccination2, df_state, 'E_pdeaths_counterfactual_resurgence_cumulative', outdir.table)
  plot_vaccine_effects_counterfactual(subset(E_pdeaths_counterfactual_resurgence_cumulative, code %in% selected_code), 
                                      subset(E_pdeaths_predict_resurgence_cumulative, code %in% selected_code), subset(resurgence_dates, code %in% selected_code), 'cumulative', 'cumulative_rperiod_selected_states', outdir.fig)
  if(any(!Code %in% selected_code)){
    plot_vaccine_effects_counterfactual(subset(E_pdeaths_counterfactual_resurgence_cumulative, !code %in% selected_code), 
                                        subset(E_pdeaths_predict_resurgence_cumulative, !code %in% selected_code), subset(resurgence_dates, !code %in% selected_code), 'cumulative', 'cumulative_rperiod', outdir.fig)
    
  }
  
  diff_E_pdeaths_counterfactual_all <- make_var_by_age_table(fit_cum, df_week2, df_age_vaccination2, 'diff_E_pdeaths_counterfactual_all', outdir.table)
  perc_E_pdeaths_counterfactual_all <- make_var_by_age_table(fit_cum, df_week2, df_age_vaccination2, 'perc_E_pdeaths_counterfactual_all', outdir.table)
  
  diff_E_pdeaths_counterfactual = make_var_by_age_by_state_table(fit_cum, df_week2, df_age_vaccination2, df_state, 'diff_E_pdeaths_counterfactual', outdir.table)
  perc_E_pdeaths_counterfactual = make_var_by_age_by_state_table(fit_cum, df_week2, df_age_vaccination2, df_state, 'perc_E_pdeaths_counterfactual', outdir.table)
  plot_vaccine_effects_counterfactual_stat(diff_E_pdeaths_counterfactual, resurgence_dates, 'difference', outdir.fig)
  plot_vaccine_effects_counterfactual_stat(perc_E_pdeaths_counterfactual, resurgence_dates, 'percentage', outdir.fig)
  find_stats_vaccine_effects(diff_E_pdeaths_counterfactual, perc_E_pdeaths_counterfactual, diff_E_pdeaths_counterfactual_all, perc_E_pdeaths_counterfactual_all, prop_vac, resurgence_dates, outdir.table)
  
  r_pdeaths = make_var_by_age_by_state_table(fit_cum, df_week2, df_age_vaccination2, df_state, 'r_pdeaths', outdir.table)
  plot_relative_resurgence_vaccine(r_pdeaths, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, outdir.fig)
  plot_relative_resurgence_vaccine2(r_pdeaths, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, T, outdir.fig)
  plot_relative_resurgence_vaccine2(r_pdeaths, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, F, outdir.fig)
  plot_relative_resurgence_vaccine_no_time(r_pdeaths, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, T, outdir.fig)
  plot_relative_resurgence_vaccine_no_time(r_pdeaths, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, F, outdir.fig)
  plot_relative_resurgence_vaccine2(subset(r_pdeaths, code %in% selected_code), prop_vac, df_age_vaccination2, df_week2, resurgence_dates, T, outdir.fig, '_selected_states')
  plot_relative_resurgence_vaccine2(subset(r_pdeaths, code %in% selected_code), prop_vac, df_age_vaccination2, df_week2, resurgence_dates, F, outdir.fig, '_selected_states')
  
  if(length(Code) > 6){
    plot_relative_resurgence_vaccine2_long(r_pdeaths, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, outdir.fig)
  }
  
  log_r_pdeaths = make_var_by_age_by_state_table(fit_cum, df_week2, df_age_vaccination2, df_state, 'log_r_pdeaths', outdir.table)
  log_r_pdeaths_predict = make_var_by_age_by_state_table(fit_cum, df_week2, df_age_vaccination2, df_state, 'log_r_pdeaths_predict', outdir.table)
  
  r_pdeaths_predict = make_var_by_age_by_state_table(fit_cum, df_week2, df_age_vaccination2, df_state, 'r_pdeaths_predict', outdir.table)
  if(length(Code) > 6){
    mid_code = round(length(Code) / 2)
    plot_PPC_relative_resurgence(subset(log_r_pdeaths, code %in% Code[1:mid_code]), subset(log_r_pdeaths_predict, code %in% Code[1:mid_code]), '_log_part_1', outdir.fig)
    plot_PPC_relative_resurgence(subset(log_r_pdeaths, code %in% Code[(1+mid_code):(mid_code*2)]), subset(log_r_pdeaths_predict, code %in% Code[(1+mid_code):(mid_code*2)]), '_log_part_2', outdir.fig)
    plot_PPC_relative_resurgence(subset(r_pdeaths, code %in% Code[1:mid_code]), subset(r_pdeaths_predict, code %in% Code[1:mid_code]), '_part_1', outdir.fig)
    plot_PPC_relative_resurgence(subset(r_pdeaths, code %in% Code[(1+mid_code):(mid_code*2)]), subset(r_pdeaths_predict, code %in% Code[(1+mid_code):(mid_code*2)]), '_part_2', outdir.fig)
  }
  else{
    plot_PPC_relative_resurgence(log_r_pdeaths, log_r_pdeaths_predict, '_log', outdir.fig)
    plot_PPC_relative_resurgence(r_pdeaths, r_pdeaths_predict, '', outdir.fig)
  }

  plot_slope_by_prop(log_r_pdeaths, prop_vac, outdir.fig)
}



# compare to DoH data
tmp <- find_cumulative_deaths_prop_givensum_state_age_multiple_states(fit_cum, date_10thcum, df_week, df_age_continuous, 
                                                                           scrapedData, 'cum.deaths', Code, outdir.table)
compare_CDCestimation_DoH_age_plot(CDC_data = copy(tmp), scraped_data = scrapedData,
                                   var.cum.deaths.CDC = c('M_abs_cum', 'CL_abs_cum', 'CU_abs_cum'), outdir = outdir.fig)
compare_CDCestimation_DoH_age_prop_plot(copy(tmp), outdir.fig)
compare_CDCestimation_DoH_age_weekly_plot(copy(tmp), outdir.fig)

# make panel figure
age_contribution_continuous_table$method = 'GP-BS-SE'
deatht$method = 'GP-BS-SE'
p2=plot_contribution_continuous_comparison_method_with_data(copy(age_contribution_continuous_table), copy(deatht), copy(data),
                                                  'GP-BS-SE', 'GP-BS-SE',
                                                  show.method = F,
                                                  heights = c(1,1), outdir.fig)



cat("\n End postprocessing_figures.R \n")



