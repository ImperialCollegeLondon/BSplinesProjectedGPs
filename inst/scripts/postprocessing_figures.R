
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
library(magick)
library(ggrepel)

indir ="~/git/BSplinesProjectedGPs/inst" # path to the repo
outdir = file.path('/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst', "results")
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
file = file.path(outdir.fit.post, paste0("posterior_samples_", run_tag,".rds"))
if(file.exists(file)){
  fit_samples <- readRDS(file=file)
}else{
  fit_samples <- rstan::extract(fit_cum)
}

# find date for the first 10th cumulative deaths
date_10thcum = subset(data, !is.na(weekly.deaths))
date_10thcum = date_10thcum[, list(weekly.deaths = sum(weekly.deaths)), by = c('date', 'code')]
date_10thcum[, cum.deaths := cumsum(weekly.deaths), by = 'code']
date_10thcum = date_10thcum[ cum.deaths >=10, list(date_10thcum = min(date)), by = 'code']

# age groups
fouragegroups = c('0-24',  '25-54', '55-84', '85+')
fiveagegroups = c('0-24', '25-54', '55-74', '75-84', '85+')

# selected states
selected_code = c('CA', 'FL', 'NY', 'TX')
selected_6_codes = c('CA', 'FL', 'NY', 'TX', 'PA')
selected_10_codes = c('CA','FL','NY','TX','PA','IL','OH','GA','NC','MI')

# table for plotting vaccine effects
min_age_index_vac = df_age_vaccination[age == '18-64']$age_index
df_age_vaccination2 = df_age_vaccination[age_index >= min_age_index_vac]
df_age_vaccination2[, age_index := age_index - min_age_index_vac + 1]

### temporary
if('Total' %in% colnames(pop_data)){
  path.to.pop.data = file.path(indir, "data", paste0("us_population.csv"))
  pop_data = as.data.table( read.csv(path.to.pop.data) )
  setnames(pop_data, 'location', 'loc_label')
}
if('age_index' %in% colnames(pop_data)){
  pop_data = select(pop_data, -age_index)
}


####

# Predictions deaths by 1-year age band
make_var_by_age_by_state_by_time_table(fit_samples, df_week, df_age_continuous, df_state, 'deaths_predict', outdir.table)

# Plots continuous and aggregated age distribution phi
age_contribution_continuous_table = make_var_by_age_by_state_by_time_table(fit_samples, df_week, df_age_continuous, df_state, 'phi', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_continuous_table, 'phi', outdir = outdir.fig)
age_contribution_discrete_table = make_var_by_age_by_state_by_time_table(fit_samples, df_week, df_age_reporting, df_state, 'phi_reduced', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_discrete_table, 'phi_reduced', outdir = outdir.fig, discrete = T)

if('phi_reduced_vac' %in% names(fit_samples)){
  make_var_by_age_by_state_by_time_table(fit_samples, df_week, df_age_vaccination2, df_state, 'phi_reduced_vac', outdir.table)
  make_var_by_age_by_state_by_time_diff_table(fit_samples, df_week, df_age_vaccination2, df_state, 'phi_reduced_vac', vaccine_data_pop, outdir.table)
  make_var_by_age_by_state_by_time_diff_over_time_table(fit_samples, df_week, df_age_vaccination2, df_state, 'phi_reduced_vac', outdir.table)
}

if('phi_predict_reduced_vac' %in% names(fit_samples)){
  make_var_by_age_by_state_by_time_table(fit_samples, df_week, df_age_vaccination2, df_state, 'phi_predict_reduced_vac', outdir.table)
  make_var_by_age_by_state_by_time_diff_table(fit_samples, df_week, df_age_vaccination2, df_state, 'phi_predict_reduced_vac', vaccine_data_pop, outdir.table)
  make_var_by_age_by_state_by_time_diff_over_time_table(fit_samples, df_week, df_age_vaccination2, df_state, 'phi_predict_reduced_vac', outdir.table)
}

# baseline contribution adjusted and non-adjusted for population composition
make_contribution_ref(fit_samples, date_10thcum, fiveagegroups, data, df_week, df_age_continuous, outdir.table)
make_contribution_ref_adj(fit_samples, date_10thcum, fiveagegroups, df_week, pop_data, outdir.table)

# contribution over time per age groups
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '0-64', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '0-74', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '0-54', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '20-64', date_10thcum, pop_data, data, outdir.table, with_empirical = F)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '55-74', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '65-79', date_10thcum, pop_data, data, outdir.table)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '65-74', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '75-84', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '65+', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '75+', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '80+', date_10thcum, pop_data, data, outdir.table)
find_contribution_one_age_group(fit_samples, df_week, df_age_continuous, df_age_reporting, '85+', date_10thcum, pop_data, data, outdir.table, with_empirical = T)


# mortality rate
mortality_rate_table = make_mortality_rate_table_discrete(fit_samples, fiveagegroups, date_10thcum, df_week, pop_data,
                                                 JHUData, df_age_continuous, 'cumulative_deaths' , nyt_data, outdir.table)
make_mortality_rate_table_discrete(fit_samples, fouragegroups, date_10thcum, df_week, pop_data, JHUData, df_age_continuous,
                                   'cumulative_deaths' , nyt_data, outdir.table, '4agegroups')
mortality_rate_table_continuous = make_mortality_rate_table_continuous(fit_samples, date_10thcum, df_week, pop_data,
                                                          JHUData, df_age_continuous, 'cumulative_deaths' , outdir.table)
plot_mortality_rate(mortality_rate_table, mortality_rate_table_continuous, outdir.fig)

# predicted weekly deaths by various age groups
deatht = make_weekly_death_rate_other_source(fit_samples, df_week, JHUData,  'alpha', df_age_continuous, outdir.table, lab = 'prediction')
tmp = make_weekly_death_rate_other_source(fit_samples, df_week, JHUData,  'alpha_reduced', df_age_reporting, outdir.table, withempirical = T, lab = 'prediction_agg')
make_weekly_death_rate_other_source(fit_samples, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                    age_groups = c('0-54', '55-74', '75+'), lab = '3agegroups', withempirical = T)
make_weekly_death_rate_other_source(fit_samples, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                    age_groups = c('0-74', '75+'), lab = '2agegroups', withempirical = T)

# compare to DoH data
if(any(Code %in% unique(scrapedData$code))){
  tmp <- find_cumulative_deaths_prop_givensum_state_age_multiple_states(fit_samples, date_10thcum, df_week, df_age_continuous,
                                                                        scrapedData, 'cum.deaths', Code, outdir.table)
  compare_CDCestimation_DoH_age_plot(CDC_data = copy(tmp), scraped_data = scrapedData,
                                     var.cum.deaths.CDC = c('M_abs_cum', 'CL_abs_cum', 'CU_abs_cum'), outdir = outdir.fig)
  compare_CDCestimation_DoH_age_prop_plot(copy(tmp), outdir.fig)
  compare_CDCestimation_DoH_age_weekly_plot(copy(tmp), outdir.fig)
}

# make panel figure
age_contribution_continuous_table$method = 'GP-BS-SE'
deatht$method = 'GP-BS-SE'
p2=plot_contribution_continuous_comparison_method_with_data(copy(age_contribution_continuous_table), copy(deatht), copy(data),
                                                            'GP-BS-SE', 'GP-BS-SE',
                                                            show.method = F,
                                                            heights = c(1,1), outdir.fig)

# resurgence deaths 
save_resurgence_dates(resurgence_dates, outdir.table)

df_week[, dummy := 1]; resurgence_dates[, dummy := 1]
df_week2 = merge(df_week, resurgence_dates, by = 'dummy', allow.cartesian=TRUE)
df_week2 = df_week2[date >= start_resurgence & date <= stop_resurgence]
df_week2[, week_index := 1:length(date), by = 'code']
df_week2 = select(df_week2, - start_resurgence, - stop_resurgence, -dummy)
df_week2 = df_week2[order(code)]
df_week = select(df_week, -dummy)

stan_data1 = add_vaccine_prop(stan_data, df_week, Code, vaccine_data, resurgence_dates)
prop_vac = prepare_prop_vac_table(stan_data1, vaccine_data, df_age_vaccination)

# plot estimate relative deaths
if('r_pdeaths' %in% names(fit_samples)){
  r_pdeaths = make_var_by_age_by_state_by_time_table(fit_samples, df_week2, df_age_vaccination2, df_state, 'r_pdeaths', outdir.table, T)
  r_pdeaths = merge(r_pdeaths, prop_vac, by = c('code', 'date'))
  for(Code in unique(r_pdeaths$code)){
    saveRDS(subset(r_pdeaths, code == Code), file = paste0(outdir.table, '-', 'r_pdeaths',  'Table_', Code, '.rds'))
  }
}

# vaccination analysis
if('intercept_resurgence0' %in% names(fit_samples)){
  
  # # p-value vaccine effects
  names_var = c('slope_resurgence0', 'vaccine_effect_slope_cross', 'intercept_resurgence0', 'vaccine_effect_intercept_cross', 'vaccine_effect_intercept_diagonal', 'vaccine_effect_slope_diagonal')
  names_var <- names_var[names_var %in% names(fit_samples)]
  save_p_value_vaccine_effects(fit_samples, names_var, outdir.table)
  
  prop_vac_counterfactual <- prepare_prop_vac_counterfactual_table(stan_data1, df_state, df_age_vaccination2, df_counterfactual)
  
  max_1864 = round(prop_vac[,list(date = min(date), prop_1 = prop_1[which.min(date)]), by = 'code'][,max(prop_1)] * 100)
  min_1864 = round(prop_vac[,list(date = min(date), prop_1 = prop_1[which.min(date)]), by = 'code'][,min(prop_1)] * 100)
  max_65p = round(prop_vac[,list(date = min(date), prop_2 = prop_2[which.min(date)]), by = 'code'][,max(prop_2)] * 100)
  min_65p = round(prop_vac[,list(date = min(date), prop_2 = prop_2[which.min(date)]), by = 'code'][,min(prop_2)] * 100)
  
  df_counterfactual = data.table(counterfactual_index = 1:6, 
                                 label_counterfactual = c(paste0('Counterfactual vaccination rate of ', max_1864, '% in 18-64'), 
                                                          paste0('Counterfactual vaccination rate of ', max_65p, '% in 65+'), 
                                                          paste0('Counterfactual vaccination rate of ', max_1864, '% in 18-64 and ', max_65p, '% in 65+'), 
                                                          paste0('Counterfactual vaccination rate of ', min_1864, '% in 18-64'), 
                                                          paste0('Counterfactual vaccination rate of ', min_65p, '% in 65+'), 
                                                          paste0('Counterfactual vaccination rate of ', min_1864, '% in 18-64 and ', min_65p, '% in 65+')))
  df_counterfactual[, label_counterfactual := factor(label_counterfactual, 
                                                     levels = df_counterfactual$label_counterfactual[c(4,1,5,2,6,3)])]

  # plot relative deaths
  # plot_relative_resurgence_vaccine(r_pdeaths, prop_vac, df_age_vaccination2, df_week2, resurgence_dates,outdir.fig)
  plot_relative_resurgence_vaccine2(r_pdeaths,  F, outdir.fig)

  # relative deaths ratio
  r_pdeaths_ratio <- make_var_by_age_by_state_by_time_table_relative(fit_samples, df_week2, df_age_vaccination2, df_state, 'r_pdeaths', 'FL', outdir.table)
  r_pdeaths_ratio[, max_week_index := max(week_index), by = 'code']
  r_pdeaths_ratio <- r_pdeaths_ratio[week_index == max_week_index]
  print(r_pdeaths_ratio)
  
  # plot predicted relative deaths
  r_pdeaths_predict = make_var_by_age_by_state_by_time_table(fit_samples, df_week2, df_age_vaccination2, df_state, 'r_pdeaths_predict', outdir.table)
  plot_PPC_relative_resurgence(r_pdeaths, r_pdeaths_predict, '', outdir.fig)
  
  # plot counterfactual analysis
  E_pdeaths_predict_resurgence_cumulative = make_var_by_age_by_state_by_time_table(fit_samples, df_week2, df_age_vaccination2, df_state, 'E_pdeaths_predict_resurgence_cumulative', outdir.table)
  E_pdeaths_counterfactual_resurgence_cumulative = make_var_by_age_by_state_by_counterfactual_table(fit_samples, df_week2, df_age_vaccination2, df_state, df_counterfactual, 'E_pdeaths_counterfactual_resurgence_cumulative', outdir.table)
  diff_E_pdeaths_counterfactual = make_inv_var_by_age_by_state_by_counterfactual_table(fit_samples, df_week2, df_age_vaccination2, df_state, df_counterfactual, 'diff_E_pdeaths_counterfactual', outdir.table)
  perc_E_pdeaths_counterfactual = make_inv_var_by_age_by_state_by_counterfactual_table(fit_samples, df_week2, df_age_vaccination2, df_state, df_counterfactual, 'perc_E_pdeaths_counterfactual', outdir.table)
  # perc_E_pdeaths_counterfactual_2 <- make_ratio_vars_by_age_state_by_counterfactual_table(fit_samples, df_week2, df_age_vaccination2, df_state, df_counterfactual, c('E_pdeaths_counterfactual_resurgence_cumulative', 'E_pdeaths_predict_resurgence_cumulative'), outdir.table)
  
  # aggregate across states
  if('diff_E_pdeaths_counterfactual_all' %in% names(fit_samples)){
    diff_E_pdeaths_counterfactual_all <- make_var_by_age_by_counterfactual_table(fit_samples, df_week2, df_age_vaccination2, df_counterfactual, 'diff_E_pdeaths_counterfactual_all', outdir.table)
    perc_E_pdeaths_counterfactual_all <- make_var_by_age_by_counterfactual_table(fit_samples, df_week2, df_age_vaccination2, df_counterfactual, 'perc_E_pdeaths_counterfactual_all', outdir.table)
    diff_E_pdeaths_counterfactual_all_2 <- make_inv_var_by_age_by_counterfactual_table(fit_samples, df_week2, df_age_vaccination2, df_counterfactual, 'diff_E_pdeaths_counterfactual_all', outdir.table)
    diff_E_pdeaths_counterfactual_allstatesages <- make_inv_var_by_counterfactual_table(fit_samples, df_counterfactual, 'diff_E_pdeaths_counterfactual_all', outdir.table)
  } else{
    diff_E_pdeaths_counterfactual_all <- find_diff_E_pdeaths_counterfactual_all(fit_samples, df_week2, df_age_vaccination2, df_counterfactual, outdir.table)
    perc_E_pdeaths_counterfactual_all <- find_perc_E_pdeaths_counterfactual_all(fit_samples, df_week2, df_age_vaccination2, df_counterfactual, outdir.table)
    diff_E_pdeaths_counterfactual_all_2 <- find_inv_diff_E_pdeaths_counterfactual_all(fit_samples, df_week2, df_age_vaccination2, df_counterfactual, outdir.table)
    diff_E_pdeaths_counterfactual_allstatesages <- find_inv_diff_E_pdeaths_counterfactual_allstatesages(fit_samples, df_week2, df_counterfactual, outdir.table)
  }
  perc_E_pdeaths_counterfactual_all_2 <- make_ratio_vars_by_age_by_counterfactual_table(fit_samples, df_age_vaccination2, df_counterfactual, c('E_pdeaths_counterfactual_resurgence_cumulative', 'E_pdeaths_predict_resurgence_cumulative'), outdir.table)
  
  # aggregate across ages
  diff_E_pdeaths_counterfactual_allages = make_var_inv_by_state_by_counterfactual_table(fit_samples, df_week2, df_state, df_counterfactual, 'diff_E_pdeaths_counterfactual', outdir.table)
  E_pdeaths_counterfactual_resurgence_cumulative_allages = make_var_by_state_by_counterfactual_table(fit_samples, df_week2, df_state, df_counterfactual, 'E_pdeaths_counterfactual_resurgence_cumulative', outdir.table)
  E_pdeaths_predict_resurgence_cumulative_allages = make_var_by_state_table(fit_samples, df_week2, df_state, 'E_pdeaths_predict_resurgence_cumulative', outdir.table)
  perc_E_pdeaths_counterfactual_allages = make_ratio_vars_by_state_by_counterfactual_table(fit_samples, df_week2, df_state, df_counterfactual, c('E_pdeaths_counterfactual_resurgence_cumulative', 'E_pdeaths_predict_resurgence_cumulative'), outdir.table)
  
  # aggregate across states and ages 
  perc_E_pdeaths_counterfactual_allstatesages <- make_ratio_vars_by_counterfactual_table(fit_samples, df_counterfactual, c('E_pdeaths_counterfactual_resurgence_cumulative', 'E_pdeaths_predict_resurgence_cumulative'), outdir.table)
  
  
  # plot
  plot_vaccine_effects_counterfactual_allages(subset(E_pdeaths_counterfactual_resurgence_cumulative_allages, code %in% selected_code), 
                                              subset(E_pdeaths_predict_resurgence_cumulative_allages, code %in% selected_code), 
                                              subset(resurgence_dates, code %in% selected_code), 'cumulative', 'cumulative_rperiod_selected_states', outdir.fig)
  plot_vaccine_effects_counterfactual_change_allages(subset(perc_E_pdeaths_counterfactual_allages, code %in% selected_10_codes), prop_vac_counterfactual, 'selected_states', 'percChange',  outdir.fig)
  plot_vaccine_effects_counterfactual_change_allages(subset(diff_E_pdeaths_counterfactual_allages, code %in% selected_10_codes), prop_vac_counterfactual, 'selected_states', 'diffChange', outdir.fig)
  
  plot_vaccine_effects_counterfactual_panel(E_pdeaths_counterfactual_resurgence_cumulative, E_pdeaths_predict_resurgence_cumulative, 
                                            perc_E_pdeaths_counterfactual, resurgence_dates, prop_vac_counterfactual, outdir.fig)
  plot_vaccine_effects_counterfactual_panel_repel(perc_E_pdeaths_counterfactual, prop_vac_counterfactual, outdir.fig)
  
  plot_vaccine_effects_counterfactual2(E_pdeaths_counterfactual_resurgence_cumulative, E_pdeaths_predict_resurgence_cumulative, resurgence_dates, 
                                       selected_code, prop_vac, '4_states', outdir.fig)
  plot_vaccine_effects_counterfactual2(E_pdeaths_counterfactual_resurgence_cumulative, E_pdeaths_predict_resurgence_cumulative, resurgence_dates, c('ID', 'CO', 'HI', 'MA'), prop_vac, 'explore', outdir.fig)
  
  find_stats_vaccine_effects(diff_E_pdeaths_counterfactual, perc_E_pdeaths_counterfactual, 
                             diff_E_pdeaths_counterfactual_all_2, perc_E_pdeaths_counterfactual_all_2,
                             diff_E_pdeaths_counterfactual_allages, perc_E_pdeaths_counterfactual_allages,
                             diff_E_pdeaths_counterfactual_allstatesages, perc_E_pdeaths_counterfactual_allstatesages,
                             prop_vac, resurgence_dates, outdir.table)
  
}

cat("\n End postprocessing_figures.R \n")


