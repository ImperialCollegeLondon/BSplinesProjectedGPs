
cat("\n Begin postprocessing_figures.R \n")

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
location.index = 5
stan_model = "210529b"
JOBID = 29051

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
  location.index <- as.numeric(args_line[[6]])
  stan_model <- args_line[[8]]
  JOBID <- as.numeric(args_line[[10]])
}


# load functions
source(file.path(indir, "functions", "postprocessing-plotting_functions.R"))
source(file.path(indir, "functions", "postprocessing-summary_functions.R"))
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
Code = locations[location.index,]$code

# load image 
load(file.path(outdir.data, paste0("stanin_", Code, "_",run_tag,".RData")))
outdir.fig = outdir.fig.post
outdir.fit = outdir.fit.post

# load fit cumulative deaths
cat("Load fits \n")
file = file.path(outdir.fit, paste0("fit_cumulative_deaths_", Code, "_", run_tag,".rds"))
fit_cum <- readRDS(file=file)

# find date for the first 10th cumulative deaths
date_10thcum = subset(data, !is.na(weekly.deaths))
date_10thcum = date_10thcum[, list(weekly.deaths = sum(weekly.deaths)), by = 'date']
date_10thcum[, cum.deaths := cumsum(weekly.deaths)]
date_10thcum = date_10thcum[ cum.deaths >=10, min(date)]
cat("The first date with >= 10th cum deaths is ", as.character(date_10thcum))
fouragegroups = c('0-24', '25-54', '55-79', '80+')
fiveagegroups = c('0-24', '25-54', '55-74', '75-84', '85+')

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
                                    reduction = c(min(vaccine_data[date %in% df_week$date, date]), max(df_week$date)))
make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                    age_groups = c('0-74', '75+'), lab = '2agegroups', withempirical = T,
                                    reduction = c(min(vaccine_data[date %in% df_week$date, date]), max(df_week$date)))
make_weekly_death_rate_other_source_posteriorsamples(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                                     age_groups = c('0-54', '55-74', '75+'), lab = '3agegroups',
                                                     reduction = c(min(vaccine_data[date %in% df_week$date, date]), max(df_week$date)))
make_weekly_death_rate_other_source_posteriorsamples(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                                     age_groups = c('0-74', '75+'), lab = '2agegroups',
                                                     reduction = c(min(vaccine_data[date %in% df_week$date, date]), max(df_week$date)))
make_weekly_death_rate_other_source_posteriorsamples(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                                     age_groups = unique(df_age_vaccination$age), lab = 'vacagegroups')


# fcompare to DoH data
if(nrow(subset(scrapedData, code == Code)) > 0 ){
  
  scrapedData = subset(scrapedData, code == Code)
  
  if(Code == 'GA')
    scrapedData = reduce_agebands_scrapedData_GA(scrapedData)
  
  tmp = find_cumulative_deaths_prop_givensum_state_age(fit_cum, date_10thcum, df_week, df_age_continuous, scrapedData, 'cum.deaths', outdir.table)
  compare_CDCestimation_DoH_age_plot(CDC_data = copy(tmp), scraped_data = scrapedData,
                                     var.cum.deaths.CDC = c('M_abs_cum', 'CL_abs_cum', 'CU_abs_cum'), outdir = outdir.fig)
  compare_CDCestimation_DoH_age_prop_plot(copy(tmp), outdir.fig)
  compare_CDCestimation_DoH_age_weekly_plot(copy(tmp), outdir.fig)
  
}


# plot vaccine effects
names <- names(fit_cum)[grepl('gamma', names(fit_cum)) & !grepl('gamma_re', names(fit_cum)) & !grepl('gamma0', names(fit_cum))]
if(length(names) != 0){
  p <- mcmc_areas(fit_cum, pars = names, prob = 0.5, prob_outer = 0.95)
  ggsave(p, file = paste0(outdir.fig, '-vaccine_effects_', Code, '.png'), h = 5, w = 6)
}

if(any(grepl('delta', names(fit_cum)))){
  p <- mcmc_areas(fit_cum, regex_pars = 'delta', prob = 0.5, prob_outer = 0.95)
  ggsave(p, file = paste0(outdir.fig, '-vaccine_effects_others_', Code, '.png'), h = 5, w = 6)
  
}

if( any(grepl('phi_wo_vaccine', names(fit_cum))) ){
  vaccine_effects = find_vaccine_effects_scaled(fit_cum, df_week, df_age_continuous, df_age_vac_effects$age, 
                                         'phi', c('', '_wo_vaccine'), outdir.table)
  phi = make_var_by_varying_age_table(fit_cum, df_week, df_age_continuous, df_age_vac_effects$age, 'phi', 'sum', outdir.table)
  phi_wo_vaccine = make_var_by_varying_age_table(fit_cum, df_week, df_age_continuous, df_age_vac_effects$age, 'phi_wo_vaccine', 'sum',outdir.table)
  plot_vaccine_effects2(vaccine_effects, phi, phi_wo_vaccine, 'phi', outdir.fig)
}

if( any(grepl('f_w_vaccine', names(fit_cum))) ){
  vaccine_effects = find_vaccine_effects_unscaled(fit_cum, df_week, df_age_continuous, df_age_vac_effects$age, 
                                         'f', c('_w_vaccine', ''), outdir.table)
  f = make_var_by_varying_age_table(fit_cum, df_week, df_age_continuous, df_age_vac_effects$age, 'f', 'mean', outdir.table)
  f_w_vaccine = make_var_by_varying_age_table(fit_cum, df_week, df_age_continuous, df_age_vac_effects$age, 'f_w_vaccine', 'mean', outdir.table)
  plot_vaccine_effects2(vaccine_effects, f_w_vaccine, f, 'f', outdir.fig)
}

if( any(grepl('f_w_vaccine_extra', names(fit_cum))) ){
  prop_tab = data.table(prop = stan_data$extra_prop_vaccinated, prop_index = 1:length(stan_data$extra_prop_vaccinated))
  vaccine_effects_extrapolated <- find_vaccine_effects_unscaled_extrapolated(fit_cum, prop_tab, df_age_continuous, 
                                                                             df_age_vac_effects$age, 'f', c('_w_vaccine_extra', ''), outdir)
  tmp = subset(vaccine_effects_extrapolated, age_index != min(age_index))
  ggplot(tmp, aes(x = prop, y = M)) +
    geom_line(aes(col = age), size = 1) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) +
    labs(x = 'Proportion of vaccinated', 
         y = paste0("Vaccine effect\n(% estimate change unscaled contribution to weekly deaths)")) +
    # scale_color_gradient2(mid = 'burlywood1', low = 'darkblue', midpoint = 1) +
    # scale_fill_gradient2(mid = 'burlywood1', low = 'darkblue', midpoint = 1) +
    scale_color_viridis_d(option = 'A', end = 0.9) + 
    scale_fill_viridis_d(option = 'A',end = 0.9) + 
    theme_bw() +
    scale_x_continuous(expand = c(0,0), labels = scales::percent_format()) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent_format()) +
    theme(legend.position = 'bottom', 
          legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white"))
  ggsave(paste0(outdir.fig, paste0('-vaccine_effects_extrapolation_', Code, '.png')), w = 5, h = 5)
  
}

weeklydv <- make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'alpha', df_age_continuous, outdir.table,
                                                age_groups = df_age_vaccination$age, lab = 'vacagegroups', withempirical = T,
                                                reduction = NULL)
weeklyf <- find_contribution_age_groups_vaccination(fit_cum, df_week, df_age_continuous, df_age_reporting, 
                                                    deathByAge, df_age_reporting$age, 
                                                    c(1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4), 'f',F, outdir.table)
weeklyphi <- find_contribution_age_groups_vaccination(fit_cum, df_week, df_age_continuous, df_age_reporting, 
                                                      deathByAge, df_age_vaccination$age,
                                                      c(1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4), 'phi',T, outdir.table)
plot_vaccine_effects(vaccine_data, weeklydv, weeklyf, weeklyphi, outdir.fig)
saveRDS(vaccine_data, paste0(outdir.table, '-', 'vaccine_data', '.rds'))


# make panel figure
age_contribution_continuous_table$method = 'GP-BS-SE'
deatht$method = 'GP-BS-SE'
p2=plot_contribution_continuous_comparison_method(copy(age_contribution_continuous_table), copy(deatht), copy(data), 
                                                  'GP-BS-SE', 'GP-BS-SE', 
                                                  show.method = F, 
                                                  heights = c(1,1)) 
ggsave(p2, file = paste0(outdir.fig, '-panel_plot_1_', Code, '.png'), w = 9, h = 7)


cat("\n End postprocessing_figures.R \n")



