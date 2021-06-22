
cat("\n Begin postprocessing_figures.R \n")

library(rstan)
library(data.table)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(cowplot)

indir = "~/git/covid19Vaccination/inst" # path to the repo
outdir = '/rds/general/user/mm3218/home/git/covid19Vaccination/inst/results/'
location.index = 47
stan_model = "210529b"
JOBID = 2117

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

# temporary
# path to population count 
path.to.pop.data = file.path(indir, "data", paste0("us_population_withnyc.rds"))
pop_data = as.data.table( reshape2::melt( readRDS(path.to.pop.data), id.vars = c('Region', 'code', 'Total')) )
setnames(pop_data, c('Region', 'variable', 'value'), c('loc_label', 'age', 'pop'))

# vaccination data 
file  = file.path(indir, 'data', 'demographic_trends_of_people_receiving_covid19_vaccinations_in_the_united_states_210520.csv')
vaccinedata = clean_vaccination_data_age(file)
file = file.path(indir, 'data', 'us_state_vaccinations_210611.csv')
vaccinedata_state = clean_vaccination_data_state(file)

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


# Plot estimated CAR covariance matrix
plot_covariance_matrix(fit_cum, outdir = outdir.fig)

# Plot estimate B-splines parameters plane 
plot_posterior_plane(fit_cum, df_week, df_age_continuous, stan_data, outdir = outdir.fig)


# Plots continuous age distribution phi
cat("\nMake continuous age distribution plots \n")
age_contribution_continuous_table = make_var_by_age_table(fit_cum, df_week, df_age_continuous, 'phi', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_continuous_table, 'phi', outdir = outdir.fig)
age_contribution_discrete_table = make_var_by_age_table(fit_cum, df_week, df_age_reporting, 'phi_reduced', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_discrete_table, 'phi_reduced', outdir = outdir.fig, discrete = T)

# baseline contribution
make_contribution_ref(fit_cum, date_10thcum, fiveagegroups, data, df_week, df_age_continuous, outdir.table)
make_contribution_ref_adj(fit_cum, date_10thcum, fiveagegroups, df_week, pop_data, outdir.table)

# contirbution over time per age groups
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '0-64', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '0-74', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '0-54', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '55-74', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '65-79', date_10thcum, pop_data, data, outdir.table)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '65-74', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '75-84', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '65+', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '75+', date_10thcum, pop_data, data, outdir.table, with_empirical = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '80+', date_10thcum, pop_data, data, outdir.table)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, df_age_reporting, '85+', date_10thcum, pop_data, data, outdir.table, with_empirical = T)

# Plot mortality rate
mortality_rate_table = make_mortality_rate_table(fit_cum, fiveagegroups, date_10thcum, df_week, pop_data, 
                                                 JHUData, df_age_continuous, 'cumulative_deaths' , outdir.table)
plot_mortality_rate(mortality_rate_table, outdir.fig)

# death ratio relative to baseline 
death_ratio_winter = make_weekly_death_rate_table(fit_cum, c('0-74', '75+'), as.Date('2020-12-21'), df_week, 
                             JHUData, data, df_age_continuous, 'cumulative_deaths' , outdir.table)
plot_death_ratio_winter(death_ratio_winter, vaccinedata, outdir.fig)

# Plot probability ratio of deaths over time
probability_ratio_table = make_probability_ratio_table(fit_cum, df_week, df_age_reporting, data, stan_data, outdir.table)
plot_probability_ratio(probability_ratio_table, df_week, stan_data, outdir.fig)


# Plot deaths ratio of deaths over time
death_ratio_table = make_death_ratio_table(fit_cum, df_week, df_age_reporting, data, outdir.table)
plot_death_ratio(death_ratio_table, outdir.fig)

# Plot imputed weekly data 
death_continuous_table = make_var_by_age_table(fit_cum, df_week, df_age_continuous, 'deaths_predict', outdir.table)
plot_imputed_deaths_by_age(death_continuous_table, 'deaths_predict', data, outdir.fig)
death_discrete_table = make_var_by_age_table(fit_cum, df_week, df_age_reporting, 'deaths_predict_state_age_strata', outdir.table)
plot_imputed_deaths_by_age(death_discrete_table, 'deaths_predict_state_age_strata', data, outdir.fig, discrete = T)

deatht = make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'phi', df_age_continuous, outdir.table)
tmp = make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'phi_reduced', df_age_reporting, outdir.table, withempirical = T)
make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'phi', df_age_continuous, outdir.table, 
                                          age_groups = c('0-54', '55-74', '75+'), lab = '3agegroups', withempirical = T,
                                    reduction = c(min(vaccinedata_state[date %in% df_week$date, date]), max(df_week$date)))
make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'phi', df_age_continuous, outdir.table, 
                                    age_groups = c('0-74', '75+'), lab = '2agegroups', withempirical = T,
                                    reduction = c(min(vaccinedata_state[date %in% df_week$date, date]), max(df_week$date)))
make_weekly_death_rate_other_source_posteriorsamples(fit_cum, df_week, JHUData,  'phi', df_age_continuous, outdir.table, 
                                                     age_groups = c('0-54', '55-74', '75+'), lab = '3agegroups', 
                                                     reduction = c(min(vaccinedata_state[date %in% df_week$date, date]), max(df_week$date)))
make_weekly_death_rate_other_source_posteriorsamples(fit_cum, df_week, JHUData,  'phi', df_age_continuous, outdir.table, 
                                                     age_groups = c('0-74', '75+'), lab = '2agegroups', 
                                                     reduction = c(min(vaccinedata_state[date %in% df_week$date, date]), max(df_week$date)))
# Plot mean age of death over time 
mean_age_death = find_mean_age_death(fit_cum, df_week, outdir.table)
plot_mean_age_death(mean_age_death, outdir.fig)


# find overall cumulative deaths (by age groups)
if(nrow(subset(scrapedData, code == Code)) > 0 ){
  
  scrapedData = subset(scrapedData, code == Code)
  
  if(Code == 'GA')
    scrapedData = reduce_agebands_scrapedData_GA(scrapedData)
  
  tmp = find_cumulative_deaths_prop_givensum_state_age(fit_cum, date_10thcum, df_week, df_age_continuous, scrapedData, 'cum.deaths', outdir.table)
  compare_CDCestimation_DoH_age_plot(CDC_data = copy(tmp), scraped_data = scrapedData, 
                                          var.cum.deaths.CDC = c('M_abs_cum', 'CL_abs_cum', 'CU_abs_cum'), outdir = outdir.fig)
  compare_CDCestimation_DoH_age_prop_plot(copy(tmp), outdir.fig)
  compare_CDCestimation_DoH_age_weekly_plot(copy(tmp), outdir.fig)
  
  tmp = make_weekly_death_rate_other_source(fit_cum, df_week, JHUData,  'phi', df_age_continuous, outdir.table, 
                                            age_groups = unique(scrapedData$age), lab = 'DoH', cumulative = T)
  compare_CDCestimation_DoH_age_plot(CDC_data = copy(tmp), scraped_data = scrapedData, 
                                     var.cum.deaths.CDC = c('M', 'CL', 'CU'), outdir = outdir.fig)

 }

# make panel figure

age_contribution_continuous_table$method = 'GP-BS-SE'
deatht$method = 'GP-BS-SE'
p2=plot_contribution_continuous_comparison_method(copy(age_contribution_continuous_table), copy(deatht), copy(data), 
                                                  'GP-BS-SE', 'GP-BS-SE', 
                                                  show.method = F, 
                                                  heights = c(1,1)) 

ggsave(p2, file = paste0(outdir.fig, '-panel_plot_1_', Code, '.png'), w = 9, h = 7)




# data = select(data, date, age, loc_label, code, weekly.deaths)
# tmp1 = data[, list(total_deaths = sum(na.omit(weekly.deaths))), by = 'date']
# data = merge(data, tmp1, by = 'date')
# data[, prop_deaths := weekly.deaths / total_deaths]
# data$method = 'observation'
# 



cat("\n End postprocessing_figures.R \n")


