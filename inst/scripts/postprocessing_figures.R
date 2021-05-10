
cat("\n Begin postprocessing_figures.R \n")

library(rstan)
library(data.table)
library(dplyr)

indir = "~/git/CDC-covid19-agespecific-mortality-data/inst" # path to the repo
outdir = '/rds/general/user/mm3218/home/git/CDC-covid19-agespecific-mortality-data/inst/results/'
location.index = 1
stan_model = "210505b1"
JOBID = 11377

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

# Plot estimated CAR covariance matrix
plot_covariance_matrix(fit_cum, outdir = outdir.fig)


# Plot estimate B-splines parameters plane 
plot_posterior_plane(fit_cum, df_week, stan_data, outdir = outdir.fig)


# Plots continuous age distribution alpha
cat("\nMake continuous age distribution plots \n")
age_contribution_continuous_table = make_var_by_age_table(fit_cum, df_week, df_age_continuous, 'phi', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_continuous_table, 'phi', "weekly COVID-19 deaths", outdir = outdir.fig)
age_contribution_discrete_table = make_var_by_age_table(fit_cum, df_week, df_age_reporting, 'phi_reduced', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_discrete_table, 'phi_reduced', "weekly COVID-19 deaths", outdir = outdir.fig, discrete = T)


# Plot imputed weekly data 
death_continuous_table = make_var_by_age_table(fit_cum, df_week, df_age_continuous, 'deaths_predict', outdir.table)
plot_var_by_age(death_continuous_table, 'deaths_predict', data, outdir.fig)
death_discrete_table = make_var_by_age_table(fit_cum, df_week, df_age_reporting, 'deaths_predict_state_age_strata', outdir.table)
plot_var_by_age(death_discrete_table, 'deaths_predict_state_age_strata', data, outdir.fig, discrete = T)


# Plot probability ratio of deaths over time
probability_ratio_table = make_probability_ratio_table(fit_cum, df_week, df_age_reporting, data, stan_data, outdir.table)
plot_probability_ratio(probability_ratio_table, df_week, stan_data, outdir.fig)


# Plot deaths ratio of deaths over time
death_ratio_table = make_death_ratio_table(fit_cum, df_week, df_age_reporting, data, outdir.table)
plot_death_ratio(death_ratio_table, outdir.fig)

# Plot mean age of death over time 
mean_age_death = find_mean_age_death(fit_cum, df_week, outdir.table)
plot_mean_age_death(mean_age_death, outdir.fig)

# plot compare to JHU and Imperial data
# find overall cumulative deaths (across age groups)
tmp = find_overall_cumulative_deaths(fit_cum, df_week, 'deaths_predict')
compare_CDCestimation_JHU_Imperial_plot(CDC_data = copy(tmp), JHU_data = JHUData, scraped_data = scrapedData,
                                        var.cum.deaths.CDC = 'M', outdir = outdir.fig)

# find overall cumulative deaths (by age groups)
if(nrow(subset(scrapedData, code == Code)) > 0 ){
  tmp = find_cumulative_deaths_state_age(fit_cum, df_week, df_age_continuous, unique(subset(scrapedData, code == Code)$age), 'deaths_predict')
  compare_CDCestimation_Imperial_age_plot(CDC_data = copy(tmp), scraped_data = scrapedData, 
                                          var.cum.deaths.CDC = 'M', outdir = outdir.fig)
  
  tmp = find_cumulative_deaths_givensum_state_age(fit_cum, df_week, df_age_continuous, scrapedData)
  compare_CDCestimation_Imperial_age_plot(CDC_data = copy(tmp), scraped_data = scrapedData, 
                                          var.cum.deaths.CDC = 'M', outdir = outdir.fig)
  
  tmp = find_phi_state_age(fit_cum, df_week, df_age_continuous, unique(subset(scrapedData, code == Code)$age))
  compare_CDCestimation_Imperial_age_prop_plot(CDC_data = copy(tmp), scraped_data = scrapedData, 
                                          var.cum.deaths.CDC = 'M', df_week, outdir = outdir.fig)
}


cat("\n End postprocessing_figures.R \n")


