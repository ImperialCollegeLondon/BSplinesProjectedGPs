
cat("\n Begin postprocessing_figures.R \n")

library(rstan)
library(data.table)
library(dplyr)
library(gridExtra)

indir = "~/git/CDC-covid19-agespecific-mortality-data/inst" # path to the repo
outdir = '/rds/general/user/mm3218/home/git/CDC-covid19-agespecific-mortality-data/inst/results/'
location.index = 1
stan_model = "210429h1"
JOBID = 10296

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

# temporary
# path to population count 
path.to.pop.data = file.path(indir, "data", paste0("us_population_withnyc.rds"))
pop_data = as.data.table( reshape2::melt( readRDS(path.to.pop.data), id.vars = c('Region', 'code', 'Total')) )
setnames(pop_data, c('Region', 'variable', 'value'), c('loc_label', 'age', 'pop'))

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
ppp = plot_posterior_plane(fit_cum, df_week, df_age_continuous, stan_data, outdir = outdir.fig)


# Plots continuous age distribution phi
cat("\nMake continuous age distribution plots \n")
age_contribution_continuous_table = make_var_by_age_table(fit_cum, df_week, df_age_continuous, 'phi', outdir.table)
pc = plot_probability_deaths_age_contribution(age_contribution_continuous_table, 'phi', outdir = outdir.fig)
age_contribution_discrete_table = make_var_by_age_table(fit_cum, df_week, df_age_reporting, 'phi_reduced', outdir.table)
plot_probability_deaths_age_contribution(age_contribution_discrete_table, 'phi_reduced', outdir = outdir.fig, discrete = T)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, '12+', outdir.table)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, '65+', outdir.table)
find_contribution_one_age_group(fit_cum, df_week, df_age_continuous, '80+', outdir.table)
make_contribution_ref(fit_cum, JHUData, data, df_week, df_age_continuous, outdir.table)
make_contribution_ref_adj(fit_cum, JHUData, data, df_week, pop_data, outdir.table)

# Plot imputed weekly data 
death_continuous_table = make_var_by_age_table(fit_cum, df_week, df_age_continuous, 'deaths_predict', outdir.table)
plot_imputed_deaths_by_age(death_continuous_table, 'deaths_predict', data, outdir.fig)
death_discrete_table = make_var_by_age_table(fit_cum, df_week, df_age_reporting, 'deaths_predict_state_age_strata', outdir.table)
plot_imputed_deaths_by_age(death_discrete_table, 'deaths_predict_state_age_strata', data, outdir.fig, discrete = T)

# Plot mortality rate
mortality_rate_table = make_mortality_rate_table(fit_cum, df_week, pop_data, JHUData, df_age_continuous, 'cumulative_deaths' , outdir.table)
plot_mortality_rate(mortality_rate_table, outdir.fig)


# Plot probability ratio of deaths over time
probability_ratio_table = make_probability_ratio_table(fit_cum, df_week, df_age_reporting, data, stan_data, outdir.table)
plot_probability_ratio(probability_ratio_table, df_week, stan_data, outdir.fig)


# Plot deaths ratio of deaths over time
death_ratio_table = make_death_ratio_table(fit_cum, df_week, df_age_reporting, data, outdir.table)
plot_death_ratio(death_ratio_table, outdir.fig)


# Plot mean age of death over time 
mean_age_death = find_mean_age_death(fit_cum, df_week, outdir.table)
plot_mean_age_death(mean_age_death, outdir.fig)


# plot compare to JHU and DoH data
# find overall cumulative deaths (across age groups)
tmp = find_overall_cumulative_deaths(fit_cum, df_week, 'deaths_predict')
compare_CDCestimation_JHU_DoH_plot(CDC_data = copy(tmp), JHU_data = JHUData, scraped_data = scrapedData,
                                        var.cum.deaths.CDC = 'M', outdir = outdir.fig)

# find overall cumulative deaths (by age groups)
if(nrow(subset(scrapedData, code == Code)) > 0 ){
  
  tmp = find_cumulative_deaths_givensum_state_age(fit_cum, df_week, df_age_continuous, scrapedData, 'cum.deaths')
  compare_CDCestimation_DoH_age_plot(CDC_data = copy(tmp), scraped_data = scrapedData, 
                                          var.cum.deaths.CDC = 'M', outdir = outdir.fig)
  
  tmp = find_phi_state_age(fit_cum, df_week, df_age_continuous, unique(subset(scrapedData, code == Code)$age))
  compare_CDCestimation_DoH_age_prop_plot(CDC_data = copy(tmp), scraped_data = scrapedData, 
                                          var.cum.deaths.CDC = 'M', df_week, outdir = outdir.fig)
}

# make panel figure
data = select(data, date, age, loc_label, code, weekly.deaths)
tmp1 = data[, list(total_deaths = sum(na.omit(weekly.deaths))), by = 'date']
data = merge(data, tmp1, by = 'date')
data[, prop_deaths := weekly.deaths / total_deaths]
data$method = 'observation'

age_contribution_discrete_table$method = 'BS-GP-SE'
p1 = plot_contribution_comparison_method(age_contribution_discrete_table, data, model_name = 'BS-GP-SE')

death_discrete_table$method = 'BS-GP-SE'
p2 = plot_death_comparison_method(death_discrete_table, data, 'BS-GP-SE')

ppp = ppp + theme(legend.position = 'left') + labs(fill = 'B-Splines\nparameters')
p = grid.arrange(p1, p2, pc[[1]], ppp, layout_matrix = rbind(c(1, 1, 1, 3), 
                                                         c(2, 2, 2, 3), 
                                                         c(NA, 4, NA, 3)), widths = c(0.2, 1, 0.2, 0.8), heights = c(1, 1, 0.9))
ggsave(p, file = paste0(outdir.fig, '-panel_plot_', Code, '.png'), w = 9, h = 9)





cat("\n End postprocessing_figures.R \n")


