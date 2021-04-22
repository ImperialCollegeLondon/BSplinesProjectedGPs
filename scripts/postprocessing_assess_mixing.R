
cat("\n Begin postprocessing_assess_mixing.R \n")

library(rstan)
library(data.table)
library(dplyr)

indir = "~/git/CDC-covid19-agespecific-mortality-data" # path to the repo
outdir = '/rds/general/user/mm3218/home/git/CDC-covid19-agespecific-mortality-data/results'
location.index = 3
stan_model = "210419"
JOBID = 18389

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
source(file.path(indir, "functions", "summary_functions.R"))
source(file.path(indir, "functions", "postprocessing-plotting_functions.R"))
source(file.path(indir, "functions", "postprocessing-summary_functions.R"))
source(file.path(indir, "functions", "stan_utility_functions.R"))

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
Code = locations[location.index,]$code

# load image 
load(file.path(outdir.data, paste0("stanin_", Code, "_",run_tag,".RData")))
outdir.fig = outdir.fig.post
outdir.fit = outdir.fit.post

# load fit cumulative deaths
cat("Load fits \n")
file = file.path(outdir.fit.post, paste0("fit_cumulative_deaths_", Code, "_", run_tag,".rds"))
fit_cum <- readRDS(file=file)

# Convergence diagnostics
cat("\nMake convergence diagnostics \n")
make_convergence_diagnostics_stats(fit_cum, outdir.table)
plot_convergence_diagnostics(fit_cum, "Daily deaths fit", 'daily', outfile = outdir.fig)

# Make predictive checks table
cat("\nMake posterior predive checks table \n")
predictive_checks_table_1 = make_predictive_checks_table(fit_cum, df_week1, df_age_reporting_1, 
                                                         data1, 'deaths_predict_state_age_strata_1', outdir.table)
predictive_checks_table_2 = make_predictive_checks_table(fit_cum, df_week2, df_age_reporting_2, 
                                                         data2, 'deaths_predict_state_age_strata_2', outdir.table)


# plot predictive checks table
cat("\nMake posterior predive checks plots \n")
plot_posterior_predictive_checks(predictive_checks_table_1, predictive_checks_table_2, 
                                 variable = "daily.deaths", 
                                 lab = "Daily COVID-19 deaths", 
                                 outdir = outdir.fig)

# plot predictive check cumulative
cat("\nMake posterior predive checks cumulative plots \n")
tmp = find_cumulative_deaths_state_age(fit_cum, df_week, df_age_continuous, unique(df_age_reporting_1$age_cat), 'deaths_predict')
setnames(tmp, 'age_state_index', 'age_index')
tmp = add_intervals_missing_data(tmp, stan_data_1)
tmp2 = find_cumulative_deaths_state_age(fit_cum, df_week, df_age_continuous, unique(df_age_reporting_2$age_cat), 'deaths_predict')
setnames(tmp2, 'age_state_index', 'age_index')
tmp2 = add_intervals_missing_data(tmp2, stan_data_2)
plot_posterior_predictive_checks_cumulative(tmp, tmp2, outdir.fig)

 
cat("\n End postprocessing_assess_mixing.R \n")

