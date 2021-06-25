
cat("\n Begin postprocessing_assess_mixing.R \n")

library(rstan)
library(data.table)
library(dplyr)

indir = "~/git/covid19Vaccination/inst" # path to the repo
outdir = '/rds/general/user/mm3218/home/git/covid19Vaccination/inst/results'
location.index = 26
stan_model = "210505b1"
JOBID = 2967

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



cat("\n End postprocessing_assess_mixing.R \n")



