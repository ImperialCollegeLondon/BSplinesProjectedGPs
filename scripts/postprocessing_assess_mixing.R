
cat("\n Begin postprocessing_assess_mixing.R \n")

library(rstan)
library(data.table)
library(dplyr)

indir = "~/git/CDC-covid19-agespecific-mortality-data" # path to the repo
outdir = file.path('~/Downloads', "results")
location.index = 2
stan_model = "210319d3"
JOBID = 782737

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
outdir.fit = file.path(outdir, run_tag, "fits")
outdir.fig = file.path(outdir, run_tag, "figure", run_tag)
outdir.table = file.path(outdir, run_tag, "table", run_tag)

if( !dir.exists( dirname(outdir.fig) ) ) dir.create( dirname(outdir.fig), recursive = T)
if( !dir.exists( dirname(outdir.table) ) ) dir.create( dirname(outdir.table), recursive = T)

# path to CDC and JHU data
path.to.CDC.data.1 = file.path(indir, "data", paste0("CDC-data-1_2021-04-11.rds"))
path.to.CDC.data.2 = file.path(indir, "data", paste0("CDC-data-2_2021-04-11.rds"))

# max age considered
age_max = 105

# Load CDC data
deathByAge_1 = readRDS(path.to.CDC.data.1) # cdc data before 2020-09-02 with first age specification 
deathByAge_2 = readRDS(path.to.CDC.data.2) # cdc data after 2020-09-02 with second age specification 

# Create age maps
create_map_age(age_max)

# locations and dates
locations = unique(deathByAge$loc_label) 
loc_name = locations[location.index]
cat("Location ", as.character(loc_name), "\n")

# reference date
ref_date = as.Date('2020-08-29')

# Prepare stan data
cat("\n Prepare stan data \n")
stan_data_1 = prepare_stan_data(deathByAge_1, loc_name, ref_date); data1 = tmp; df_week1 = df_week
stan_data_2 = prepare_stan_data(deathByAge_2, loc_name, ref_date); data2 = tmp; df_week2 = df_week
stan_data = merge_stan_data(stan_data_1, stan_data_2)

# load fit cumulative deaths
cat("Load fits \n")
file = file.path(outdir.fit, paste0("fit_cumulative_deaths_", Code, "_", run_tag,".rds"))
fit_cum <- readRDS(file=file)

# Convergence diagnostics
cat("\nMake convergence diagnostics \n")
make_convergence_diagnostics_stats(fit_cum, outdir.table)
plot_convergence_diagnostics(fit_cum, "Cumulative deaths fit", 'cum', 
                             outfile = paste0(outdir.fig, "-convergence_diagnostics-"))

# Make predictive checks table
cat("\nMake posterior predive checks table \n")
predictive_checks_table_1 = make_predictive_checks_table(fit_cum, df_week1, df_age_reporting_1, 
                                                         data1, 'deaths_predict_state_age_strata_1', outdir.table)
predictive_checks_table_2 = make_predictive_checks_table(fit_cum, df_week2, df_age_reporting_2, 
                                                         data2, 'deaths_predict_state_age_strata_2', outdir.table)

# plot predictive checks table
cat("\nMake posterior predive checks plots \n")
plot_posterior_predictive_checks(predictive_checks_table_1, predictive_checks_table_2, 
                                 variable = "COVID.19.Deaths", 
                                 lab = "Cumulative COVID-19 deaths", 
                                 outdir = outdir.fig)
 

cat("\n End postprocessing_assess_mixing.R \n")

