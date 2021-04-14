
cat("\n Begin postprocessing_figures.R \n")


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

# path to CDC and JHU data
path.to.CDC.data = file.path(indir, "data", paste0("CDC-data_2021-04-11.rds"))

# max age considered
age_max = 105

# Load CDC data
deathByAge = readRDS(path.to.CDC.data)

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

df_week = rbind(df_week1, df_week2)
df_week$week_index = 1:nrow(df_week)

# load fit cumulative deaths
cat("Load fits \n")
file = file.path(outdir.fit, paste0("fit_cumulative_deaths_", Code, "_", run_tag,".rds"))
fit_cum <- readRDS(file=file)

# Plots continuous age distribution alpha
cat("\nMake continuous age distribution plots \n")
plot_continuous_age_contribution(fit_cum, df_age_continuous, df_week, "cumulative COVID-19 deaths", 
                                 outdir = outdir.fig)

# Plot estimated CAR covariance matrix
plot_covariance_matrix(fit_cum, outdir = outdir.fig)

# Plot estimate plane with CAR of ICAR
plot_posterior_plane(fit_cum, df_week, outdir = outdir.fig)

# Plot probability ratio of deaths over time
probability_ratio_table = make_probability_ratio_table(fit_cum, df_week, df_age_reporting_2, data1, data2, outdir.table)
plot_probability_ratio(probability_ratio_table, outdir.fig)


cat("\n End postprocessing_figures.R \n")


