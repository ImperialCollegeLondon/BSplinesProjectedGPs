library(rstan)
library(data.table)
library(dplyr)
library(tidyverse)
library(viridis)

indir = "/rds/general/user/mm3218/home/git/CDC-covid19-agespecific-mortality-data/inst/" # path to the repo
outdir = file.path(indir, "results")
stan_model = "210505b1"
JOBID = 5791

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-stan_model')
  stopifnot(args_line[[7]]=='-JOBID')
  indir <- args_line[[2]]
  outdir <- args_line[[4]]
  stan_model <- args_line[[6]]
  JOBID <- as.numeric(args_line[[8]])
}

# load functions
source(file.path(indir, "functions", "postprocessing-plotting_functions.R"))
source(file.path(indir, "functions", "postprocessing-summary_functions.R"))
source(file.path(indir, "functions", "postprocessing-utils.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir.fit = file.path(outdir, run_tag, "fits")
outdir.data = file.path(outdir, run_tag, "data")
outdir.table.post = file.path(outdir, run_tag, "table", run_tag)
outdir.fig.post = file.path(outdir, run_tag, "figure", run_tag)

# find locs
files = list.files(path = outdir.fit)
locs = unique(gsub(paste0("fit_cumulative_deaths_(.+)_.*"), "\\1", files))
locs = locs[!grepl('location', locs)]

# load image 
load(file.path(outdir.data, paste0("stanin_", locs[1], "_",run_tag,".RData")))
outdir.fig = outdir.fig.post
outdir.tab = outdir.table.post

# plot mean age of deatg over time
mean_age_death = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mean_age_death[[i]] = readRDS(paste0(outdir.tab, '-MeanAgeOfDeath_', locs[i], '.rds'))
}
mean_age_death = do.call('rbind', mean_age_death)
mean_age_death = merge(mean_age_death, unique(select(deathByAge, code, loc_label)), by = 'code')
plot_mean_age_death(copy(mean_age_death), outdir.fig)

# plot contribution over time
# files = list.files(path = dirname(outdir.tab))
# files = files[grepl('ProbabilityRatioTable', files)]
# locs = unique(gsub(paste0(run_tag, '-ProbabilityRatioTable_(.+).rds'), "\\1", files))

ratio_contribution = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  ratio_contribution[[i]] = readRDS(paste0(outdir.tab, '-ProbabilityRatioTable_', locs[i], '.rds'))
}
ratio_contribution = do.call('rbind', ratio_contribution)
ratio_contribution = merge(ratio_contribution, unique(select(deathByAge, code, loc_label)), by = 'code')
plot_ratio_contribution(copy(ratio_contribution), outdir.fig)





