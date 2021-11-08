
cat("\n Begin postprocessing_union.R \n")


library(rstan)
library(data.table)
library(dplyr)
library(tidyverse)
library(viridis)
library(gridExtra)
library(ggpubr)
library(jcolors)

indir ="~/git/BSplinesProjectedGPs/inst" # path to the repo
outdir = file.path('/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst', "results")
states = strsplit('CA,FL',',')[[1]]
stan_model = "211014b"
JOBID = 7259

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

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir.data = file.path(outdir, run_tag, "data")
outdir.fit.post = file.path(outdir, run_tag, "fits")
outdir.fig.post = file.path(outdir, run_tag, "figure", run_tag)
outdir.table = file.path(outdir, run_tag, "table", run_tag)

# find locations
locations = readRDS( file.path(outdir.fit.post, paste0("location_", run_tag,".rds")) )
loc_name = locations[code %in% states,]$loc_label
locs = locations[code %in% states, ]$code
region_name = data.table(loc_label = loc_name, code = locs, state_index = 1:length(locs))

# load image 
load(file.path(outdir.data, paste0("stanin_",run_tag,".RData")))
outdir.fig = outdir.fig.post

# date 
start_vaccine = vaccine_data[prop > 0 & date %in% df_week$date, min(date)]

#
# mortality rate over time
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir.table, '-MortalityRateTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)

plot_mortality_rate_all_states(mortality_rate, outdir.fig)
find_statistics_mortality_rate(mortality_rate, outdir.table)


#
# Plot weekly deaths over time
death3 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  death3[[i]] = readRDS(paste0(outdir.table, '-DeathByAgeTable_3agegroups_', locs[i], '.rds'))
}
death3 = do.call('rbind', death3)
plot_mortality_all_states(death3, resurgence_dates, outdir.fig)


#
# proportion weekly deaths after vaccine
propdeath3 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  propdeath3[[i]] = readRDS(paste0(outdir.table, '-DeathByAgeprop_Table_3agegroups_', locs[i], '.rds'))
}
propdeath3 = do.call('rbind', propdeath3)
find_prop_deaths_vaccine_statistics(propdeath3, start_vaccine, min(resurgence_dates$start_resurgence), outdir.table)


# 
# vaccination effect on contribution
if(!is.null(stan_data$prop_vac)){
  contribution = vector(mode = 'list', length = length(locs))
  for(i in seq_along(locs)){
    contribution[[i]] = readRDS(paste0(outdir.table, '-phi_reduced_vacTable_', locs[i], '.rds'))
  }
  contribution = do.call('rbind', contribution)
  plot_contribution_vaccine(contribution, vaccine_data, resurgence_dates, outdir.fig)
  find_regime_state(contribution, vaccine_data, resurgence_dates, start_vaccine, outdir.table)
}


#
#  contribution baseline
contribution_ref = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution_ref[[i]] = readRDS(paste0(outdir.table, '-contribution_refTable_', locs[i], '.rds'))
}
contribution_ref = do.call('rbind', contribution_ref)

contribution_ref_adj = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution_ref_adj[[i]] = readRDS(paste0(outdir.table, '-contribution_ref_adjTable_', locs[i], '.rds'))
}
contribution_ref_adj = do.call('rbind', contribution_ref_adj)

plot_contribution_ref_all_states(contribution_ref, contribution_ref_adj, outdir.fig)

contribution_baseline = statistics_contributionref_all_states(contribution_ref_adj, outdir.table)


cat("\n End postprocessing_union.R \n")



