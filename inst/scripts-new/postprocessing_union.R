library(rstan)
library(data.table)
library(dplyr)
library(tidyverse)
library(viridis)
library(grid)
library(ggpubr)
library(jcolors)

indir ="~/git/covid19Vaccination/inst" # path to the repo
outdir = file.path('~/Downloads/', "results")
states = strsplit('CA,TX',',')[[1]]
stan_model = "211014b"
JOBID = 3541

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
outdir.fig.post = file.path(outdir, run_tag, "figure", run_tag)
outdir.table = file.path(outdir, run_tag, "table", run_tag)

# find locations
locations = readRDS( file.path(outdir.fit.post, paste0("location_", run_tag,".rds")) )
loc_name = locations[code %in% states,]$loc_label
locs = locations[code %in% states, ]$code
region_name = data.table(loc_label = loc_name, code = Code, state_index = 1:length(Code))

# load image 
load(file.path(outdir.data, paste0("stanin_",run_tag,".RData")))
outdir.fig = outdir.fig.post


#
# mortality rate over time
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir, '-MortalityRateTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = merge(mortality_rate, region_name, by = 'code')
mortality_rate = subset(mortality_rate, code %in% selected_states)
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)

plot_mortality_rate_all_states(mortality_rate, outdir)
find_statistics_mortality_rate(mortality_rate, outdir)


#
# Plot weekly deaths over time
death3 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  death3[[i]] = readRDS(paste0(outdir, '-DeathByAgeTable_3agegroups_', locs[i], '.rds'))
}
death3 = do.call('rbind', death3)
death3 = merge(death3, region_name, by = 'code')
death3 = subset(death3, code %in% selected_states)
plot_mortality_all_states(death3, start_resurgence, outdir)


#
# proportion weekly deaths after vaccine
propdeath3 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  propdeath3[[i]] = readRDS(paste0(outdir, '-DeathByAgeprop_Table_3agegroups_', locs[i], '.rds'))
}
propdeath3 = do.call('rbind', propdeath3)
propdeath3 = merge(propdeath3, region_name, by = 'code')
propdeath3 = subset(propdeath3, code %in% selected_states)
find_prop_deaths_vaccine_statistics(propdeath3, start_vaccine, start_resurgence, outdir)


# 
# vaccination effect on contribution
contribution = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution[[i]] = readRDS(paste0(outdir, '-posterior_table_phi_', locs[i], '.rds'))
}
contribution = do.call('rbind', contribution)
contribution = subset(contribution, code %in% selected_states)
plot_contribution_vaccine(contribution, vaccine_data, outdir)
find_regime_state(contribution, vaccine_data, start_resurgence, start_vaccine, outdir)



#
#  contribution baseline
contribution_ref = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution_ref[[i]] = readRDS(paste0(outdir, '-contribution_refTable_', locs[i], '.rds'))
}
contribution_ref = do.call('rbind', contribution_ref)
contribution_ref = merge(contribution_ref, region_name, by = 'code')
contribution_ref = subset(contribution_ref, code %in% selected_states)

contribution_ref_adj = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution_ref_adj[[i]] = readRDS(paste0(outdir, '-contribution_ref_adjTable_', locs[i], '.rds'))
}
contribution_ref_adj = do.call('rbind', contribution_ref_adj)
contribution_ref_adj = merge(contribution_ref_adj, region_name, by = 'code')
contribution_ref_adj = subset(contribution_ref_adj, code %in% selected_states)

plot_contribution_ref_all_states(contribution_ref, contribution_ref_adj, outdir)

contribution_baseline = statistics_contributionref_all_states(contribution_ref_adj, outdir)






