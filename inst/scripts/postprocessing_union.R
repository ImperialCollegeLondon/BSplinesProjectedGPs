library(rstan)
library(data.table)
library(dplyr)
library(tidyverse)
library(viridis)
library(grid)
library(ggpubr)

indir = "~/git/covid19Vaccination/inst/" # path to the repo
outdir = file.path(indir, "results")
stan_model = "210529b"
JOBID = 29051

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
source(file.path(indir, "functions", "postprocessing-statistics_functions.R"))
source(file.path(indir, "functions", "summary_functions.R"))
source(file.path(indir, "functions", "postprocessing-utils.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir = file.path(outdir, run_tag, run_tag)

# states with too few deaths
rm_states = c('HI', 'VT', 'AK', 'WY')
selected_states = c('CA', 'FL', 'TX', 'WA', 'NY')

# find locs
files = list.files(path = dirname(outdir))
files = files[grepl('predictive_checks_table_', files)]
locs = unique(gsub(paste0(run_tag, '-predictive_checks_table_(.+).rds'), "\\1", files))

# find region name
region_name = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  region_name[[i]] = readRDS(paste0(outdir, '-predictive_checks_table_', locs[i], '.rds'))
}
region_name = do.call('rbind', region_name)
region_name = unique(select(region_name, code, loc_label))



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
# vaccination effect on resurgences
weekly_deaths = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  weekly_deaths[[i]] = readRDS(paste0(outdir, '-PosteriorsamplesDeathByAgeTable_vacagegroups_', locs[i], '.rds'))
}
weekly_deaths = do.call('rbind', weekly_deaths)
weekly_deaths = merge(weekly_deaths, region_name, by = 'code')
weekly_deaths = subset(weekly_deaths, !code %in% rm_states)

vaccine_data = readRDS(paste0(outdir, '-vaccine_data.rds'))

start_resurgence <- as.Date('2021-07-03')
pick_resurgence <- as.Date('2021-08-28')

unique(weekly_deaths$date)
weekly_deaths1 = copy(weekly_deaths)
weekly_deaths = subset(weekly_deaths1, code %in% selected_states)
data_res = find_vaccine_effects(weekly_deaths, vaccine_data, start_resurgence, pick_resurgence)

variable= 'y_counterfactual'
data_res1 = subset(data_res[[1]], var == variable)
data_res2 = subset(data_res[[2]], date >= as.Date('2021-01-01'))
data_res4 <- copy(data_res[[4]])
data_res5 = copy(data_res[[5]])

plot_vaccine_effects_counterfactual(data_res1, data_res2, outdir)

plot_relative_resurgence_vaccine(data_res2, data_res4, outdir)

plot_estimate_vaccine(data_res5, outdir)

find_stats_vaccine_effects(start_resurgence, pick_resurgence, data_res1, data_res4, outdir)
  



#
# Deaths

# Plot weekly deaths over time
death3 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  death3[[i]] = readRDS(paste0(outdir, '-DeathByAgeTable_3agegroups_', locs[i], '.rds'))
}
death3 = do.call('rbind', death3)
death3 = merge(death3, region_name, by = 'code')
death3 = subset(death3, code %in% selected_states)

#subset(death3, date >= as.Date('2021-01-01'))
plot_mortality_all_states(death3, outdir)



# 
# vaccination effect on contribution
contribution = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution[[i]] = readRDS(paste0(outdir, '-posterior_table_phi_', locs[i], '.rds'))
}
contribution = do.call('rbind', contribution)
contribution = subset(contribution, code %in% selected_states)
plot_contribution_vaccine(contribution, vaccine_data, outdir)
find_regime_state(contribution, vaccine_data, start_resurgence,  outdir)

# contribution75 = copy(contribution)


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








###################################
##################################
################################


#
# contribution over time
contribution054 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution054[[i]] = readRDS(paste0(outdir, '-Contribution_Age_0-54_', locs[i], '.rds'))
}
contribution054 = do.call('rbind', contribution054)
contribution054 = merge(contribution054, region_name, by = 'code')
contribution054 = subset(contribution054, !code %in% rm_states)

contribution5574 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution5574[[i]] = readRDS(paste0(outdir, '-Contribution_Age_55-74_', locs[i], '.rds'))
}
contribution5574 = do.call('rbind', contribution5574)
contribution5574 = merge(contribution5574, region_name, by = 'code')
contribution5574 = subset(contribution5574, !code %in% rm_states)

contribution75 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution75[[i]] = readRDS(paste0(outdir, '-Contribution_Age_75+_', locs[i], '.rds'))
}
contribution75 = do.call('rbind', contribution75)
contribution75 = merge(contribution75, region_name, by = 'code')
contribution75 = subset(contribution75, !code %in% rm_states)

contribution = rbind(contribution75, rbind(contribution5574, contribution054))
plot_contribution_all_states(contribution, vaccinedata, outdir)
find_regime_state(contribution75, vaccinedata_state, rm_states,  outdir)



#
# Plot weekly deaths over time
death3 = vector(mode = 'list', length = length(locs))
death2 = vector(mode = 'list', length = length(locs))
deathpost2 = vector(mode = 'list', length = length(locs))
deathpost3 = vector(mode = 'list', length = length(locs))
propdeath = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  death3[[i]] = readRDS(paste0(outdir, '-DeathByAgeTable_3agegroups_', locs[i], '.rds'))
  death2[[i]] = readRDS(paste0(outdir, '-DeathByAgeTable_2agegroups_', locs[i], '.rds'))
  deathpost2[[i]] = readRDS(paste0(outdir, '-PosteriorsamplesDeathByAgeprop_Table_2agegroups_', locs[i], '.rds'))
  deathpost3[[i]] = readRDS(paste0(outdir, '-PosteriorsamplesDeathByAgeprop_Table_3agegroups_', locs[i], '.rds'))
  propdeath[[i]] = readRDS(paste0(outdir, '-DeathByAgeprop_Table_2agegroups_', locs[i], '.rds'))
  
}
death3 = do.call('rbind', death3)
death3 = merge(death3, region_name, by = 'code')
death3 =  subset(death3, !code %in% rm_states)
plot_mortality_all_states(death3, data, vaccinedata_state[,min(date)], outdir)

deathpost2 = do.call('rbind', deathpost2)
deathpost2 = summary_death_all_states(deathpost2, rm_states)

deathpost3 = do.call('rbind', deathpost3)
deathpost3 = summary_death_all_states(deathpost3, rm_states)

death2 = do.call('rbind', death2)
death2 = merge(death2, region_name, by = 'code')
death2 =  subset(death2, !code %in% rm_states)

propdeath = do.call('rbind', propdeath)
propdeath = merge(propdeath, region_name, by = 'code')
propdeath =  subset(propdeath, !code %in% rm_states)

tmp = find_statistics_weekly_deaths(death2, propdeath, deathpost2, deathpost3, vaccinedata_state, outdir)




