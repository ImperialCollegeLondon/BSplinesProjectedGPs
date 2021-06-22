library(rstan)
library(data.table)
library(dplyr)
library(tidyverse)
library(viridis)

indir = "~/git/covid19Vaccination/inst/" # path to the repo
outdir = file.path(indir, "results")
stan_model = "210529b"
JOBID = 2117

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

#load vaccination data
file  = file.path(indir, 'data', 'demographic_trends_of_people_receiving_covid19_vaccinations_in_the_united_states_210520.csv')
vaccinedata_age = clean_vaccination_data_age(file)

file = file.path(indir, 'data', 'us_state_vaccinations_210611.csv')
vaccinedata_state = clean_vaccination_data_state(file)

# create loc division map
loc_div = data.table(code = c(c("CT", "ME", "MA", "NH", "RI", "VT"), c("NJ", "NY", "PA", "NYC"), c("IL", "IN", "MI", "OH", "WI"), c("IA", "KS", "MN", "MO", "NE", "ND", "SD"),
                              c("DE", "DC", "FL", "GA", "MD", "NC", "SC", "VA", "WV"), c("AL", "KY", "MS", "TN", "AR", "LA", "OK", "TX"),c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY"),
                              c("AK", "CA", "HI", "OR", "WA")),
                     division = c(rep("New\nEngland", 6), rep("Middle\nAtlantic", 4), rep("East North\nCentral", 5), rep("West North\nCentral", 7), 
                                  rep("South\nAtlantic", 9), rep("South\nCentral", 8), rep("Mountain", 8), rep("Pacific", 5)))

# states with too few deaths
rm_states = c('HI', 'VT', 'AK', 'WY')

#
# find locs
files = list.files(path = dirname(outdir))
files = files[grepl('ProbabilityRatioTable', files)]
locs = unique(gsub(paste0(run_tag, '-ProbabilityRatioTable_(.+).rds'), "\\1", files))

# find region name
region_name = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  region_name[[i]] = readRDS(paste0(outdir, '-predictive_checks_table_', locs[i], '.rds'))
}
region_name = do.call('rbind', region_name)
region_name = unique(select(region_name, code, loc_label))



#
# plot mortality rate over time
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir, '-MortalityRateTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = merge(mortality_rate, region_name, by = 'code')
mortality_rate = subset(mortality_rate, !code %in% rm_states)
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)
plot_mortality_rate_all_states(mortality_rate, outdir)

max_date = format(unique(mortality_rate$date), '%B %d, %Y')
dold2p = mortality_rate[age == '85+' & M > 0.025, loc_label]
dold2p_n = paste0(paste0(dold2p[-length(dold2p)], collapse = ', '), ' and ', dold2p[length(dold2p)])
state_max = mortality_rate[age == '85+',]
state_max = state_max[order(M, decreasing = T)]
state_max = state_max[,loc_label][1:5]
state_max = paste0(paste0(state_max[-length(state_max)], collapse = ', '), ' and ', state_max[length(state_max)])

d5574 = mortality_rate[age == '55-74',paste0(round((median(M)*100),2),'\\%')]
d7584 = mortality_rate[age == '75-84',paste0(round((median(M)*100),0),'\\%')]

mortalityCAWA = subset(mortality_rate, code %in% c('CA', 'WA') & age == '85+')[, list(loc_label = loc_label,
                                                                                      M = paste0(format(round((M*100),2),nsmall = 2),'\\%'),
                                                                                      CL = paste0(format(round((CL*100),2),nsmall = 2),'\\%'),
                                                                                      CU = paste0(format(round((CU*100),2),nsmall = 2),'\\%')) ] 

mortality_stats = list(max_date = max_date, nstates2p = length(dold2p), dold2p_n,
                       state_max, d5574, d7584, mortalityCAWA)
saveRDS(mortality_stats, file = paste0(outdir, '-mortality_stats.rds'))


#
# trend in contribution
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
plot_contribution_all_states(contribution, vaccinedata_state, outdir)
find_regime_state(contribution75, vaccinedata_state, rm_states,  outdir)

#
# Plot magnitude of mortality
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
  propdeath[[i]] = readRDS(paste0(outdir, '-DeathByAgeprop_Table_phi_', locs[i], '.rds'))
  
}
death3 = do.call('rbind', death3)
death3 = merge(death3, region_name, by = 'code')
death3 =  subset(death3, !code %in% rm_states)
plot_mortality_all_states(death3, data, vaccinedata_state, outdir)

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

tmp = find_statistics_weekly_deaths(death2, propdeath, deathpost2, vaccinedata_state, outdir)


#
# plot contribution baseline
contribution_ref = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution_ref[[i]] = readRDS(paste0(outdir, '-contribution_refTable_', locs[i], '.rds'))
}
contribution_ref = do.call('rbind', contribution_ref)
contribution_ref = merge(contribution_ref, region_name, by = 'code')
contribution_ref = merge(contribution_ref, loc_div, by = 'code')
contribution_ref = subset(contribution_ref, !code %in% rm_states)

contribution_ref_adj = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution_ref_adj[[i]] = readRDS(paste0(outdir, '-contribution_ref_adjTable_', locs[i], '.rds'))
}
contribution_ref_adj = do.call('rbind', contribution_ref_adj)
contribution_ref_adj = merge(contribution_ref_adj, region_name, by = 'code')
contribution_ref_adj = merge(contribution_ref_adj, loc_div, by = 'code')
contribution_ref_adj = subset(contribution_ref_adj, !code %in% rm_states)

plot_contribution_ref_all_states(contribution_ref, contribution_ref_adj, outdir)

# stats
contribution_baseline = statistics_contributionref_all_states(contribution_ref_adj)
saveRDS(contribution_baseline, file = paste0(outdir, '-contribution_ref_adj_stats.rds'))








