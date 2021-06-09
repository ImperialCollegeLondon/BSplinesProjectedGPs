library(rstan)
library(data.table)
library(dplyr)
library(tidyverse)
library(viridis)

indir = "~/git/CDC-covid19-agespecific-mortality-data/inst/" # path to the repo
outdir = file.path(indir, "results", 'new')
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
vaccinedata = as.data.table(read.csv(file))
vaccinedata = subset(vaccinedata, Demographic.Group %in% c("Ages_65-74_yrs", "Ages_75+_yrs"), 
                     select = c('Date', 'Demographic.Group', 'People.with.at.least.one.dose', 'People.who.are.fully.vaccinated', 'Census'))
setnames(vaccinedata, 1:4, c('date', 'age', 'count_vaccinated_1dosep', 'count_vaccinated_fully'))
vaccinedata[, age := gsub('Ages_(.+)_yrs', '\\1', age)]
tmp = vaccinedata[, list(count_vaccinated_1dosep= sum(count_vaccinated_1dosep), 
                         count_vaccinated_fully = sum(count_vaccinated_fully), 
                         Census = sum(Census)), by = c('date')]
tmp[, age := '65+']
vaccinedata = rbind(vaccinedata, tmp)
vaccinedata[, prop_vaccinated_1dosep := count_vaccinated_1dosep / Census]
vaccinedata[, prop_vaccinated_fully := count_vaccinated_fully / Census]

# create loc division map
loc_div = data.table(code = c(c("CT", "ME", "MA", "NH", "RI", "VT"), c("NJ", "NY", "PA", "NYC"), c("IL", "IN", "MI", "OH", "WI"), c("IA", "KS", "MN", "MO", "NE", "ND", "SD"),
                              c("DE", "DC", "FL", "GA", "MD", "NC", "SC", "VA", "WV"), c("AL", "KY", "MS", "TN", "AR", "LA", "OK", "TX"),c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY"),
                              c("AK", "CA", "HI", "OR", "WA")),
                     division = c(rep("New\nEngland", 6), rep("Middle\nAtlantic", 4), rep("East North\nCentral", 5), rep("West North\nCentral", 7), 
                                  rep("South\nAtlantic", 9), rep("South\nCentral", 8), rep("Mountain", 8), rep("Pacific", 5)))

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
# plot contribution over time
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir, '-MortalityRateTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = merge(mortality_rate, region_name, by = 'code')
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date))
plot_mortality_rate_all_states(mortality_rate, outdir)

max_date = format(unique(mortality_rate$date), '%B %d, %Y')
dold2p = mortality_rate[age == '80+' & M > 0.02, loc_label]
dold2p_n = paste0(paste0(dold2p[-length(dold2p)], collapse = ', '), ' and ', dold2p[length(dold2p)])
dold3p = mortality_rate[age == '80+' & M > 0.03, loc_label]
dold3p_n = paste0(paste0(dold3p[-length(dold3p)], collapse = ', '), ' and ', dold3p[length(dold3p)])

mortality_stats = list(max_date = max_date, nstates2p = length(dold2p), nstates3p = length(dold3p),
     dold3p_n = dold3p_n, dold2p_n = dold2p_n)
saveRDS(mortality_stats, file = paste0(outdir, '-mortality_stats.rds'))


#
# trend in contribution
contribution074 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution074[[i]] = readRDS(paste0(outdir, '-Contribution_Age_0-74_', locs[i], '.rds'))
}
contribution074 = do.call('rbind', contribution074)
contribution074 = merge(contribution074, region_name, by = 'code')

contribution75 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution75[[i]] = readRDS(paste0(outdir, '-Contribution_Age_75+_', locs[i], '.rds'))
}
contribution75 = do.call('rbind', contribution75)
contribution75 = merge(contribution75, region_name, by = 'code')

contribution = rbind(contribution75, contribution074)

death2agegroups = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  death2agegroups[[i]] = readRDS(paste0(outdir, '-DeathByAgeTable_2agegroups_', locs[i], '.rds'))
}
death2agegroups = do.call('rbind', death2agegroups)
death2agegroups = merge(death2agegroups, region_name, by = 'code')

# find states slow, fast, plateau
beta = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  y = subset(contribution75, code == locs[i] & date >= as.Date('2020-12-01') & date < as.Date('2021-05-01'))$M
  x = 1:length(y)
  fit1 = lm(y ~ x)
  
  y = subset(contribution75, code == locs[i] & date >= as.Date('2021-05-01'))$M
  x = 1:length(y)
  fit2 = lm(y ~ x)
  
  beta[[i]] = data.table(code = locs[i], loc_label = region_name[code == locs[i], loc_label],
                         betatot = fit1$coefficients[2],  betalast = fit2$coefficients[2])
}
beta = do.call('rbind', beta)
beta =  subset(beta, !code %in% c('HI', 'VT', 'AK'))

plateau = beta[betalast > 0, list(code = code, loc_label = loc_label)]
slowd = beta[betatot > summary(beta$betatot)[5],  list(code = code, loc_label = loc_label)]
fastd = beta[betatot < summary(beta$betatot)[2],  list(code = code, loc_label = loc_label)]

contribution_stats = statistics_contribution_all_states(contribution75)
contribution_stats[['plateau']] = paste0(paste0(plateau$loc_label[-nrow(plateau)], collapse = ', '), ' and ', plateau$loc_label[nrow(plateau)])
contribution_stats[['slowd']] = paste0(paste0(slowd$loc_label[-nrow(slowd)], collapse = ', '), ' and ', slowd$loc_label[nrow(slowd)])
contribution_stats[['fastd']] = paste0(paste0(fastd$loc_label[-nrow(fastd)], collapse = ', '), ' and ', fastd$loc_label[nrow(fastd)])
contribution_stats[['date']] = format(as.Date('2021-05-01'),  '%B %d, %Y')

saveRDS(contribution_stats, file = paste0(outdir, '-contribution_rel_adj_stats.rds'))


selected_states = c('CA', 'ID', 'AR')
stopifnot(selected_states[1] %in% plateau$code) 
stopifnot(selected_states[2] %in% slowd$code) 
stopifnot(selected_states[3] %in% fastd$code) 

contribution75 =  subset(contribution75, !code %in% c('HI', 'VT', 'AK'))
death2agegroups =  subset(death2agegroups, !code %in% c('HI', 'VT', 'AK'))


plot_contribution_magnitude_all_states(contribution, death2agegroups, selected_states, outdir)
  


#
# plot contribution baseline
contribution_ref = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution_ref[[i]] = readRDS(paste0(outdir, '-contribution_refTable_', locs[i], '.rds'))
}
contribution_ref = do.call('rbind', contribution_ref)
contribution_ref = merge(contribution_ref, region_name, by = 'code')
contribution_ref = merge(contribution_ref, loc_div, by = 'code')

contribution_ref_adj = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution_ref_adj[[i]] = readRDS(paste0(outdir, '-contribution_ref_adjTable_', locs[i], '.rds'))
}
contribution_ref_adj = do.call('rbind', contribution_ref_adj)
contribution_ref_adj = merge(contribution_ref_adj, region_name, by = 'code')
contribution_ref_adj = merge(contribution_ref_adj, loc_div, by = 'code')

contribution_ref = subset(contribution_ref, !code %in% c('AK', 'HI', 'VT'))
contribution_ref_adj = subset(contribution_ref_adj, !code %in% c('AK', 'HI', 'VT'))

plot_contribution_ref_all_states(contribution_ref, contribution_ref_adj, outdir)

# stats
contribution_baseline = statistics_contributionref_all_states(contribution_ref_adj)
saveRDS(contribution_baseline, file = paste0(outdir, '-contribution_ref_adj_stats.rds'))




# 
# plot_contribution_all_states(contribution, vaccinedata, outdir)
# 
# # stats

# 

# 






