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

# laod CDC death by age data 
deathByAge = readRDS( file.path(indir, "data", paste0("CDC-data_2021-06-02.rds")))
age_max = 105
create_map_age(age_max)
fouragegroups = c('0-24', '25-54', '55-74', '75-84', '85+')
df_age = data.table(age = fouragegroups)
df_age[, age_index := 1:nrow(df_age)]
df_age[, age_from := gsub('(.+)-.*', '\\1', age)]
df_age[, age_to := gsub('.*-(.+)', '\\1', age)]
df_age[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
df_age[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
df_age[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
df_age[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
set(df_age, NULL, 'age_from', df_age[,as.numeric(age_from)])
set(df_age, NULL, 'age_to', df_age[,as.numeric(age_to)])
set(df_age, NULL, 'age_to_index', df_age[,as.numeric(age_to_index)])
set(df_age, NULL, 'age_from_index', df_age[,as.numeric(age_from_index)])
df_age_reporting[, age_state_index := which(df_age$age_from <= age_from & df_age$age_to >= age_to), by = 'age_index']

date_10thcum = subset(deathByAge, !is.na(weekly.deaths))
date_10thcum = date_10thcum[, list(weekly.deaths = sum(weekly.deaths)), by = c('date', 'code')]
date_10thcum[, cum.deaths := cumsum(weekly.deaths), by = 'code']
date_10thcum = date_10thcum[ cum.deaths >=10 , list(min_date = min(date)), by = 'code']

# ref week
df_week = unique(select(deathByAge, date, code))
df_week[, month := as.numeric(format(date, '%m'))]
df_week[grepl('2021', date), month := month + 12]
df_week = merge(df_week, date_10thcum, by = 'code')
tmp1 = df_week[date >= min_date, list(min_month = min(month)), by = 'code']
df_week = merge(df_week, tmp1, by = 'code')
df_week1= df_week[month >= min_month & month <= min_month + 2, , by = 'code']

deathByAge = merge(deathByAge, df_age_reporting, by = 'age_from')
deathByAge = merge(deathByAge, df_week1, by = c('date', 'code'))
deathByAge = deathByAge[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c('code', 'age_state_index')]
tmp1 = deathByAge[, list(total_deaths = sum(na.omit(weekly.deaths))), by = c( 'code')]
deathByAge = merge(deathByAge, tmp1, by = c('code'))
deathByAge[, prop_deaths := weekly.deaths / total_deaths]
deathByAge[, age := fouragegroups[age_state_index]]


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
contribution064 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution064[[i]] = readRDS(paste0(outdir, '-Contribution_Age_0-64_', locs[i], '.rds'))
}
contribution064 = do.call('rbind', contribution064)
contribution064 = merge(contribution064, region_name, by = 'code')

contribution6574 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution6574[[i]] = readRDS(paste0(outdir, '-Contribution_Age_65-74_', locs[i], '.rds'))
}
contribution6574 = do.call('rbind', contribution6574)
contribution6574 = merge(contribution6574, region_name, by = 'code')

contribution75 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution75[[i]] = readRDS(paste0(outdir, '-Contribution_Age_75+_', locs[i], '.rds'))
}
contribution75 = do.call('rbind', contribution75)
contribution75 = merge(contribution75, region_name, by = 'code')

contribution = rbind(contribution75, rbind(contribution064, contribution6574))

plot_contribution_all_states(contribution, vaccinedata, outdir)

# stats
contribution_stats = statistics_contribution_all_states(contribution064, contribution6574, contribution75, vaccinedata)
saveRDS(contribution_stats, file = paste0(outdir, '-contribution_rel_adj_stats.rds'))


# plot contribution over time
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

deathByAge = merge(deathByAge, region_name, by = 'code')
deathByAge = merge(deathByAge, loc_div, by = 'code')

plot_contribution_ref_all_states(contribution_ref, contribution_ref_adj, deathByAge, outdir)
  
# stats
contribution_baseline = statistics_contributionref_all_states(contribution_ref_adj)
saveRDS(contribution_baseline, file = paste0(outdir, '-contribution_ref_adj_stats.rds'))










