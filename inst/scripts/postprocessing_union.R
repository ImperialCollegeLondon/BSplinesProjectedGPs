library(rstan)
library(data.table)
library(dplyr)
library(tidyverse)
library(viridis)

indir = "~/git/CDC-covid19-agespecific-mortality-data/inst/" # path to the repo
outdir = file.path(indir, "results")
stan_model = "210505c1"
JOBID = 8843

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

# plot mean age of deatg over time
# mean_age_death = vector(mode = 'list', length = length(locs))
# for(i in seq_along(locs)){
#   mean_age_death[[i]] = readRDS(paste0(outdir, '-MeanAgeOfDeath_', locs[i], '.rds'))
# }
# mean_age_death = do.call('rbind', mean_age_death)
# mean_age_death = merge(mean_age_death, region_name, by = 'code')
# plot_mean_age_death(copy(mean_age_death), outdir)
# 
# # plot contribution over time
# ratio_contribution = vector(mode = 'list', length = length(locs))
# for(i in seq_along(locs)){
#   ratio_contribution[[i]] = readRDS(paste0(outdir, '-ProbabilityRatioTable_', locs[i], '.rds'))
# }
# ratio_contribution = do.call('rbind', ratio_contribution)
# ratio_contribution = merge(ratio_contribution, region_name, by = 'code')
# plot_ratio_contribution(copy(ratio_contribution), outdir)

# plot contribution over time
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir, '-MortalityRateTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = merge(mortality_rate, region_name, by = 'code')
plot_mortality_rate_all_states(mortality_rate, outdir)


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


plot_contribution_ref_all_states(contribution_ref, contribution_ref_adj, outdir)
  



