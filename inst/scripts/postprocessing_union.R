
cat("\n Begin postprocessing_union.R \n")

library(rstan)
library(data.table)
library(dplyr)
library(tidyverse)
library(viridis)
library(gridExtra)
library(ggpubr)
library(jcolors)
library(facetscales)
library(usmap)

indir ="~/git/BSplinesProjectedGPs/inst" # path to the repo
outdir = file.path('/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst', "results")
# states = strsplit('FL',',')[[1]]
# states = strsplit('CA,FL,NY,TX,PA,IL,OH,GA,NC,MI',',')[[1]]
states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME",
           "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
           "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")

states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", 
            "MI", "MN", "MO", "MS", "MT", "NC", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
            "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
paste0(states,collapse =  ',')
stan_model = "220209a"
JOBID = 8004

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

# 4 states
selected_codes = c('CA', 'FL', 'NY', 'TX')
selected_10_codes = c('CA','FL','NY','TX','PA','IL','OH','GA','NC','MI')

#
# Mortality rate

# mortality rate over time continuous
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir.table, '-MortalityRateContinuousTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)
plot_mortality_rate_continuous_all_states(mortality_rate, outdir.fig)

# mortality rate over time discrete
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir.table, '-MortalityRateTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbindlist', list(l = mortality_rate, fill=TRUE))
mortality_rateJ21 <- subset(mortality_rate, date == '2021-06-05')
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)
# mortality_rate[is.na(M), M := c(2.056321e-03, 1.684421e-02)]
crude_mortality_rate = find_crude_mortality_rate(mortality_rate, df_age_continuous, df_age_reporting, pop_data)
plot_mortality_rate_all_states(mortality_rate, crude_mortality_rate, outdir.fig)
plot_mortality_rate_all_states2(mortality_rate, outdir.fig)
plot_mortality_rate_all_states_map(mortality_rate, outdir.fig)
plot_relative_mortality_all_states(mortality_rateJ21, nyt_data, outdir.fig)

# mortality rate correlation with longtermdeaths in care home facilities
mortality_rateJ21 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rateJ21[[i]] = readRDS(paste0(outdir.table, '-PosteriorSampleMortalityRateTable_', locs[i], '.rds'))
  mortality_rateJ21[[i]] <- mortality_rateJ21[[i]][date == '2021-06-05']
}
mortality_rateJ21 = do.call('rbindlist', list(l = mortality_rateJ21, fill=TRUE))
save_mortality_rate_correlation_longtermdeaths(mortality_rateJ21, nyt_data, region_name, outdir.table)


# aggregate across states
mortality_rate_posterior_samples = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate_posterior_samples[[i]] = readRDS(paste0(outdir.table, '-PosteriorSampleMortalityRateTable_', locs[i], '.rds'))
}
mortality_rate_posterior_samples = do.call('rbind', mortality_rate_posterior_samples)
mortality_rate_across_states <- find_mortality_rate_across_states(mortality_rate_posterior_samples)


# statistics
find_statistics_mortality_rate(mortality_rate, mortality_rate_across_states, outdir.table)


#
# Plot weekly deaths over time
death3 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  death3[[i]] = readRDS(paste0(outdir.table, '-DeathByAgeTable_3agegroups_', locs[i], '.rds'))
}
death3 = do.call('rbind', death3)
if(!'prop_vac_start_counterfactual' %in% names(stan_data)) resurgence_dates = NULL
plot_mortality_all_states(subset(death3, code %in% selected_codes), resurgence_dates,'selectedStates', outdir.fig)
if(any(!locs %in% selected_codes))
  plot_mortality_all_states(subset(death3, !code %in% selected_codes), resurgence_dates,'otherStates', outdir.fig)

# if(length(locs) > 6){
#   mid_locs = floor(length(locs) / 2)
# 
#   plot_mortality_all_states(subset(death3, code %in% locs[1:mid_locs]), resurgence_dates, 'allStates_part1', outdir.fig)
#   plot_mortality_all_states(subset(death3, code %in% locs[(mid_locs+1):(mid_locs*2)]), resurgence_dates, 'allStates_part2', outdir.fig)
# 
# } else{
#   plot_mortality_all_states(death3, resurgence_dates, 'allStates', outdir.fig)
# }


#
# predictions
predictions = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  predictions[[i]] = readRDS(paste0(outdir.table, '-predictive_checks_table_', locs[i], '.rds'))
}
predictions = do.call('rbind', predictions)
predictions <- select(predictions, - min.sum.weekly.deaths, - max.sum.weekly.deaths, - sum.weekly.deaths, - weekly.deaths, - inside.CI, 
                      -state_index, - week_index, - age_index)
predictions <- predictions[order(loc_label, date, age)]

dir = file.path(gsub('(.+)\\/results.*', '\\1', outdir.table), 'results', 'predictions')
dir.create(dir)
saveRDS(predictions, file = file.path(dir, 'predicted_weekly_deaths.rds'))


# 
# vaccination effect on contribution predict
contribution = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution[[i]] = readRDS(paste0(outdir.table, '-phi_predict_reduced_vacTable_', locs[i], '.rds'))
}
contribution = do.call('rbind', contribution)

mid_code = round(length(locs) / 2)
contribution[, M_median := median(M), by = c('date', 'age')]
plot_contribution_vaccine(contribution, vaccine_data, 'predict_all', outdir.fig)
plot_contribution_vaccine(subset(contribution, code %in% selected_codes), vaccine_data,  'predict_selected_codes',outdir.fig)

## contribution 
contribution = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  contribution[[i]] = readRDS(paste0(outdir.table, '-phi_reduced_vacTable_', locs[i], '.rds'))
}
contribution = do.call('rbind', contribution)
contribution[, M_median := median(M), by = c('date', 'age')]
mid_code = round(length(locs) / 2)
plot_contribution_vaccine(contribution, vaccine_data, 'all', outdir.fig)
plot_contribution_vaccine(subset(contribution, code %in% selected_codes), vaccine_data,  'selected_codes',outdir.fig)
plot_contribution_vaccine(subset(contribution, code %in% 'FL'), vaccine_data,  'FL',outdir.fig)

## contribution shift
locs_plus_US <- c(locs, 'US')
contributiondiff = vector(mode = 'list', length = length(locs_plus_US))
for(i in seq_along(locs_plus_US)){
  contributiondiff[[i]] = readRDS(paste0(outdir.table, '-phi_reduced_vacDiffTable_', locs_plus_US[i], '.rds'))
}
contributiondiff = do.call('rbind', contributiondiff)
save_statistics_contributiondiff(contributiondiff, outdir.table)

# Resurgence
r_pdeaths  = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  r_pdeaths[[i]] = readRDS(paste0(outdir.table, '-r_pdeathsTable_', locs[i], '.rds'))
}
r_pdeaths = do.call('rbind', r_pdeaths)

p4 <- plot_relative_resurgence_vaccine2(r_pdeaths, F, outdir.fig, '_selected_states', selected_codes)
p_all <- plot_relative_resurgence_vaccine_end_3(r_pdeaths, T, outdir.fig, '_all_states')
p <- ggarrange(p4, p_all,  ncol = 1, heights = c(0.4,0.6))
ggsave(p, file =paste0(outdir.fig, '-relative_deaths_vaccine_coverage_panel.png'), w = 8, h = 9)


cat("\n End postprocessing_union.R \n")



