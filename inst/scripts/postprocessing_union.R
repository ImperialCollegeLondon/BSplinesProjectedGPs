
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
states = strsplit('CA,FL,NY,TX,PA,IL,OH,GA,NC,MI',',')[[1]]
stan_model = "220131a"
JOBID = 30502

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

# date 
start_vaccine = vaccine_data[prop > 0 & date %in% df_week$date, min(date)]

if(is.null(stan_data[['week_indices_resurgence']])){
  resurgence_dates = NULL
}
  
#
# mortality rate over time discrete
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir.table, '-MortalityRateTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)

p <- plot_mortality_rate_all_states(mortality_rate, outdir.fig)
find_statistics_mortality_rate(mortality_rate, outdir.table)
limits_mortality_rate <- range(c(mortality_rate$CL, mortality_rate$CU))


#
# mortality rate over time continuous
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir.table, '-MortalityRateContinuousTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)
p_NY <- plot_mortality_rate_continuous_all_states(mortality_rate, limits_mortality_rate, outdir.fig)

# make panel mortaliry
p_NY <- ggarrange(p_NY , labels = 'A', label.y = 1) #, label.x = 0.03
p <- ggarrange(p, labels = 'B')

panel <- grid.arrange(grobs = list(p_NY, p), widths = c(1, 0.7), heights = c(0.8, 0.1, 1.05), 
                      layout_matrix = rbind(c(1, NA),
                                            c(NA, NA),
                                            c(2, 2)),
                      left = paste0('Predicted COVID-19 attributable mortality rates\nas of ', format(unique(mortality_rate$date), '%B %Y')))
ggsave(panel, file = paste0(outdir.fig, '-MortalityRate_panel_plot.png'), w = 7, h = 6)


#
# Plot weekly deaths over time
death3 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  death3[[i]] = readRDS(paste0(outdir.table, '-DeathByAgeTable_3agegroups_', locs[i], '.rds'))
}
death3 = do.call('rbind', death3)
plot_mortality_all_states(subset(death3, code %in% selected_codes), resurgence_dates,'selectedStates', outdir.fig)
if(any(!locs %in% selected_codes))
  plot_mortality_all_states(subset(death3, !code %in% selected_codes), resurgence_dates,'otherStates', outdir.fig)

if(length(locs) > 6){
  mid_locs = floor(length(locs) / 2)
  
  plot_mortality_all_states(subset(death3, code %in% locs[1:mid_locs]), resurgence_dates, 'allStates_part1', outdir.fig)
  plot_mortality_all_states(subset(death3, code %in% locs[(mid_locs+1):(mid_locs*2)]), resurgence_dates, 'allStates_part2', outdir.fig)
  
} else{
  plot_mortality_all_states(death3, resurgence_dates, 'allStates', outdir.fig)
}


#
# proportion weekly deaths after vaccine
if(!is.null(resurgence_dates)){
  propdeath3 = vector(mode = 'list', length = length(locs))
  for(i in seq_along(locs)){
    propdeath3[[i]] = readRDS(paste0(outdir.table, '-DeathByAgeprop_Table_3agegroups_', locs[i], '.rds'))
  }
  propdeath3 = do.call('rbind', propdeath3)
  find_prop_deaths_vaccine_statistics(propdeath3, start_vaccine, min(resurgence_dates$start_resurgence), outdir.table)

# 
# vaccination effect on contribution
  contribution = vector(mode = 'list', length = length(locs))
  for(i in seq_along(locs)){
    contribution[[i]] = readRDS(paste0(outdir.table, '-phi_reduced_vacTable_', locs[i], '.rds'))
  }
  contribution = do.call('rbind', contribution)
  
  mid_code = round(length(locs) / 2)
  plot_contribution_vaccine(subset(contribution, code %in% locs[1:mid_code]), vaccine_data, resurgence_dates, 'part_1', outdir.fig)
  plot_contribution_vaccine(subset(contribution, code %in% locs[(mid_code+1):(mid_code*2)]), vaccine_data, resurgence_dates,  'part_2',outdir.fig)
  
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



