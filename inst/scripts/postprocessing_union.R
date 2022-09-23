
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
library(ggrepel)

indir ="~/git/BSplinesProjectedGPs/inst" # path to the repo
outdir = file.path('/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst', "results")
states <- 'CA'
stan_model = "220209a"
JOBID = 3541

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
selected_15_codes <- c('CA','FL','NY','TX','PA','IL','OH','GA','NC','MI','NJ','VA','WA','AZ','MA')
selected_20_codes <- c('CA','FL','NY','TX','PA','IL','OH','GA','NC','MI','NJ','VA','WA','AZ','MA','TN','IN','MD','MO','WI')
selected_30_codes <- c('CA','FL','NY','TX','PA','IL','OH','GA','NC','MI','NJ','VA','WA','AZ','MA','TN','IN','MD','MO','WI','CO','MN','SC','AL','LA','KY','OR','UT','IA','NV')


#
# Comparison to DoH data
if(length(locs) == 50){
  tab_doh = list(); k = 1
  for(i in seq_along(locs)){
    file <- paste0(outdir.table, '-CumDeathsComp_ScrapedData_', locs[i], '.rds')
    if(file.exists(file)){
      tab_doh[[k]]  = readRDS(file)[[2]]
      k = k + 1
    }
  }
  tab_doh = as.data.table(do.call('rbind', tab_doh))
  tab_doh[, method := 'BSGP']; model_name = 'BSGP'
  save_doh_comparison(tab_doh, region_name, outdir.table)
}


#
# Mortality rate

# mortality rate over time continuous
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir.table, '-MortalityRateContinuousTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)
if(all(selected_codes %in% mortality_rate[, unique(code)])){
  plot_mortality_rate_continuous_all_states(mortality_rate, selected_codes, outdir.fig)
}

# mortality rate over time discrete
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir.table, '-MortalityRateTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbindlist', list(l = mortality_rate, fill=TRUE))
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)
crude_mortality_rate = find_crude_mortality_rate(mortality_rate, df_age_continuous, df_age_reporting, pop_data)
plot_mortality_rate_all_states(mortality_rate, crude_mortality_rate, outdir.fig)
plot_mortality_rate_all_states2(mortality_rate, outdir.fig)

# aggregate across states
mortality_rate_posterior_samples = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate_posterior_samples[[i]] = readRDS(paste0(outdir.table, '-PosteriorSampleMortalityRateTable_', locs[i], '.rds'))
}
mortality_rate_posterior_samples = as.data.table(do.call('rbind', mortality_rate_posterior_samples))
mortality_rate_across_states <- find_mortality_rate_across_states(mortality_rate_posterior_samples)

# statistics
find_statistics_mortality_rate(mortality_rate, mortality_rate_across_states, outdir.table)

# mortality rate over time discrete 4 gae groups
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir.table, '-MortalityRate4agegroupsTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbindlist', list(l = mortality_rate, fill=TRUE))
mortality_rateJ21 <- subset(mortality_rate, date == '2021-06-05')
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)
p <- plot_mortality_rate_all_states_map(mortality_rate, mortality_rateJ21, outdir.fig)
p2 <- plot_relative_mortality_all_states(mortality_rateJ21, nyt_data, outdir.fig)
p <- ggarrange(p, labels = 'A')
p2 <- ggarrange(p2, labels = 'B')
pp <- grid.arrange(p, p2, layout_matrix = rbind(c(1, 2), c(1, NA)), widths = c(0.45, 0.55), heights = c(0.6, 0.4))
ggsave(pp, file = paste0(outdir.fig, paste0('-MortalityRateRelative_CareHome_panel.png')), w = 9.5, h = 5.5)
ggsave(pp, file = paste0(outdir.fig, paste0('-MortalityRateRelative_CareHome_panel.pdf')), w = 9.5, h = 5.5)

mortality_rate_posterior_samples = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate_posterior_sample = readRDS(paste0(outdir.table, '-PosteriorSampleMortalityRate4agegroupsTable_', locs[i], '.rds'))
  mortality_rate_posterior_samples[[i]] <- mortality_rate_posterior_sample[date == '2021-06-05']
}
mortality_rate_posterior_samples = do.call('rbindlist', list(l = mortality_rate_posterior_samples, fill=TRUE))
save_mortality_rate_correlation_longtermdeaths(mortality_rate_posterior_samples, mortality_rateJ21, nyt_data, region_name, outdir.table)


#
# Plot weekly deaths over time
death3 = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  death3[[i]] = readRDS(paste0(outdir.table, '-DeathByAgeTable_3agegroups_', locs[i], '.rds'))
}
death3 = do.call('rbind', death3)
if(!'prop_vac_start_counterfactual' %in% names(stan_data)) resurgence_dates <- find_resurgence_dates(JHUData, deathByAge, locations$code)
if(all(selected_codes %in% death3[, unique(code)])){
  plot_mortality_all_states(subset(death3, code %in% selected_codes), resurgence_dates, 'selectedStates', outdir.fig)
}else{
  plot_mortality_all_states(death3, resurgence_dates, 'selectedStates', outdir.fig)
}
if(any(!locs %in% selected_codes & locs %in% selected_10_codes))
  plot_mortality_all_states(subset(death3, !code %in% selected_codes & code %in% selected_10_codes), resurgence_dates, 'otherStates', outdir.fig)
if(length(locs) > 10){
  locs_not_selected = locs[!locs %in% selected_10_codes]
  mid_locs = floor(length(locs_not_selected) / 2)
  plot_mortality_all_states(subset(death3, code %in% locs_not_selected[1:mid_locs]), resurgence_dates, 'allStates_part1', outdir.fig)
  plot_mortality_all_states(subset(death3, code %in% locs_not_selected[(mid_locs+1):(mid_locs*2)]), resurgence_dates, 'allStates_part2', outdir.fig)

}


#
# predictions
predictions = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  predictions[[i]] = readRDS(paste0(outdir.table, '-DeathByAge', 'Table_', 'prediction_agg', '_', locs[i], '.rds'))
}

predictions = do.call('rbind', predictions)
predictions <- select(predictions, - week_index, - state_index, - mean, - emp_JHU, -age_index, -emp, 
                      -emp_adj)
predictions <- predictions[order(loc_label, date, age)]

dir = file.path(gsub('(.+)\\/results.*', '\\1', outdir.table), 'results', 'predictions')
dir.create(dir)
saveRDS(predictions, file = file.path(dir, 'predicted_weekly_deaths.rds'))


# 
# Predicted contibution 
file <- paste0(outdir.table, '-phi_predict_reduced_vacTable_', locs[1], '.rds')
if(file.exists(file)){
  contribution = vector(mode = 'list', length = length(locs))
  for(i in seq_along(locs)){
    contribution[[i]] = readRDS(paste0(outdir.table, '-phi_predict_reduced_vacTable_', locs[i], '.rds'))
  }
  contribution = do.call('rbind', contribution)
  
  mid_code = round(length(locs) / 2)
  contribution[, M_median := median(M), by = c('date', 'age')]
  plot_contribution_vaccine(contribution, vaccine_data_pop, 'predict_all', outdir.fig)
  if(all(selected_codes %in% contribution[, unique(code)])){
    plot_contribution_vaccine(subset(contribution, code %in% selected_codes), vaccine_data_pop,  'predict_selected_codes',outdir.fig)
  }
  
  ## shift over time
  locs_plus_US <- c(locs, 'US')
  contributiondiff = vector(mode = 'list', length = length(locs_plus_US))
  for(i in seq_along(locs_plus_US)){
    contributiondiff[[i]] = readRDS(paste0(outdir.table, '-phi_predict_reduced_vacDiffTable_', locs_plus_US[i], '.rds'))
  }
  contributiondiff = do.call('rbind', contributiondiff)
  plot_contributiondiff_map(contributiondiff, 'diff1', outdir.fig, 'predict')
  plot_contributiondiff_map(contributiondiff, 'diff2', outdir.fig, 'predict')
  save_statistics_contributiondiff(contributiondiff, outdir.table, 'predict')
}

# 
# Estimated contibution 
file <- paste0(outdir.table, '-phi_reduced_vacTable_', locs[1], '.rds')
if(file.exists(file)){
  contribution = vector(mode = 'list', length = length(locs))
  for(i in seq_along(locs)){
    contribution[[i]] = readRDS(paste0(outdir.table, '-phi_reduced_vacTable_', locs[i], '.rds'))
  }
  contribution = do.call('rbind', contribution)
  contribution[, M_median := median(M), by = c('date', 'age')]
  mid_code = round(length(locs) / 2)
  plot_contribution_vaccine(contribution, vaccine_data_pop, 'all', outdir.fig)
  if(all(selected_codes %in% contribution[, unique(code)])){
    plot_contribution_vaccine(subset(contribution, code %in% selected_codes), vaccine_data_pop,  'selected_codes',outdir.fig)
  }else{
    plot_contribution_vaccine(contribution, vaccine_data_pop,  'selected_codes',outdir.fig)
  }
  if(all(c('TX', 'NH') %in% locs)) plot_contribution_vaccine(subset(contribution, code %in% c('TX', 'NH')), vaccine_data_pop,  'TXNH',outdir.fig)
  
  ## shift over time
  locs_plus_US <- c(locs, 'US')
  contributiondiff = vector(mode = 'list', length = length(locs_plus_US))
  for(i in seq_along(locs_plus_US)){
    contributiondiff[[i]] = readRDS(paste0(outdir.table, '-phi_reduced_vacDiffTable_', locs_plus_US[i], '.rds'))
  }
  contributiondiff = do.call('rbind', contributiondiff)
  plot_contributiondiff_map(contributiondiff, 'diff1', outdir.fig)
  plot_contributiondiff_map(contributiondiff, 'diff2', outdir.fig)
  save_statistics_contributiondiff(contributiondiff, outdir.table)
}


#
# Resurgence during the summer 2021
file <- paste0(outdir.table, '-r_pdeathsTable_', locs[1], '.rds')
if(file.exists(file)){
  r_pdeaths  = vector(mode = 'list', length = length(locs))
  for(i in seq_along(locs)){
    r_pdeaths[[i]] = readRDS(paste0(outdir.table, '-r_pdeathsTable_', locs[i], '.rds'))
  }
  r_pdeaths = do.call('rbind', r_pdeaths)
  
  if(all(selected_codes %in% r_pdeaths[, unique(code)])){
    p4 <- plot_relative_resurgence_vaccine2(r_pdeaths, F, outdir.fig, '_selected_states', selected_codes)
  }else{
    p4 <- plot_relative_resurgence_vaccine2(r_pdeaths, F, outdir.fig, '_selected_states', r_pdeaths[, unique(code)])
  }
  
  if(all(selected_10_codes %in% locs)){
    plot_relative_resurgence_vaccine2_long(r_pdeaths, outdir.fig, '_10_states', selected_10_codes)
    plot_relative_resurgence_vaccine2_long(r_pdeaths, outdir.fig, '_other_states', locs[!locs %in% selected_10_codes])
    
    p_all <- plot_relative_resurgence_vaccine_end_3(r_pdeaths, T, outdir.fig, '_all_states')
    plot_relative_resurgence_vaccine_panel(p4, p_all, '50_log', outdir.fig)
    p_all2 <- plot_relative_resurgence_vaccine_end_3(r_pdeaths, T, outdir.fig, '_10_states', selected_10_codes)
    plot_relative_resurgence_vaccine_panel(p4, p_all2, '10_log', outdir.fig)
    p_all2 <- plot_relative_resurgence_vaccine_end_3(r_pdeaths, T, outdir.fig, '_15_states', selected_15_codes)
    plot_relative_resurgence_vaccine_panel(p4, p_all2, '15_log', outdir.fig)
    p_all2 <- plot_relative_resurgence_vaccine_end_3(r_pdeaths, T, outdir.fig, '_20_states', selected_20_codes)
    plot_relative_resurgence_vaccine_panel(p4, p_all2, '20_log', outdir.fig)
    p_all2 <- plot_relative_resurgence_vaccine_end_3(r_pdeaths, T, outdir.fig, '_30_states', selected_30_codes)
    plot_relative_resurgence_vaccine_panel(p4, p_all2, '30_log', outdir.fig)
    p_all2 <- plot_relative_resurgence_vaccine_end_3(r_pdeaths, F, outdir.fig, '_30_states', selected_30_codes)
    plot_relative_resurgence_vaccine_panel(p4, p_all2, '30', outdir.fig)
    
  }
}


cat("\n End postprocessing_union.R \n")



