
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
states = strsplit('CA,FL,NY,TX',',')[[1]]
# states = strsplit('CA,FL,NY,TX,PA,IL,OH,GA,NC,MI',',')[[1]]
# states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME",
#            "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
#            "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
stan_model = "220209a"
JOBID = 1811791

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
selected_codes_10 <- c('CA','FL','NY','TX','PA','IL','OH','GA','NC','MI')
  
#
# mortality rate over time discrete
mortality_rate = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  mortality_rate[[i]] = readRDS(paste0(outdir.table, '-MortalityRateTable_', locs[i], '.rds'))
}
mortality_rate = do.call('rbind', mortality_rate)
mortality_rate = subset(mortality_rate, date == max(mortality_rate$date)-7)

plot_mortality_rate_all_states(mortality_rate, outdir.fig)

mortality_rate_10 <- mortality_rate[code %in% selected_codes_10]
p <- plot_mortality_rate_all_states(mortality_rate_10, outdir.fig)
find_statistics_mortality_rate(mortality_rate_10, outdir.table)
limits_mortality_rate <- range(c(mortality_rate_10$CL, mortality_rate_10$CU))


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
if(any(locs %in% selected_codes)){
  plot_mortality_all_states(subset(death3, code %in% selected_codes),'selectedStates', outdir.fig)
}
if(any(locs %in% selected_codes_10 & !locs %in% selected_codes)){
  plot_mortality_all_states(subset(death3, code %in% selected_codes_10 & !code %in% selected_codes),'otherStates', outdir.fig)
}


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



cat("\n End postprocessing_union.R \n")



