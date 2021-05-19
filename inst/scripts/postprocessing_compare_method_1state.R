library(rstan)
library(data.table)
library(dplyr)
library(gridExtra)

indir = "~/git/CDC-covid19-agespecific-mortality-data/inst" # path to the repo
outdir = '~/git/CDC-covid19-agespecific-mortality-data/inst/results/testing_knots'
location.index = 1
stan_model = c("210429b2", '210513b', '210505c1', '210513b')
JOBID = c(4116, 13254, 12392, 6187)
model_name = c('GP-BS-CAR', 'GP-SE', 'GP-BS-SE', 'GP-IN')

# load functions
source(file.path(indir, "functions", "postprocessing-plotting_functions.R"))
source(file.path(indir, "functions", "plotting_functions.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir.table = file.path(outdir, run_tag,  run_tag)

# load data
data = readRDS( paste0(outdir.table[1], '-predictive_checks_table_AL.rds') )
data = select(data, date, age, loc_label, code, weekly.deaths)
tmp1 = data[, list(total_deaths = sum(na.omit(weekly.deaths))), by = 'date']
data = merge(data, tmp1, by = 'date')
data[, prop_deaths := weekly.deaths / total_deaths]
data$method = 'observation'

# load contribution table discrete
tab_cd = list()
for(i in seq_along(JOBID)){
  tab = readRDS( paste0(outdir.table[i], '-phi_reducedTable_AL.rds') )
  tab$method = model_name[i]
  tab_cd[[i]] = tab
}
tab_cd = do.call('rbind', tab_cd)
p_cd = plot_contribution_comparison_method(tab_cd,  data, model_name)


# death by age figure 
tab_d = list()
for(i in seq_along(JOBID)){
  tab = readRDS( paste0(outdir.table[i], '-deaths_predict_state_age_strataTable_AL.rds') )
  tab$method = model_name[i]
  tab_d[[i]] = tab
}
tab_d = do.call('rbind', tab_d)
p_d = plot_death_comparison_method(tab_d,  data, model_name)

p = grid.arrange(p_cd, p_d, nrow = 2)
# ggsave()


# load contribution table continuous
tab_cc = list()
for(i in seq_along(JOBID)){
  tab = readRDS( paste0(outdir.table[i], '-phiTable_AL.rds') )
  tab$method = model_name[i]
  tab_cc[[i]] = tab
}
tab_cc = do.call('rbind', tab_cc)
p = plot_contribution_continuous_comparison_method(tab_cc, 'GP-BS-SE')
# ggsave()


# surface 

