library(rstan)
library(data.table)
library(dplyr)
library(gridExtra)
library(loo)
library(ggpubr)

indir = "~/git/CDC-covid19-agespecific-mortality-data/inst" # path to the repo
outdir = '~/git/CDC-covid19-agespecific-mortality-data/inst/results'
location.index = 1
stan_model = c("210429h1", '210513b', '210429b2', '210505c1')
JOBID = c(17173, 4197, 3708, 3934)
model_name = c('GP-SE', 'GP-BS-I', 'GP-BS-CAR', 'GP-BS-SE')

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

p = ggarrange(p_d, p_cd, nrow = 2, labels = c('A', 'B'), label.x = 0.1)
ggsave(p, file= paste0(outdir.table[1], '-deaths_predict_state_age_strata_compmethod_AL.png'), w = 10, h = 5)


# load contribution table continuous
tab_cc = list()
for(i in seq_along(JOBID)){
  tab = readRDS( paste0(outdir.table[i], '-phiTable_AL.rds') )
  tab$method = model_name[i]
  tab_cc[[i]] = tab
}
tab_cc = do.call('rbind', tab_cc)
p = plot_contribution_continuous_comparison_method(tab_cc, 'GP-BS-SE', model_name)
ggsave(p, file= paste0(outdir.table[1], '-phi_short_compmethod_AL.png'), w = 6, h = 4)



# LOO
LOO = list()
for(i in seq_along(JOBID)){
  LOO[[i]] = readRDS( paste0(outdir.table[i], "-LOO_AL.rds") )
}

loo_compare(LOO[[2]], LOO[[3]], LOO[[4]])

loo_compare(LOO[[1]], LOO[[4]])

