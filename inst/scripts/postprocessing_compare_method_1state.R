library(rstan)
library(data.table)
library(dplyr)
library(gridExtra)
library(loo)
library(ggpubr)

indir = "~/git/BSplinesProjectedGPs/inst" # path to the repo
outdir = '~/git/BSplinesProjectedGPs/inst/results'
stan_model = c('211031c2', '211031d2', '211031e2', '211031b2')
JOBID = c(20426, 20145, 19986, 4715)
model_name = c('Standard 2D GP', 'Standard B-splines\nsurface', 'P-splines surface', 'Regularised B-splines\nprojected 2D GP')

# load functions
source(file.path(indir, "functions", "postprocessing-plotting_functions.R"))
source(file.path(indir, "functions", "plotting_functions.R"))
source(file.path(indir, "functions", "postprocessing-utils.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir = file.path(outdir, run_tag,  run_tag)

#
# get region name
files = list.files(path = dirname(outdir[1]))
files = files[grepl('CumDeathsComp_ScrapedData', files)]
locs = unique(gsub(paste0(run_tag[1], '-CumDeathsComp_ScrapedData_(.+).rds'), "\\1", files))

region_name = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  region_name[[i]] = readRDS(paste0(outdir[1], '-phiTable_', locs[i], '.rds'))
}
region_name = do.call('rbind', region_name)
region_name = unique(select(region_name, code, loc_label))


#
# load predictions
tab_doh = list(); k = 1
for(i in seq_along(JOBID)){
  for(j in seq_along(locs)){
    tab = readRDS( paste0(outdir[i], '-CumDeathsComp_ScrapedData_', locs[j], '.rds') )
    tab[[2]]$method = model_name[i]
    tab_doh[[k]] = tab[[2]]
    k = k + 1
  }
}
tab_doh = do.call('rbind', tab_doh)

tab = tab_doh[, list(prop_correct_predictions = mean(prop.death.inside.CI)), by = c('code', 'method')]
tab = merge(tab, region_name, by = 'code')
tab[, propT_n := paste0(format(round(prop_correct_predictions * 100, 2), nsmall = 2), '\\%')]
tab[, method := factor(method, levels = model_name)]

tab2 = tab_doh[, list(avgpropT = mean(prop.death.inside.CI)),  by = c('method')]
tab2 = tab2[, avg_propT_n := paste0(format(round(avgpropT * 100, 2), nsmall = 2), '\\%')]

tab = reshape2::dcast(tab, loc_label~ method, value.var = 'propT_n')
tab = rbind(tab,  data.table(loc_label = 'Average', reshape2::dcast(tab2, .~ method, value.var = 'avg_propT_n')[,-1]) )

stat = paste0(format(round(range(tab2$avgpropT[-4])* 100, 2), nsmall = 2), '\\%')

saveRDS(list(tab, stat), paste0(outdir[4], '-compDoH.rds'))


#
# Compare contribution 

# load predicted contribution 
tab_cc = list()
for(i in seq_along(JOBID)){
  tab = readRDS( paste0(outdir[i], '-phiTable_FL.rds') )
  tab$method = model_name[i]
  tab_cc[[i]] = tab
}
tab_cc = do.call('rbind', tab_cc)

# load predicted weekly deaths 
tab_d = readRDS( paste0(outdir[4], '-DeathByAgeTable_alpha_FL.rds') )
tab_d$method = 'Regularised B-splines\nprojected 2D GP'

# load CDC data
data = readRDS( paste0(outdir[1], '-predictive_checks_table_FL.rds') )
data = select(data, date, age, loc_label, code, weekly.deaths)
tmp1 = data[, list(total_deaths = sum(na.omit(weekly.deaths))), by = 'date']
data = merge(data, tmp1, by = 'date')
data[, prop_deaths := weekly.deaths / total_deaths]
data$method = 'observation'

p = plot_contribution_continuous_comparison_method_with_data(tab_cc, tab_d, data, 'Regularised B-splines\nprojected 2D GP', model_name)
ggsave(p, file= paste0(outdir[4], '-phi_short_compmethod_with_data_FL.png'), w = 10, h = 12)
p = plot_contribution_continuous_comparison_method(tab_cc, 'Regularised B-splines\nprojected 2D GP', model_name)
ggsave(p, file= paste0(outdir[4], '-phi_short_compmethod_FL.png'), w = 8, h = 6)


#
# LOO
LOO = list()
for(i in seq_along(JOBID)){
  LOO[[i]] = readRDS( paste0(outdir[i], "-LOO.rds") )
}

loo_compare(LOO)

tmp = loo_compare(LOO[[1]], LOO[[4]])
tmp = gsub(" ", "", format( round(tmp, digits = 2), nsmall = 2) )
saveRDS(tmp, file = paste0(outdir[4], "-LOO_comp.rds"))

#
# time of execution
time = list()
for(i in seq_along(JOBID)){
  time[[i]] = readRDS( paste0(outdir[i], "-time_elapsed.rds") )
}
time = do.call('rbind', time)
time = data.table(time = round(time[,1] / time[4,1], 2), method = model_name)
saveRDS(time, file = paste0(outdir[4], "-time_execution_com.rds"))

time[,1] / 60 / 60 / 8

