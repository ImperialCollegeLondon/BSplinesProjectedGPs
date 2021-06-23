library(rstan)
library(data.table)
library(dplyr)
library(gridExtra)
library(loo)
library(ggpubr)

indir = "~/git/covid19Vaccination/inst" # path to the repo
outdir = '~/git/covid19Vaccination/inst/results'
location.index = 1
stan_model = c('210429h1', '210529c', '210529b')
JOBID = c(1954, 499, 2117)
model_name = c('Standard 2D GP', 'Standard B-splines\nsurface', 'Regularised B-splines\nprojected 2D GP')

# load functions
source(file.path(indir, "functions", "postprocessing-plotting_functions.R"))
source(file.path(indir, "functions", "plotting_functions.R"))
source(file.path(indir, "functions", "postprocessing-utils.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir.table = file.path(outdir, run_tag,  run_tag)

# load data
data = readRDS( paste0(outdir.table[1], '-predictive_checks_table_FL.rds') )
data = select(data, date, age, loc_label, code, weekly.deaths)
tmp1 = data[, list(total_deaths = sum(na.omit(weekly.deaths))), by = 'date']
data = merge(data, tmp1, by = 'date')
data[, prop_deaths := weekly.deaths / total_deaths]
data$method = 'observation'

# states with too few deaths
rm_states = c('HI', 'VT', 'AK', 'WY')


files = list.files(path = dirname(outdir.table[3]))
files = files[grepl('CumDeathsComp_ScrapedData', files)]
locs = unique(gsub(paste0(run_tag[3], '-CumDeathsComp_ScrapedData_(.+).rds'), "\\1", files))
region_name = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  region_name[[i]] = readRDS(paste0(outdir.table[3], '-predictive_checks_table_', locs[i], '.rds'))
}
region_name = do.call('rbind', region_name)
region_name = unique(select(region_name, code, loc_label))

tab_doh = list(); k = 1
for(i in seq_along(JOBID)){
  for(j in seq_along(locs)){
    tab = readRDS( paste0(outdir.table[i], '-CumDeathsComp_ScrapedData_', locs[j], '.rds') )
    tab[[2]]$method = model_name[i]
    tab_doh[[k]] = tab[[2]]
    k = k + 1
  }
}
tab_doh = do.call('rbind', tab_doh)
tab_doh = subset(tab_doh, !code %in% rm_states)
tab_doh[weekly.deaths_total == 0, prop.weekly.deaths := 0]
tab_doh[, prop.death.inside.CI := (prop.weekly.deaths <= CU_prop & prop.weekly.deaths >= CL_prop)]

tab = tab_doh[, list(propT = mean(prop.death.inside.CI)), by = c('code', 'method')]
tab = merge(tab, region_name, by = 'code')
tab[, propT_n := paste0(format(round(propT * 100, 2), nsmall = 2), '\\%')]
tab[, method := factor(method, levels = model_name)]
tab2 = tab[, list(avgpropT = mean(propT)),  by = c('method')]
tab2 = tab2[, avg_propT_n := paste0(format(round(avgpropT * 100, 2), nsmall = 2), '\\%')]

tab = reshape2::dcast(tab, loc_label~ method, value.var = 'propT_n')
tab = rbind(tab,  data.table(loc_label = 'Average', reshape2::dcast(tab2, .~ method, value.var = 'avg_propT_n')[,-1]) )
stat = paste0(format(round(range(tab2$avgpropT)* 100, 2), nsmall = 2), '\\%')
saveRDS(list(tab, stat), paste0(outdir.table[1], '-compDoH.rds'))

# load contribution table continuous
tab_cc = list()
for(i in seq_along(JOBID)){
  tab = readRDS( paste0(outdir.table[i], '-phiTable_FL.rds') )
  tab$method = model_name[i]
  tab_cc[[i]] = tab
}
tab_cc = do.call('rbind', tab_cc)
tab_d = readRDS( paste0(outdir.table[3], '-DeathByAgeTable_alpha_FL.rds') )
tab_d$method = 'Regularised B-splines\nprojected 2D GP'

p = plot_contribution_continuous_comparison_method(tab_cc, tab_d, data, 
                                                   'Regularised B-splines\nprojected 2D GP', model_name)
ggsave(p, file= paste0(outdir.table[1], '-phi_short_compmethod_FL.png'), w = 9, h = 10)


# LOO
LOO = list()
for(i in seq_along(JOBID)){
  LOO[[i]] = readRDS( paste0(outdir.table[i], "-LOO_FL.rds") )
}

loo_compare(LOO[[2]], LOO[[3]])

tmp = loo_compare(LOO[[1]], LOO[[3]])
tmp = gsub(" ", "", format( round(tmp, digits = 2), nsmall = 2) )
saveRDS(tmp, file = paste0(outdir.table[1], "-LOO_comp.rds"))

# 
# 
# # load contribution table discrete
# tab_cd = list()
# for(i in seq_along(JOBID)){
#   tab = readRDS( paste0(outdir.table[i], '-phi_reducedTable_AL.rds') )
#   tab$method = model_name[i]
#   tab_cd[[i]] = tab
# }
# tab_cd = do.call('rbind', tab_cd)
# p_cd = plot_contribution_comparison_method(tab_cd,  data, model_name)
# 
# 
# # death by age figure 
# tab_d = list()
# for(i in seq_along(JOBID)){
#   tab = readRDS( paste0(outdir.table[i], '-deaths_predict_state_age_strataTable_AL.rds') )
#   tab$method = model_name[i]
#   tab_d[[i]] = tab
# }
# tab_d = do.call('rbind', tab_d)
# p_d = plot_death_comparison_method(tab_d,  data, model_name)
# 
# p = ggarrange(p_d, p_cd, nrow = 2, labels = c('A', 'B'), label.x = 0.11)
# ggsave(p, file= paste0(outdir.table[1], '-deaths_predict_state_age_strata_compmethod_AL.png'), w = 10, h = 5)


# load comparison to DoH
# tab_doh = list()
# for(i in seq_along(JOBID)){
#   tab = readRDS( paste0(outdir.table[i], '-CumDeathsComp_ScrapedData_FL.rds') )
#   scraped_data = tab[[1]]
#   tab[[2]]$method = model_name[i]
#   tab_doh[[i]] = tab[[2]]
# }
# tab_doh = do.call('rbind', tab_doh)
# 
# p = compare_CDCestimation_DoH_age_plot_compmethod(tab_doh, scraped_data, model_name, 
#                                                   selected_method = 'Regularised B-splines\nprojected GP')
# ggsave(p, file = paste0(outdir.table[1], '-comparison_DoH_CDC_uncertainty_', unique(tab_doh$code), '_commethod.png'), w = 7, h = 9, limitsize = F)
# 
