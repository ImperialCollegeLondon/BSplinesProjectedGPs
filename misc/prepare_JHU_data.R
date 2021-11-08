library(data.table)
library(dplyr)

indir ="~/git/BSplinesProjectedGPs/inst" # path to the repo

path.to.CDC.data = file.path(indir, "data", paste0("CDC-data_2021-09-29.rds"))
deathByAge = readRDS(path.to.CDC.data) # cdc data 

tmp = read.csv(url('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv'))
tmp = as.data.table(reshape2::melt(tmp, id.vars = c('UID', "iso2", 'iso3', 'code3', 'FIPS', 'Admin2', 'Province_State', 'Country_Region', 'Lat', 'Long_', 'Combined_Key', 'Population')))
tmp[, date := as.Date(variable, format = 'X%m.%d.%y')]
tmp = select(tmp, Province_State, date, value)
setnames(tmp, c('Province_State', 'value'), c('loc_label', 'cumulative_deaths'))
tmp = tmp[loc_label %in% deathByAge$loc_label]
tmp = tmp[, list(cumulative_deaths = sum(cumulative_deaths)), by = c('date', 'loc_label')]
tmp[, daily_deaths := c(diff(cumulative_deaths), NA), by ='loc_label' ]
tmp[daily_deaths < 0, daily_deaths := 0]
tmp = tmp[!is.na(daily_deaths)]


tmp = merge(tmp, unique(select(deathByAge, code, loc_label)), by = 'loc_label')
tmp = select(tmp, -loc_label)

tmp = subset(tmp, date <= max(deathByAge$date) + 7)

path.to.JHU.data = file.path(indir, "data", paste0("jhu_data_2021-09-30.rds"))
saveRDS(tmp, path.to.JHU.data) 

