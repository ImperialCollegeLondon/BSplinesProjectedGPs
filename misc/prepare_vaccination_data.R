library(data.table)
library(dplyr)

indir ="~/git/BSplinesProjectedGPs/inst" # path to the repo
indir2 ="~/git/BSplinesProjectedGPs/misc" # path to the repo

# data for location name
path.to.CDC.data = file.path(indir, "data", paste0("CDC-data_2021-09-29.rds"))
locname_data = unique(select(readRDS(path.to.CDC.data), loc_label, code)) # cdc data 

# population data
path.to.popdata = file.path(indir, "data", paste0("us_population_withnyc.rds"))
pop_data = readRDS(path.to.popdata)
pop_data = as.data.table(select(reshape2::melt(pop_data, id.vars = c('Region', 'code')), -Region))

# vaccination data
path.to.data = file.path(indir2, "data-vaccination", paste0("COVID-19_Vaccinations_in_the_United_States_Jurisdiction_210930.csv"))
data = as.data.table(read.csv(path.to.data)) 

data[, date := as.Date(Date, '%m/%d/%Y')]
setnames(data, 'Location', 'code')
data = subset(data, code %in% unique(locname_data$code))
data = select(data, date, code, Series_Complete_12Plus, Series_Complete_12PlusPop_Pct,
              Series_Complete_18Plus, Series_Complete_18PlusPop_Pct,
              Series_Complete_65Plus, Series_Complete_65PlusPop_Pct)

data = as.data.table(reshape2::melt(data, id.vars = c('date', 'code')))
data[, age.group := gsub('Series_Complete_(.+)Plus(.*)', '\\1', variable)]
data[, variable := ifelse(grepl('Pct', variable), 'prop', 'abs')]
data = as.data.table( reshape2::dcast(data, date + code + age.group ~ variable, value.var = 'value') )

# ensure that the series is strictly increasing
data = data[order(code, age.group, date)]
data[, abs := cummax(abs), by = c('code', 'age.group')]

# find unique population size (change over time in the raw data)
tmp = subset(data, date == max(date))
tmp[, pop := abs / (prop / 100)]

# before this date, 12+ was not reported
date.min.12p = data[age.group == '12' & abs == 0, list(max_date = max(date)), by = 'code']

# find absolute by age group (from 12+, 18+, 65+, to 12-17, 18-64, 65+)
data[, abs := abs( c( diff(abs), abs[length(abs)] ) ), by = c('date', 'code')]

# find proportion
data = merge(select(data, -prop), 
             subset(tmp, select = c('code', 'age.group', 'pop')), by = c('code', 'age.group'))
data[, prop := abs / pop]
data = select(data, -abs, -pop)

# replace prop for 12+ by 0 before date.min.12p
data = merge(data, date.min.12p, by = 'code')
data[date <= max_date & age.group == 12, prop := 0]
data = select(data, -max_date)

# new age group 
data[age.group == '12', age.group := '12-17']
data[age.group == '18', age.group := '18-64']
data[age.group == '65', age.group := '65+']

# take rolling average
ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}
data[, prop := ma(prop, 7), by = c('code', 'age.group')]
data = data[!is.na(prop)]

# add young age group
tmp = unique(select(data, date, code))
tmp[, age.group := '0-11']
tmp[, prop := 0]
data = rbind(data, tmp)

# find pop
pop_data = subset(pop_data, variable != 'Total')
pop_data[, variable := gsub('\\+', '-105', variable)]
tmp = unique(select(pop_data, variable))
tmp = tmp[, age.min := gsub('(.+)\\-(.*)', '\\1', variable), by = 'variable']
tmp = tmp[, age.max := gsub(paste0(age.min, '\\-(.+)'), '\\1', variable), by = 'age.min']
tmp = tmp[, list(age = age.min:age.max), by = 'variable']
pop_data = merge(pop_data, tmp, by = 'variable', allow.cartesian=TRUE)
pop_data[, pop := value / length(value), by = c('variable', 'code')]
pop_data = select(pop_data, code, age, pop)

# enlarge by 1-y age band
tmp = unique(select(data, age.group))
tmp[, age.group := gsub('\\+', '-105', age.group)]
tmp = tmp[, age.min := gsub('(.+)\\-(.*)', '\\1', age.group), by = 'age.group']
tmp = tmp[, age.max := gsub(paste0(age.min, '\\-(.+)'), '\\1', age.group), by = 'age.min']
tmp = tmp[, list(age = age.min:age.max), by = 'age.group']
data[, age.group := gsub('\\+', '-105', age.group)]
data = merge(tmp, data, by = 'age.group', allow.cartesian=TRUE)
data = select(data, code, date, age, prop)

# merge with pop
data = merge(data, pop_data, by = c('code', 'age'))

# find name of location 
data = merge(data, locname_data, by = 'code')

# sanity checks
stopifnot(all(data$prop <= 1 & data$prop >= 0))
stopifnot(all(!is.na(data)))
nrow(data) == length(unique(data$date)) * length(unique(data$code)) * length(unique(data$age))

# keep first date of non zera
start_date = data[prop > 0, min(date) ]
data = data[date >= start_date]

# save
saveRDS(data, file.path(indir, "data", paste0("vaccination-prop-2021-10-14.rds")))

        