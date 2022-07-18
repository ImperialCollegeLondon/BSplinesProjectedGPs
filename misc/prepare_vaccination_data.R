library(data.table)
library(dplyr)

indir ="~/git/BSplinesProjectedGPs/inst" # path to the repo
indir2 ="~/git/BSplinesProjectedGPs/misc" # path to the repo

# data for location name
path.to.CDC.data = file.path(indir, "data", paste0("CDC-data_2022-02-06.rds"))
locname_data = unique(select(readRDS(path.to.CDC.data), loc_label, code)) # cdc data 
deathByAge = readRDS(path.to.CDC.data) # cdc data 

# population data
path.to.popdata = file.path(indir, "data", paste0("us_population.csv"))
pop_data = as.data.table( read.csv(path.to.popdata) )
pop_data = select(pop_data, code, age, pop)

# vaccination data
path.to.data = file.path(indir2, "data-vaccination", paste0("COVID-19_Vaccination_Trends_in_the_United_States_National_and_Jurisdictional.csv"))
data = as.data.table(read.csv(path.to.data)) 

data[Date == '12/13/2020' & Location =='AK']
data[, date := as.Date(Date, '%m/%d/%Y')]
data <- data[date_type == 'Admin']
setnames(data, 'Location', 'code')
data = subset(data, code %in% unique(locname_data$code))

if(0){
  library(ggplot2)
  ggplot(data, aes(x= date , y = Series_Complete_Pop_Pct)) + geom_line()
}

# ensure that the series is strictly increasing
data = data[order(code, date)]
data[, abs := cummax(Series_Complete_Pop_Pct), by = c('code')]

data[, prop := Series_Complete_Pop_Pct / 100]
data = select(data, -Series_Complete_Pop_Pct)

# find name of location 
data = merge(data, locname_data, by = 'code')

# sanity checks
stopifnot(all(data$prop <= 1 & data$prop >= 0))
stopifnot(all(!is.na(data)))
nrow(data) == length(unique(data$date)) * length(unique(data$code)) 
data[,range(date)]

# keep first date of non zera
start_date = data[prop > 0, min(date) ]
data = data[date >= start_date]

# keep date before max date of deathByAge
data <- data[date <= max(deathByAge$date) + 7]

# save
saveRDS(data, file.path(indir, "data", paste0("vaccination-prop-pop-2022-02-06.rds")))
