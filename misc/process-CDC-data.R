library(rstan)
library(data.table)
library(dplyr)

indir = "~/git/CDC-covid19-agespecific-mortality-data" # path to the repo
outdir = file.path(indir, 'inst', "data")

# load functions
source(file.path(indir, 'misc', "functions.R"))
source(file.path(indir, 'misc', "utils.R"))

# max age considered
age_max = 105

# Gather CDC data
last.day = Sys.Date() - 1 # yesterday

# first part only available with Male and Female separation
deathByAge_Male = prepare_CDC_data(last.day, age_max, sex = 'Male', indir)
deathByAge_Male = find_daily_deaths(deathByAge_Male)
deathByAge_Female = prepare_CDC_data(last.day, age_max, sex = 'Female', indir)
deathByAge_Female = find_daily_deaths(deathByAge_Female)
deathByAge = merge_deathByAge_over_Sex(copy(deathByAge_Male), copy(deathByAge_Female))

deathByAge_AllSexes = prepare_CDC_data(last.day, age_max, sex = 'All Sexes', indir)
deathByAge_AllSexes = find_daily_deaths(deathByAge_AllSexes, rm.COVID.19.Deaths = F)
deathByAge_AllSexes = select(deathByAge_AllSexes, -c(date_idx, min_date_idx, max_date_idx))
deathByAge_res = incorporate_AllSexes_information(deathByAge, deathByAge_AllSexes)

deathByAge_res[!is.na(daily.deaths), min.sum.daily.deaths := NA]
deathByAge_res[!is.na(daily.deaths), max.sum.daily.deaths := NA]
deathByAge_res[!is.na(daily.deaths), sum.daily.deaths := NA]

# last change
# on 2020-07-04 repeated update so the daily death on 2020-07-04 and 2020-07-11 are not obtainable
all(na.omit(subset(deathByAge_res, date == "2020-07-04")$daily.deaths == 0))
deathByAge_res = subset(deathByAge_res, !date %in% c(as.Date("2020-07-04"), as.Date("2020-07-11")))

saveRDS(deathByAge_res, file.path(outdir, paste0('CDC-data_', last.day, '.rds')))



