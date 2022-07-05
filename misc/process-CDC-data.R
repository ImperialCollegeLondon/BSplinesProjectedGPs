library(data.table)
library(dplyr)

indir = "~/git/BSplinesProjectedGPs" # path to the repo
outdir = file.path(indir, 'inst', "data")

# load functions
source(file.path(indir, 'misc', "functions.R"))
source(file.path(indir, 'misc', "utils.R"))

# max age considered
age_max = 105

# Gather CDC data
last.week = Sys.Date() - 1 # yesterweek

# first part only available with Male and Female separation
deathByAge_Male = prepare_CDC_data(last.week, age_max, sex = 'Male', indir)
deathByAge_Male = find_weekly_deaths(deathByAge_Male)
deathByAge_Female = prepare_CDC_data(last.week, age_max, sex = 'Female', indir)
deathByAge_Female = find_weekly_deaths(deathByAge_Female)
deathByAge = merge_deathByAge_over_Sex(copy(deathByAge_Male), copy(deathByAge_Female))

# second part available for all sexes 
deathByAge_AllSexes = prepare_CDC_data(last.week, age_max, sex = 'All Sexes', indir)
deathByAge_AllSexes = find_weekly_deaths(deathByAge_AllSexes, rm.COVID.19.Deaths = F)
deathByAge_AllSexes = select(deathByAge_AllSexes, -c(date_idx, min_date_idx, max_date_idx))
deathByAge_res = incorporate_AllSexes_information(deathByAge, deathByAge_AllSexes)

# sum NY and NYC 
deathByAge_final <- merge_deathByAge_over_NY_NYC(deathByAge_res)

# rm boudaries is weekly deaths is not NA
deathByAge_final[!is.na(weekly.deaths), min.sum.weekly.deaths := NA]
deathByAge_final[!is.na(weekly.deaths), max.sum.weekly.deaths := NA]
deathByAge_final[!is.na(weekly.deaths), sum.weekly.deaths := NA]

# all weekly deaths are 0 on 2020-06-27 --> repeated update
# 2020-07-04 and 2020-06-27 are missing
# deathByAge_res[, list(sum(na.omit(weekly.deaths))), by = 'date']
all(na.omit(subset(deathByAge_final, date == "2020-06-27")$weekly.deaths == 0))
deathByAge_final = subset(deathByAge_final, !date %in% c(as.Date("2020-07-04"), as.Date("2020-06-27")))

saveRDS(deathByAge_final, file.path(outdir, paste0('CDC-data_', last.week, '.rds')))

