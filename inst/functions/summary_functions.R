create_map_age = function(age_max){
  # create map by 5-year age bands
  df_age_continuous <<- data.table(age_from = 0:age_max,
                                   age_to = 0:age_max,
                                   age_index = 1:(age_max+1),
                                   age = 0:age_max)
  
  # create map for reporting age groups before 2020-09-02
  df_age_reporting <<- data.table(age_from = c(0,1, 5,15,25,35,45,55,65,75,85),
                                  age_to = c(0,4,14,24,34,44,54,64,74,84,age_max),
                                  age_index = 1:11,
                                  age = c('0-0', '1-4', '5-14', '15-24', '25-34', '35-44', '45-54', '55-64', '65-74', '75-84', '85+'))
  df_age_reporting[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_age_reporting[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  
}

clean_vaccination_data = function(file){
  vaccinedata = as.data.table(read.csv(file))
  vaccinedata = subset(vaccinedata, Demographic.Group %in% c("Ages_65-74_yrs", "Ages_75+_yrs"), 
                       select = c('Date', 'Demographic.Group', 'People.with.at.least.one.dose', 'People.who.are.fully.vaccinated', 'Census'))
  setnames(vaccinedata, 1:4, c('date', 'age', 'count_vaccinated_1dosep', 'count_vaccinated_fully'))
  vaccinedata[, age := gsub('Ages_(.+)_yrs', '\\1', age)]
  tmp = vaccinedata[, list(count_vaccinated_1dosep= sum(count_vaccinated_1dosep), 
                           count_vaccinated_fully = sum(count_vaccinated_fully), 
                           Census = sum(Census)), by = c('date')]
  tmp[, age := '65+']
  vaccinedata = rbind(vaccinedata, tmp)
  vaccinedata[, prop_vaccinated_1dosep := count_vaccinated_1dosep / Census]
  vaccinedata[, prop_vaccinated_fully := count_vaccinated_fully / Census]
  vaccinedata[, date := as.Date(date)]
  return(vaccinedata)
}

reduce_agebands_scrapedData_GA = function(tmp)
{
  old_age = data.table(age = c("0-4","5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", 
                               "65-69", "70-74", "75-79", "80-84", "85-89", "90+"), 
                       age_index = c(rep(1:8, each = 2), rep(9, 3)))
  new_age = data.table(age_index = 1:9, 
                       age = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"))
  
  tmp = as.data.table(tmp)
  tmp = merge(tmp, old_age, by = 'age')
  tmp = tmp[, list(cum.deaths = sum(cum.deaths), 
                   daily.deaths = sum(daily.deaths)), by= c('code', 'date', 'age_index')]
  tmp = merge(tmp, new_age, by = 'age_index')
  tmp = select(tmp, -age_index)
  return(tmp)
}
