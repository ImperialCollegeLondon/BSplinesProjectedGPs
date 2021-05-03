prepare_CDC_data = function(last.day,age_max,sex,indir, check_increasing_cumulative_deaths =0)
  {
  
  path_to_data = file.path(indir, 'misc', 'data')
  
  dates = seq.Date(as.Date("2020-03-01"), last.day, by = "day")
  
  file_format = ".csv"
  
  # find dates with data
  data_files = list.files(path_to_data, full.names = T)
  dates = as.Date(gsub( ".*\\/CDC_data_(.+).csv*", "\\1", data_files))
  
  # we remove the first date because it the data set is cut after OH
  dates = dates[-1]
  
  tmp = vector(mode = 'list', length = length(dates))
  idx.rm = c()
  for(t in 1:length(dates)){
    
    csv_file = file.path(path_to_data, paste0('CDC_data_', dates[t], '.csv'))
    tmp[[t]] = as.data.table( read.csv(csv_file) ) 
    
    # account for the change of variables over time
    if('state' %in% names(tmp[[t]])) setnames(tmp[[t]], 'state', 'State')
    if('age_group' %in% names(tmp[[t]])) setnames(tmp[[t]], 'age_group', 'Age.group')
    if('sex' %in% names(tmp[[t]])) setnames(tmp[[t]], 'sex', 'Sex')
    if('covid_19_deaths' %in% names(tmp[[t]])) setnames(tmp[[t]], 'covid_19_deaths', 'COVID.19.Deaths')
    
    if('end_week' %in% names(tmp[[t]])) setnames(tmp[[t]], 'end_week', 'End.Week')
    if('End.Date' %in% names(tmp[[t]])) setnames(tmp[[t]], 'End.Date', 'End.Week')
    if('Age.Group' %in% names(tmp[[t]])) setnames(tmp[[t]], 'Age.Group', 'Age.group')
    
    if('Group'%in% names(tmp[[t]])) tmp[[t]] = subset(tmp[[t]], Group == 'By Total')
    
    stopifnot(length(unique(tmp[[t]]$End.Week)) == 1)
    
    if(dates[t] == "2020-07-24") # fix bug in the data
      tmp[[t]][, End.Week := '07/18/2020']
    
    if(t > 1) 
      if(unique(tmp[[t]]$End.Week) == unique(tmp[[t-1]]$End.Week) )
        {
        idx.rm = c(idx.rm, t)
        next
      }
    
    tmp[[t]] = select(tmp[[t]], State, 'End.Week', Sex, Age.group, COVID.19.Deaths)
  }
  tmp = do.call('rbind', tmp[-idx.rm])
  
  # set date variable
  tmp[, date := as.Date(End.Week, format = '%m/%d/%Y')]
  tmp = select(tmp, -End.Week)

  # check that all the weeks are recorded
  real_dates = seq.Date(as.Date("2020-05-02"), as.Date("2021-04-24"), by = 'week')
  stopifnot( length( real_dates[!real_dates %in% unique(tmp$date)] ) == 0)
  
  # boundaries if deaths is missing
  tmp[, min_COVID.19.Deaths := 1]
  tmp[, max_COVID.19.Deaths := 9]
  
  # choose sex
  tmp = subset(tmp, Sex == sex)
  
  # bugfix if there is a 0 before and after NA
  tmp = fix_inconsistent_NA_between0(tmp)
  tmp = fix_inconsistent_NA_betweenpos(tmp)
  
  # ensure increasing cumulative deaths
  if(check_increasing_cumulative_deaths){ # costly computationally, run when new data
    tmp = ensure_increasing_cumulative_deaths_origin(tmp)
  } else{ # bugfix already noticed
    tmp = bugfix_nonincreasing_cumulative_deaths(tmp)
  }
  
  # rename age groups
  tmp[, Age.group := ifelse(Age.group == "Under 1 year", "0-0", 
                            ifelse(Age.group == "85 years and over", "85+", gsub("(.+) years", "\\1", Age.group)))]
  
  # group age groups
  tmp = group.age.specification.1(tmp)
  
  # rm US and add code
  setnames(tmp, c('Age.group', 'State'), c("age", "loc_label"))
  tmp = subset(tmp, loc_label != 'United States')
  tmp = merge(tmp, map_statename_code, by.x = 'loc_label', by.y = 'State')
  
  # find age from and age to
  tmp[, age_from := as.numeric(ifelse(grepl("\\+", age), gsub("(.+)\\+", "\\1", age), gsub("(.+)-.*", "\\1", age)))]
  tmp[, age_to := as.numeric(ifelse(grepl("\\+", age), age_max, gsub(".*-(.+)", "\\1", age)))]
  
  # order
  tmp = tmp[order(loc_label, date, age)]
  
  # ensure that it is cumulative
  tmp1 = tmp[, list(noncum = na.omit(COVID.19.Deaths) <= cummax(na.omit(COVID.19.Deaths))), by = c('loc_label', 'date', 'age')]
  stopifnot(all(tmp1$noncum))
  
  # plot
  if(0){
    ggplot(tmp, aes(x = date, y = COVID.19.Deaths, col = age)) + 
      geom_line() + 
      facet_grid(loc_label~.,  scales = 'free') + 
      theme_bw()
  }
  
  return(tmp)
}

group.age.specification.1 = function(tmp)
{

  # factor age
  tmp = subset(tmp, Age.group %in% c('0-0', '1-4', '5-14', '15-24', '25-34', '35-44', '45-54', '55-64', '65-74', '75-84', '85+'))
  tmp[, Age.group := factor(Age.group, c('0-0', '1-4', '5-14', '15-24', '25-34', '35-44', '45-54', '55-64', '65-74', '75-84', '85+'))]
  
  # rm overall
  tmp = subset(tmp, !is.na(Age.group))
  
  # check that the number of age group is the same for every stata/date combinations
  tmp1 = tmp[, list(N = .N), by = c('State', 'date')]
  stopifnot(all(tmp1$N == 11))
  
  # sanity checks
  tmp1 = tmp[, list(idx_last_NA = which(is.na(COVID.19.Deaths))), by = c("State", 'Age.group')]
  tmp2 = tmp1[, list(min_idx_last_NA = min(idx_last_NA), max_idx_last_NA = max(idx_last_NA)), by = c("State", 'Age.group')]
  tmp1 = merge(tmp1, tmp2, by = c("State", 'Age.group'))
  tmp1 = tmp1[, list(all.inside = all( c(unique(min_idx_last_NA):unique(max_idx_last_NA)) %in% idx_last_NA )), by = c("State", 'Age.group')]
  stopifnot(all(tmp1$all.inside == T))
  
  return(tmp)
}

find_daily_deaths = function(tmp, rm.COVID.19.Deaths= T)
  {

  # ensure incrinsing cumulative deaths
  tmp = ensure_increasing_cumulative_deaths(tmp)
  
  # find date index with NA
  tmp1 = tmp[, list(idx_NA = which(is.na(COVID.19.Deaths))), by = c("loc_label", 'age')]
  tmp1[, min_idx_NA := min(idx_NA), by = c("loc_label", 'age')]
  tmp1[, max_idx_NA := max(idx_NA), by = c("loc_label", 'age')]
  tmp2 = tmp1[, list(is.inside = all( seq(unique(min_idx_NA), unique(max_idx_NA), 1) %in% idx_NA)), by = c("loc_label", 'age')]
  stopifnot(all(tmp2$is.inside == T))
  
  tmp2 = tmp1[, list(min_idx_NA = min(idx_NA), max_idx_NA = max(idx_NA)), by = c("loc_label", 'age')]
  tmp1 = unique(select(tmp, "loc_label", 'age'))
  tmp2 = merge(tmp1, tmp2, by = c("loc_label", 'age'), all.x = T)
  tmp = merge(tmp, tmp2, by = c("loc_label", 'age'))
  tmp[, date_idx := 1:length(date), by = c("loc_label", 'age')]
  tmp1 = tmp[, list(min_date_idx = min(date_idx), max_date_idx = max(date_idx)), by = c("loc_label", 'age')]
  tmp = merge(tmp, tmp1, by = c("loc_label", 'age'))
  
  # find daily deaths
  tmp[, daily.deaths := c(NA, diff(COVID.19.Deaths))]
  tmp[, daily.deaths := as.numeric(daily.deaths)]
  tmp[, min.sum.daily.deaths := NA]
  tmp[, min.sum.daily.deaths := as.numeric(min.sum.daily.deaths)]
  tmp[, max.sum.daily.deaths := NA]
  tmp[, max.sum.daily.deaths := as.numeric(max.sum.daily.deaths)]
  tmp[, sum.daily.deaths := NA]
  tmp[, sum.daily.deaths := as.numeric(sum.daily.deaths)]
  
  # first day, daily death is NA
  tmp[date_idx == 1, daily.deaths := NA]

  # boundaries deaths if NA
  # contained within the period
  tmp[min_idx_NA != min_date_idx & max_idx_NA != max_date_idx, sum.daily.deaths := COVID.19.Deaths[max_idx_NA + 1], by = c("loc_label", 'age')]
  
  # beginning of the period
  tmp[min_idx_NA == min_date_idx & max_idx_NA != max_date_idx, min.sum.daily.deaths := COVID.19.Deaths[max_idx_NA + 1] - max_COVID.19.Deaths, by = c("loc_label", 'age')]
  tmp[min_idx_NA == min_date_idx & max_idx_NA != max_date_idx, max.sum.daily.deaths := COVID.19.Deaths[max_idx_NA + 1] - min_COVID.19.Deaths, by = c("loc_label", 'age')]
  
  # end of the period 
  tmp[min_idx_NA != min_date_idx & max_idx_NA == max_date_idx, min.sum.daily.deaths := min_COVID.19.Deaths, by = c("loc_label", 'age')]
  tmp[min_idx_NA != min_date_idx & max_idx_NA == max_date_idx, max.sum.daily.deaths := max_COVID.19.Deaths, by = c("loc_label", 'age')]
  
  # entire period
  tmp[min_idx_NA == min_date_idx & max_idx_NA == max_date_idx, min.sum.daily.deaths := 0, by = c("loc_label", 'age')]
  tmp[min_idx_NA == min_date_idx & max_idx_NA == max_date_idx, max.sum.daily.deaths := max_COVID.19.Deaths - 1, by = c("loc_label", 'age')]
  
  # remove first date 
  dates = sort(unique(tmp$date))
  tmp = subset(tmp, date != dates[1])
  
  # checks
  stopifnot(nrow(tmp[daily.deaths < 0]) == 0)
  tmp1 = subset(tmp, is.na(daily.deaths))
  stopifnot(all( !is.na(tmp1$min.sum.daily.deaths) | !is.na(tmp1$sum.daily.deaths) ))
  stopifnot(all( !is.na(tmp1$max.sum.daily.deaths) | !is.na(tmp1$sum.daily.deaths) ))
  tmp2 = subset(tmp1, !is.na(max.sum.daily.deaths))
  stopifnot(all(!is.na(tmp2$min.sum.daily.deaths)))
  stopifnot(all(tmp2$min.sum.daily.deaths < tmp2$max.sum.daily.deaths))
  tmp2 = subset(tmp1, !is.na(sum.daily.deaths))
  stopifnot(all(tmp2$sum.daily.deaths >= 0))
  
  tmp = select(tmp, -min_COVID.19.Deaths, -max_COVID.19.Deaths, -min_idx_NA, -max_idx_NA)
  
  if(rm.COVID.19.Deaths)
    tmp = select(tmp, -COVID.19.Deaths)
  
  return(tmp)
}

merge_deathByAge_over_Sex = function(tmp1, tmp2)
  {
  
  setnames(tmp1, c('daily.deaths', 'sum.daily.deaths', 'min.sum.daily.deaths', 'max.sum.daily.deaths'), 
           c('daily.deaths.sex_1', 'sum.daily.deaths.sex_1', 'min.sum.daily.deaths.sex_1', 'max.sum.daily.deaths.sex_1'))
  setnames(tmp2, c('daily.deaths', 'sum.daily.deaths', 'min.sum.daily.deaths', 'max.sum.daily.deaths'), 
           c('daily.deaths.sex_2', 'sum.daily.deaths.sex_2', 'min.sum.daily.deaths.sex_2', 'max.sum.daily.deaths.sex_2'))
  
  # find date index with NA
  tmp3 = tmp1[, list(idx_NA = which(is.na(daily.deaths.sex_1))), by = c("loc_label", 'age')]
  tmp3 = tmp3[, list(min_idx_NA.sex_1 = min(idx_NA), max_idx_NA.sex_1 = max(idx_NA)), by = c("loc_label", 'age')]
  tmp1 = merge(tmp1, tmp3, by = c("loc_label", 'age'), all.x = T)
  tmp3 = tmp2[, list(idx_NA = which(is.na(daily.deaths.sex_2))), by = c("loc_label", 'age')]
  tmp3 = tmp3[, list(min_idx_NA.sex_2 = min(idx_NA), max_idx_NA.sex_2 = max(idx_NA)), by = c("loc_label", 'age')]
  tmp2 = merge(tmp2, tmp3, by = c("loc_label", 'age'), all.x = T)
  
  tmp = unique(select(tmp1, loc_label, code, age, age_from, age_to))
  
  tmp1 = merge(tmp1, tmp2, by = c('age', 'date', 'loc_label', 'date_idx', 'min_date_idx', 'max_date_idx'))
  
  # reajust date idx after emoval of the first day
  tmp1 = select(tmp1, -date_idx, -min_date_idx, -max_date_idx)
  tmp1[, date_idx := 1:length(date), by = c("loc_label", 'age')]
  tmp3 = tmp1[, list(min_date_idx = min(date_idx), max_date_idx = max(date_idx)), by = c("loc_label", 'age')]
  tmp1 = merge(tmp1, tmp3, by = c("loc_label", 'age'))
  
  tmp1[, daily.deaths := as.numeric(NA)]
  tmp1[, daily.deaths := as.numeric(daily.deaths)]
  tmp1[, min.sum.daily.deaths := as.numeric(NA)]
  tmp1[, min.sum.daily.deaths := as.numeric(min.sum.daily.deaths)]
  tmp1[, max.sum.daily.deaths := as.numeric(NA)]
  tmp1[, max.sum.daily.deaths := as.numeric(max.sum.daily.deaths)]
  tmp1[, sum.daily.deaths := as.numeric(NA)]
  tmp1[, sum.daily.deaths := as.numeric(sum.daily.deaths)]

  # both are at after the beginning of the period and both finish before the interval: 1
  # both are at after the beginning of the period and both finish after the interval: 1
  # both are at the beginning of the period and both finish after the interval: 1
  # both are at the beginning of the period and both finish before the interval: 1
  
  # both are at the beginning of the period and one finish after the interval: 2
  # both are after the beginning of the period and one finish after the interval: 2
  # one is on the beginning of the period and  both finish before the period: 2
  # one is on the beginning of the period and both finish after the period: 2
  
  # one is on the beginning of the period and one does not finish before the period: 4
  
  # both are after the beginning of the period and both finish before the interval
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
        max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
        max_idx_NA.sex_1 < max_idx_NA.sex_2, 
      sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
                          sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1, 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1, 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 , by = c('age', 'loc_label')]
  
  # both are after the beginning of the period and both finish after the interval
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx & 
         (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         !is.na(sum.daily.deaths.sex_1) & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2, by = c('age', 'loc_label')]

  # both are at the beginning of the period and both finish after the interval
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         !is.na(sum.daily.deaths.sex_1) & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2, by = c('age', 'loc_label')]
  
  # both are at the beginning of the period and both finish before the interval
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2  & (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2  & !is.na(sum.daily.deaths.sex_1) & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2)+ 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & !is.na(sum.daily.deaths.sex_1) & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_1 +  sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & !is.na(sum.daily.deaths.sex_1) & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2, by = c('age', 'loc_label')]
  
  # both are at the beginning of the period and one finish after the interval
  # sex 2 finished after the interval
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         !is.na(sum.daily.deaths.sex_1) & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  
  # sex 1 finished after the interval
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         (is.na(sum.daily.deaths.sex_1) |is.na(sum.daily.deaths.sex_2)), 
       max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         !is.na(sum.daily.deaths.sex_1) & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  
  # both are after the beginning of the period and one finish after the interval
  # sex 2 finish after
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         is.na(sum.daily.deaths.sex_2), 
       min.sum.daily.deaths := sum.daily.deaths.sex_1 + min.sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         is.na(sum.daily.deaths.sex_2), 
       max.sum.daily.deaths := sum.daily.deaths.sex_1 + max.sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         !is.na(sum.daily.deaths.sex_1), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  
  
  # sex 1 finish after
   tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
          is.na(sum.daily.deaths.sex_1), 
       min.sum.daily.deaths := min.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
   tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
          max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
          is.na(sum.daily.deaths.sex_1), 
        max.sum.daily.deaths := max.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
          sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
   tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
          max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
          !is.na(sum.daily.deaths.sex_1), 
        sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
          sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
   
  # one is at the beginning of the period and both finish before the period
  # sex 1 is at the beginning
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & is.na(sum.daily.deaths.sex_1), 
       min.sum.daily.deaths := min.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & is.na(sum.daily.deaths.sex_1), 
       max.sum.daily.deaths := max.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & !is.na(sum.daily.deaths.sex_1), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & is.na(sum.daily.deaths.sex_1), 
       min.sum.daily.deaths := min.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & is.na(sum.daily.deaths.sex_1), 
       max.sum.daily.deaths := max.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & !is.na(sum.daily.deaths.sex_1), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & is.na(sum.daily.deaths.sex_1), 
       min.sum.daily.deaths := min.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & is.na(sum.daily.deaths.sex_1), 
       max.sum.daily.deaths := max.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2, by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & !is.na(sum.daily.deaths.sex_1), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2, by = c('age', 'loc_label')]
  
  # sex 2 is at the beginning
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & is.na(sum.daily.deaths.sex_2), 
       min.sum.daily.deaths := min.sum.daily.deaths.sex_2 + sum.daily.deaths.sex_1 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & is.na(sum.daily.deaths.sex_2), 
       max.sum.daily.deaths := max.sum.daily.deaths.sex_2 + sum.daily.deaths.sex_1 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_2 + sum.daily.deaths.sex_1 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & is.na(sum.daily.deaths.sex_2), 
       min.sum.daily.deaths := min.sum.daily.deaths.sex_2 + sum.daily.deaths.sex_1 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & is.na(sum.daily.deaths.sex_2), 
       max.sum.daily.deaths := max.sum.daily.deaths.sex_2 + sum.daily.deaths.sex_1 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_2 + sum.daily.deaths.sex_1 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & is.na(sum.daily.deaths.sex_2), 
       min.sum.daily.deaths := min.sum.daily.deaths.sex_2 + sum.daily.deaths.sex_1 , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & is.na(sum.daily.deaths.sex_2), 
       max.sum.daily.deaths := max.sum.daily.deaths.sex_2 + sum.daily.deaths.sex_1, by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & !is.na(sum.daily.deaths.sex_2), 
       sum.daily.deaths := sum.daily.deaths.sex_2 + sum.daily.deaths.sex_1, by = c('age', 'loc_label')]
  
  # one is at the beginning of the period and both finish after the period
  # sex 1 is at the beginning
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx, 
       min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1)  + 
         min.sum.daily.deaths.sex_2, by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx, 
       max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         max.sum.daily.deaths.sex_2, by = c('age', 'loc_label')]

  # sex 2 is at the beginning
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx , 
       min.sum.daily.deaths := min.sum.daily.deaths.sex_1 + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx, 
       max.sum.daily.deaths := max.sum.daily.deaths.sex_1 + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2), by = c('age', 'loc_label')]
  
  # one is at the beginning of the period and one finish after the period
  # sex 1 is at the beginning and sex 1 finish after the period
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         is.na(sum.daily.deaths.sex_1), 
       min.sum.daily.deaths := min.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx&
         is.na(sum.daily.deaths.sex_1), 
         max.sum.daily.deaths := max.sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx&
         !is.na(sum.daily.deaths.sex_1), 
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  
  # sex 1 is at the beginning and sex 2 does not finish  before the period
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.daily.deaths.sex_1) | is.na(sum.daily.deaths.sex_2)),
         min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         (is.na(sum.daily.deaths.sex_1) | is.na(sum.daily.deaths.sex_2)),
         max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         !is.na(sum.daily.deaths.sex_1) & !is.na(sum.daily.deaths.sex_2),
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  
  # sex 2 is at the beginning and sex 1 does not finish before the period
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         (is.na(sum.daily.deaths.sex_1) | is.na(sum.daily.deaths.sex_2)),
       min.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), min.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), min.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         (is.na(sum.daily.deaths.sex_1) | is.na(sum.daily.deaths.sex_2)),
       max.sum.daily.deaths := ifelse(is.na(unique(sum.daily.deaths.sex_1)), max.sum.daily.deaths.sex_1, sum.daily.deaths.sex_1) + 
         ifelse(is.na(unique(sum.daily.deaths.sex_2)), max.sum.daily.deaths.sex_2, sum.daily.deaths.sex_2) + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         !is.na(sum.daily.deaths.sex_1) & !is.na(sum.daily.deaths.sex_2),
       sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  
  # sex 2 is at the beginning and sex 2 does not finish before the period
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         is.na(sum.daily.deaths.sex_2),
         min.sum.daily.deaths := sum.daily.deaths.sex_1 + min.sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         is.na(sum.daily.deaths.sex_2),
         max.sum.daily.deaths := sum.daily.deaths.sex_1 + max.sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         !is.na(sum.daily.deaths.sex_2),
          sum.daily.deaths := sum.daily.deaths.sex_1 + sum.daily.deaths.sex_2 + 
         sum(daily.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  
  # sex 1 doesn't have any NA
  tmp1[is.na(min_idx_NA.sex_1) & !is.na(min_idx_NA.sex_2) & !is.na(sum.daily.deaths.sex_2), 
       daily.deaths := sum(daily.deaths.sex_1[unique(min_idx_NA.sex_2):unique(max_idx_NA.sex_2)]) + 
         sum.daily.deaths.sex_2, by = c('age', 'loc_label')]
  tmp1[is.na(min_idx_NA.sex_1) & !is.na(min_idx_NA.sex_2) & is.na(sum.daily.deaths.sex_2), 
       min.sum.daily.deaths := sum(daily.deaths.sex_1[unique(min_idx_NA.sex_2):unique(max_idx_NA.sex_2)]) + 
         min.sum.daily.deaths.sex_2, by = c('age', 'loc_label')]
  tmp1[is.na(min_idx_NA.sex_1) & !is.na(min_idx_NA.sex_2) & is.na(sum.daily.deaths.sex_2), 
       max.sum.daily.deaths := sum(daily.deaths.sex_1[unique(min_idx_NA.sex_2):unique(max_idx_NA.sex_2)]) + 
         max.sum.daily.deaths.sex_2, by = c('age', 'loc_label')]
  
  # sex 2 doesn't have any NA
  tmp1[!is.na(min_idx_NA.sex_1) & is.na(min_idx_NA.sex_2) & !is.na(sum.daily.deaths.sex_1), 
       daily.deaths := sum(daily.deaths.sex_2[unique(min_idx_NA.sex_1):unique(max_idx_NA.sex_1)]) + 
         sum.daily.deaths.sex_1, by = c('age', 'loc_label')]
  tmp1[!is.na(min_idx_NA.sex_1) & is.na(min_idx_NA.sex_2) & is.na(sum.daily.deaths.sex_1), 
       min.sum.daily.deaths := sum(daily.deaths.sex_2[unique(min_idx_NA.sex_1):unique(max_idx_NA.sex_1)]) + 
         min.sum.daily.deaths.sex_1, by = c('age', 'loc_label')]
  tmp1[!is.na(min_idx_NA.sex_1) & is.na(min_idx_NA.sex_2) & is.na(sum.daily.deaths.sex_1), 
       max.sum.daily.deaths := sum(daily.deaths.sex_2[unique(min_idx_NA.sex_1):unique(max_idx_NA.sex_1)]) + 
         max.sum.daily.deaths.sex_1, by = c('age', 'loc_label')]

  
  # non missing daily deaths
  tmp1[daily.deaths.sex_1 >= 0 & daily.deaths.sex_2 >= 0, daily.deaths := daily.deaths.sex_1 + daily.deaths.sex_2, by = c('age', 'loc_label')]
  
  # checks 
  tmp2 = subset(tmp1, is.na(daily.deaths))
  which(!is.na(tmp2$min.sum.daily.deaths) | !is.na(tmp2$sum.daily.deaths))
  tmp2[ !(!is.na(tmp2$min.sum.daily.deaths) | !is.na(tmp2$sum.daily.deaths)) ]
  stopifnot(all( !is.na(tmp2$min.sum.daily.deaths) | !is.na(tmp2$sum.daily.deaths) ))
  stopifnot(all( !is.na(tmp2$max.sum.daily.deaths) | !is.na(tmp2$sum.daily.deaths)))
  stopifnot(!any( !is.na(tmp2$min.sum.daily.deaths) & !is.na(tmp2$sum.daily.deaths) ) )
  stopifnot(!any( !is.na(tmp2$max.sum.daily.deaths) & !is.na(tmp2$sum.daily.deaths)))
  tmp3 = subset(tmp1, !is.na(min.sum.daily.deaths))
  stopifnot(all(!is.na(tmp3$max.sum.daily.deaths)))
  stopifnot(all(tmp3$min.sum.daily.deaths < tmp3$max.sum.daily.deaths))
  tmp3 = subset(tmp1, !is.na(max.sum.daily.deaths))
  stopifnot(all(!is.na(tmp3$min.sum.daily.deaths)))
  stopifnot(all(tmp3$min.sum.daily.deaths < tmp3$max.sum.daily.deaths))
  tmp3 = subset(tmp1, !is.na(sum.daily.deaths))
  stopifnot(all(tmp3$sum.daily.deaths >= 0))
  
  for(Loc in unique(tmp1$loc_label)){
    for(Age in unique(tmp1$age)){
      tmp3 = subset(tmp1, loc_label == Loc & age == Age)
      .idx_missing = which(is.na(tmp3$daily.deaths))
      
      if(length(.idx_missing) == 0)
        next
      stopifnot( length(unique(tmp3$min.sum.daily.deaths[.idx_missing])) == 1)
      stopifnot( length(unique(tmp3$max.sum.daily.deaths[.idx_missing])) == 1)
      stopifnot( length(unique(tmp3$sum.daily.deaths[.idx_missing])) == 1)
    }
  }
  
  # final
  tmp1 = select(tmp1, 'date', 'loc_label', 'age', 'min.sum.daily.deaths', 'max.sum.daily.deaths', 'sum.daily.deaths', 'daily.deaths')
  tmp = merge(tmp, tmp1, by = c('age', 'loc_label'))
  
  return(tmp)
}

incorporate_AllSexes_information = function(tmp1, tmp2)
{
  tmp = incorporate_AllSexes_boundary_information(tmp1, tmp2)
  tmp = incorporate_AllSexes_totalsum_information(tmp)
  
  return(tmp)
}

incorporate_AllSexes_totalsum_information = function(tmp)
{
  
  for(Loc in unique(tmp$loc_label)){
    for(Age in unique(tmp$age)){
      
      tmp3 = subset(tmp, loc_label == Loc & age == Age)
      .idx_missing = which(is.na(tmp3$daily.deaths))
      
      if(length(.idx_missing) == 0)
        next
      
      if(all(!is.na(tmp3$sum.daily.deaths[.idx_missing])))
        next
      
      stopifnot( length(unique(tmp3$min.sum.daily.deaths[.idx_missing])) <= 2)
      stopifnot( length(unique(tmp3$max.sum.daily.deaths[.idx_missing])) <= 2)
      
      stopifnot(tmp3[.idx_missing]$min.sum.daily.deaths[1] >= tmp3[.idx_missing]$min.sum.daily.deaths[length(.idx_missing)])
      stopifnot(tmp3[.idx_missing]$max.sum.daily.deaths[1] >= tmp3[.idx_missing]$max.sum.daily.deaths[length(.idx_missing)])
      
      if(tmp3[.idx_missing]$min.sum.daily.deaths[1] > tmp3[.idx_missing]$min.sum.daily.deaths[length(.idx_missing)])
        tmp3[.idx_missing]$min.sum.daily.deaths = tmp3[.idx_missing]$min.sum.daily.deaths[length(.idx_missing)]
      
      if(tmp3[.idx_missing]$max.sum.daily.deaths[1] > tmp3[.idx_missing]$max.sum.daily.deaths[length(.idx_missing)])
        tmp3[.idx_missing]$max.sum.daily.deaths =  tmp3[.idx_missing]$max.sum.daily.deaths[length(.idx_missing)]
      
      if(max(.idx_missing) == nrow(tmp3) & min(.idx_missing) != 1){
        tmp3[.idx_missing]$min.sum.daily.deaths = 1
        tmp3[.idx_missing]$max.sum.daily.deaths = 9
      }
      
      tmp = anti_join(tmp, tmp3, by = c('loc_label', 'age'))
      tmp = rbind(tmp, tmp3)
    }
  }
  
  tmp = tmp[order(loc_label, age, date)]
  return(tmp)
}

incorporate_AllSexes_boundary_information = function(tmp1, tmp2)
{
  
  locations = unique(tmp1$loc_label)
  ages = unique(tmp1$age)
  
  for(Loc in locations){
    for(Age in ages){
      
      tmp3 = subset(tmp1, loc_label == Loc & age == Age)
      tmp4 = subset(tmp2, loc_label == Loc & age == Age)
      tmp3 = subset(tmp3, date < min(tmp4$date))
      
      .idx_missing.daily = which(is.na(tmp3$daily.deaths))
      .idx_missing.daily2 = which(is.na(tmp4$daily.deaths))
      .idx_non_missing.cum = which(!is.na(tmp4$COVID.19.Deaths))
      
      if(length(.idx_missing.daily) == 0)
        next
      if(length(.idx_non_missing.cum) == 0)
        next
      if(!nrow(tmp3) %in% .idx_missing.daily)
        next
      
      # verify that sum daily death is equal to the first value of COVID.19.Deaths
      if(!is.na(tmp3[nrow(tmp3)]$sum.daily.deaths) & !is.na(tmp4[1]$COVID.19.Deaths))
      {
        if(is.na(tmp4[1]$daily.deaths)){
          tmp4[1]$daily.deaths = tmp3[nrow(tmp3)]$sum.daily.deaths
          tmp2 = anti_join(tmp2, tmp4, by = c('loc_label', 'age'))
          tmp2 = rbind(tmp2, tmp4)
          
        }
        if(tmp3[nrow(tmp3)]$sum.daily.deaths == tmp4[1]$COVID.19.Deaths)
          next
        if(tmp3[nrow(tmp3)]$sum.daily.deaths != tmp4[1]$COVID.19.Deaths){
          tmp3[.idx_missing.daily]$sum.daily.deaths = tmp4[1]$COVID.19.Deaths
          tmp1 = anti_join(tmp1, tmp3, by = c('loc_label', 'age'))
          tmp1 = rbind(tmp1, tmp3)
          
        }
      }

      # fill sum daily death with first value of COVID.19.Deaths if NA doesn't start at the beginning of the period
      first_non_missing_cum = tmp4$COVID.19.Deaths[.idx_non_missing.cum[1]]
      if(1 %in% .idx_missing.daily & is.na(tmp3[nrow(tmp3)]$sum.daily.deaths))
      {
        stopifnot(all(is.na(tmp3[.idx_missing.daily]$sum.daily.deaths )))
        tmp3[.idx_missing.daily]$max.sum.daily.deaths = first_non_missing_cum
        
        .idx_missing.daily2 = c(.idx_missing.daily2, .idx_non_missing.cum[1])
        stopifnot(all(is.na(tmp4[.idx_missing.daily2]$sum.daily.deaths )))
        tmp4[.idx_missing.daily2]$min.sum.daily.deaths = unique(tmp3[.idx_missing.daily]$min.sum.daily.deaths)
        tmp4[.idx_missing.daily2]$max.sum.daily.deaths = first_non_missing_cum
        tmp4[.idx_non_missing.cum[1]]$daily.deaths = NA
      } else {
        
        tmp3[.idx_missing.daily]$sum.daily.deaths = first_non_missing_cum
        tmp3[.idx_missing.daily]$min.sum.daily.deaths = NA
        tmp3[.idx_missing.daily]$max.sum.daily.deaths = NA
        .idx_missing.daily2 = c(.idx_missing.daily2, .idx_non_missing.cum[1])
        tmp4[.idx_missing.daily2]$sum.daily.deaths = first_non_missing_cum
        tmp4[.idx_missing.daily2]$min.sum.daily.deaths = NA
        tmp4[.idx_missing.daily2]$max.sum.daily.deaths = NA
        tmp4[.idx_non_missing.cum[1]]$daily.deaths = NA
      }
      
      tmp1 = anti_join(tmp1, tmp3, by = c('loc_label', 'age'))
      tmp1 = rbind(tmp1, tmp3)
      
      tmp2 = anti_join(tmp2, tmp4, by = c('loc_label', 'age'))
      tmp2 = rbind(tmp2, tmp4)

    }
  }
  
  tmp2 = select(tmp2,  - COVID.19.Deaths )
  tmp1 = anti_join(tmp1, tmp2, by = c('loc_label', 'age', 'date'))
  tmp = rbind(tmp1, select(tmp2, -Sex))
  tmp = tmp[order(loc_label, age, date)]
  
  return(tmp)
}

