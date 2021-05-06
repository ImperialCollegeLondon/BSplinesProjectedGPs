prepare_CDC_data = function(last.week,age_max,sex,indir, check_increasing_cumulative_deaths =0)
  {
  
  path_to_data = file.path(indir, 'misc', 'data')
  
  dates = seq.Date(as.Date("2020-03-01"), last.week, by = "week")
  
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
  
  # rm missing
  tmp = subset(tmp, !is.na(Age.group))
  
  # check that the number of age group is the same for every state/date combinations
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

find_weekly_deaths = function(tmp, rm.COVID.19.Deaths= T)
  {

  # ensure increasing cumulative deaths
  tmp = ensure_increasing_cumulative_deaths(tmp)
  
  # find date index with NA
  tmp1 = tmp[, list(idx_NA = which(is.na(COVID.19.Deaths))), by = c("loc_label", 'age')]
  tmp1[, min_idx_NA := min(idx_NA), by = c("loc_label", 'age')]
  tmp1[, max_idx_NA := max(idx_NA), by = c("loc_label", 'age')]
  
  # check that na appear series
  tmp2 = tmp1[, list(is.inside = all( seq(unique(min_idx_NA), unique(max_idx_NA), 1) %in% idx_NA)), by = c("loc_label", 'age')]
  stopifnot(all(tmp2$is.inside == T))
  
  # merge
  tmp2 = tmp1[, list(min_idx_NA = min(idx_NA), max_idx_NA = max(idx_NA)), by = c("loc_label", 'age')]
  tmp1 = unique(select(tmp, "loc_label", 'age'))
  tmp2 = merge(tmp1, tmp2, by = c("loc_label", 'age'), all.x = T)
  tmp = merge(tmp, tmp2, by = c("loc_label", 'age'))
  
  # find date index
  tmp[, date_idx := 1:length(date), by = c("loc_label", 'age')]
  tmp1 = tmp[, list(min_date_idx = min(date_idx), max_date_idx = max(date_idx)), by = c("loc_label", 'age')]
  tmp = merge(tmp, tmp1, by = c("loc_label", 'age'))
  
  # find weekly deaths
  tmp[, weekly.deaths := c(NA_integer_, diff(COVID.19.Deaths))]
  tmp[, min.sum.weekly.deaths := NA_integer_]
  tmp[, max.sum.weekly.deaths := NA_integer_]
  tmp[, sum.weekly.deaths := NA_integer_]
  
  # first week, weekly death is NA
  tmp[date_idx == 1, weekly.deaths := NA_integer_]

  # boundaries deaths if NA
  # contained within the period
  tmp[min_idx_NA != min_date_idx & max_idx_NA != max_date_idx, sum.weekly.deaths := COVID.19.Deaths[max_idx_NA + 1], by = c("loc_label", 'age')]
  
  # beginning of the period
  tmp[min_idx_NA == min_date_idx & max_idx_NA != max_date_idx, min.sum.weekly.deaths := COVID.19.Deaths[max_idx_NA + 1] - max_COVID.19.Deaths, by = c("loc_label", 'age')]
  tmp[min_idx_NA == min_date_idx & max_idx_NA != max_date_idx, max.sum.weekly.deaths := COVID.19.Deaths[max_idx_NA + 1] - min_COVID.19.Deaths, by = c("loc_label", 'age')]
  
  # end of the period 
  tmp[min_idx_NA != min_date_idx & max_idx_NA == max_date_idx, min.sum.weekly.deaths := min_COVID.19.Deaths, by = c("loc_label", 'age')]
  tmp[min_idx_NA != min_date_idx & max_idx_NA == max_date_idx, max.sum.weekly.deaths := max_COVID.19.Deaths, by = c("loc_label", 'age')]
  
  # entire period
  tmp[min_idx_NA == min_date_idx & max_idx_NA == max_date_idx, min.sum.weekly.deaths := 0, by = c("loc_label", 'age')]
  tmp[min_idx_NA == min_date_idx & max_idx_NA == max_date_idx, max.sum.weekly.deaths := max_COVID.19.Deaths - 1, by = c("loc_label", 'age')]
  
  # remove first date 
  dates = sort(unique(tmp$date))
  tmp = subset(tmp, date != dates[1])
  
  # sanity checks
  stopifnot(nrow(tmp[weekly.deaths < 0]) == 0)
  tmp1 = subset(tmp, is.na(weekly.deaths))
  stopifnot(all( !is.na(tmp1$min.sum.weekly.deaths) | !is.na(tmp1$sum.weekly.deaths) ))
  stopifnot(all( !is.na(tmp1$max.sum.weekly.deaths) | !is.na(tmp1$sum.weekly.deaths) ))
  tmp2 = subset(tmp1, !is.na(max.sum.weekly.deaths))
  stopifnot(all(!is.na(tmp2$min.sum.weekly.deaths)))
  stopifnot(all(tmp2$min.sum.weekly.deaths < tmp2$max.sum.weekly.deaths))
  tmp2 = subset(tmp1, !is.na(sum.weekly.deaths))
  stopifnot(all(tmp2$sum.weekly.deaths >= 0))
  
  tmp = select(tmp, -min_COVID.19.Deaths, -max_COVID.19.Deaths, -min_idx_NA, -max_idx_NA)
  
  if(rm.COVID.19.Deaths)
    tmp = select(tmp, -COVID.19.Deaths)
  
  return(tmp)
}

merge_deathByAge_over_Sex = function(tmp1, tmp2)
{
  
  setnames(tmp1, c('weekly.deaths', 'sum.weekly.deaths', 'min.sum.weekly.deaths', 'max.sum.weekly.deaths'), 
           c('weekly.deaths.sex_1', 'sum.weekly.deaths.sex_1', 'min.sum.weekly.deaths.sex_1', 'max.sum.weekly.deaths.sex_1'))
  setnames(tmp2, c('weekly.deaths', 'sum.weekly.deaths', 'min.sum.weekly.deaths', 'max.sum.weekly.deaths'), 
           c('weekly.deaths.sex_2', 'sum.weekly.deaths.sex_2', 'min.sum.weekly.deaths.sex_2', 'max.sum.weekly.deaths.sex_2'))
  
  # find date index with NA
  tmp3 = tmp1[, list(idx_NA = which(is.na(weekly.deaths.sex_1))), by = c("loc_label", 'age')]
  tmp3 = tmp3[, list(min_idx_NA.sex_1 = min(idx_NA), max_idx_NA.sex_1 = max(idx_NA)), by = c("loc_label", 'age')]
  tmp1 = merge(tmp1, tmp3, by = c("loc_label", 'age'), all.x = T)
  tmp3 = tmp2[, list(idx_NA = which(is.na(weekly.deaths.sex_2))), by = c("loc_label", 'age')]
  tmp3 = tmp3[, list(min_idx_NA.sex_2 = min(idx_NA), max_idx_NA.sex_2 = max(idx_NA)), by = c("loc_label", 'age')]
  tmp2 = merge(tmp2, tmp3, by = c("loc_label", 'age'), all.x = T)
  
  tmp = unique(select(tmp1, loc_label, code, age, age_from, age_to))
  
  tmp1 = merge(tmp1, tmp2, by = c('age', 'date', 'loc_label', 'date_idx', 'min_date_idx', 'max_date_idx'))
  
  # reajust date idx after removal of the first week
  tmp1 = select(tmp1, -date_idx, -min_date_idx, -max_date_idx)
  tmp1[, date_idx := 1:length(date), by = c("loc_label", 'age')]
  tmp3 = tmp1[, list(min_date_idx = min(date_idx), max_date_idx = max(date_idx)), by = c("loc_label", 'age')]
  tmp1 = merge(tmp1, tmp3, by = c("loc_label", 'age'))
  
  tmp1[, weekly.deaths := NA_integer_]
  tmp1[, min.sum.weekly.deaths := NA_integer_]
  tmp1[, max.sum.weekly.deaths := NA_integer_]
  tmp1[, sum.weekly.deaths := NA_integer_]
  
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
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1, 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1, 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 , by = c('age', 'loc_label')]
  
  # both are after the beginning of the period and both finish after the interval
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx & 
         (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         !is.na(sum.weekly.deaths.sex_1) & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  
  # both are at the beginning of the period and both finish after the interval
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         !is.na(sum.weekly.deaths.sex_1) & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  
  # both are at the beginning of the period and both finish before the interval
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2  & (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2  & !is.na(sum.weekly.deaths.sex_1) & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2)+ 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & !is.na(sum.weekly.deaths.sex_1) & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 +  sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & !is.na(sum.weekly.deaths.sex_1) & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  
  # both are at the beginning of the period and one finish after the interval
  # sex 2 finished after the interval
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         !is.na(sum.weekly.deaths.sex_1) & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  
  # sex 1 finished after the interval
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         (is.na(sum.weekly.deaths.sex_1) |is.na(sum.weekly.deaths.sex_2)), 
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         !is.na(sum.weekly.deaths.sex_1) & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  
  # both are after the beginning of the period and one finish after the interval
  # sex 2 finish after
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         is.na(sum.weekly.deaths.sex_2), 
       min.sum.weekly.deaths := sum.weekly.deaths.sex_1 + min.sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         is.na(sum.weekly.deaths.sex_2), 
       max.sum.weekly.deaths := sum.weekly.deaths.sex_1 + max.sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         !is.na(sum.weekly.deaths.sex_1), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  
  
  # sex 1 finish after
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         is.na(sum.weekly.deaths.sex_1), 
       min.sum.weekly.deaths := min.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         is.na(sum.weekly.deaths.sex_1), 
       max.sum.weekly.deaths := max.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         !is.na(sum.weekly.deaths.sex_1), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  
  # one is at the beginning of the period and both finish before the period
  # sex 1 is at the beginning
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & is.na(sum.weekly.deaths.sex_1), 
       min.sum.weekly.deaths := min.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & is.na(sum.weekly.deaths.sex_1), 
       max.sum.weekly.deaths := max.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & !is.na(sum.weekly.deaths.sex_1), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & is.na(sum.weekly.deaths.sex_1), 
       min.sum.weekly.deaths := min.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & is.na(sum.weekly.deaths.sex_1), 
       max.sum.weekly.deaths := max.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & !is.na(sum.weekly.deaths.sex_1), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & is.na(sum.weekly.deaths.sex_1), 
       min.sum.weekly.deaths := min.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & is.na(sum.weekly.deaths.sex_1), 
       max.sum.weekly.deaths := max.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & !is.na(sum.weekly.deaths.sex_1), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  
  # sex 2 is at the beginning
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & is.na(sum.weekly.deaths.sex_2), 
       min.sum.weekly.deaths := min.sum.weekly.deaths.sex_2 + sum.weekly.deaths.sex_1 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & is.na(sum.weekly.deaths.sex_2), 
       max.sum.weekly.deaths := max.sum.weekly.deaths.sex_2 + sum.weekly.deaths.sex_1 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_1 < max_idx_NA.sex_2 & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_2 + sum.weekly.deaths.sex_1 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):unique(max_idx_NA.sex_2)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & is.na(sum.weekly.deaths.sex_2), 
       min.sum.weekly.deaths := min.sum.weekly.deaths.sex_2 + sum.weekly.deaths.sex_1 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & is.na(sum.weekly.deaths.sex_2), 
       max.sum.weekly.deaths := max.sum.weekly.deaths.sex_2 + sum.weekly.deaths.sex_1 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 < max_idx_NA.sex_1 & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_2 + sum.weekly.deaths.sex_1 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):unique(max_idx_NA.sex_1)]), by = c('age', 'loc_label')]
  
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & is.na(sum.weekly.deaths.sex_2), 
       min.sum.weekly.deaths := min.sum.weekly.deaths.sex_2 + sum.weekly.deaths.sex_1 , by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & is.na(sum.weekly.deaths.sex_2), 
       max.sum.weekly.deaths := max.sum.weekly.deaths.sex_2 + sum.weekly.deaths.sex_1, by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 != max_date_idx & 
         max_idx_NA.sex_2 == max_idx_NA.sex_1 & !is.na(sum.weekly.deaths.sex_2), 
       sum.weekly.deaths := sum.weekly.deaths.sex_2 + sum.weekly.deaths.sex_1, by = c('age', 'loc_label')]
  
  # one is at the beginning of the period and both finish after the period
  # sex 1 is at the beginning
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx, 
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1)  + 
         min.sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx, 
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         max.sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  
  # sex 2 is at the beginning
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx , 
       min.sum.weekly.deaths := min.sum.weekly.deaths.sex_1 + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 == max_date_idx, 
       max.sum.weekly.deaths := max.sum.weekly.deaths.sex_1 + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2), by = c('age', 'loc_label')]
  
  # one is at the beginning of the period and one finish after the period
  # sex 1 is at the beginning and sex 1 finish after the period
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         is.na(sum.weekly.deaths.sex_1), 
       min.sum.weekly.deaths := min.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx&
         is.na(sum.weekly.deaths.sex_1), 
       max.sum.weekly.deaths := max.sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx&
         !is.na(sum.weekly.deaths.sex_1), 
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  
  # sex 1 is at the beginning and sex 2 does not finish  before the period
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         (is.na(sum.weekly.deaths.sex_1) | is.na(sum.weekly.deaths.sex_2)),
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         (is.na(sum.weekly.deaths.sex_1) | is.na(sum.weekly.deaths.sex_2)),
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 == min_date_idx & min_idx_NA.sex_2 != min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx&
         !is.na(sum.weekly.deaths.sex_1) & !is.na(sum.weekly.deaths.sex_2),
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  
  # sex 2 is at the beginning and sex 1 does not finish before the period
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         (is.na(sum.weekly.deaths.sex_1) | is.na(sum.weekly.deaths.sex_2)),
       min.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), min.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), min.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         (is.na(sum.weekly.deaths.sex_1) | is.na(sum.weekly.deaths.sex_2)),
       max.sum.weekly.deaths := ifelse(is.na(unique(sum.weekly.deaths.sex_1)), max.sum.weekly.deaths.sex_1, sum.weekly.deaths.sex_1) + 
         ifelse(is.na(unique(sum.weekly.deaths.sex_2)), max.sum.weekly.deaths.sex_2, sum.weekly.deaths.sex_2) + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 == max_date_idx & max_idx_NA.sex_2 != max_date_idx &
         !is.na(sum.weekly.deaths.sex_1) & !is.na(sum.weekly.deaths.sex_2),
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_2[(unique(max_idx_NA.sex_2)+1):(unique(max_idx_NA.sex_1)-1)]), by = c('age', 'loc_label')]
  
  # sex 2 is at the beginning and sex 2 does not finish before the period
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         is.na(sum.weekly.deaths.sex_2),
       min.sum.weekly.deaths := sum.weekly.deaths.sex_1 + min.sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         is.na(sum.weekly.deaths.sex_2),
       max.sum.weekly.deaths := sum.weekly.deaths.sex_1 + max.sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  tmp1[min_idx_NA.sex_1 != min_date_idx & min_idx_NA.sex_2 == min_date_idx & 
         max_idx_NA.sex_1 != max_date_idx & max_idx_NA.sex_2 == max_date_idx &
         !is.na(sum.weekly.deaths.sex_2),
       sum.weekly.deaths := sum.weekly.deaths.sex_1 + sum.weekly.deaths.sex_2 + 
         sum(weekly.deaths.sex_1[(unique(max_idx_NA.sex_1)+1):(unique(max_idx_NA.sex_2)-1)]), by = c('age', 'loc_label')]
  
  # sex 1 doesn't have any NA
  tmp1[is.na(min_idx_NA.sex_1) & !is.na(min_idx_NA.sex_2) & !is.na(sum.weekly.deaths.sex_2), 
       weekly.deaths := sum(weekly.deaths.sex_1[unique(min_idx_NA.sex_2):unique(max_idx_NA.sex_2)]) + 
         sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  tmp1[is.na(min_idx_NA.sex_1) & !is.na(min_idx_NA.sex_2) & is.na(sum.weekly.deaths.sex_2), 
       min.sum.weekly.deaths := sum(weekly.deaths.sex_1[unique(min_idx_NA.sex_2):unique(max_idx_NA.sex_2)]) + 
         min.sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  tmp1[is.na(min_idx_NA.sex_1) & !is.na(min_idx_NA.sex_2) & is.na(sum.weekly.deaths.sex_2), 
       max.sum.weekly.deaths := sum(weekly.deaths.sex_1[unique(min_idx_NA.sex_2):unique(max_idx_NA.sex_2)]) + 
         max.sum.weekly.deaths.sex_2, by = c('age', 'loc_label')]
  
  # sex 2 doesn't have any NA
  tmp1[!is.na(min_idx_NA.sex_1) & is.na(min_idx_NA.sex_2) & !is.na(sum.weekly.deaths.sex_1), 
       weekly.deaths := sum(weekly.deaths.sex_2[unique(min_idx_NA.sex_1):unique(max_idx_NA.sex_1)]) + 
         sum.weekly.deaths.sex_1, by = c('age', 'loc_label')]
  tmp1[!is.na(min_idx_NA.sex_1) & is.na(min_idx_NA.sex_2) & is.na(sum.weekly.deaths.sex_1), 
       min.sum.weekly.deaths := sum(weekly.deaths.sex_2[unique(min_idx_NA.sex_1):unique(max_idx_NA.sex_1)]) + 
         min.sum.weekly.deaths.sex_1, by = c('age', 'loc_label')]
  tmp1[!is.na(min_idx_NA.sex_1) & is.na(min_idx_NA.sex_2) & is.na(sum.weekly.deaths.sex_1), 
       max.sum.weekly.deaths := sum(weekly.deaths.sex_2[unique(min_idx_NA.sex_1):unique(max_idx_NA.sex_1)]) + 
         max.sum.weekly.deaths.sex_1, by = c('age', 'loc_label')]
  
  
  # non missing weekly deaths
  tmp1[weekly.deaths.sex_1 >= 0 & weekly.deaths.sex_2 >= 0, weekly.deaths := weekly.deaths.sex_1 + weekly.deaths.sex_2, by = c('age', 'loc_label')]
  
  # sanity checks 
  tmp2 = subset(tmp1, is.na(weekly.deaths))
  which(!is.na(tmp2$min.sum.weekly.deaths) | !is.na(tmp2$sum.weekly.deaths))
  tmp2[ !(!is.na(tmp2$min.sum.weekly.deaths) | !is.na(tmp2$sum.weekly.deaths)) ]
  stopifnot(all( !is.na(tmp2$min.sum.weekly.deaths) | !is.na(tmp2$sum.weekly.deaths) ))
  stopifnot(all( !is.na(tmp2$max.sum.weekly.deaths) | !is.na(tmp2$sum.weekly.deaths)))
  stopifnot(!any( !is.na(tmp2$min.sum.weekly.deaths) & !is.na(tmp2$sum.weekly.deaths) ) )
  stopifnot(!any( !is.na(tmp2$max.sum.weekly.deaths) & !is.na(tmp2$sum.weekly.deaths)))
  tmp3 = subset(tmp1, !is.na(min.sum.weekly.deaths))
  stopifnot(all(!is.na(tmp3$max.sum.weekly.deaths)))
  stopifnot(all(tmp3$min.sum.weekly.deaths < tmp3$max.sum.weekly.deaths))
  tmp3 = subset(tmp1, !is.na(max.sum.weekly.deaths))
  stopifnot(all(!is.na(tmp3$min.sum.weekly.deaths)))
  stopifnot(all(tmp3$min.sum.weekly.deaths < tmp3$max.sum.weekly.deaths))
  tmp3 = subset(tmp1, !is.na(sum.weekly.deaths))
  stopifnot(all(tmp3$sum.weekly.deaths >= 0))
  
  for(Loc in unique(tmp1$loc_label)){
    for(Age in unique(tmp1$age)){
      tmp3 = subset(tmp1, loc_label == Loc & age == Age)
      .idx_missing = which(is.na(tmp3$weekly.deaths))
      
      if(length(.idx_missing) == 0)
        next
      stopifnot( length(unique(tmp3$min.sum.weekly.deaths[.idx_missing])) == 1)
      stopifnot( length(unique(tmp3$max.sum.weekly.deaths[.idx_missing])) == 1)
      stopifnot( length(unique(tmp3$sum.weekly.deaths[.idx_missing])) == 1)
    }
  }
  
  # final
  tmp1 = select(tmp1, 'date', 'loc_label', 'age', 'min.sum.weekly.deaths', 'max.sum.weekly.deaths', 'sum.weekly.deaths', 'weekly.deaths')
  tmp = merge(tmp, tmp1, by = c('age', 'loc_label'))
  
  return(tmp)
}

incorporate_AllSexes_information = function(tmp1, tmp2)
{
  # check that the total sum exact value is the same btw the Male and Female merged data set and the all sexe
  tmp = incorporate_AllSexes_totalsum_information(tmp1, tmp2)
  
  # check that the total sum boundaries are the same btw the Male and Female merged data set and the all sexe
  tmp = incorporate_AllSexes_boundary_information(tmp)
  
  return(tmp)
}

incorporate_AllSexes_boundary_information = function(tmp)
{
  
  for(Loc in unique(tmp$loc_label)){
    for(Age in unique(tmp$age)){
      
      tmp3 = subset(tmp, loc_label == Loc & age == Age)
      .idx_missing = which(is.na(tmp3$weekly.deaths))
      
      if(length(.idx_missing) == 0)
        next
      
      if(all(!is.na(tmp3$sum.weekly.deaths[.idx_missing])))
        next
      
      stopifnot( length(unique(tmp3$min.sum.weekly.deaths[.idx_missing])) <= 2)
      stopifnot( length(unique(tmp3$max.sum.weekly.deaths[.idx_missing])) <= 2)
      
      stopifnot(tmp3[.idx_missing]$min.sum.weekly.deaths[1] >= tmp3[.idx_missing]$min.sum.weekly.deaths[length(.idx_missing)])
      stopifnot(tmp3[.idx_missing]$max.sum.weekly.deaths[1] >= tmp3[.idx_missing]$max.sum.weekly.deaths[length(.idx_missing)])
      
      # check that the boundaries of weekly death sum are unique
      # it is possible that they are greater in the Male and Female merged data set 
      # so we can use the All sexes to find tighter boundaries
      if(tmp3[.idx_missing]$min.sum.weekly.deaths[1] > tmp3[.idx_missing]$min.sum.weekly.deaths[length(.idx_missing)])
        tmp3[.idx_missing]$min.sum.weekly.deaths = tmp3[.idx_missing]$min.sum.weekly.deaths[length(.idx_missing)]
      
      if(tmp3[.idx_missing]$max.sum.weekly.deaths[1] > tmp3[.idx_missing]$max.sum.weekly.deaths[length(.idx_missing)])
        tmp3[.idx_missing]$max.sum.weekly.deaths =  tmp3[.idx_missing]$max.sum.weekly.deaths[length(.idx_missing)]
      
      if(max(.idx_missing) == nrow(tmp3) & min(.idx_missing) != 1){
        tmp3[.idx_missing]$min.sum.weekly.deaths = 1
        tmp3[.idx_missing]$max.sum.weekly.deaths = 9
      }
      
      tmp = anti_join(tmp, tmp3, by = c('loc_label', 'age'))
      tmp = rbind(tmp, tmp3)
    }
  }
  
  tmp = tmp[order(loc_label, age, date)]
  return(tmp)
}

incorporate_AllSexes_totalsum_information = function(tmp1, tmp2)
{
  
  locations = unique(tmp1$loc_label)
  ages = unique(tmp1$age)
  
  for(Loc in locations){
    for(Age in ages){
      
      # male and female merged data set
      tmp3 = subset(tmp1, loc_label == Loc & age == Age)
      # all sexes dataset
      tmp4 = subset(tmp2, loc_label == Loc & age == Age)
      tmp3 = subset(tmp3, date < min(tmp4$date))
      
      .idx_missing.weekly = which(is.na(tmp3$weekly.deaths))
      .idx_missing.weekly2 = which(is.na(tmp4$weekly.deaths))
      .idx_non_missing.cum = which(!is.na(tmp4$COVID.19.Deaths))
      
      if(length(.idx_missing.weekly) == 0)
        next
      if(length(.idx_non_missing.cum) == 0)
        next
      if(!nrow(tmp3) %in% .idx_missing.weekly)
        next
      
      # verify that sum weekly death is equal to the first value of COVID.19.Deaths
      if(!is.na(tmp3[nrow(tmp3)]$sum.weekly.deaths) & !is.na(tmp4[1]$COVID.19.Deaths))
      {
        if(is.na(tmp4[1]$weekly.deaths)){
          tmp4[1]$weekly.deaths = tmp3[nrow(tmp3)]$sum.weekly.deaths
          tmp2 = anti_join(tmp2, tmp4, by = c('loc_label', 'age'))
          tmp2 = rbind(tmp2, tmp4)
          
        }
        if(tmp3[nrow(tmp3)]$sum.weekly.deaths == tmp4[1]$COVID.19.Deaths)
          next
        if(tmp3[nrow(tmp3)]$sum.weekly.deaths != tmp4[1]$COVID.19.Deaths){
          tmp3[.idx_missing.weekly]$sum.weekly.deaths = tmp4[1]$COVID.19.Deaths
          tmp1 = anti_join(tmp1, tmp3, by = c('loc_label', 'age'))
          tmp1 = rbind(tmp1, tmp3)
          
        }
      }

      # fill sum weekly death with first value of COVID.19.Deaths if NA doesn't start at the beginning of the period
      first_non_missing_cum = tmp4$COVID.19.Deaths[.idx_non_missing.cum[1]]
      if(1 %in% .idx_missing.weekly & is.na(tmp3[nrow(tmp3)]$sum.weekly.deaths))
      {
        stopifnot(all(is.na(tmp3[.idx_missing.weekly]$sum.weekly.deaths )))
        tmp3[.idx_missing.weekly]$max.sum.weekly.deaths = first_non_missing_cum
        
        .idx_missing.weekly2 = c(.idx_missing.weekly2, .idx_non_missing.cum[1])
        stopifnot(all(is.na(tmp4[.idx_missing.weekly2]$sum.weekly.deaths )))
        tmp4[.idx_missing.weekly2]$min.sum.weekly.deaths = unique(tmp3[.idx_missing.weekly]$min.sum.weekly.deaths)
        tmp4[.idx_missing.weekly2]$max.sum.weekly.deaths = first_non_missing_cum
        tmp4[.idx_non_missing.cum[1]]$weekly.deaths = NA_integer_
      } else {
        
        tmp3[.idx_missing.weekly]$sum.weekly.deaths = first_non_missing_cum
        tmp3[.idx_missing.weekly]$min.sum.weekly.deaths = NA_integer_
        tmp3[.idx_missing.weekly]$max.sum.weekly.deaths = NA_integer_
        .idx_missing.weekly2 = c(.idx_missing.weekly2, .idx_non_missing.cum[1])
        tmp4[.idx_missing.weekly2]$sum.weekly.deaths = first_non_missing_cum
        tmp4[.idx_missing.weekly2]$min.sum.weekly.deaths = NA_integer_
        tmp4[.idx_missing.weekly2]$max.sum.weekly.deaths = NA_integer_
        tmp4[.idx_non_missing.cum[1]]$weekly.deaths = NA_integer_
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

