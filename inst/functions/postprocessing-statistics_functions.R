statistics_contributionref_all_states = function(contribution_ref_adj, outdir){
  

  tmp = copy( subset(contribution_ref_adj, !code %in% c('HI', 'VT', 'AK')))
  tmp = tmp[order(M, decreasing = T)]
  statemax80 = tmp[age == '85+', list(loc_label = loc_label, 
                                      paste0(round(M*100, 2), '\\%'), 
                                      paste0(round(CL*100, 2), '\\%'),
                                      paste0(round(CU*100, 2), '\\%'))][1,]
  statemax5574 = tmp[age == '55-74', list(loc_label = loc_label, 
                                          paste0(round(M*100, 2), '\\%'), 
                                          paste0(round(CL*100, 2), '\\%'),
                                          paste0(round(CU*100, 2), '\\%'))][1,]
  tmp = tmp[order(M, decreasing = F)]
  statemin80 = tmp[age == '85+',list(loc_label = loc_label, 
                                     paste0(round(M*100, 2), '\\%'), 
                                     paste0(round(CL*100, 2), '\\%'),
                                     paste0(round(CU*100, 2), '\\%'))][1]
  statemin5574 = tmp[age == '55-74', list(loc_label = loc_label, 
                                          paste0(round(M*100, 2), '\\%'), 
                                          paste0(round(CL*100, 2), '\\%'),
                                          paste0(round(CU*100, 2), '\\%'))][1,]
  
  average80 = tmp[age == '85+',list(paste0(round(mean(M)*100, 2), '\\%'), 
                                    paste0(round(mean(CL)*100, 2), '\\%'),
                                    paste0(round(mean(CU)*100, 2), '\\%'))]
  average5574 = tmp[age == '55-74',list(paste0(round(mean(M)*100, 2), '\\%'), 
                                        paste0(round(mean(CL)*100, 2), '\\%'),
                                        paste0(round(mean(CU)*100, 2), '\\%'))]
  average7584 = tmp[age == '75-84',list(paste0(round(mean(M)*100, 2), '\\%'), 
                                        paste0(round(mean(CL)*100, 2), '\\%'),
                                        paste0(round(mean(CU)*100, 2), '\\%'))]
  average024 = tmp[age == '0-24',list(paste0(round(mean(M)*100, 2), '\\%'), 
                                      paste0(round(mean(CL)*100, 2), '\\%'),
                                      paste0(round(mean(CU)*100, 2), '\\%'))]
  average2554 = tmp[age == '25-54',list(paste0(round(mean(M)*100, 2), '\\%'), 
                                        paste0(round(mean(CL)*100, 2), '\\%'),
                                        paste0(round(mean(CU)*100, 2), '\\%'))]


  
  contribution_baseline = list(statemax85 = statemax80, statemin85 = statemin80, 
                               statemax5574 = statemax5574, statemin5574 = statemin5574, 
                               average85 = average80, average7584 = average7584, average5574 = average5574, 
                               average2554 = average2554, average024 = average024)
  
  saveRDS(contribution_baseline, file = paste0(outdir, '-contribution_ref_adj_stats.rds'))
  
  # return(contribution_baseline)
}

find_statistics_mortality_rate <- function(mortality_rate, mortality_rate_across_states, outdir){
  
  mortality_stats = list()
  
  max_date = format(unique(mortality_rate$date), '%B %d, %Y')
  mortality_stats[[1]] = max_date
  
  prop.85 <- 0.04
  tmp = mortality_rate[age == '85+' & M > prop.85, .(loc_label, M)]
  tmp = tmp[order(M, decreasing = T)]
  tmp = tmp[,as.character(loc_label)]
  if(length(tmp) > 0){
    tmp <- paste0(paste0(tmp[1:(length(tmp) - 1)], collapse = ', '), ' and ', tmp[length(tmp)])
    mortality_stats[[2]] = list(prop.85 * 100, tmp)
  } 

  prop.75 <- 0.015
  tmp = mortality_rate[age == '75-84' & M > prop.75, .(loc_label, M)]
  tmp = tmp[order(M, decreasing = T)]
  tmp = tmp[,as.character(loc_label)]
  if(length(tmp) > 0){
    tmp <- paste0(paste0(tmp[1:(length(tmp) - 1)], collapse = ', '), ' and ', tmp[length(tmp)])
    mortality_stats[[3]] = list(prop.75 * 100, tmp)
  }

  prop.55 <- 0.005
  tmp = mortality_rate[age == '55-74' & M > prop.55, .(loc_label, M)]
  tmp = tmp[order(M, decreasing = T)]
  tmp = tmp[,as.character(loc_label)]
  if(length(tmp) > 0){
    tmp <- paste0(paste0(tmp[1:(length(tmp) - 1)], collapse = ', '), ' and ', tmp[length(tmp)])
    mortality_stats[[4]] = list(prop.55 * 100, tmp)
  }
  
  prop.25 <- 0.001
  tmp = mortality_rate[age == '25-54' & M > prop.25, .(loc_label, M)]
  tmp = tmp[order(M, decreasing = T)]
  tmp = tmp[,as.character(loc_label)]
  if(length(tmp) > 0){
    tmp <- paste0(paste0(tmp[1:(length(tmp) - 1)], collapse = ', '), ' and ', tmp[length(tmp)])
    mortality_stats[[5]] = list(prop.25 * 100, tmp)
  }
  
  mortality_rate_across_states[, M := round(M * 100, digits = 3)]
  mortality_rate_across_states[, CL := round(CL * 100, digits = 3)]
  mortality_rate_across_states[, CU := round(CU * 100, digits = 3)]
  
  mortality_stats[[6]] = mortality_rate_across_states
  
  saveRDS(mortality_stats, file = paste0(outdir, '-mortality_stats.rds'))
  
  return(mortality_stats)
}

save_p_value_vaccine_effects <- function(samples, names, outdir){
  
  names_fit <- names(samples)
  names <- names_fit[ grepl(paste(paste0('^',names),collapse = '|'),names_fit) ]
  
  tmp <- lapply(names, function(x) find_p_value_vaccine_effect(samples[[x]], x))
  tmp <- do.call('rbind', tmp)
  
  saveRDS(tmp, paste0(outdir.table, '-p_value_vaccine_effects.rds'))
  
}

save_mortality_rate_correlation_longtermdeaths <- function(mortality_rateJ21, nyt_data, region_name, outdir.table){
 
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # find relative value to 85+
  mortality_rateJ21[, value_maxage := value[age_index == max(age_index)], by = c('iterations', 'loc_label', 'week_index')]
  mortality_rateJ21[, value_rel := value_maxage / value]
  
  # find correlation coefficients 
  tmp3 <- merge(mortality_rateJ21, nyt_data, by.y = 'STATE', by.x = 'loc_label')
  tmp3 <- tmp3[!is.na(value_rel) & !is.na(SHARE_DEATHS)]
  tmp3[, value_corr := cor(value_rel, SHARE_DEATHS), by = c('week_index', 'age_index', 'iterations')]
  
  # quantile
  tmp2 = tmp3[, list(q= quantile(value_corr, prob=ps, na.rm = T), q_label=p_labs), by=c('week_index', 'age_index')]	
  tmp2 = dcast(tmp2, week_index + age_index ~ q_label, value.var = "q")
  
  saveRDS(tmp2, paste0(outdir.table, '-mortality_rateJ21.rds'))
}

save_resurgence_dates <- function(resurgence_dates, outdir){
  resurgence_dates[,start_resurgence_name := format(start_resurgence, '%d %b, %Y')]
  resurgence_dates[,stop_resurgence_name := format(stop_resurgence, '%d %b, %Y')]
  
  tmp <- merge(resurgence_dates, df_state, by = 'code')[, .(loc_label, start_resurgence_name, stop_resurgence_name)]
  
  saveRDS(tmp, paste0(outdir, '-resurgence_dates.rds'))
}

find_stats_vaccine_effects <- function(data_res1, data_res2, data_res3, data_res4,
                                       data_res5, data_res6, data_res7, data_res8, prop_vac, resurgence_dates, outdir){
  
  data_res1 = merge(data_res1, resurgence_dates, by = 'code')
  prop_vac = merge(prop_vac, resurgence_dates, by = 'code')
  data_res2 = merge(data_res2, resurgence_dates, by = 'code')
  data_res5 = merge(data_res5, resurgence_dates, by = 'code')
  data_res6 = merge(data_res6, resurgence_dates, by = 'code')
  
  stat = list(format(c(min(resurgence_dates$start_resurgence), max(resurgence_dates$stop_resurgence)),  '%B %d, %Y'),
              subset(data_res1, date == stop_resurgence)[, list(M = round(M), 
                                                                CL = round(CL), 
                                                                CU = round(CU)), by = c('counterfactual_index', 'age', 'loc_label')],
              subset(prop_vac, date == start_resurgence)[, list(min_3 = paste0(round(min(prop_1*100), 2), '\\%'),
                                                                max_3 = paste0(round(max(prop_1*100), 2), '\\%'),
                                                                min_4 = paste0(round(min(prop_2*100), 2), '\\%'),
                                                                max_4 = paste0(round(max(prop_2*100), 2), '\\%'))],
              subset(data_res2, date == stop_resurgence)[, list(M = format(round(M*100, digits = 2), nsmall = 2), 
                                                                CL = format(round(CL*100, digits = 2), nsmall = 2), 
                                                                CU = format(round(CU*100, digits = 2), nsmall = 2)), by = c('counterfactual_index', 'age', 'loc_label')],
              subset(data_res3, week_index == max(week_index))[, list(M = round(M), 
                                                                      CL = round(CL), 
                                                                      CU = round(CU)), by = c('counterfactual_index', 'age')],
              subset(data_res4, week_index == max(week_index))[, list(M = format(round(M*100, digits = 2), nsmall = 2), 
                                                                      CL = format(round(CL*100, digits = 2), nsmall = 2), 
                                                                      CU = format(round(CU*100, digits = 2), nsmall = 2)), by = c('counterfactual_index', 'age')],
              subset(data_res5, date == stop_resurgence)[, list(M = round(M), 
                                                                CL = round(CL), 
                                                                CU = round(CU)), by = c('counterfactual_index', 'loc_label')],
              subset(data_res6, date == stop_resurgence)[, list(M = format(round((M)*100, digits = 2), nsmall = 2), 
                                                                CL = format(round((CL)*100, digits = 2), nsmall = 2), 
                                                                CU = format(round((CU)*100, digits = 2), nsmall = 2)), by = c('counterfactual_index', 'loc_label')],
              subset(data_res7, week_index == max(week_index))[, list(M = round(M), 
                                                                      CL = round(CL), 
                                                                      CU = round(CU)), by = c('counterfactual_index')],
              subset(data_res8, week_index == max(week_index))[, list(M = format(round((M)*100, digits = 2), nsmall = 2), 
                                                                      CL = format(round((CL)*100, digits = 2), nsmall = 2), 
                                                                      CU = format(round((CU)*100, digits = 2), nsmall = 2)), by = c('counterfactual_index')]
              
              
  )
  saveRDS(stat, file = paste0(outdir, paste0('-Mortality_counterfactual.rds')))
  
  return(stat)
}

save_statistics_contributiondiff <- function(contributiondiff, outdir.table, lab = NULL){
  l = list()
  tmp <- contributiondiff[, list(code_min = code[which.min(M)], 
                                 code_max = code[which.max(M)]), by = c('variable', 'age')]
  contributiondiff <- merge(contributiondiff, tmp, by = c('variable', 'age'))
  contributiondiff[, `:=`(M = round(M, 4), CL = round(CL, 4), CU = round(CU, 4))]
  l[['code_min_1']] = contributiondiff[age == '65+' & code == code_min & variable == 'diff1', .(loc_label, M, CL, CU)]
  l[['code_max_1']] = contributiondiff[age == '65+' & code == code_max & variable == 'diff1', .(loc_label, M, CL, CU)]
  l[['code_US_1']] = contributiondiff[age == '65+' & code == 'US' & variable == 'diff1', .(loc_label, M, CL, CU)]
  l[['code_min_2']] = contributiondiff[age == '65+' & code == code_min & variable == 'diff2', .(loc_label, M, CL, CU)]
  l[['code_max_2']] = contributiondiff[age == '65+' & code == code_max & variable == 'diff2', .(loc_label, M, CL, CU)]
  l[['code_US_2']] = contributiondiff[age == '65+' & code == 'US' & variable == 'diff2', .(loc_label, M, CL, CU)]
  
  saveRDS(l, file = paste0(outdir.table, '-contributiondiff', lab, '.rds'))
  
}


