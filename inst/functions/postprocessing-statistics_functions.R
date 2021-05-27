statistics_contribution_all_states = function(contribution064, contribution6574, contribution75, vaccinedata){
  
  tmp75  = subset(contribution75, !code %in% c('HI', 'VT', 'AK'))
  tmp6574 = subset(contribution6574, !code %in% c('HI', 'VT', 'AK'))
  tmp064  = subset(contribution064, !code %in% c('HI', 'VT', 'AK'))
  
  date_red = tmp75[date > as.Date('2020-09-01'), list(MM = mean(M_rel_adj)), by = 'date']
  date_red = date_red[MM < 1, min(date)]
  red_baseline = tmp75[date==date_red,list(M_adj = paste0(round(mean(M_adj)*100, 2), '\\%'), 
                                           CL_adj = paste0(round(mean(CL_adj)*100, 2), '\\%'),
                                           CU_adj = paste0(round(mean(CU_adj)*100, 2), '\\%'), 
                                           M_rel_adj = paste0(round((1 - mean(M_rel_adj))*100, 2), '\\%'), 
                                           CL_rel_adj = paste0(round((1 - mean(CU_rel_adj))*100, 2), '\\%'),
                                           CU_rel_adj = paste0(round((1 - mean(CL_rel_adj))*100, 2), '\\%'))]
  
  # longest time
  date_d = tmp75[M_rel_adj > 1 & date > "2020-07-18" , list(min_date = max(date)), by = 'loc_label']
  date_d = merge(date_d, tmp75, by = 'loc_label')
  date_d = date_d[M_rel_adj < 1 & date >= min_date, list(min_date = min(date)), by = 'loc_label']
  datemin = date_d[min_date == min(min_date),list(loc_label = loc_label, date = format(min_date,  '%d %B, %Y'))]
  datemax = date_d[min_date == max(min_date),list(loc_label = loc_label, date = format(min_date,  '%d %B, %Y'))]
  
  # fastest decrease
  speed = merge(date_d, tmp75, by = 'loc_label')
  speed[, dummy := (min_date + 8*7) %in% date , by = 'loc_label']
  speed = speed[dummy == 1 & (date == min_date + 8*7|date == min_date)]
  speed[, dummy := ifelse(date == min_date, 'date_1', 'date_2')]
  speed = as.data.table(reshape2::dcast(speed, loc_label ~ dummy, value.var = 'M_rel_adj'))
  speed[, ratio := date_2 / date_1]
  state_fast=speed[order(ratio), loc_label][1:3]
  state_slow=speed[order(ratio, decreasing = T),loc_label][1:3]
  
  # greatest decrease
  dec = merge(date_d, tmp75[date == max(date)], by = 'loc_label')
  decmax = dec[order(M_rel_adj), list(M_rel_adj = paste0(round((1-(M_rel_adj))*100, 2), '\\%'), 
                                      CL_rel_adj = paste0(round((1-(CU_rel_adj))*100, 2), '\\%'),
                                      CU_rel_adj = paste0(round((1-(CL_rel_adj))*100, 2), '\\%'),
                                      M_adj = paste0(round(M_adj*100, 2), '\\%'), 
                                      CL_adj = paste0(round(CL_adj*100, 2), '\\%'),
                                      CU_adj = paste0(round(CU_adj*100, 2), '\\%')), by = 'loc_label'][1:3]
  decmin = dec[order(M_rel_adj, decreasing = T), list(M_rel_adj = paste0(round((1-(M_rel_adj))*100, 2), '\\%'), 
                                                      CL_rel_adj = paste0(round((1-(CU_rel_adj))*100, 2), '\\%'),
                                                      CU_rel_adj = paste0(round((1-(CL_rel_adj))*100, 2), '\\%'),
                                                      M_adj = paste0(round(M_adj*100, 2), '\\%'), 
                                                      CL_adj = paste0(round(CL_adj*100, 2), '\\%'),
                                                      CU_adj = paste0(round(CU_adj*100, 2), '\\%')), by = 'loc_label'][1:3]
  
  # compensate
  comp6574 = merge(tmp6574, date_d, by = 'loc_label')
  comp6574 = comp6574[date >= min_date, list(MM_6574 = mean(M_rel_adj)), by = 'loc_label']
  comp064 = merge(tmp064, date_d, by = 'loc_label')
  comp064= comp064[date >= min_date, list(MM_064 = mean(M_rel_adj)), by = 'loc_label']
  comp = merge(comp6574, comp064, by = 'loc_label')
  comp[, ratio := MM_6574 /  MM_064]
  stat_comp = comp[ratio > 1, loc_label]
  stat_comp = paste0(paste0(stat_comp[-length(stat_comp)], collapse = ', '), ' and ', stat_comp[length(stat_comp)])
  
  # vaccine
  vac = vaccinedata[prop_vaccinated_fully >= 0.25, list(date = min(date), prop = 0.25), by = 'age']
  tmp1 = vaccinedata[prop_vaccinated_fully >= 0.5, list(date = min(date), prop = 0.5), by = 'age']
  vac = rbind(vac, tmp1)
  vac[, date := as.Date(date)]
  
  date75025 = vac[prop == 0.25 & age == '75+',date]
  tmp75025 = tmp75[date >= date75025 & date <= date75025 + 6, ]
  vac75025av = tmp75025[, list(M_adj = paste0(round(mean(M_adj)*100, 2), '\\%'), 
                               CL_adj = paste0(round(mean(CL_adj)*100, 2), '\\%'),
                               CU_adj = paste0(round(mean(CU_adj)*100, 2), '\\%'), 
                               M_rel_adj = paste0(round((1 - mean(M_rel_adj))*100, 2), '\\%'), 
                               CL_rel_adj = paste0(round((1 - mean(CU_rel_adj))*100, 2), '\\%'),
                               CU_rel_adj = paste0(round((1 - mean(CL_rel_adj))*100, 2), '\\%'))]
  vac75025max = tmp75025[order(M_rel_adj), list(M_adj = paste0(round((M_adj)*100, 2), '\\%'), 
                                                CL_adj = paste0(round((CL_adj)*100, 2), '\\%'),
                                                CU_adj = paste0(round((CU_adj)*100, 2), '\\%'), 
                                                M_rel_adj = paste0(round((1 - (M_rel_adj))*100, 2), '\\%'), 
                                                CL_rel_adj = paste0(round((1 - (CL_rel_adj))*100, 2), '\\%'),
                                                CU_rel_adj = paste0(round((1 - (CU_rel_adj))*100, 2), '\\%'),
                                                loc_label = loc_label)][1:3,]
  vac75025min = tmp75025[order(M_rel_adj,decreasing = T), list(M_adj = paste0(round((M_adj)*100, 2), '\\%'), 
                                                               CL_adj = paste0(round((CL_adj)*100, 2), '\\%'),
                                                               CU_adj = paste0(round((CU_adj)*100, 2), '\\%'), 
                                                               M_rel_adj = paste0(round((1 - (M_rel_adj))*100, 2), '\\%'), 
                                                               CL_rel_adj = paste0(round((1 - (CL_rel_adj))*100, 2), '\\%'),
                                                               CU_rel_adj = paste0(round((1 - (CU_rel_adj))*100, 2), '\\%'),
                                                               loc_label = loc_label)][1:3,]
  
  date7505 = vac[prop == 0.5 & age == '75+',date]
  tmp7505 = tmp75[date >= date7505 & date <= date7505 + 6, ]
  vac7505av= tmp7505[, list(M_adj = paste0(round(mean(M_adj)*100, 2), '\\%'), 
                            CL_adj = paste0(round(mean(CL_adj)*100, 2), '\\%'),
                            CU_adj = paste0(round(mean(CU_adj)*100, 2), '\\%'), 
                            M_rel_adj = paste0(round((1 - mean(M_rel_adj))*100, 2), '\\%'), 
                            CL_rel_adj = paste0(round((1 - mean(CU_rel_adj))*100, 2), '\\%'),
                            CU_rel_adj = paste0(round((1 - mean(CL_rel_adj))*100, 2), '\\%'))]
  vac7505max = tmp7505[order(M_rel_adj), list(M_adj = paste0(round((M_adj)*100, 2), '\\%'), 
                                              CL_adj = paste0(round((CL_adj)*100, 2), '\\%'),
                                              CU_adj = paste0(round((CU_adj)*100, 2), '\\%'), 
                                              M_rel_adj = paste0(round((1 - (M_rel_adj))*100, 2), '\\%'), 
                                              CL_rel_adj = paste0(round((1 - (CU_rel_adj))*100, 2), '\\%'),
                                              CU_rel_adj = paste0(round((1 - (CL_rel_adj))*100, 2), '\\%'),
                                              loc_label = loc_label)][1:3,]
  vac7505min = tmp7505[order(M_rel_adj,decreasing = T), list(M_adj = paste0(round((M_adj)*100, 2), '\\%'), 
                                                             CL_adj = paste0(round((CL_adj)*100, 2), '\\%'),
                                                             CU_adj = paste0(round((CU_adj)*100, 2), '\\%'), 
                                                             M_rel_adj = paste0(round((1 - (M_rel_adj))*100, 2), '\\%'), 
                                                             CL_rel_adj = paste0(round((1 - (CU_rel_adj))*100, 2), '\\%'),
                                                             CU_rel_adj = paste0(round((1 - (CL_rel_adj))*100, 2), '\\%'),
                                                             loc_label = loc_label)][1:3,]
  
  contribution_stats = list(baseline_date =format(date_red, '%d %B, %Y'), end_period = format(max(tmp75$date), '%d %B, %Y'),
                            baseline_stat = red_baseline,
                            decreasemin = datemin, decreasemax = datemax, 
                            fastestdecrease = state_fast, slowestdecrease = state_slow, 
                            maxdecrease = decmax, mindecrease = decmin,
                            compensation = stat_comp,
                            list(date75025 = format(date75025,  '%d %B, %Y'), vac75025av = vac75025av, vac75025max = vac75025max, vac75025min = vac75025min),
                            list(date7505 = format(date7505,  '%d %B, %Y'), vac7505av = vac7505av, vac7505max = vac7505max,vac7505min = vac7505min))
  
  return(contribution_stats)
}



statistics_contributionref_all_states = function(contribution_ref_adj){
  tmp = copy( subset(contribution_ref_adj, !code %in% c('HI', 'VT', 'AK')))
  tmp = tmp[order(M, decreasing = T)]
  statemax80 = tmp[age == '85+', list(loc_label = loc_label, 
                                      paste0(round(M*100, 2), '\\%'), 
                                      paste0(round(mean(CL)*100, 2), '\\%'),
                                      paste0(round(mean(CU)*100, 2), '\\%'))][1,]
  statemax5574 = tmp[age == '55-74', list(loc_label = loc_label, 
                                          paste0(round(M*100, 2), '\\%'), 
                                          paste0(round(mean(CL)*100, 2), '\\%'),
                                          paste0(round(mean(CU)*100, 2), '\\%'))][1,]
  tmp = tmp[order(M, decreasing = F)]
  statemin80 = tmp[age == '85+',list(loc_label = loc_label, 
                                     paste0(round(M*100, 2), '\\%'), 
                                     paste0(round(mean(CL)*100, 2), '\\%'),
                                     paste0(round(mean(CU)*100, 2), '\\%'))][1]
  statemin5574 = tmp[age == '55-74', list(loc_label = loc_label, 
                                          paste0(round(M*100, 2), '\\%'), 
                                          paste0(round(mean(CL)*100, 2), '\\%'),
                                          paste0(round(mean(CU)*100, 2), '\\%'))][1,]
  
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
  
  return(contribution_baseline)
}

