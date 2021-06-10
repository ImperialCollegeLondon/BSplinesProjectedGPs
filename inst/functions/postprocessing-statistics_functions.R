statistics_contribution_all_states = function(contribution75){
  
  tmp75  = subset(contribution75, !code %in% c('HI', 'VT', 'AK'))

  # longest time
  date_d = tmp75[M_rel_adj > 1 & date > "2020-12-01" , list(min_date = max(date)), by = 'loc_label']
  date_d = merge(date_d, tmp75, by = 'loc_label')
  date_d = date_d[M_rel_adj < 1 & date >= min_date, list(min_date = min(date)), by = 'loc_label']
  datemin = date_d[min_date == min(min_date),list(loc_label = loc_label, date = format(min_date,  '%B %d, %Y'))]
  datemax = date_d[min_date == max(min_date),list(loc_label = loc_label, date = format(min_date,  '%B %d, %Y'))]
  
  contribution_stats = list(date_d = date_d, datemin = datemin, datemax = datemax)
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
  
  tmp = contribution_ref_adj[, list(M = mean(M), CL = mean(CL), CU = mean(CU)), by = c('age', 'division')]
  
  tmp = tmp[order(M, decreasing = T)]
  divisionmax80 = tmp[age == '85+',division][1:2]
  divisionmax5574 = tmp[age == '55-74',division][1:2]
  
  tmp = tmp[order(M, decreasing = F)]
  divisionmin80 = tmp[age == '85+',division][1:2]
  divisionmin5574 = tmp[age == '55-74',division][1:2]

  
  contribution_baseline = list(statemax85 = statemax80, statemin85 = statemin80, 
                               statemax5574 = statemax5574, statemin5574 = statemin5574, 
                               average85 = average80, average7584 = average7584, average5574 = average5574, 
                               average2554 = average2554, average024 = average024)
  
  return(contribution_baseline)
}

