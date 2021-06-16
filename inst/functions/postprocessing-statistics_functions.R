

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


find_regime_state = function(contribution75, vaccinedata_state, outdir){
  
  locs = unique(contribution75$code)
  
  date_break= as.Date('2021-05-01')
  date_before_vaccine = as.Date('2020-10-31')
  
  contribution_stats = list()
  contribution_stats[['date_before']] = format(date_before_vaccine,  '%B %d, %Y')
  contribution_stats[['date_break']] = format(date_break,  '%B %d, %Y')
  
  # find states slow, fast, plateau
  beta = vector(mode = 'list', length = length(locs))
  for(i in seq_along(locs)){
    y = subset(contribution75, code == locs[i] & date >= date_before_vaccine & date < date_break)$M
    x = 1:length(y)
    fit1 = lm(y ~ x)
    
    y = subset(contribution75, code == locs[i] & date >= date_break)$M
    x = 1:length(y)
    fit2 = lm(y ~ x)
    
    beta[[i]] = data.table(code = locs[i], loc_label = region_name[code == locs[i], loc_label],
                           betatot = fit1$coefficients[2],  betalast = fit2$coefficients[2])
  }
  beta = do.call('rbind', beta)
  beta =  subset(beta, !code %in% c('HI', 'VT', 'AK'))
  
  plateau = beta[betalast > 0, loc_label]
  beta = beta[order(betatot, decreasing = T)]
  slowd = beta[betatot > summary(beta$betatot)[5],  loc_label][1:5]
  beta = beta[order(betatot)]
  fastd = beta[betatot < summary(beta$betatot)[2],  loc_label][1:5]
  
  contribution_stats[['plateau']] = paste0(paste0(plateau[-length(plateau)], collapse = ', '), ' and ', plateau[length(plateau)])
  contribution_stats[['slowd']] = paste0(paste0(slowd[-length(slowd)], collapse = ', '), ' and ', slowd[length(slowd)])
  contribution_stats[['fastd']] = paste0(paste0(fastd[-length(fastd)], collapse = ', '), ' and ', fastd[length(fastd)])
  
  con_bv = contribution75[date == date_before_vaccine, list(paste0(round(mean(M)*100, 2), '\\%'), 
                                                   paste0(round(mean(CL)*100, 2), '\\%'),
                                                   paste0(round(mean(CU)*100, 2), '\\%'))]
  
  con_a = contribution75[date == date_break, list(paste0(round(mean(M)*100, 2), '\\%'), 
                                                   paste0(round(mean(CL)*100, 2), '\\%'),
                                                   paste0(round(mean(CU)*100, 2), '\\%'))]
  contribution_stats[['con_bv']] = con_bv
  contribution_stats[['con_a']] = con_a
  
  tmp = merge(vaccinedata_state, contribution75, by = c('date', 'loc_label'))
  tmp = tmp[date <= date_break]
  fit = lm(M ~ prop_vaccinated_1dosep, data = tmp)
  coefficients = fit$coefficients
  
  contribution_stats[['date_vacs']] = format(min(tmp$date),  '%B %d, %Y')
  contribution_stats[['beta']] = -round(coefficients, digits =2)
  
  saveRDS(contribution_stats, file = paste0(outdir, '-contribution_vaccination_stats.rds'))
  
  return(contribution_stats)
}
