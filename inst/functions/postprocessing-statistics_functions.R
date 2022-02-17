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

find_statistics_mortality_rate <- function(mortality_rate, outdir){
  max_date = format(unique(mortality_rate$date), '%B %d, %Y')
  dold2p = mortality_rate[age == '85+' & M > 0.025, loc_label]
  dold2p_n = paste0(paste0(dold2p[-length(dold2p)], collapse = ', '), ' and ', dold2p[length(dold2p)])
  state_max = mortality_rate[age == '85+',]
  state_max = state_max[order(M, decreasing = T)]
  state_max = state_max[,loc_label][1]
  # state_max = paste0(paste0(state_max[-length(state_max)], collapse = ', '), ' and ', state_max[length(state_max)])
  
  d5574 = mortality_rate[age == '55-74',paste0(round((median(M)*100),2),'\\%')]
  d7584 = mortality_rate[age == '75-84',paste0(round((median(M)*100),0),'\\%')]
  
  mortality_stats = list(max_date = max_date, nstates2p = length(dold2p), dold2p_n,
                         state_max, d5574, d7584)
  saveRDS(mortality_stats, file = paste0(outdir, '-mortality_stats.rds'))
}
