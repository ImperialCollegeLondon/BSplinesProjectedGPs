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
  
  mortality_rate_across_states[, M := round(M * 100, digits = 2)]
  mortality_rate_across_states[, CL := round(CL * 100, digits = 2)]
  mortality_rate_across_states[, CU := round(CU * 100, digits = 2)]
  
  mortality_stats[[6]] = mortality_rate_across_states
  
  saveRDS(mortality_stats, file = paste0(outdir, '-mortality_stats.rds'))
  
  return(mortality_stats)
}
