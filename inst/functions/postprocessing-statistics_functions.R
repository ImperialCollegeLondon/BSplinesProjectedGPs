
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

find_regime_state = function(contribution75, vaccine_data, resurgence_dates, date_before_vaccine, outdir){
  
  delay = 2*7
  
  contribution_stats = list()
  
  contribution_stats[['date_before']] = format(date_before_vaccine,  '%B %d, %Y')
  contribution_stats[['start_resurgence']] = format(min(resurgence_dates$start_resurgence),  '%B %d, %Y')
  
  contribution75 = merge(contribution75, resurgence_dates, by = 'code')
  
  con_bv = contribution75[date == date_before_vaccine, list(paste0(round(mean(M)*100, 2), '\\%'), 
                                                   paste0(round(mean(CL)*100, 2), '\\%'),
                                                   paste0(round(mean(CU)*100, 2), '\\%')), by = 'age']
  
  con_a = contribution75[date == start_resurgence-delay, list(paste0(round(mean(M)*100, 2), '\\%'), 
                                                   paste0(round(mean(CL)*100, 2), '\\%'),
                                                   paste0(round(mean(CU)*100, 2), '\\%')), by = 'age']
  contribution_stats[['con_bv']] = con_bv
  contribution_stats[['con_a']] = con_a
  
  
  # find vaccine effect
  vaccine_data_sum = vaccine_data[age >= 18, list(prop = mean(prop)), by = c('code', 'date')]
  vaccine_data_sum[, date := date + delay]
  contribution75 = merge(contribution75, vaccine_data_sum, by = c('code', 'date'))
  
  beta = contribution75[, {
    M100 = M * 100
    prop100 = prop* 100
    fit1 <- lm(M100 ~ prop100)
    list(betap = fit1$coefficients[2])}, by = c('code', 'loc_label', 'age')]
  
  beta1 = beta[, list(beta_M = median(betap), 
                     beta_CL = quantile(betap, probs = 0.025), 
                     beta_CU = quantile(betap, probs = 0.975)), by = 'age']
  
  beta2 = beta[, list(min_state = loc_label[betap == min(betap)],
                      max_state = loc_label[betap == max(betap)]), by = 'age']
  
  contribution_stats[['beta']] = beta1[age == '65+', list(beta_M = -round(beta_M*10,2),
                                                         beta_CL = -round(beta_CU*10,2),
                                                         beta_CU = -round(beta_CL*10, 2))]
  contribution_stats[['state_r']] = beta2
  
  saveRDS(contribution_stats, file = paste0(outdir, '-contribution_vaccination_stats.rds'))
  
  return(contribution_stats)
}

find_statistics_weekly_deaths = function(death, propdeath, deathpost, deathpost2, vaccinedata_state, outdir)
{
  
  locs = unique(propdeath$code)
  
  date_before_vaccine = min(vaccinedata_state$date)
  date_end = max(death$date)
  
  death_75 = subset(deathpost, age== '75+')
  death_074 = subset(deathpost, age== '0-74')
  death_054 = subset(deathpost2, age == '0-54')
  death_5574 = subset(deathpost2, age == '55-74')
  
  propdeath75 = subset(propdeath, age== '75+')
  propdeath074 = subset(propdeath, age== '0-74')

  # in the united states
  b75 = format(round(death_75[, list(M = M_week1, CL = CL_week1, CU = CU_week1)]), big.mark=",") 
  a75 = format(round(death_75[, list(M = M_week2, CL = CL_week2, CU = CU_week2)]), big.mark=",") 
  p75 = paste0(format(round((1 - death_75[, list(M = M, CL = CL, CU = CU)])*100, 2), nsmall = 2), '\\%')
  
  b074 = format(round(death_074[, list(M = M_week1, CL = CL_week1, CU = CU_week1)]), big.mark=",") 
  a074 = format(round(death_074[, list(M = M_week2, CL = CL_week2, CU = CU_week2)]), big.mark=",") 
  p074 = paste0(format(round((1 - death_074[, list(M = M, CL = CL, CU = CU)])*100, 2), nsmall = 2), '\\%')

  b054 = format(round(death_054[, list(M = M_week1, CL = CL_week1, CU = CU_week1)]), big.mark=",") 
  a054 = format(round(death_054[, list(M = M_week2, CL = CL_week2, CU = CU_week2)]), big.mark=",") 
  p054 = paste0(format(round((1 - death_054[, list(M = M, CL = CL, CU = CU)])*100, 2), nsmall = 2), '\\%')
  
  b5574 = format(round(death_5574[, list(M = M_week1, CL = CL_week1, CU = CU_week1)]), big.mark=",") 
  a5574 = format(round(death_5574[, list(M = M_week2, CL = CL_week2, CU = CU_week2)]), big.mark=",") 
  p5574 = paste0(format(round((1 - death_5574[, list(M = M, CL = CL, CU = CU)])*100, 2), nsmall = 2), '\\%')
  
  
  # empirical 
  subset(death, date ==  date_end & age == '75+')
  emp = subset(death,date %in% c(date_before_vaccine, date_end))[, list(emp = sum(na.omit(emp))), by = c('date', 'age')]
  emp = as.data.table( reshape2::dcast(emp, age~ date, value.var = 'emp') )
  setnames(emp, 2:3, c('week1', 'week2'))
  emp[, prop := paste0(format(round((1 - (week2 / week1))*100, 2), nsmall = 2), '\\%')]
  emp[, week1 := format(week1, big.mark=",") ]
  emp[, week2 := format(week2, big.mark=",") ]
  
  # state fastest
  sf75 = propdeath75[order(M, decreasing = T), list(loc_label = loc_label, 
                                                    M = paste0(format(round((1 - M)*100, 2), nsmall = 2), '\\%'), 
                                                    CL = paste0(format(round((1-CL)*100, 2), nsmall = 2), '\\%'),
                                                    CU = paste0(format(round((1-CU)*100, 2), nsmall = 2), '\\%')), by = 'loc_label'][1:3]
  sl75 = propdeath75[order(M, decreasing = F), list(loc_label = loc_label, 
                                                    M = paste0(format(round((1-M)*100, 2), nsmall = 2), '\\%'), 
                                                    CL = paste0(format(round((1-CL)*100, 2), nsmall = 2), '\\%'),
                                                    CU = paste0(format(round((1-CU)*100, 2), nsmall = 2), '\\%')), by = 'loc_label'][1:3]
  
  # state slowest
  sf074 = propdeath074[order(M, decreasing = T), list(loc_label = loc_label, 
                                                      M = gsub(" ", '', paste0(format(round( (1 - M)*100, 2), nsmall = 2), '\\%')), 
                                                      CL = gsub(" ", '',paste0(format(round((1-CL)*100, 2), nsmall = 2), '\\%')),
                                                      CU = gsub(" ", '',paste0(format(round((1-CU)*100, 2), nsmall = 2), '\\%'))), by = 'loc_label'][1:3]
  sl074 = propdeath074[order(M, decreasing = F), list(loc_label = loc_label, 
                                                      M = gsub(" ", '',paste0(format(round((1 - M)*100, 2), nsmall = 2), '\\%')), 
                                                      CL = gsub(" ", '',paste0(format(round((1-CL)*100, 2), nsmall = 2), '\\%')),
                                                      CU = gsub(" ", '',paste0(format(round((1-CU)*100, 2), nsmall = 2), '\\%'))), by = 'loc_label'][1:3]
  
  # increase in death
  unique(death2$date)
  datemarch = as.Date('2021-02-27')
  dateapril = as.Date('2021-03-27')
  
  death_75 = subset(death2, age== '75+' & date %in% c(datemarch, '2021-05-15'))
  death_75 = as.data.table( reshape2::dcast(death_75, loc_label~ date, value.var = 'M') )
  setnames(death_75, 2:3, c('week1', 'week2'))
  tmp = death_75[, list(d = week2- week1, r = week2/week1), by = 'loc_label']
  imarch = tmp[r > 1 & d > 20]$loc_label
  
  death_75 = subset(death2, age== '75+' & date %in% c(dateapril, '2021-05-15'))
  death_75 = as.data.table( reshape2::dcast(death_75, loc_label~ date, value.var = 'M') )
  setnames(death_75, 2:3, c('week1', 'week2'))
  tmp = death_75[, list(d = week2- week1, r = week2/week1), by = 'loc_label']
  iapril = tmp[r > 1 & d > 5]$loc_label
  iapril = iapril[!iapril %in% imarch]
  iapril = paste0(paste0(iapril[-length(iapril)], collapse = ', '), ' and ', iapril[length(iapril)])
  
  death_stats = list(date_before= format(date_before_vaccine,  '%B %d, %Y'),
                     date_end = format(date_end,  '%B %d, %Y'), 
                     length = as.numeric( (date_end - date_before_vaccine) / 7 - 1 ), 
                     us75 = list(b75, a75, p75), us074 = list(b074, a074, p074), 
                     s75 = list(sf75, sl75), s074 = list(sf074, sl074), emp = emp,
                     u054 = list(b054, a054, p054), u5574 = list(b5574, a5574, p5574),
                     id = list(imarch, format(datemarch, '%B %d, %Y'), iapril, format(dateapril, '%B %d, %Y') ))

  saveRDS(death_stats, file = paste0(outdir, '-absolutedeaths.rds'))
  
  return(death_stats)
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

save_p_value_vaccine_effects <- function(samples, names, outdir){
  
  names_fit <- names(samples)
  names <- names_fit[ grepl(paste(paste0('^',names),collapse = '|'),names_fit) ]
  
  tmp <- lapply(names, function(x) find_p_value_vaccine_effect(samples[[x]], x))
  tmp <- do.call('rbind', tmp)
  
  saveRDS(tmp, paste0(outdir.table, '-p_value_vaccine_effects.rds'))
  
}

save_resurgence_dates <- function(resurgence_dates, outdir){
  resurgence_dates[,start_resurgence_name := format(start_resurgence, '%d %b, %Y')]
  resurgence_dates[,stop_resurgence_name := format(stop_resurgence, '%d %b, %Y')]
  
  tmp <- merge(resurgence_dates, df_state, by = 'code')[, .(loc_label, start_resurgence, stop_resurgence)]
  
  saveRDS(tmp, paste0(outdir, '-resurgence_dates.rds'))
}


find_stats_vaccine_effects <- function(data_res1, data_res2, data_res3, data_res4, prop_vac, resurgence_dates, outdir){
  
  data_res1 = merge(data_res1, resurgence_dates, by = 'code')
  prop_vac = merge(prop_vac, resurgence_dates, by = 'code')
  data_res2 = merge(data_res2, resurgence_dates, by = 'code')
  
  stat = list(format(c(min(resurgence_dates$start_resurgence), max(resurgence_dates$stop_resurgence)),  '%B %d, %Y'),
              subset(data_res1, date == stop_resurgence)[, list(M = round(M), 
                                                                CL = round(CU), 
                                                                CU = round(CL)), by = c('counterfactual_index', 'age', 'loc_label')],
              subset(prop_vac, date == start_resurgence)[, list(min_3 = paste0(round(min(prop_1*100), 2), '\\%'),
                                                                                                   max_3 = paste0(round(max(prop_1*100), 2), '\\%'),
                                                                                                   min_4 = paste0(round(min(prop_2*100), 2), '\\%'),
                                                                                                   max_4 = paste0(round(max(prop_2*100), 2), '\\%'))],
              subset(data_res2, date == stop_resurgence)[, list(M = format(round(M*100, digits = 2), nsmall = 2), 
                                                                CL = format(round(CU*100, digits = 2), nsmall = 2), 
                                                                CU = format(round(CL*100, digits = 2), nsmall = 2)), by = c('counterfactual_index', 'age', 'loc_label')],
              subset(data_res3, week_index == max(week_index))[, list(M = round(M), 
                                                                CL = round(CU), 
                                                                CU = round(CL)), by = c('counterfactual_index', 'age')],
              subset(data_res4, week_index == max(week_index))[, list(M = format(round(M*100, digits = 2), nsmall = 2), 
                                                                CL = format(round(CU*100, digits = 2), nsmall = 2), 
                                                                CU = format(round(CL*100, digits = 2), nsmall = 2)), by = c('counterfactual_index', 'age')]
              
              )
  saveRDS(stat, file = paste0(outdir, paste0('-Mortality_counterfactual.rds')))
  
  return(stat)
}

find_stats_vaccine_effects_old <- function(data_res1, data_res2, prop_vac, resurgence_dates, outdir){
  
  data_res1 = merge(data_res1, resurgence_dates, by = 'code')
  prop_vac = merge(prop_vac, resurgence_dates, by = 'code')
  data_res2 = merge(data_res2, resurgence_dates, by = 'code')
  
  stat = list(format(c(min(resurgence_dates$start_resurgence), max(resurgence_dates$stop_resurgence)),  '%B %d, %Y'),
              subset(data_res1, date == stop_resurgence)[, list(M = round(M), 
                                                                CL = round(CU), 
                                                                CU = round(CL)), by = c('age', 'loc_label')],
              subset(prop_vac, date == start_resurgence)[, list(min_3 = paste0(round(min(prop_1*100), 2), '\\%'),
                                                                max_3 = paste0(round(max(prop_1*100), 2), '\\%'),
                                                                min_4 = paste0(round(min(prop_2*100), 2), '\\%'),
                                                                max_4 = paste0(round(max(prop_2*100), 2), '\\%'))],
              subset(data_res2, date == stop_resurgence)[, list(M = format(round(M*100, digits = 2), nsmall = 2), 
                                                                CL = format(round(CU*100, digits = 2), nsmall = 2), 
                                                                CU = format(round(CL*100, digits = 2), nsmall = 2)), by = c('age', 'loc_label')]
              
  )
  saveRDS(stat, file = paste0(outdir, paste0('-Mortality_counterfactual.rds')))
  
  return(stat)
}


find_prop_deaths_vaccine_statistics <- function(propdeath3, start_vaccine, start_resurgence, outdir){
  
  dmean75 = propdeath3[age == '75+', list(M = paste0(format(round((- mean(M))*100, 2), nsmall = 2), '\\%'), 
                                          CL = paste0(format(round((-mean(CU))*100, 2), nsmall = 2), '\\%'),
                                          CU = paste0(format(round((-mean(CL))*100, 2), nsmall = 2), '\\%'))]
  dmean5574 = propdeath3[age == '55-74', list(M = paste0(format(round((- mean(M))*100, 2), nsmall = 2), '\\%'), 
                                              CL = paste0(format(round((-mean(CU))*100, 2), nsmall = 2), '\\%'),
                                              CU = paste0(format(round((-mean(CL))*100, 2), nsmall = 2), '\\%'))]
  dmean054 = propdeath3[age == '0-54', list(M = paste0(format(round((- mean(M))*100, 2), nsmall = 2), '\\%'), 
                                            CL = paste0(format(round((-mean(CU))*100, 2), nsmall = 2), '\\%'),
                                            CU = paste0(format(round((-mean(CL))*100, 2), nsmall = 2), '\\%'))]
  # state fastest and slowest
  sf75 = propdeath3[order(M, decreasing = T) & age == '75+', list(loc_label = loc_label, 
                                                                  M = paste0(format(round((- M)*100, 2), nsmall = 2), '\\%'), 
                                                                  CL = paste0(format(round((-CU)*100, 2), nsmall = 2), '\\%'),
                                                                  CU = paste0(format(round((-CL)*100, 2), nsmall = 2), '\\%')), by = 'loc_label'][c(1,length(locs))]
  
  
  tmp <- list(format(c(start_vaccine, start_resurgence-7),  '%B %d, %Y'), 
              as.numeric( (start_resurgence-7 - start_vaccine) / 7 - 1 ),
              list(dmean75, dmean5574, dmean054), 
              sf75)
  
  saveRDS(tmp, paste0(outdir, '-prop_red_deaths_vaccine.rds'))
}

