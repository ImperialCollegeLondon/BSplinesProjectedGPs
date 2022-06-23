make_predictive_checks_table = function(fit, fit_samples, df_week, df_state_age, data, deaths_predict_var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
  setnames(tmp1, 2:4, c('state_index', 'age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_state_age$age[age_index]]
  
  tmp1 = merge(tmp1, data, by = c('date', 'age', 'state_index'))
  tmp1[, age := factor(age, levels = levels(data$age))]
  
  # stat
  tmp1[, inside.CI := weekly.deaths <= CU & weekly.deaths >= CL]
  
  # save
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-predictive_checks_table_', Code, '.rds'))
  }
  
  return(tmp1)
}

make_convergence_diagnostics_stats = function(fit, re, outdir)
  {
  
  stopifnot(!is.null(fit))
  
  summary = rstan::summary(fit)$summary
  eff_sample_size_cum = summary[,9][!is.na(summary[,9])]
  Rhat_cum = summary[,10][!is.na(summary[,10])]
  cat("the minimum and maximum effective sample size are ", range(eff_sample_size_cum), "\n")
  cat("the minimum and maximum Rhat are ", range(Rhat_cum), "\n")
  if(min(eff_sample_size_cum) < 500) cat('\nEffective sample size smaller than 500 \n')
  
  tryCatch({
    sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
    sampler_diagnostics <- data.table()
    for (i in colnames(sampler_params[[1]])) {
      tmp <- data.table(t(sapply(sampler_params, function(x) quantile(x[, i],probs = c(0.025,0.5,0.975)))))
      tmp[, diagnostics:=i ]
      tmp[, chain:= seq_len(length(sampler_params))]
      sampler_diagnostics <- rbind(sampler_diagnostics, tmp)
    }
    print(sampler_diagnostics)
  }, error = function(e) e)
  
  
  
  check_all_diagnostics(fit, outdir)
  
  # compute WAIC and LOO
  tryCatch({
    
    if('log_lik' %in% names(re)){
      log_lik <- loo::extract_log_lik(fit)
      log_lik = log_lik[!is.na(log_lik[,1]),]
      .WAIC = loo::waic(log_lik)
      .LOO = loo::loo(log_lik)
      print(.WAIC); print(.LOO)
      WAIC = .WAIC$pointwise
      LOO = .LOO$pointwise
    }} , error = function(e) e)
  
  # time of execution
  time = sum(rstan::get_elapsed_time(fit))
  
  # save
  saveRDS(eff_sample_size_cum, file = paste0(outdir, "-eff_sample_size_cum.rds"))
  saveRDS(Rhat_cum, file = paste0(outdir, "-Rhat_cum.rds"))
  saveRDS(.WAIC, file = paste0(outdir, "-WAIC.rds"))
  saveRDS(.LOO, file = paste0(outdir, "-LOO.rds"))
  saveRDS(sampler_diagnostics, file = paste0(outdir, "-sampler_diagnostics.rds"))
  saveRDS(time, file = paste0(outdir, "-time_elapsed.rds"))
  
  for(i in Code){
    saveRDS(log_lik, file = paste0(outdir, "-log_lik_", i, ".rds"))
  }
  
  return(summary)
}

make_var_by_age_by_state_by_time_table = function(fit_samples, df_week, df_state_age, df_state, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:4, c('state_index', 'age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]

  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'Table_', Code, '.rds'))
  }

  return(tmp1)
}

find_diff_E_pdeaths_counterfactual_all = function(fit_samples, df_week, df_state_age, df_counterfactual, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[['diff_E_pdeaths_counterfactual']]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'age_index','week_index', 'counterfactual_index')]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index', 'week_index', 'counterfactual_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + week_index + age_index ~ q_label, value.var = "q")
  
  # tmp1 = merge(tmp1, df_week, by = 'week_index')
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  tmp1 <- merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  saveRDS(tmp1, file = paste0(outdir, '-', 'diff_E_pdeaths_counterfactual_all',  'AllStatesTable.rds'))
  
  return(tmp1)
}

find_inv_diff_E_pdeaths_counterfactual_all = function(fit_samples, df_week, df_state_age, df_counterfactual, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[['diff_E_pdeaths_counterfactual']]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'age_index','week_index', 'counterfactual_index')]
  
  tmp1[, value := - value]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index', 'week_index', 'counterfactual_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + week_index + age_index ~ q_label, value.var = "q")
  
  # tmp1 = merge(tmp1, df_week, by = 'week_index')
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  tmp1 <- merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  saveRDS(tmp1, file = paste0(outdir, '-', 'diff_E_pdeaths_counterfactual_all',  'AllStatesTable.rds'))
  
  return(tmp1)
}

find_inv_diff_E_pdeaths_counterfactual_allstatesages = function(fit_samples, df_week, df_counterfactual, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[['diff_E_pdeaths_counterfactual']]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations','week_index', 'counterfactual_index')]
  
  tmp1[, value := - value]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'counterfactual_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + week_index ~ q_label, value.var = "q")
  
  # tmp1 = merge(tmp1, df_week, by = 'week_index')
  
  tmp1 <- merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  saveRDS(tmp1, file = paste0(outdir, '-', 'diff_E_pdeaths_counterfactual_all',  'AllStatesAllAgesTable.rds'))
  
  return(tmp1)
}

find_perc_E_pdeaths_counterfactual_all = function(fit_samples, df_week, df_state_age, df_counterfactual, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[['diff_E_pdeaths_counterfactual']]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'age_index','week_index', 'counterfactual_index')]
  
  tmp2 = as.data.table( reshape2::melt(fit_samples[['E_pdeaths_predict_resurgence_cumulative']]) )
  setnames(tmp2, 2:4, c('state_index', 'age_index','week_index'))
  tmp2 <- tmp2[, list(value_d = sum(value)), by = c('iterations', 'age_index','week_index')]
  
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'age_index','week_index'))
  tmp1[, value := value / value_d]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index', 'week_index', 'counterfactual_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + week_index + age_index ~ q_label, value.var = "q")
  
  # tmp1 = merge(tmp1, df_week, by = 'week_index')
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  tmp1 <- merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  saveRDS(tmp1, file = paste0(outdir, '-', 'perc_E_pdeaths_counterfactual_all',  'AllStatesTable.rds'))
  
  return(tmp1)
}

make_ratio_vars_by_age_state_by_counterfactual_table = function(fit_samples, df_week, df_state_age, df_state, df_counterfactual, vars_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[vars_name[1]]]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  
  tmp2 = as.data.table( reshape2::melt(fit_samples[[vars_name[2]]]) )
  setnames(tmp2, 2:5, c('state_index', 'age_index','week_index', 'value_denominator'))
  
  tmp1 <- merge(tmp1, tmp2,  by = c('iterations', 'age_index', 'state_index','week_index'))
  tmp1[, value := value / value_denominator]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'age_index', 'state_index', 'week_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + age_index + state_index + week_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  tmp1 = merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', vars_name[1], 'RatioTable_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

make_var_by_age_by_state_by_time_table_relative = function(fit_samples, df_week, df_state_age, df_state, var_name, Code, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:4, c('state_index', 'age_index','week_index'))
  
  
  tmp1[, value_rel := value[state_index == df_state[, which(code == Code)] ], by = c('age_index', 'week_index', 'iterations')]
  tmp1[, value := value / value_rel]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'Table_', Code, '.rds'))
    
  }
  
  return(tmp1)
}


make_var_by_age_by_state_table = function(fit_samples, df_state_age, df_state, var_name, outdir, vaccine_data = NULL){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:3, c('state_index', 'age_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'age_index')]	
  tmp1 = dcast(tmp1, state_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  if(!is.null(vaccine_data)){
    tmp1 <- merge(tmp1, vaccine_data, by= c('code', 'loc_label', 'age'))
  }
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'Table_', Code, '.rds'))
    
  }
  
  return(tmp1)
}


make_var_by_state_table = function(fit_samples, df_week, df_state, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:4, c('state_index', 'age_index','week_index'))
  
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'state_index','week_index')]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'week_index')]	
  tmp1 = dcast(tmp1, state_index + week_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'AllAgesTable_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

make_var_by_age_table = function(fit_samples, df_week, df_state_age, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:3, c('age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c( 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1,  week_index + age_index ~ q_label, value.var = "q")
  # tmp1 = merge(tmp1, df_week, by = 'week_index')
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  saveRDS(tmp1, file = paste0(outdir, '-', var_name,  'AllStatesTable.rds'))

  
  return(tmp1)
}

make_var_cum_by_age_table = function(fit_samples, df_week, df_state_age, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:4, c('state_index', 'age_index','week_index'))
  tmp1[, value := cumsum(value), by = c('state_index', 'age_index', 'iterations')]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'CumTable_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

make_contribution_ref = function(fit_samples, data_10thdeaths, fiveagegroups, data, df_week, df_age_continuous, outdir){
  
  # df age
  df_age = data.table(age = fiveagegroups)
  df_age[, age_index := 1:nrow(df_age)]
  df_age[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  df_age[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_age[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  set(df_age, NULL, 'age_from', df_age[,as.numeric(age_from)])
  set(df_age, NULL, 'age_to', df_age[,as.numeric(age_to)])
  set(df_age, NULL, 'age_to_index', df_age[,as.numeric(age_to_index)])
  set(df_age, NULL, 'age_from_index', df_age[,as.numeric(age_from_index)])
  df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
  df_age_continuous[, age_state_index := which(df_age$age_from_index <= age_index & df_age$age_to_index >= age_index), by = 'age_index']
  
  # week data table
  df_week[, dummy := 1]
  data_10thdeaths[, dummy := 1]
  df_week1 = merge(df_week, data_10thdeaths, by = 'dummy', allow.cartesian = T)
  df_week1 = df_week1[date >= date_10thcum]
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1[, keep := month %in% min(df_week1$month):(min(df_week1$month) + 2), by = 'code']
  df_week1 = df_week1[keep == T]
  
  # empirical estimate
  df_age_reporting[, age_state_index := which(df_age$age_from <= age_from & df_age$age_to >= age_to), by = 'age_index']
  empirical = merge(data, df_age_reporting, by = 'age_from')
  empirical = merge(empirical, df_week1, by = c('date', 'code'))
  empirical = empirical[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c( 'age_state_index', 'code')]
  empirical[, emp_est := weekly.deaths / sum(na.omit(weekly.deaths)), by = c('code')]
  setnames(empirical, 'age_state_index', 'age_index')
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[['phi']]) )
  setnames(tmp1,2:4, c('state_index', 'age_index','week_index'))
  
  # by selected age groups
  tmp1 = merge(tmp1, df_age_continuous, by = 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('iterations','state_index',  'week_index', 'age_state_index')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # take reference
  df_week1 = merge(df_week1, df_state, by = 'code')
  tmp1 = merge(tmp1, select(df_week1, week_index, state_index), by = c('state_index', 'week_index'))
  tmp1 = tmp1[, list(value = mean(value)), by = c('iterations', 'state_index', 'age_index')]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'age_index')]	
  tmp1 = dcast(tmp1, state_index + age_index ~ q_label, value.var = "q")
  
  tmp1[, age := fiveagegroups[age_index]]
  tmp1[, age := factor(age, levels = fiveagegroups)]
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  tmp1 = merge(tmp1, select(empirical, age_index, emp_est, code), by = c('age_index', 'code'))
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', 'contribution_ref',  'Table_', Code, '.rds'))
    
  }

  return(tmp1)
}

make_contribution_ref_adj = function(fit_samples, data_10thdeaths, fouragegroups, df_week, pop_data, outdir){
  
  # df age
  df_age = data.table(age = fouragegroups)
  df_age[, age_index := 1:nrow(df_age)]
  df_age[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  df_age[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_age[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  set(df_age, NULL, 'age_from', df_age[,as.numeric(age_from)])
  set(df_age, NULL, 'age_to', df_age[,as.numeric(age_to)])
  set(df_age, NULL, 'age_to_index', df_age[,as.numeric(age_to_index)])
  set(df_age, NULL, 'age_from_index', df_age[,as.numeric(age_from_index)])
  df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
  df_age_continuous[, age_state_index := which(df_age$age_from_index <= age_index & df_age$age_to_index >= age_index), by = 'age_index']
  
  # find population proportion
  pop_data1 = copy(pop_data)
  pop_data1[, age_state_index := which(df_age$age_from <= age & df_age$age_to >= age), by = 'age']
  pop_data1 = pop_data1[, list(pop = sum(pop)), by = c('age_state_index', 'code')]
  setnames(pop_data1, 'age_state_index', 'age_index')
  pop_data1[, age := fouragegroups[age_index]]
  tmp <- pop_data1[, list(Total = sum(pop)), by = 'code']
  pop_data1 <- merge(tmp, pop_data1, by = 'code')
  pop_data1[, pop_prop := pop / Total]
  tmp = pop_data1[, list(pop = sum(pop)), by = 'age']
  tmp[, pop_prop_US := pop / sum(pop_data$pop)]
  pop_data1 = merge(pop_data1, select(tmp, age, pop_prop_US), by = 'age')
  pop_data1 = subset(pop_data1, code %in% df_state$code)
  
  # ref week
  df_week[, dummy := 1]
  data_10thdeaths[, dummy := 1]
  df_week1 = merge(df_week, data_10thdeaths, by = 'dummy', allow.cartesian = T)
  df_week1 = df_week1[date >= date_10thcum]
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1[, keep := month %in% min(df_week1$month):(min(df_week1$month) + 2), by = 'code']
  df_week1 = df_week1[keep == T]
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[['phi']]) )
  setnames(tmp1,2:4, c('state_index', 'age_index','week_index'))
  
  # group by age group
  tmp1 = merge(tmp1, df_age_continuous, by = 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('iterations','state_index',  'week_index', 'age_state_index')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # take reference
  df_week1 = merge(df_week1, df_state, by = 'code')
  tmp1 = merge(tmp1, select(df_week1, week_index, state_index), by = c('state_index', 'week_index'))
  tmp1 = tmp1[, list(value = mean(value)), by = c('iterations', 'state_index', 'age_index')]
  
  # adjust by population proportion
  pop_data1 = merge(pop_data1, df_state, by = 'code')
  tmp1 = merge(tmp1, select(pop_data1, age_index, pop_prop, pop_prop_US, state_index), by = c('age_index', 'state_index'))
  tmp1[, value := value / pop_prop * pop_prop_US, by = 'state_index']
  tmp1[, value := value / sum(value), by = c('iterations', 'state_index')]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index', 'state_index')]	
  tmp1 = dcast(tmp1, state_index + age_index ~ q_label, value.var = "q")
  
  tmp1[, age := df_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_age$age)]
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', 'contribution_ref_adj',  'Table_', Code, '.rds'))
  }

  return(tmp1)
}

find_contribution_one_age_group = function(fit_samples, df_week, df_age_continuous, df_age_reporting, age_group, 
                                           date_10thcum, pop_data, data, outdir, with_empirical = F){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # find ref week 
  df_week[, dummy := 1]
  date_10thcum[, dummy := 1]
  df_week1 = merge(df_week, date_10thcum, by = 'dummy', allow.cartesian = T)
  df_week1 = df_week1[date >= date_10thcum]
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1[, keep := month %in% min(df_week1$month):(min(df_week1$month) + 2), by = 'code']
  df_week1 = df_week1[keep == T]
  
  # df age
  if(grepl('\\+', age_group))
  {
    age_groups = c(paste0('0-', as.numeric(gsub('(.+)\\+', '\\1', age_group))-1), age_group)
  } else if( grepl('0-.*',age_group) & gsub('(.+)0-.*', '\\1', age_group) == age_group )
  {
    
    age_groups = c(age_group, paste0(as.numeric(gsub('.*-(.+)', '\\1', age_group))+1, '+'))
  } else{
    age_groups = c(paste0('0-', as.numeric(gsub('(.+)-.*', '\\1', age_group))-1), age_group,
                   paste0(as.numeric(gsub('.*-(.+)', '\\1', age_group))+1, '+'))
    
  }

  df_age_state = data.table(age = age_groups)
  df_age_state[, age_index := 1:nrow(df_age_state)]
  df_age_state[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age_state[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age_state[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age_state[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  df_age_state[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_age_state[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  set(df_age_state, NULL, 'age_from', df_age_state[,as.numeric(age_from)])
  set(df_age_state, NULL, 'age_to', df_age_state[,as.numeric(age_to)])
  set(df_age_state, NULL, 'age_to_index', df_age_state[,as.numeric(age_to_index)])
  set(df_age_state, NULL, 'age_from_index', df_age_state[,as.numeric(age_from_index)])
  df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
  df_age_continuous[, age_state_index := which(df_age_state$age_from_index <= age_index & df_age_state$age_to_index >= age_index), by = 'age_index']
  
  # find population proportion
  pop_data1 = copy(pop_data)
  pop_data1[, age_state_index := which(df_age_state$age_from <= age & df_age_state$age_to >= age), by = c('age')]
  
  tmp <- pop_data1[, list(pop = sum(pop)), by = 'age']
  tmp[, pop_prop_US := pop / sum(tmp$pop)]
  pop_data1 = merge(pop_data1, select(tmp, age, pop_prop_US), by = 'age')

  tmp <- pop_data1[, list(Total = sum(pop)), by = 'code']
  pop_data1 = merge(pop_data1, tmp, by = 'code')
  
  pop_data1[, pop_prop := pop / Total]
  pop_data1 = pop_data1[, list(pop = sum(pop), pop_prop = sum(pop_prop), pop_prop_US = sum(pop_prop_US)), by = c('Total', 'code', 'age_state_index')]
  pop_data1[, multiplier := pop_prop / pop_prop_US]
  pop_data1 = subset(pop_data1, code %in% df_state$code)
  stopifnot(nrow(subset(pop_data1, code == df_state$code[1])) == length(age_groups))
  
  # empirical estimate
  if(with_empirical){
    df_age_reporting[, age_state_index := which(df_age_state$age_from <= age_from & df_age_state$age_to >= age_to), by = 'age_index']
    empirical = merge(data, df_age_reporting, by = 'age_from')
    empirical = empirical[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c('date', 'age_state_index', 'code')]
    baseline_emp = merge(empirical, df_week1, by = c('date', 'code'))
    baseline_emp = baseline_emp[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c( 'age_state_index', 'code')]
    baseline_emp[, emp_ref := weekly.deaths / sum(na.omit(weekly.deaths)), by = 'code']
    empirical1 = empirical[, list(weekly.deaths_t = sum(na.omit(weekly.deaths))), by = c( 'date', 'code')]
    empirical = merge(empirical, empirical1, by = c( 'date', 'code'))
    empirical = merge(empirical, select(baseline_emp, age_state_index, emp_ref, code), by = c( 'age_state_index', 'code'))
    empirical = merge(empirical, pop_data1, by = c('age_state_index', 'code'))
    empirical[, emp := weekly.deaths / weekly.deaths_t]
    empirical[, emp_rel := emp / emp_ref]
    empirical[, emp_adj := emp * multiplier]
    empirical[, emp_adj := emp_adj / sum(emp_adj), by = c('date', 'code')]
    empirical[, emp_ref_adj :=  emp_ref * multiplier]
    empirical[, emp_ref_adj :=  emp_ref_adj / sum(emp_ref_adj), by = c('date', 'code')]
    empirical[, emp_rel_adj :=  emp_adj / emp_ref_adj]
    
    empirical = subset(empirical, age_state_index == which(age_groups == age_group))
    setnames(empirical, 'age_state_index', 'age_index')
    
    empirical = select(empirical, date, code, age_index, emp, emp_rel, emp_adj, emp_rel_adj)
    
  }

  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[['phi']]) )
  setnames(tmp1,2:4, c('state_index', 'age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('state_index', 'week_index', 'iterations', 'age_state_index')]
  
  # # find reference
  df_week1 = merge(df_week1, df_state, by = 'code')
  tmp1 = merge(tmp1, select(df_week1, week_index, state_index), by = c('state_index', 'week_index'))
  tmp2 = tmp1[, list(value_ref = mean(value)), by = c('iterations', 'state_index', 'age_state_index')]
  tmp1 = merge(tmp1, tmp2, by = c('iterations', 'state_index', 'age_state_index'))
  tmp1[, value_rel := value / value_ref]
  
  # adjust by population proportion
  pop_data1 = merge(pop_data1, df_state, by = 'code')
  tmp1 = merge(tmp1, pop_data1, by = c('age_state_index', 'state_index'))
  tmp1[, value_adj := value * multiplier]
  tmp1[, value_adj := value_adj / sum(value_adj), by = c('iterations', 'week_index', 'state_index')]
  tmp1[, value_ref_adj := value_ref * multiplier]
  tmp1[, value_ref_adj := value_ref_adj / sum(value_ref_adj), by = c('iterations', 'week_index', 'state_index')]
  tmp1[, value_rel_adj := value_adj / value_ref_adj ]
  
  # keep only interested age group 
  tmp1 = subset(tmp1, age_state_index == which(age_groups == age_group))
  
  # take quantiles
  tmp2 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'week_index')]	
  tmp2 = dcast(tmp2, state_index + week_index ~ q_label, value.var = "q")
  
  tmp3 = tmp1[, list( 	q= quantile(value_rel, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'week_index')]	
  tmp3 = dcast(tmp3, state_index + week_index ~ q_label, value.var = "q")
  setnames(tmp3, p_labs, paste0(p_labs, '_rel'))
  tmp2 = merge(tmp2, tmp3, by =c('state_index', 'week_index'))
  
  tmp3 = tmp1[, list( 	q= quantile(value_adj, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'week_index')]	
  tmp3 = dcast(tmp3, state_index + week_index ~ q_label, value.var = "q")
  setnames(tmp3, p_labs, paste0(p_labs, '_adj'))
  tmp2 = merge(tmp2, tmp3, by = c('state_index', 'week_index'))
  
  tmp1 = tmp1[, list( 	q= quantile(value_rel_adj, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'week_index')]	
  tmp1 = dcast(tmp1, state_index + week_index ~ q_label, value.var = "q")
  setnames(tmp1, p_labs, paste0(p_labs, '_rel_adj'))
  tmp1 = merge(tmp2, tmp1, by = c('state_index', 'week_index'))
  
  tmp1[, age := age_group]
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  
  tmp1 = merge(tmp1, date_10thcum, by = 'code')
  tmp1[, after.10thcumdeaths := date >= date_10thcum]
  
  if(with_empirical){
    tmp1= merge(tmp1, empirical, by = c('date', 'code'), all.x = T)

  }

  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-Contribution_Age_', age_group,'_', Code, '.rds'))
  }

  return(tmp1)
}


find_cumsum_nonr_deaths_state_age = function(fit_samples, df_week, df_age_continuous, state_age_groups, stan_data, deaths_predict_var){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # df age
  df_age_state = data.table(age = state_age_groups)
  df_age_state[, age_index := 1:nrow(df_age_state)]
  df_age_state[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age_state[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age_state[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age_state[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  df_age_state[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_age_state[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  set(df_age_state, NULL, 'age_from', df_age_state[,as.numeric(age_from)])
  set(df_age_state, NULL, 'age_to', df_age_state[,as.numeric(age_to)])
  set(df_age_state, NULL, 'age_to_index', df_age_state[,as.numeric(age_to_index)])
  set(df_age_state, NULL, 'age_from_index', df_age_state[,as.numeric(age_from_index)])
  df_age_continuous[, age_index := as.numeric(1:nrow(df_age_continuous))]
  df_age_continuous[, age_state_index := which(df_age_state$age_from_index <= age_index & df_age_state$age_to_index >= age_index), by = 'age_index']
  
  # include missing information
  df_missing = vector(mode = 'list', length = stan_data$M)
  for(m in 1:stan_data$M){
    df_missing_m = NULL
    for(x in 1:stan_data$N_missing[[m]]){
      df_missing_m = rbind(df_missing_m, data.table(idx_serie_missing = x, week_index = stan_data$idx_weeks_missing_min[[m]][x]:stan_data$idx_weeks_missing_max[[m]][x])
      )
    }
    
    tmp = data.table(age_index = stan_data$age_missing[[m]],
                     min_count_censored = stan_data$min_count_censored[[m]],
                     max_count_censored = stan_data$max_count_censored[[m]],
                     sum_count_censored = stan_data$sum_count_censored[[m]],
                     idx_serie_missing = 1:length(stan_data$age_missing[[m]]))
    tmp[min_count_censored == -1, min_count_censored := NA]
    tmp[max_count_censored == -1, max_count_censored := NA]
    tmp[sum_count_censored == -1, sum_count_censored := NA]
    df_missing[[m]] = merge(df_missing_m, tmp, by = 'idx_serie_missing')
    df_missing[[m]]$state_index = m
  }
  df_missing = do.call('rbind', df_missing)
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
  setnames(tmp1, 2:4, c('state_index', 'age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('state_index', 'week_index', 'age_state_index', 'iterations')]
  tmp1 = merge(tmp1, df_missing, by.x = c('age_state_index', 'state_index', 'week_index'), by.y = c('age_index','state_index', 'week_index'))
  
  # take the cum sum
  tmp1 = tmp1[, value := cumsum(value), by = c('iterations', 'age_state_index', 'state_index')]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index', 'week_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_state_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1 = merge(tmp1, df_age_state, by.x = 'age_state_index', by.y = 'age_index')
  tmp1 = merge(tmp1, df_missing, by.x = c('age_state_index','week_index', 'state_index'), by.y = c('age_index','week_index', 'state_index'))
  
  setnames(tmp1, 'age_state_index', 'age_index')
  
  return(tmp1)
}

find_sum_nonr_deaths_state_age = function(fit_samples, df_age_continuous, state_age_groups, stan_data, deaths_predict_var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # df age
  df_age_state = data.table(age = state_age_groups)
  df_age_state[, age_index := 1:nrow(df_age_state)]
  df_age_state[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age_state[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age_state[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age_state[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  df_age_state[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_age_state[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  set(df_age_state, NULL, 'age_from', df_age_state[,as.numeric(age_from)])
  set(df_age_state, NULL, 'age_to', df_age_state[,as.numeric(age_to)])
  set(df_age_state, NULL, 'age_to_index', df_age_state[,as.numeric(age_to_index)])
  set(df_age_state, NULL, 'age_from_index', df_age_state[,as.numeric(age_from_index)])
  df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
  df_age_continuous[, age_state_index := which(df_age_state$age_from_index <= age_index & df_age_state$age_to_index >= age_index), by = 'age_index']
  
  # include missing information
  df_missing = vector(mode = 'list', length = stan_data$M)
  for(m in 1:stan_data$M){
    df_missing_m = NULL
    for(x in 1:stan_data$N_missing[[m]]){
      df_missing_m = rbind(df_missing_m, data.table(idx_serie_missing = x, week_index = stan_data$idx_weeks_missing_min[[m]][x]:stan_data$idx_weeks_missing_max[[m]][x])
      )
    }
    
    tmp = data.table(age_index = stan_data$age_missing[[m]],
                     min_count_censored = stan_data$min_count_censored[[m]],
                     max_count_censored = stan_data$max_count_censored[[m]],
                     sum_count_censored = stan_data$sum_count_censored[[m]],
                     idx_serie_missing = 1:length(stan_data$age_missing[[m]]))
    tmp[min_count_censored == -1, min_count_censored := NA]
    tmp[max_count_censored == -1, max_count_censored := NA]
    tmp[sum_count_censored == -1, sum_count_censored := NA]
    df_missing[[m]] = merge(df_missing_m, tmp, by = 'idx_serie_missing')
    df_missing[[m]]$state_index = m
  }
  df_missing = do.call('rbind', df_missing)
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
  setnames(tmp1, 2:4, c('state_index', 'age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('state_index','week_index', 'age_state_index', 'iterations')]
  tmp1 = merge(tmp1, df_missing, by.x = c('state_index','age_state_index','week_index'), by.y = c('state_index','age_index','week_index'))
  
  # take the cum sum
  tmp1 = tmp1[, list(value = sum(value)), by = c('state_index','age_state_index', 'iterations')]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('state_index','age_state_index')]	
  tmp1 = dcast(tmp1,state_index + age_state_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  tmp1 = merge(tmp1, df_age_state, by.x = 'age_state_index', by.y = 'age_index')
  tmp1 = merge(tmp1, unique(select(df_missing, -week_index)),by.x = c('age_state_index', 'state_index'), by.y = c('age_index', 'state_index'))
  
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # inside CI?
  tmp1[is.na(sum_count_censored), inside.CI := (CL >= min_count_censored & CL <= max_count_censored) | (CU >= min_count_censored & CU <= max_count_censored) | (CL <= min_count_censored & CU >= max_count_censored)]
  tmp1[!is.na(sum_count_censored), inside.CI := sum_count_censored >= CL & sum_count_censored <= CU]
  
  # save
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-predictive_checks_table2_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

find_cumulative_deaths_prop_givensum_state_age_multiple_states <- function(fit_samples, date_10thcum, df_week, df_age_continuous, scrapedData, cum.death.var, 
                                                                           Code, outdir){
  if(nrow(subset(scrapedData, code %in% Code)) > 0 ){
    
    tmp = list(); i = 1
    for(c in Code){

      scrapedData_c = subset(scrapedData, code == c)
      
      if(c == 'GA'){
        scrapedData_c = reduce_agebands_scrapedData_GA(scrapedData_c)
      }
      
      if(nrow(scrapedData_c) > 0){
        tmp[[i]] = find_cumulative_deaths_prop_givensum_state_age(fit_samples, date_10thcum, df_week, df_age_continuous, scrapedData_c, cum.death.var, c, outdir)
        i = i + 1
      }
      
    }

    tmp = do.call('rbind',tmp)
    tmp = merge(tmp, df_state, by = 'code')
    
    return(tmp)
  }
  
  
}

find_cumulative_deaths_prop_givensum_state_age = function(fit_samples, date_10thcum, df_week, df_age_continuous, data_comp, cum.death.var, c, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # df age
  df_age_state = data.table(age = unique(subset(data_comp, code == c)$age))
  df_age_state[, age_index := 1:nrow(df_age_state)]
  df_age_state[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age_state[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age_state[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age_state[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  df_age_state[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_age_state[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  set(df_age_state, NULL, 'age_from', df_age_state[,as.numeric(age_from)])
  set(df_age_state, NULL, 'age_to', df_age_state[,as.numeric(age_to)])
  set(df_age_state, NULL, 'age_to_index', df_age_state[,as.numeric(age_to_index)])
  set(df_age_state, NULL, 'age_from_index', df_age_state[,as.numeric(age_from_index)])
  df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
  df_age_continuous[, age_state_index := which(df_age_state$age_from_index <= age_index & df_age_state$age_to_index >= age_index), by = 'age_index']
  
  # ref week
  df_week1 = subset(df_week, date >= date_10thcum[code == c]$date_10thcum)
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1 = subset(df_week1, month %in% min(df_week1$month):(min(df_week1$month) + 2))
  
  # data comp
  df_week_10thcum = subset(df_week, date >= date_10thcum[code == c]$date_10thcum)
  tmp2 = as.data.table( subset(data_comp, code == c))
  tmp2[, date := as.Date(date)]
  tmp2 = tmp2[, list(cum.death = sum(get(cum.death.var))), by = c('date', 'age')]
  tmp2 = tmp2[order( date)]
  dummy = min(df_week_10thcum$date) %in% tmp2$date 
  last.cum.deaths = sum(subset(tmp2, date == ifelse(dummy, min(df_week_10thcum$date), 
                                                    (min(tmp2$date):(min(tmp2$date)+6))[min(tmp2$date):(min(tmp2$date)+6) %in% df_week_10thcum$date] ) )$cum.death )
  tmp2 = unique(subset(tmp2, date %in% c(max(df_week_10thcum$date) + 7, df_week_10thcum$date)))
  tmp2[, weekly.deaths := c(diff(cum.death), NA), by = 'age']
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = 'date')
  tmp2 = select(tmp2, week_index, age, weekly.deaths)
  tmp3 = tmp2[, list(weekly.deaths_total = sum(weekly.deaths)), by = 'week_index']
  tmp3 = subset(tmp3, !is.na(weekly.deaths_total))
  tmp2 = merge(tmp2, tmp3, by = 'week_index')
  tmp2[, prop.weekly.deaths :=weekly.deaths / weekly.deaths_total]
  
  # adjust dfweek for cum
  df_week_adj = df_week[, list(date = date + 7), by = 'week_index']
  
  # find pi noisy
  tmp1 = as.data.table( reshape2::melt(fit_samples['deaths_predict']) )
  setnames(tmp1,2:4, c('state_index', 'age_index','week_index'))
  tmp1 = subset(tmp1, state_index == df_state[code == c]$state_index)
  
  # find proportion with overdispersion noise
  tmpt = tmp1[, list(deaths_predict_t = sum(value)), by = c('week_index', 'iterations')]
  tmp1 = merge(tmp1, tmpt, by = c('week_index', 'iterations'))
  tmp1[, value := value / deaths_predict_t]
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # predict weekly deaths
  tmp5 = as.data.table( reshape2::melt(fit_samples['alpha']) )
  setnames(tmp5,2:4, c('state_index', 'age_index','week_index'))
  tmp5 = subset(tmp5, state_index == df_state[code == c]$state_index)
  
  tmp5 = merge(tmp5, df_age_continuous, 'age_index')
  tmp5 = tmp5[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
  setnames(tmp5, 'age_state_index', 'age_index')
  
  # add cumsum from weeks before 10th date
  tmp4 = subset(tmp5, week_index %in% df_week1$week_index)
  tmp4 = tmp4[, list(value_ref = mean(value)), by = c('age_index', 'iterations')]
  tmp4[, value_ref := rdirmnom(1, last.cum.deaths, value_ref), by = c('iterations')]
  
  # predict  weekly deaths
  tmp5 = merge(tmp5, tmp3, by = 'week_index')
  tmp5[, value_abs := rdirmnom(1, weekly.deaths_total, value), by = c('iterations', 'week_index')]
  
  # find cumulative death
  tmp5 = merge(tmp5, tmp4, by = c('age_index', 'iterations'))
  tmp5[, value_abs_cum := cumsum(value_abs), by = c('iterations', 'age_index')]
  tmp5[, value_abs_cum := value_abs_cum + value_ref]
  
  #merge
  tmp5 = select(tmp5, -value)
  tmp1 = merge(tmp1, tmp5, by = c('iterations', 'age_index', 'week_index'))
  
  # take quantiles
  tmp4 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_index')]	
  tmp4 = dcast(tmp4, week_index + age_index ~ q_label, value.var = "q")
  setnames(tmp4, p_labs, paste0(p_labs, '_prop'))
  
  tmp5 = tmp1[, list(  q= quantile(value_abs, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_index')]	
  tmp5 = dcast(tmp5, week_index + age_index ~ q_label, value.var = "q")
  setnames(tmp5, p_labs, paste0(p_labs, '_abs_weekly'))
  tmp4 = merge(tmp4, tmp5, by = c('week_index', 'age_index'))
  
  tmp5 = tmp1[, list( 	q= quantile(value_abs_cum, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_index')]	
  tmp5 = dcast(tmp5, week_index + age_index ~ q_label, value.var = "q")
  setnames(tmp5, p_labs, paste0(p_labs, '_abs_cum'))
  tmp4 = merge(tmp4, tmp5, by = c('week_index', 'age_index'))
  
  tmp1 = merge(tmp4, select(df_age_state, age,age_index), by = 'age_index')
  tmp1 = merge(tmp1, tmp2, by = c('week_index', 'age'))
  
  tmp1 = merge(tmp1, df_week_adj, by = 'week_index')
  tmp4 = as.data.table( select(subset(data_comp, code == c), age, date, cum.deaths) )
  tmp4[, date := as.Date(date)]
  tmp1 = merge(tmp1, tmp4, by = c('date', 'age'))
  
  tmp1[weekly.deaths_total == 0, prop.weekly.deaths := 0]
  tmp1[, prop.death.inside.CI := (prop.weekly.deaths <= CU_prop & prop.weekly.deaths >= CL_prop)]
  tmp1[, weekly.death.inside.CI := ( weekly.deaths <= CU_abs_weekly & weekly.deaths  >= CL_abs_weekly)] 
  tmp1[, cum.death.inside.CI := (cum.deaths  <= CU_abs_cum & cum.deaths  >= CL_abs_cum)] 
  
  tmp1[, code := c]
  tmp1[, age := factor(age, levels = as.character(df_age_state[order(age_from)]$age))]
  
  saveRDS(list(data_comp, tmp1), file = paste0(outdir, '-CumDeathsComp_', 'ScrapedData', '_', c, '.rds'))
  
  return(tmp1)
}

make_mortality_rate_table_discrete = function(fit_samples, fouragegroups, date_10thcum, df_week, pop_data, data_comp, df_age_continuous, cum.death.var, nyt_data, outdir, lab = NULL){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # df age
  df_age = data.table(age = fouragegroups)
  df_age[, age_index := 1:nrow(df_age)]
  df_age[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  df_age[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_age[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  set(df_age, NULL, 'age_from', df_age[,as.numeric(age_from)])
  set(df_age, NULL, 'age_to', df_age[,as.numeric(age_to)])
  set(df_age, NULL, 'age_to_index', df_age[,as.numeric(age_to_index)])
  set(df_age, NULL, 'age_from_index', df_age[,as.numeric(age_from_index)])
  df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
  df_age_continuous[, age_state_index := which(df_age$age_from_index <= age_index & df_age$age_to_index >= age_index), by = 'age_index']
  
  # ref week
  df_week[, dummy := 1]
  date_10thcum[, dummy := 1]
  df_week1 = merge(df_week, date_10thcum, by = 'dummy', allow.cartesian = T)
  df_week1 = df_week1[date >= date_10thcum]
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1[, keep := month %in% min(df_week1$month):(min(df_week1$month) + 2), by = 'code']
  df_week1 = df_week1[keep == T]
  
  # find population proportion
  pop_data1 = copy(pop_data)
  pop_data1[, age_state_index := which(df_age$age_from <= age & df_age$age_to >= age), by = 'age']
  tmp <- pop_data1[, list(Total = sum(pop)), by = 'code']
  pop_data1 <- merge(pop_data1, tmp, by = 'code')
  pop_data1 = pop_data1[, list(pop = sum(pop)), by = c('Total', 'age_state_index', 'code')]
  setnames(pop_data1, 'age_state_index', 'age_index')
  pop_data1[, age := fouragegroups[age_index]]
  pop_data1[, pop_prop := pop / Total]
  tmp = pop_data1[, list(pop = sum(pop)), by = 'age']
  tmp[, pop_prop_US := pop / sum(pop_data$pop)]
  pop_data1 = merge(pop_data1, select(tmp, age, pop_prop_US), by = 'age')
  pop_data1 = subset(pop_data1, code %in% df_state$code)
  
  # week index
  df_week[, dummy := 1]
  date_10thcum[, dummy := 1]
  df_week = merge(df_week, date_10thcum, by = 'dummy', allow.cartesian=TRUE)
  df_week_10thcum = subset(df_week, date >= date_10thcum)
  data_comp = as.data.table( subset(data_comp, code %in% df_state$code))
  data_comp[, date := as.Date(date)]
  tmp2 = data_comp[, list(cum.death = sum(get(cum.death.var))), by = c('date', 'code')]
  tmp2 = tmp2[order(code, date)]
  last.cum.deaths = list(); tmp3 = list()
  for(Code in df_state$code){
    df_week_10thcum_m = subset(df_week_10thcum, code == Code)
    tmp2_m = subset(tmp2, code == Code)
    dummy = min(df_week_10thcum_m$date)  %in% tmp2_m$date
    last.cum.deaths[[Code]] = max(subset(tmp2_m, date == ifelse(dummy, min(df_week_10thcum_m$date) , 
                                                      (min(tmp2_m$date):(min(tmp2_m$date)+6))[min(tmp2_m$date):(min(tmp2_m$date)+6) %in% df_week_10thcum_m$date] ) )$cum.death )
    tmp3[[Code]] = unique(subset(tmp2, code == Code & date %in% c(df_week_10thcum_m$date, max(df_week_10thcum_m$date)+7)))
  }
  tmp2 = do.call('rbind', tmp3)
  tmp2[, weekly.deaths := c(diff(cum.death),NA), by = 'code']
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = c('date', 'code'))
  tmp2 = select(tmp2, code, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))
  tmp2[weekly.deaths <0, weekly.deaths := 0]
  
  # adjust dfweek for cum
  df_week_adj = unique(df_week[, list(date = date + 7), by = 'week_index'])
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples['alpha']) )
  setnames(tmp1,2:4, c('state_index', 'age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('state_index', 'week_index', 'age_state_index', 'iterations')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # find value from weeks before 10th date
  df_week1 = merge(df_week1, df_state, by = 'code')
  tmp3 = merge(tmp1, select(df_week1, week_index, state_index), by = c('state_index', 'week_index'))
  tmp3 = tmp3[, list(value_ref = mean(value)), by = c('state_index', 'age_index', 'iterations')]
  tmp3 = merge(tmp3, data.table(state_index = 1:length(last.cum.deaths), 
                         last.cum.deaths = unlist(last.cum.deaths) ), by = 'state_index')
  tmp3[, value_ref := rdirmnom(1, last.cum.deaths, value_ref), by =  c('state_index', 'iterations')]
  
  #find value from weeks after 10th date
  tmp2 = merge(tmp2, df_state, by = 'code')
  tmp1 = merge(tmp1, tmp2, by = c('week_index', 'state_index'))
  tmp1[, value := rdirmnom(1, weekly.deaths, value), by =  c('iterations', 'state_index', 'week_index')]
  
  # find cumulative death
  tmp1 = merge(tmp1, tmp3, by = c('age_index', 'iterations', 'state_index'))
  tmp1[, value := cumsum(value), by = c('iterations', 'state_index', 'age_index')]
  tmp1[, value := value + value_ref]
  
  # find mortality rate
  tmp1 = merge(tmp1, pop_data1, by = c('age_index', 'code'))
  tmp1[, value := value / pop]
  
  # save
  tmp1 = merge(tmp1, df_week_adj, by = 'week_index')
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', 'PosteriorSampleMortalityRate',lab, 'Table_', Code, '.rds'))
  }
  
  # find relative value to 85+
  tmp1[, value_maxage := value[age_index == max(age_index)], by = c('iterations', 'state_index', 'week_index')]
  tmp1[, value_rel := value_maxage / value]
  
  # find correlation coefficients 
  tmp3 <- merge(nyt_data, df_state[, .(loc_label, state_index)], by.x = 'STATE', by.y = 'loc_label')
  tmp3 <- merge(tmp1, tmp3[, .(state_index, SHARE_DEATHS)], by = 'state_index')
  tmp3 <- tmp3[!is.na(value_rel) & !is.na(SHARE_DEATHS)]
  tmp3[, value_corr := cor(value_rel, SHARE_DEATHS), by = c('week_index', 'age_index', 'iterations')]
  
  # quantile
  tmp2 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('state_index', 'week_index', 'age_index')]	
  tmp2 = dcast(tmp2, state_index + week_index + age_index ~ q_label, value.var = "q")
  
  if(nrow(tmp3) > 1){
    tmp3 = tmp3[, list(q= quantile(na.omit(value_corr), prob=ps, na.rm = T), q_label=paste0(p_labs, '_corr')), by=c('week_index', 'age_index')]	
    tmp3 = dcast.data.table(tmp3,  week_index + age_index ~ q_label, value.var = "q")
    tmp2 = merge(tmp2, tmp3, by=c('week_index', 'age_index'))
  }

  tmp1 = tmp1[, list(q= quantile(value_rel, prob=ps, na.rm = T), q_label=paste0(p_labs, '_rel')), by=c('state_index', 'week_index', 'age_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_index ~ q_label, value.var = "q")
  tmp1 = merge(tmp1, tmp2, by=c('state_index', 'week_index', 'age_index'))
  
  tmp1 = merge(tmp1, df_week_adj, by = 'week_index')
  tmp1[, age := df_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_age$age)]
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', 'MortalityRate', lab, 'Table_', Code, '.rds'))
    
  }

  return(tmp1)
}

make_mortality_rate_table_continuous = function(fit_samples, date_10thcum, df_week, pop_data, 
                                                data_comp, df_age_continuous, cum.death.var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # ref week
  df_week[, dummy := 1]
  date_10thcum[, dummy := 1]
  df_week1 = merge(df_week, date_10thcum, by = 'dummy', allow.cartesian = T)
  df_week1 = df_week1[date >= date_10thcum]
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1[, keep := month %in% min(df_week1$month):(min(df_week1$month) + 2), by = 'code']
  df_week1 = df_week1[keep == T]
  
  # group 85+ together
  df_age <- copy(df_age_continuous)
  df_age[age >= 85, age_index := age_index[age == 85]]
  
  # find population proportion
  pop_data1 = copy(pop_data)
  tmp <- pop_data1[, list(Total = sum(pop)), by = 'code']
  pop_data1 <- merge(pop_data1, tmp, by = 'code')
  pop_data1[, pop_prop := pop / Total]
  pop_data1 <- merge(pop_data1, df_age, by = 'age')
  pop_data1 <- pop_data1[, list(pop = sum(pop), 
                                age = min(age)), by = c('age_index', 'code')]
  tmp = pop_data1[, list(pop = sum(pop)), by = 'age']
  tmp[, pop_prop_US := pop / sum(pop_data$pop)]
  pop_data1 = merge(pop_data1, select(tmp, age, pop_prop_US), by = 'age')
  pop_data1 = subset(pop_data1, code %in% df_state$code)
  df_age <- df_age[, list(age_from = min(age_from), age_to = max(age_to), age = min(age)), by = 'age_index']
  
  # week index
  df_week[, dummy := 1]
  date_10thcum[, dummy := 1]
  df_week = merge(df_week, date_10thcum, by = 'dummy', allow.cartesian=TRUE)
  df_week_10thcum = subset(df_week, date >= date_10thcum)
  data_comp = as.data.table( subset(data_comp, code %in% df_state$code))
  data_comp[, date := as.Date(date)]
  tmp2 = data_comp[, list(cum.death = sum(get(cum.death.var))), by = c('date', 'code')]
  tmp2 = tmp2[order(code, date)]
  last.cum.deaths = list(); tmp3 = list()
  for(Code in df_state$code){
    df_week_10thcum_m = subset(df_week_10thcum, code == Code)
    tmp2_m = subset(tmp2, code == Code)
    dummy = min(df_week_10thcum_m$date)  %in% tmp2_m$date
    last.cum.deaths[[Code]] = max(subset(tmp2_m, date == ifelse(dummy, min(df_week_10thcum_m$date) , 
                                                                (min(tmp2_m$date):(min(tmp2_m$date)+6))[min(tmp2_m$date):(min(tmp2_m$date)+6) %in% df_week_10thcum_m$date] ) )$cum.death )
    tmp3[[Code]] = unique(subset(tmp2, code == Code & date %in% c(df_week_10thcum_m$date, max(df_week_10thcum_m$date)+7)))
  }
  tmp2 = do.call('rbind', tmp3)
  tmp2[, weekly.deaths := c(diff(cum.death),NA), by = 'code']
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = c('date', 'code'))
  tmp2 = select(tmp2, code, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))
  tmp2[weekly.deaths <0, weekly.deaths := 0]
  
  # adjust dfweek for cum
  df_week_adj = unique(df_week[, list(date = date + 7), by = 'week_index'])
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples['alpha']) )
  setnames(tmp1,2:4, c('state_index', 'age_index','week_index'))

  # aggregate 85+
  tmp1[age_index > 85, age_index := 86]
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'state_index', 'age_index','week_index')]
  
  # find value from weeks before 10th date
  df_week1 = merge(df_week1, df_state, by = 'code')
  tmp3 = merge(tmp1, select(df_week1, week_index, state_index), by = c('state_index', 'week_index'))
  tmp3 = tmp3[, list(value_ref = mean(value)), by = c('state_index', 'age_index', 'iterations')]
  tmp3 = merge(tmp3, data.table(state_index = 1:length(last.cum.deaths), 
                                last.cum.deaths = unlist(last.cum.deaths) ), by = 'state_index')
  tmp3[, value_ref := rdirmnom(1, last.cum.deaths, value_ref), by =  c('state_index', 'iterations')]
  
  #find value from weeks after 10th date
  tmp2 = merge(tmp2, df_state, by = 'code')
  tmp1 = merge(tmp1, tmp2, by = c('week_index', 'state_index'))
  tmp1[, value := rdirmnom(1, weekly.deaths, value), by =  c('iterations', 'state_index', 'week_index')]
  
  # find cumulative death
  tmp1 = merge(tmp1, tmp3, by = c('age_index', 'iterations', 'state_index'))
  tmp1[, value := cumsum(value), by = c('iterations', 'state_index', 'age_index')]
  tmp1[, value := value + value_ref]
  
  # find mortality rate
  tmp1 = merge(tmp1, pop_data1, by = c('age_index', 'code'))
  tmp1[, value := value / pop]
  
  # quantiles
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('state_index', 'week_index', 'age_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week_adj, by = 'week_index')

  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  tmp1 = merge(tmp1, df_age, by = 'age_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', 'MortalityRateContinuous', 'Table_', Code, '.rds'))
  }
  
  return(tmp1)
}

make_weekly_death_rate_other_source = function(fit_samples, df_week, data_comp, var.alpha, df_age, outdir, age_groups = NULL, 
                                               lab = NULL, withempirical = F){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(!is.null(age_groups)){
    # df age
    df_age = data.table(age = age_groups)
    df_age[, age_index := 1:nrow(df_age)]
    df_age[, age_from := gsub('(.+)-.*', '\\1', age)]
    df_age[, age_to := gsub('.*-(.+)', '\\1', age)]
    df_age[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
    df_age[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
    df_age[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
    df_age[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
    set(df_age, NULL, 'age_from', df_age[,as.numeric(age_from)])
    set(df_age, NULL, 'age_to', df_age[,as.numeric(age_to)])
    set(df_age, NULL, 'age_to_index', df_age[,as.numeric(age_to_index)])
    set(df_age, NULL, 'age_from_index', df_age[,as.numeric(age_from_index)])
    df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
    df_age_continuous[, age_state_index := which(df_age$age_from_index <= age_index & df_age$age_to_index >= age_index), by = 'age_index']
    
  }
  
  # ref data
  data_comp = as.data.table(data_comp)
  data_comp = subset(data_comp, code %in% df_state$code)
  data_comp[, date := as.Date(date)]
  tmp2 = data_comp[, list(daily_deaths = sum(daily_deaths)), by = c('date', 'code')]
  tmp2 = tmp2[order(code, date)]
  tmp2[, cum.death := cumsum(daily_deaths), by = 'code']
  tmp2 = unique(subset(tmp2, date %in% c(df_week$date, max(df_week$date)+7)))
  tmp2[, weekly.deaths := c(diff(cum.death),NA), by = 'code']
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = 'date')
  tmp2 = select(tmp2, code, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))
  tmp2[weekly.deaths<0, weekly.deaths := 0]
  tmp2 = merge(tmp2, df_state, by = 'code')
  
  if(withempirical){
    
    if(!is.null(age_groups)){
      df_age_reporting[, age_state_index := which(df_age$age_from <= age_from & df_age$age_to >= age_to), by = 'age_index']
      empirical = merge(data, df_age_reporting, by = 'age_from')
      empirical = merge(empirical, df_week, by = 'date')
      empirical = empirical[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c('code', 'age_state_index', 'week_index')]
      setnames(empirical, 'age_state_index', 'age_index')
    } else{
      empirical = merge(data, df_age_reporting, by = 'age_from')
      empirical = merge(empirical, df_week, by = 'date')
    }
    setnames(empirical, 'weekly.deaths', 'emp')
    tmp_emp = empirical[, list(total_deaths = sum(na.omit(emp))), by = c('code', 'week_index')]
    empirical = merge(empirical, tmp_emp, by = c('code', 'week_index'))
    empirical[, prop_deaths := emp / total_deaths]
    empirical = merge(empirical, tmp2, by = c('code','week_index'))
    empirical[, emp_adj := weekly.deaths * prop_deaths]
    empirical = select(empirical, emp, emp_adj, code, week_index, age_index)
  }
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var.alpha]]) )
  setnames(tmp1,2:4, c('state_index', 'age_index','week_index'))
  
  # sum by state age group
  if(!is.null(age_groups)){
    tmp1 = merge(tmp1, df_age_continuous, 'age_index')
    tmp1 = tmp1[, list(value = sum(value)), by = c('state_index', 'week_index', 'age_state_index', 'iterations')]
    setnames(tmp1, 'age_state_index', 'age_index')
  }
  
  # prediction weekly deaths
  tmp1 = merge(tmp1, tmp2, by = c('week_index', 'state_index'))
  tmp1[, value := rdirmnom(1, weekly.deaths, value), by = c('state_index', 'week_index', 'iterations')]
  
  # quantiles
  tmp1 = tmp1[, list(q= c(quantile(value, prob=ps, na.rm = T), mean(value)), q_label=c(p_labs, 'mean')), 
              by=c('state_index', 'week_index', 'age_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')

  tmp1[, age := df_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_age$age)]

  tmp1 = merge(tmp1, tmp2, by = c('week_index', 'state_index'), all.x = T)
  setnames(tmp1, 'weekly.deaths', 'emp_JHU')
  
  if(withempirical){
    tmp1 = merge(tmp1, empirical, by = c('code', 'week_index', 'age_index'), all.x = T)
  }
  
  for(Code in unique(tmp1$code)){
    if(!is.null(lab)){
      file = paste0(outdir, '-', 'DeathByAge', 'Table_', lab, '_', Code, '.rds')
      
    }else{
      file =  paste0(outdir, '-', 'DeathByAge', 'Table_', var.alpha, '_', Code, '.rds')
    }
    
    saveRDS(subset(tmp1, code == Code), file =file)
    
  }

  
  return(tmp1)
}

summary_death_all_states = function(death2, rm_states){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  death2 =  subset(death2, !code %in% rm_states)
  
  # sum across codes
  death2 = death2[, list(week1 = sum(week1), week2 = sum(week2)), by = c('age_index', 'age', 'iterations')]
  
  # find proportion
  death2[, value := week2 / week1]
  
  # quantiles
  tmp1 = death2[, list( 	q= quantile(value, prob=ps, na.rm = T),
                         q_label=p_labs), 
                by=c('age', 'age_index')]	
  tmp1 = dcast(tmp1, age + age_index ~ q_label, value.var = "q")
  
  tmp2 = death2[, list( 	q= quantile(week1, prob=ps, na.rm = T),
                         q_label=p_labs), 
                by=c('age', 'age_index')]	
  tmp2 = dcast(tmp2, age + age_index ~ q_label, value.var = "q")
  setnames(tmp2, p_labs, paste0(p_labs, '_week1'))
  tmp1= merge(tmp1, tmp2, by=c('age', 'age_index'))
  
  death2 = death2[, list( 	q= quantile(week2, prob=ps, na.rm = T),
                           q_label=p_labs), 
                  by=c('age', 'age_index')]	
  death2 = dcast(death2, age + age_index ~ q_label, value.var = "q")
  setnames(death2, p_labs, paste0(p_labs, '_week2'))
  tmp1= merge(tmp1, death2, by=c('age', 'age_index'))
  
  return(tmp1)
}

make_lambda_table <- function(fit_samples, stan_data, df_week, df_state, outdir){
  
  ps <- c(0.5, 0.1, 0.9)
  p_labs <- c('M','CL','CU')
  
  tmp <- as.data.table(reshape2::melt(fit_samples[['lambda_raw']]))
  setnames(tmp, 2:3, c('state_index', 'week_index_womissing'))
  tmp = tmp[, list(q= quantile(value, prob=ps, na.rm = T),
                           q_label=p_labs), 
                  by=c('state_index', 'week_index_womissing')]	
  tmp = dcast(tmp, state_index + week_index_womissing ~ q_label, value.var = "q")
  tmp[, type := 'posterior']
  
  tmp1 <- lapply(1:length(stan_data$lambda_prior_parameters), function(i){
    tmp <- as.data.table(stan_data$lambda_prior_parameters[[i]])
    setnames(tmp, 'V1', 'value')
    tmp[, week_index_womissing := 1:nrow(tmp)]
    tmp[, state_index := i]
    tmp
  })
  tmp1 <- do.call('rbind', tmp1)
  
  
  if(with_cmdstan){
    model_code = 'sensitivity analysis 1 on lambda'
  }else{
    model_code <- model@model_code
  }
  
  if(!grepl('sensitivity analysis 1 on lambda|sensitivity analysis 2 on lambda', model_code)){

    tmp1 <- tmp1[, list( 	q= qexp(ps, 1/value),
                          q_label=p_labs), 
                 by=c('state_index', 'week_index_womissing')]	
  }
  
  if(grepl('sensitivity analysis 1 on lambda', model_code)){

    tmp1 <- tmp1[, list( 	q= extraDistr::qtnorm(ps, mean = value, sd = 100, a = 0),
                          q_label=p_labs), 
                 by=c('state_index', 'week_index_womissing')]	
  }

  tmp1 = dcast(tmp1, state_index + week_index_womissing ~ q_label, value.var = "q")
  tmp1[, type := 'prior']
  
  tmp <- rbind(tmp, tmp1)
  tmp[, type :=factor(type, levels = c('prior', 'posterior'))]
  
  df_week_womissing <- copy(df_week[stan_data$IDX_WEEKS_OBSERVED,])
  df_week_womissing[, week_index_womissing:= 1:nrow(df_week_womissing)]
  tmp <- merge(tmp, df_week_womissing, by = 'week_index_womissing')
  tmp <- merge(tmp, df_state, by = 'state_index')
  
  for(Code in unique(tmp$code)){
    saveRDS(subset(tmp, code == Code), paste0(outdir, '-lambda_prior_posterior_', Code, '.rds'))
  }

  saveRDS(tmp, paste0(outdir, '-lambda_prior_posterior.rds'))
  return(tmp)
}

make_var_base_model_table <- function(fit_samples, stan_data, df_state, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  .math_name = c('nu', 'gamma[1]', 'gamma[2]', 'zeta')
  variable_name <- c('nu', 'gamma_gp1', 'gamma_gp2', 'zeta_gp')
  df_variable_name <- data.table(variable_name = variable_name, math_name = .math_name)

  tmp <- lapply(variable_name, function(x){
    tmp <- as.data.table(reshape2::melt(fit_samples[[x]]))
    setnames(tmp, 2, c('state_index'))
    tmp[, variable_name := x]
    })
  tmp <- do.call('rbind', tmp)
  tmp = tmp[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('state_index', 'variable_name')]
  tmp = dcast(tmp, state_index + variable_name ~ q_label, value.var = "q")
  tmp[, type := 'posterior']
  
  if(with_cmdstan){
    model_code <- 'sensitivity analysis 1 on gamma,sensitivity analysis 1 on zeta,sensitivity analysis 1 on nu'
  }else{
    model_code <- model@model_code
  }
  
  # prior gamma
  tmp1 <- data.table(expand.grid(state_index = df_state$state_index, variable_name = c('gamma_gp1', 'gamma_gp2')))
  if(!grepl('sensitivity analysis 1 on gamma|sensitivity analysis 2 on gamma', model_code)){
    tmp1 <- tmp1[, list(q= invgamma::qinvgamma(ps, 2, 2), q_label=p_labs),by=c('state_index', 'variable_name')]	
  }
  if(grepl('sensitivity analysis 1 on gamma', model_code)){
    tmp1 <- tmp1[, list(q= invgamma::qinvgamma(ps, 5, 5), q_label=p_labs),by=c('state_index', 'variable_name')]	
  }
  if(grepl('sensitivity analysis 2 on gamma', model_code)){
    tmp1 <- tmp1[, list(q= extraDistr::qtnorm(ps, mean = 0, sd = 5, a = 0), q_label=p_labs),by=c('state_index', 'variable_name')]	
  }
  tmp1 = dcast(tmp1, state_index + variable_name ~ q_label, value.var = "q")
  
  # prior zeta
  tmp2 <- data.table(state_index = df_state$state_index, variable_name = 'zeta_gp')
  if(!grepl('sensitivity analysis 1 on zeta|sensitivity analysis 2 on zeta', model_code)){
    tmp2 <- tmp2[, list(q= extraDistr::qhcauchy(ps,  1), q_label=p_labs), by=c('state_index', 'variable_name')]	
  }
  if(grepl('sensitivity analysis 1 on zeta', model_code)){
    tmp2 <- tmp2[, list(q= extraDistr::qhcauchy(ps,  10), q_label=p_labs),by=c('state_index', 'variable_name')]	
  }
  if(grepl('sensitivity analysis 2 on zeta', model_code)){
    tmp2 <- tmp2[, list(q= extraDistr::qinvchisq(ps,  1), q_label=p_labs),by=c('state_index', 'variable_name')]	
  }
  tmp2 = dcast(tmp2, state_index + variable_name ~ q_label, value.var = "q")
  tmp1 <- rbind(tmp1, tmp2)
  
  # prior nu 
  tmp2 <- data.table(state_index = df_state$state_index, variable_name = 'nu')
  if(!grepl('sensitivity analysis 1 on nu|sensitivity analysis 2 on nu', model_code)){
    tmp2 <- tmp2[, list(q= qexp(ps, 1), q_label=p_labs),by=c('state_index', 'variable_name')]	
  } 
  if(grepl('sensitivity analysis 1 on nu', model_code)){
    tmp2 <- data.table(expand.grid(state_index = df_state$state_index, variable_name = 'nu',
                                   samples_nu_inverse = truncnorm::rtruncnorm(10000, a=0,  mean = 0, sd = 5)))
    tmp2[, samples_nu := (1/samples_nu_inverse)]
    tmp2 <- tmp2[, list(q= quantile(samples_nu, prob=ps, na.rm = T), q_label=p_labs), by=c('state_index', 'variable_name')]	
  }
  if(grepl('sensitivity analysis 2 on nu', model_code)){
    tmp2 <- tmp2[, list(q= qexp(ps, 10), q_label=p_labs),by=c('state_index', 'variable_name')]	
  }
  if(grepl('sensitivity analysis 3 on nu', model_code)){
    tmp2 <- tmp2[, list(q=  extraDistr::qhcauchy(ps,  1), q_label=p_labs),by=c('state_index', 'variable_name')]	
  }
  tmp2 = dcast(tmp2, state_index + variable_name ~ q_label, value.var = "q")
  tmp1 <- rbind(tmp1, tmp2)
  tmp1[, type := 'prior']
  
  tmp <- rbind(tmp, tmp1)
  tmp[, type :=factor(type, levels = c('prior', 'posterior'))]
  
  tmp <- merge(tmp, df_variable_name, by = 'variable_name')
  tmp <- merge(tmp, df_state, by = 'state_index')
  
  tmp[, math_name := factor(math_name, levels = .math_name)]
  # tmp[, math_name := paste0(math_name, '\\["', loc_label, '"\\]')]
  
  for(Code in unique(tmp$code)){
    saveRDS(subset(tmp, code == Code), paste0(outdir, '-var_base_model_prior_posterior_', Code, '.rds'))
  }
  
  saveRDS(tmp, paste0(outdir, '-var_base_model_prior_posterior.rds'))
  
  return(tmp)
  
}

find_mortality_rate_across_states <- function(mortality_rate_posterior_samples){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  df_pop <- unique(mortality_rate_posterior_samples[, .(code, pop, age)])
  df_pop[, prop_state := pop / sum(pop), by = c('age')]
  
  tmp = subset(mortality_rate_posterior_samples, date == max(mortality_rate_posterior_samples$date)-7)
  tmp[, death := value * pop]
  tmp <- merge(tmp, df_pop, by = c('code', 'pop', 'age'))
  tmp[, mortality_rate := death / pop]
  
  tmp <- tmp[, list(value = sum(mortality_rate * prop_state)), by = c('iterations', 'age')]
  
  tmp = tmp[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), 
            by=c('age')]	
  tmp = dcast(tmp, age ~ q_label, value.var = "q")
  
  return(tmp)
}

prepare_prop_vac_table <- function(stan_data, vaccine_data, df_age_vaccination){
  
  delay = 2*7
  age_index_min = df_age_vaccination[age == '18-64']$age_index
  vaccine_data[, age_index := which(df_age_vaccination$age_from <= age & df_age_vaccination$age_to >= age), by = 'age']
  tmp = vaccine_data[age_index >= age_index_min, list(prop = unique(prop)), by = c('code', 'date', 'loc_label', 'age_index')]
  tmp[, date := date + delay]
  tmp = merge(tmp, resurgence_dates, by = 'code')
  tmp = tmp[code %in% Code]
  tmp = tmp[date %in% df_week$date]
  # resurgence period
  tmp = tmp[date >= start_resurgence & date <= stop_resurgence]
  tmp[, idx.resurgence.date := 1:length(date), by = c('age_index', 'loc_label')]
  # set proportion
  tmp1 = tmp[, list(pre_prop = prop[date == start_resurgence]), by = c('loc_label', 'age_index')]
  tmp = merge(tmp, tmp1,  by = c('loc_label', 'age_index'))
  
  tmp[, cat := paste0('prop_', age_index - 2)]
  tmp = data.table(reshape2::dcast(tmp, code + date ~ cat, value.var = 'prop'))
  
  return(tmp)
}

prepare_prop_vac_counterfactual_table <- function(stan_data, df_state, df_age_vaccination2, df_counterfactual){
  tmp <- lapply(1:length(stan_data$prop_vac_start_counterfactual), function(x){
    tmp <- as.data.table( reshape2::melt(stan_data$prop_vac_start_counterfactual[[x]]) )
    tmp$age_index <- x
    tmp
  })
  tmp <- do.call('rbind', tmp)
  setnames(tmp, 1:2, c('counterfactual_index', 'state_index'))
  
  tmp1 <- lapply(1:length(stan_data$prop_vac_start), function(x){
    tmp <- as.data.table( reshape2::melt(stan_data$prop_vac_start[[x]]) )
    tmp$state_index <- 1:nrow(tmp)
    tmp$age_index <- x
    tmp
  })
  tmp1 <- do.call('rbind', tmp1)
  setnames(tmp1, 'value', 'value_true')
  
  tmp <- merge(tmp, tmp1, by = c('age_index', 'state_index'))
  tmp[, diff_value := value - value_true]
  
  # tmp <- merge(tmp, df_state, by = 'state_index')
  # tmp <- merge(tmp, df_counterfactual, by = 'counterfactual_index')
  # tmp <- merge(tmp, df_age_vaccination2, by = 'age_index')
  
  return(tmp)
}

make_var_by_age_by_state_by_counterfactual_table = function(fit_samples, df_week, df_state_age, df_state, df_counterfactual, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'state_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  tmp1 = merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'Table_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

make_inv_var_by_age_by_state_by_counterfactual_table = function(fit_samples, df_week, df_state_age, df_state, df_counterfactual, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  tmp1[, value := -value]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'state_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  tmp1 = merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-inv_', var_name,  'Table_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

make_var_by_age_by_counterfactual_table = function(fit_samples, df_week, df_state_age, df_counterfactual, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:4, c('counterfactual_index', 'age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1,counterfactual_index +week_index + age_index ~ q_label, value.var = "q")
  # tmp1 = merge(tmp1, df_week, by = 'week_index')
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  tmp1 <- merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  saveRDS(tmp1, file = paste0(outdir, '-', var_name,  'AllStatesTable.rds'))
  
  
  return(tmp1)
}

make_inv_var_by_age_by_counterfactual_table = function(fit_samples, df_week, df_state_age, df_counterfactual, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  
  tmp1[, value := - value]
  
  setnames(tmp1, 2:4, c('counterfactual_index', 'age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1,counterfactual_index +week_index + age_index ~ q_label, value.var = "q")
  # tmp1 = merge(tmp1, df_week, by = 'week_index')
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  tmp1 <- merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  saveRDS(tmp1, file = paste0(outdir, '-', var_name,  'AllStatesTable.rds'))
  
  
  return(tmp1)
}

make_ratio_vars_by_age_by_counterfactual_table = function(fit_samples, df_state_age, df_counterfactual, vars_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[vars_name[1]]]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'counterfactual_index', 'age_index','week_index')]
  
  tmp2 = as.data.table( reshape2::melt(fit_samples[[vars_name[2]]]) )
  setnames(tmp2, 2:5, c('state_index', 'age_index','week_index', 'value_denominator'))
  tmp2 <- tmp2[, list(value_denominator = sum(value_denominator)), by = c('iterations', 'age_index','week_index')]
  
  tmp1 <- merge(tmp1, tmp2,  by = c('iterations', 'age_index','week_index'))
  tmp1[, value := value / value_denominator]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + age_index + week_index ~ q_label, value.var = "q")
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  tmp1 = merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', vars_name[1], 'AllStatesRatioTable_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

make_var_inv_by_state_by_counterfactual_table = function(fit_samples, df_week, df_state, df_counterfactual, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'counterfactual_index', 'state_index','week_index')]
  
  tmp1[, value := - value]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'state_index', 'week_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + state_index + week_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  tmp1 = merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'AllAgesTable_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

make_var_by_state_by_counterfactual_table = function(fit_samples, df_week, df_state, df_counterfactual, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'counterfactual_index', 'state_index','week_index')]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'state_index', 'week_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + state_index + week_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  tmp1 = merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'AllAgesTable_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

make_ratio_vars_by_state_by_counterfactual_table = function(fit_samples, df_week, df_state, df_counterfactual, vars_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[vars_name[1]]]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'counterfactual_index', 'state_index','week_index')]
  
  tmp2 = as.data.table( reshape2::melt(fit_samples[[vars_name[2]]]) )
  setnames(tmp2, 2:5, c('state_index', 'age_index','week_index', 'value_denominator'))
  tmp2 <- tmp2[, list(value_denominator = sum(value_denominator)), by = c('iterations', 'state_index','week_index')]
  
  tmp1 <- merge(tmp1, tmp2,  by = c('iterations', 'state_index','week_index'))
  tmp1[, value := (value - value_denominator) / value_denominator]
  # tmp1[, value := value  / value_denominator]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'state_index', 'week_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + state_index + week_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  tmp1 = merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', vars_name[1], 'RatioAllAgesTable_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

make_inv_var_by_counterfactual_table = function(fit_samples, df_counterfactual, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:4, c('counterfactual_index', 'age_index','week_index'))
  
  tmp1 <- tmp1[, list(value = sum(value)), by = c('counterfactual_index', 'iterations','week_index')]
  
  tmp1[, value := - value]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'week_index')]	
  tmp1 = dcast(tmp1,counterfactual_index +week_index  ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  saveRDS(tmp1, file = paste0(outdir, '-', var_name,  'AllStatesAllAgesTable.rds'))
  
  
  return(tmp1)
}

make_ratio_vars_by_counterfactual_table = function(fit_samples, df_counterfactual, vars_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[vars_name[1]]]) )
  setnames(tmp1, 2:5, c('counterfactual_index', 'state_index', 'age_index','week_index'))
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'counterfactual_index','week_index')]
  
  tmp2 = as.data.table( reshape2::melt(fit_samples[[vars_name[2]]]) )
  setnames(tmp2, 2:5, c('state_index', 'age_index','week_index', 'value_denominator'))
  tmp2 <- tmp2[, list(value_denominator = sum(value_denominator)), by = c('iterations','week_index')]
  
  tmp1 <- merge(tmp1, tmp2,  by = c('iterations','week_index'))
  tmp1[, value := value / value_denominator]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('counterfactual_index', 'week_index')]	
  tmp1 = dcast(tmp1, counterfactual_index + week_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_counterfactual, by = 'counterfactual_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', vars_name[1], 'RatioAllAgesAllStatesTable_', Code, '.rds'))
    
  }
  
  return(tmp1)
}

find_p_value_vaccine_effect <- function(sample, name){
  
  tmp <- as.data.table( reshape2::melt(sample) )
  tmp[, value := value >= 0]
  
  if(ncol(tmp) == 2){
    tmp = tmp[, list(M= mean(value)*100)]	
    tmp[, name := name]
  } else if(ncol(tmp) == 3){
    tmp = tmp[, list(M= mean(value)*100),
              by=c('Var2')]	
    tmp[, name := paste0(name, '_', Var2)]
    
  }
  
  tmp <- tmp[, .(name, M)]
  return(tmp)
}

find_crude_mortality_rate <- function(mortality_rate, df_age_continuous, df_age_reporting, pop_data){
  
  fouragegroups <- unique(mortality_rate$age)
  
  # df age
  df_age = data.table(age = fouragegroups)
  df_age[, age_index := 1:nrow(df_age)]
  df_age[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  df_age[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_age[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  set(df_age, NULL, 'age_from', df_age[,as.numeric(age_from)])
  set(df_age, NULL, 'age_to', df_age[,as.numeric(age_to)])
  set(df_age, NULL, 'age_to_index', df_age[,as.numeric(age_to_index)])
  set(df_age, NULL, 'age_from_index', df_age[,as.numeric(age_from_index)])
  df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
  df_age_continuous[, age_state_index := which(df_age$age_from_index <= age_index & df_age$age_to_index >= age_index), by = 'age_index']
  df_age_reporting[, age_state_index := unique(df_age_continuous$age_state_index[which(df_age_continuous$age_from == age_from)]), by = 'age_index']
  
  # find population proportion
  pop_data1 = copy(pop_data)
  pop_data1[, age_state_index := which(df_age$age_from <= age & df_age$age_to >= age), by = 'age']
  tmp <- pop_data1[, list(Total = sum(pop)), by = 'code']
  pop_data1 <- merge(pop_data1, tmp, by = 'code')
  pop_data1 = pop_data1[, list(pop = sum(pop)), by = c('Total', 'age_state_index', 'code')]
  setnames(pop_data1, 'age_state_index', 'age_index')
  pop_data1[, age := fouragegroups[age_index]]
  pop_data1[, pop_prop := pop / Total]
  tmp = pop_data1[, list(pop = sum(pop)), by = 'age']
  tmp[, pop_prop_US := pop / sum(pop_data$pop)]
  pop_data1 = merge(pop_data1, select(tmp, age, pop_prop_US), by = 'age')
  
  # find mortality
  file = file.path(dirname(indir), paste0('misc/data/CDC_data_', unique(mortality_rate$date) + 13,'.csv'))
  tmp <- as.data.table(read.csv(file))
  tmp[, date := as.Date(End.Date, format = '%m/%d/%Y')]
  tmp <- tmp[date == unique(mortality_rate$date) + 7]
  tmp <- tmp[Sex == 'All Sexes']
  tmp <- tmp[Group == 'By Total']
  tmp[, age := ifelse(Age.Group == "Under 1 year", "0-0", 
                      ifelse(Age.Group == "85 years and over", "85+", gsub("(.+) years", "\\1", Age.Group)))]
  
  # state
  setnames(tmp, 'State', 'loc_label')
  tmp <- merge(tmp, unique(mortality_rate[, .(loc_label, code)]))

  # merge age
  tmp <- merge(tmp, df_age_reporting, by = 'age')
  tmp <- tmp[, list(total_deaths = sum(na.omit(COVID.19.Deaths))), by = c('loc_label', 'age_state_index', 'code')]
  setnames(tmp, 'age_state_index', 'age_index')
  
  # merge pop and find mortality rate
  tmp <- merge(tmp, pop_data1, by = c('code', 'age_index'))
  tmp[, crude_mortality_rate := total_deaths / pop]
  
  return(tmp)
}


make_forest_plot_table <- function(summary, df_age_vaccination2, df_state, names, math_name, groups, groups_levels){
  
  tmp <- summary[ grepl(paste(paste0('^',names),collapse = '|'),rownames(summary)) ,]
  # tmp1 <- tmp[grepl('diagonal', rownames(tmp)),]
  # tmp <- tmp[!grepl('diagonal', rownames(tmp)),]
  
  variables <- rownames(tmp); group = c()
  for(x in 1:length(variables)) group[x] = groups[which(grepl(gsub('(.+)\\[.*', '\\1', variables[x]), names))]
  for(x in 1:length(names)) variables = gsub(names[x], math_name[x], variables)
  for(x in df_age_vaccination2$age_index) variables = gsub(paste0('\\[',x), paste0('\\["', df_age_vaccination2$age[x], '"'), variables) 
  # for(x in df_age_vaccination2$age_index) variables[grepl('vac', variables) & !grepl(', ', variables)] = sub(paste0('\\[\"',df_age_vaccination2$age[x]), paste0('\\[\"', df_age_vaccination2$age[-x],', ',df_age_vaccination2$age[x]), variables[grepl('vac', variables) & !grepl(', ', variables)])
  # for(x in df_age_vaccination2$age_index) variables[grepl('vac', variables)] = gsub('(.+),.*', '\\1"\\]', variables[grepl('vac', variables)])
  
  for(x in df_state$state_index) variables[!grepl('vac', variables)] = gsub(paste0('\",', x, '\\]'), paste0(', ',df_state$loc_label[x], '"\\]'), variables[!grepl('vac', variables)]) 
  
  variables = gsub('\\+\\+', '\\+', variables)
  variables = gsub('\\+', 'plus', variables)
  
  age1 = gsub('(.+),.*', '\\1', gsub('.*\\["(.+)','\\1', variables))
  loc1 = gsub('.*, (.+)"\\]','\\1', variables)
  idx.swap = which(!grepl('\\]', age1) & !grepl('[0-9]', loc1))
  
  for(x in 1:length(idx.swap)){
    
    variables[idx.swap[x]] = gsub(age1[idx.swap[x]],paste0(age1[idx.swap[x]], '1'), variables[idx.swap[x]])
    variables[idx.swap[x]] = gsub(loc1[idx.swap[x]],age1[idx.swap[x]], variables[idx.swap[x]])
    variables[idx.swap[x]] = gsub(paste0(age1[idx.swap[x]], '1'),loc1[idx.swap[x]], variables[idx.swap[x]])
    
  }
  
  variables = gsub('plus', '+', variables)
  
  tmp <- as.data.table(tmp)
  tmp[, variable := variables]
  tmp[, group := factor(group, levels = groups_levels)]
  
  setnames(tmp, c('50%', '2.5%', '97.5%'), c('M', 'CL', "CU"))
  
}

make_var_by_age_by_state_by_time_diff_table = function(fit_samples, df_week, df_state_age, df_state, var_name, vaccine_data, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # delay = 7*2
  vaccine_data_pop = vaccine_data_pop[, list(prop = unique(prop)), by = c('code', 'date', 'loc_label')]
  vaccine_data_pop <- vaccine_data_pop[date %in% df_week[, unique(date)]]
  tmp <- vaccine_data_pop[, list(mindate = min(date[prop > 0.05])), by = 'code']
  tmp1 <- merge(tmp, df_state, by = 'code')
  
  tmp = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp, 2:4, c('state_index', 'age_index','week_index'))
  
  tmp1 <- merge(tmp, tmp1, by = 'state_index')
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  # difference value 1
  tmp1 <- tmp1[, list(diff1 = value[date == mindate] - value[date == min(date)], 
                      diff2 =value[date == (mindate + 7*8)] -  value[date == mindate]), by = c('state_index', 'iterations', 'age_index')]
  tmp1 <- melt.data.table(tmp1, id.vars = c('state_index', 'iterations', 'age_index'))
  
  tmp2 <- tmp1[, list(value = mean(na.omit(value))), by = c('iterations', 'age_index', 'variable')]

  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), 
              by=c('state_index', 'age_index', 'variable')]	
  tmp1 = dcast(tmp1, variable + state_index + age_index ~ q_label, value.var = "q")
  tmp1 = merge(tmp1, df_state, by = 'state_index')

  tmp2 = tmp2[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), 
              by=c('age_index', 'variable')]	
  tmp2 = dcast(tmp2, variable + age_index ~ q_label, value.var = "q")
  tmp2[, state_index := 0]
  tmp2[, code := 'US']
  tmp2[, loc_label := 'United States']
  
  tmp1 <- rbind(tmp1, tmp2)
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'DiffTable_', Code, '.rds'))
  }
  
  return(tmp1)
}

make_var_by_age_by_state_by_time_diff_over_time_table = function(fit_samples, df_week, df_state_age, df_state, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:4, c('state_index', 'age_index','week_index'))
  
  # difference value 1
  tmp1[, baselinevalue := value[week_index == min(week_index)], by = c('iterations', 'age_index', 'state_index')]
  tmp1[, value := value - baselinevalue]
  
  tmp2 <- tmp1[, list(value = mean(na.omit(value))), by = c('iterations', 'age_index', 'week_index')]
  
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), 
              by=c('state_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1, week_index + state_index + age_index ~ q_label, value.var = "q")
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  tmp2 = tmp2[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), 
              by=c('age_index', 'week_index')]	
  tmp2 = dcast(tmp2, week_index + age_index ~ q_label, value.var = "q")
  tmp2[, state_index := 0]
  tmp2[, code := 'US']
  tmp2[, loc_label := 'United States']
  
  tmp1 <- rbind(tmp1, tmp2)
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, df_week, by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, df_week, by = 'week_index')
  }
  
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name,  'DiffOverTimeTable_', Code, '.rds'))
  }
  
  return(tmp1)
}

