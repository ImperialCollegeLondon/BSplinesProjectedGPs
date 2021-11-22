make_predictive_checks_table = function(fit, df_week, df_state_age, data, deaths_predict_var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
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

  dir = gsub('(.+)/results/.*', '\\1', outdir.table), 'results/', 'predictions')
  dir.create(dir)
  saveRDS(tmp1, file = file.path(dir, 'predicted_weekly_deaths.rds'))
  
  return(tmp1)
}

make_convergence_diagnostics_stats = function(fit, outdir)
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
  re = rstan::extract(fit)
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
  
  return(summary)
}

make_probability_ratio_table = function(fit, df_week, df_state_age, data, stan_data, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples$probability_ratio_age_strata) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index', 'week_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  tmp1[, code := Code]
  
  # find significant shift
  tmp1[, shift := F]
  tmp1[CL > 1, shift := T]
  tmp1[CU < 1, shift := T]
  tmp2 = tmp1[, list(shift.tot = (sum(shift) > 0)), by = 'age']
  tmp1 = merge(tmp1, tmp2, by = 'age')
  
  # find empirical estimate
  tmp = select(data, weekly.deaths, date, age)
  tmp2 = tmp[, list(total.deaths = sum(na.omit(weekly.deaths))), by = 'date']
  tmp = merge(tmp, tmp2, by = 'date')
  tmp[, emp.prob := weekly.deaths / total.deaths]
  ref_dates = subset(df_week, week_index %in% stan_data$w_ref_index)$date
  tmp2 = subset(tmp, date %in% ref_dates)
  tmp2 = tmp2[, list(emp.prob = mean(emp.prob)), by = c('age')]
  setnames(tmp2, 'emp.prob', 'emp.prob.ref')
  tmp = merge(tmp, tmp2, by = c('age'))
  tmp[, emp.prob.ratio := emp.prob / emp.prob.ref]
  subset(tmp, age %in% unique(data$age))
  
  tmp1 = merge(tmp1, select(tmp, -weekly.deaths), by = c('age', 'date'), all.x = T)
  
  # save
  saveRDS(tmp1, file = paste0(outdir, '-ProbabilityRatioTable_', Code, '.rds'))

  return(tmp1)
}

make_var_by_age_by_state_table = function(fit, df_week, df_state_age, df_state, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
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

make_var_by_age_table = function(fit, df_week, df_state_age, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
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


make_var_cum_by_age_table = function(fit, df_week, df_state_age, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
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

make_var_cum_by_age_table_counterfactual = function(fit, df_week, df_week_counterfactual, resurgence_dates, df_state_age, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # find start resurgence date
  df_week[, dummy := 1]; 
  df_week3 <- merge(df_week, resurgence_dates, by = 'dummy', allow.cartesian=TRUE)
  df_week3[, after_resurgence := date >= start_resurgence]
  df_week3 = merge(df_week3, df_state, by = 'code')
  
  df_week_counterfactual = df_week3[(after_resurgence)]
  df_week_counterfactual[, week_index_resurgence := 1:length(date), by = 'code']

  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name[1]]]) )
  setnames(tmp1, 2:4, c('state_index', 'age_index','week_index'))
  tmp1 = merge(tmp1, select(df_week3, state_index, week_index, after_resurgence), by = c('state_index', 'week_index'))
  tmp1 = tmp1[(!after_resurgence)]
  tmp1 = select(tmp1, -after_resurgence)
  
  tmp2 = as.data.table( reshape2::melt(fit_samples[[var_name[2]]]) )
  setnames(tmp2, 2:4, c('state_index', 'age_index','week_index_resurgence'))
  tmp2 = merge(tmp2, select(df_week_counterfactual, state_index, week_index, week_index_resurgence), by = c('state_index', 'week_index_resurgence'))
  tmp2 = select(tmp2,  -week_index_resurgence)

  tmp1 <- rbind(tmp1, tmp2)
  tmp1[, value := cumsum(value), by = c('state_index', 'age_index', 'iterations')]
  
  tmp1 = tmp1[, list(q=quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('state_index', 'age_index', 'week_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, select(df_state, loc_label, code, state_index), by = 'state_index')
  
  if('code' %in% names(df_week)){
    tmp1 = merge(tmp1, select(df_week, week_index, code, date), by = c('week_index', 'code'))
  }else{
    tmp1 = merge(tmp1, select(df_week, week_index, date), by = 'week_index')
  }
  
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', var_name[2],  'CumTable_', Code, '.rds'))
    
  }
  
  return(tmp1)
}


make_var_by_age_age_table = function(fit, df_state_age, prop_prediction, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, 2:4, c('age_index_recipient', 'age_index_source', 'prop_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index_recipient', 'age_index_source', 'prop_index')]	
  tmp1 = dcast(tmp1, age_index_recipient + age_index_source + prop_index ~ q_label, value.var = "q")
  
  tmp1[, age_recipient := df_state_age$age[age_index_recipient]]
  tmp1[, age_source := df_state_age$age[age_index_source]]
  tmp1 = merge(tmp1, prop_prediction, by = 'prop_index')
  
  return(tmp1)
}


make_var_by_varying_age_table = function(fit, df_week, df_age_continuous, age_groups, var_name, operation, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  # df age
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
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  if(operation == 'sum')
    tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'iterations', 'age_state_index')]
  if(operation == 'mean')
    tmp1 = tmp1[, list(value = mean(value)), by = c('week_index', 'iterations', 'age_state_index')]
  
  # take quantiles
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, week_index + age_state_index ~ q_label, value.var = "q")
  setnames(tmp1, 'age_state_index', 'age_index')
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1 = merge(tmp1, data.table(age_index = 1:length(age_groups), age = age_groups), by = 'age_index')
  
  return(tmp1)
}


make_contribution_ref = function(fit, data_10thdeaths, fiveagegroups, data, df_week, df_age_continuous, outdir){
  
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
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
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

make_contribution_ref_adj = function(fit, data_10thdeaths, fouragegroups, df_week, pop_data, outdir){
  
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
  pop_data1[, age_index := 1:nrow(pop_data1)]
  pop_data1[, age_from := gsub('(.+)-.*', '\\1', age)]
  pop_data1[, age_to := gsub('.*-(.+)', '\\1', age)]
  pop_data1[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  pop_data1[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  pop_data1[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  pop_data1[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  set(pop_data1, NULL, 'age_from', pop_data1[,as.numeric(age_from)])
  set(pop_data1, NULL, 'age_to', pop_data1[,as.numeric(age_to)])
  set(pop_data1, NULL, 'age_to_index', pop_data1[,as.numeric(age_to_index)])
  set(pop_data1, NULL, 'age_from_index', pop_data1[,as.numeric(age_from_index)])
  pop_data1[, age_state_index := which(df_age$age_from_index <= age_from_index & df_age$age_to_index >= age_to_index), by = 'age_index']
  pop_data1 = pop_data1[, list(pop = sum(pop)), by = c('Total', 'age_state_index', 'code')]
  setnames(pop_data1, 'age_state_index', 'age_index')
  pop_data1[, age := fouragegroups[age_index]]
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
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
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

find_contribution_one_age_group = function(fit, df_week, df_age_continuous, df_age_reporting, age_group, 
                                           date_10thcum, pop_data, data, outdir, with_empirical = F){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
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
  pop_data1[, age_from := gsub('(.+)-.*', '\\1', age)]
  pop_data1[, age_to := gsub('.*-(.+)', '\\1', age)]
  pop_data1[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  pop_data1[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  set(pop_data1, NULL, 'age_from', pop_data1[,as.numeric(age_from)])
  set(pop_data1, NULL, 'age_to', pop_data1[,as.numeric(age_to)])
  tmp = pop_data1[, list(pop = sum(pop)), by = 'age']
  tmp[, pop_prop_US := pop / sum(pop_data$pop)]
  pop_data1 = merge(pop_data1, select(tmp, age, pop_prop_US), by = 'age')
  pop_data1[, age_state_index := which(df_age_state$age_from <= age_from & df_age_state$age_to >= age_to), by = c('loc_label', 'age_from')]
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


find_contribution_age_groups_vaccination = function(fit, df_week, df_age_continuous, df_age_reporting, deathByAge, 
                                                    age_groups, age_indices, var, empirical, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
  # df age
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
  
  # empirical
  if(empirical){
    df_age_close_vaccination = copy(df_age_reporting)
    df_age_close_vaccination[, age_index := age_indices]
    emp = merge(deathByAge, df_age_close_vaccination, by = c('age_from', 'age_to', 'age'))
    emp = emp[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c('code', 'date', 'loc_label', 'age_index')]
    tmp2 = emp[, list(total_deaths = sum(na.omit(weekly.deaths))), by = c('date', 'code')]
    emp = merge(emp, tmp2, by = c('date', 'code'))
    emp[, prop_deaths := weekly.deaths / total_deaths]
  }

  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var]]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'iterations', 'age_state_index')]
  
  
  # take quantiles
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, week_index + age_state_index ~ q_label, value.var = "q")
  setnames(tmp1, 'age_state_index', 'age_index')

  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  if(empirical)
    tmp1 = merge(tmp1, emp, by = c('date', 'age_index', 'code'))
  tmp1 = merge(tmp1, data.table(age_index = 1:length(age_groups), age = age_groups), by = 'age_index')
  
  # save
  saveRDS(tmp1, file = paste0(outdir, '-posterior_table_', var, '_', Code, '.rds'))
  
  return(tmp1)
}



find_cumsum_nonr_deaths_state_age = function(fit, df_week, df_age_continuous, state_age_groups, stan_data, deaths_predict_var){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
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

find_sum_nonr_deaths_state_age = function(fit, df_age_continuous, state_age_groups, stan_data, deaths_predict_var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
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

find_mean_age_death = function(fit, df_week, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[['phi']]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = tmp1[, list( 	value= sum(age_index*value)), 
              by=c('iterations', 'week_index')]	
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index')]	
  tmp1 = dcast(tmp1, week_index  ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, code := Code]
  
  # save
  saveRDS(tmp1, file = paste0(outdir, '-MeanAgeOfDeath_', Code, '.rds'))
  
  return(tmp1)
}

find_cumulative_deaths_prop_givensum_state_age_multiple_states <- function(fit, date_10thcum, df_week, df_age_continuous, scrapedData, cum.death.var, 
                                                                           Code, outdir){
  if(nrow(subset(scrapedData, code %in% Code)) > 0 ){
    
    tmp = list(); i = 1
    for(c in Code){

      scrapedData_c = subset(scrapedData, code == c)
      
      if(c == 'GA'){
        scrapedData_c = reduce_agebands_scrapedData_GA(scrapedData_c)
      }
      
      if(nrow(scrapedData_c) > 0){
        tmp[[i]] = find_cumulative_deaths_prop_givensum_state_age(fit, date_10thcum, df_week, df_age_continuous, scrapedData_c, cum.death.var, c, outdir)
        i = i + 1
      }
      
    }

    tmp = do.call('rbind',tmp)
    tmp = merge(tmp, df_state, by = 'code')
    
    return(tmp)
  }
  
  
}


find_cumulative_deaths_prop_givensum_state_age = function(fit, date_10thcum, df_week, df_age_continuous, data_comp, cum.death.var, c, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
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

make_mortality_rate_table = function(fit_cum, fouragegroups, date_10thcum, df_week, pop_data, data_comp, df_age_continuous, cum.death.var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit_cum)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
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
  pop_data1[, age_index := 1:nrow(pop_data1)]
  pop_data1[, age_from := gsub('(.+)-.*', '\\1', age)]
  pop_data1[, age_to := gsub('.*-(.+)', '\\1', age)]
  pop_data1[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  pop_data1[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
  pop_data1[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  pop_data1[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  pop_data1[, age_state_index := which(df_age$age_from_index <= age_from_index & df_age$age_to_index >= age_to_index), by = 'age_index']
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
  
  # quantiles
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('state_index', 'week_index', 'age_index')]	
  tmp1 = dcast(tmp1, state_index + week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week_adj, by = 'week_index')
  tmp1[, age := df_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_age$age)]
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  for(Code in unique(tmp1$code)){
    saveRDS(subset(tmp1, code == Code), file = paste0(outdir, '-', 'MortalityRate', 'Table_', Code, '.rds'))
    
  }

  return(tmp1)
}

make_weekly_death_rate_other_source = function(fit_cum, df_week, data_comp, var.alpha, df_age, outdir, age_groups = NULL, 
                                               lab = NULL, withempirical = F, reduction = NULL){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit_cum)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
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
  
  if(!is.null(reduction)){
    tmp3 = tmp1[week_index %in% df_week[date %in%reduction, week_index]]
    tmp3 = select(tmp3, - weekly.deaths)
    tmp3 = as.data.table( dcast.data.table(tmp3, state_index + age_index + iterations ~ week_index, value.var = 'value') )
    setnames(tmp3, 4:5, c('week1', 'week2'))
    tmp3[, value := (week2 - week1) / week1, by = c('iterations', 'state_index', 'age_index')]
    tmp3 = tmp3[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('state_index', 'age_index')]	
    tmp3 = dcast(tmp3, state_index + age_index ~ q_label, value.var = "q")
    tmp3 = merge(tmp3, df_state, by = 'state_index')
    tmp3[, age := df_age$age[age_index]]
    tmp3[, age := factor(age, levels = df_age$age)]
  }
  
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
    
    if(!is.null(reduction)){
      if(!is.null(lab)){
        file = paste0(outdir, '-', 'DeathByAge', 'prop_Table_', lab, '_', Code, '.rds')
        
      }else{
        file =  paste0(outdir, '-', 'DeathByAge', 'prop_Table_', var.alpha, '_', Code, '.rds')
      }
      
      saveRDS(subset(tmp3, code == Code), file= file)
    }
  }

  
  return(tmp1)
}

make_weekly_death_rate_other_source_posteriorsamples = function(fit_cum, df_week, data_comp, var.alpha, df_age, outdir, age_groups = NULL, 
                                               lab = NULL,reduction = NULL){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit_cum)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
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
  data_comp = subset(data_comp, code == Code)
  data_comp[, date := as.Date(date)]
  tmp2 = data_comp[, list(daily_deaths = sum(daily_deaths)), by = 'date']
  tmp2 = tmp2[order( date)]
  tmp2[, cum.death := cumsum(daily_deaths)]
  tmp2 = unique(subset(tmp2, date %in% c(df_week$date, max(df_week$date)+7)))
  tmp2[, weekly.deaths := c(diff(cum.death),NA)]
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = 'date')
  tmp2 = select(tmp2, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))
  tmp2[weekly.deaths <0, weekly.deaths := 0]
  
  # adjust dfweek for cum
  df_week_adj = df_week[, list(date = date + 7), by = 'week_index']
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var.alpha]]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  if(!is.null(age_groups)){
    tmp1 = merge(tmp1, df_age_continuous, 'age_index')
    tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
    setnames(tmp1, 'age_state_index', 'age_index')
  }
  
  # predict weekly death
  tmp1 = merge(tmp1, tmp2, by = 'week_index')
  tmp1[, value := rdirmnom(1, weekly.deaths, value), by = c('week_index', 'iterations')]
  
  if(!is.null(reduction)){
    tmp3 = tmp1[week_index %in% df_week[date %in%reduction, week_index]]
    tmp3 = select(tmp3, - weekly.deaths)
    tmp3 = as.data.table( dcast.data.table(tmp3, age_index + iterations ~ week_index, value.var = 'value') )
    setnames(tmp3, 3:4, c('week1', 'week2'))
    tmp3[, value := week2 / week1, by = c('iterations', 'age_index')]
    tmp3[, code := Code]
    tmp3[, age := df_age$age[age_index]]
    tmp3[, age := factor(age, levels = df_age$age)]
  }
  
  tmp1[, code := Code]
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_age$age)]
  

  if(!is.null(lab)){
    file = paste0(outdir, '-', 'PosteriorsamplesDeathByAge', 'Table_', lab, '_', Code, '.rds')
    
  }else{
    file =  paste0(outdir, '-', 'PosteriorsamplesDeathByAge', 'Table_', var.alpha, '_', Code, '.rds')
  }
  
  saveRDS(tmp1, file =file)
  
  if(!is.null(reduction)){
    if(!is.null(lab)){
      file = paste0(outdir, '-', 'PosteriorsamplesDeathByAge', 'prop_Table_', lab, '_', Code, '.rds')
      
    }else{
      file =  paste0(outdir, '-', 'PosteriorsamplesDeathByAge', 'prop_Table_', var.alpha, '_', Code, '.rds')
    }
    
    saveRDS(tmp3, file= file)
  }
  
  return(tmp1)
}

find_vaccine_effects_scaled <- function(fit, df_week, df_age_continuous, age_groups, var, suffix_var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  # df age
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
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[[paste0(var, suffix_var[1])]]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp2 = as.data.table( reshape2::melt(fit_samples[[paste0(var,  suffix_var[2])]]) )
  setnames(tmp2, c('Var2', 'Var3', 'value'), c('age_index','week_index', 'value_wo_vaccine'))
  tmp1 = merge(tmp1, tmp2, by = c('iterations', 'age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value), 
                     value_wo_vaccine = sum(value_wo_vaccine)), by = c('week_index', 'iterations', 'age_state_index')]
  tmp1 = tmp1[, value := value - value_wo_vaccine, by = c('week_index', 'iterations', 'age_state_index')]

  # take quantiles
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, week_index + age_state_index ~ q_label, value.var = "q")
  setnames(tmp1, 'age_state_index', 'age_index')
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1 = merge(tmp1, data.table(age_index = 1:length(age_groups), age = age_groups), by = 'age_index')
  
  file =  paste0(outdir, '-', 'VaccineEffects_', var, '_', Code, '.rds')
  saveRDS(tmp1, file)
  
  return(tmp1)
  
}

find_vaccine_effects_old <- function(weekly_deaths, vaccine_data, start_resurgence, pick_resurgence)
{
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  df_age_vaccination = unique(select(weekly_deaths, age_index, age))
  df_age_vaccination[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age_vaccination[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age_vaccination[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age_vaccination[grepl('\\+', age_to), age_to := max(vaccine_data$age)]
  set(df_age_vaccination, NULL, 'age_from', df_age_vaccination[,as.numeric(age_from)])
  set(df_age_vaccination, NULL, 'age_to', df_age_vaccination[,as.numeric(age_to)])
  vaccine_data[, age_index := which(df_age_vaccination$age_from <= age & df_age_vaccination$age_to >= age), by = 'age']
  
  delay = 2*7
  tmp = vaccine_data[, list(prop = unique(prop)), by = c('code', 'date', 'loc_label', 'age_index')]
  tmp[, age_index := paste0('prop_', age_index)]
  tmp = as.data.table( reshape2::dcast(tmp, code + date + loc_label ~ age_index, value.var = 'prop') )
  tmp[, date := date + delay]
  # tmp = subset(tmp, date == start_resurgence-delay)

  # population
  tmp1 = vaccine_data[, list(pop = sum(pop)), by = c('code', 'loc_label', 'age_index')]
  tmp2 = tmp1[, list(pop_avg = mean(pop)), by = c('age_index')]
  tmp1 = merge(tmp1, tmp2, by = c('age_index'))
  
  # weekly deaths before resurgence
  tmp2 = subset(weekly_deaths, date %in% c(start_resurgence - 7*1:2))
  tmp2 = tmp2[, list(weekly_deaths_before_resurgence = mean(value)), by = c('code', 'age_index', 'iterations')]
  
  # adjusted weekly deaths
  tmp2 = merge(weekly_deaths, tmp2, by = c('code', 'age_index', 'iterations'))
  # tmp2 = merge(tmp2, tmp1, by = c('code', 'loc_label', 'age_index'))
  tmp2[, prop_increase := value / (weekly_deaths_before_resurgence)]
  tmp2[value == 0 & weekly_deaths_before_resurgence == 0, prop_increase := 1]
  tmp2[value > 0 & weekly_deaths_before_resurgence == 0, prop_increase := value / 0.1]
  stopifnot(sum(is.na(tmp2$prop_increase)) == 0 )
  stopifnot(sum(tmp2$prop_increase == Inf) == 0 )
  tmp2 = subset(tmp2, age %in% c('18-64', '65+'))
  
  # find date index
  tmp3 = merge(tmp2, tmp, by = c('code', 'loc_label', 'date'))
  tmp3 = subset(tmp3, date >= start_resurgence)
  tmp4 = unique(select(tmp3, date))
  tmp4[, min_date := min(date)]
  tmp4[, max_date := max(date)]
  tmp4[, date_index := which(seq.Date(min_date, max_date, 7) == date), by = 'date']
  tmp3 = merge(tmp3, tmp4, by = 'date')
  
  # regression to find the vaccine effects
  tmp4 = unique(select(tmp3, loc_label, date, age, prop_3, prop_4))
  tmp4[, prop_3_start_resurgence := prop_3[date == start_resurgence], by = c('loc_label', 'age')]
  tmp4[, prop_4_start_resurgence := prop_4[date == start_resurgence], by = c('loc_label', 'age')]
  # tmp4[, prop_start_resurgence_cat := ifelse(prop_3_start_resurgence < 0.45 & prop_4_start_resurgence < 0.8, 1, 
  #                                             ifelse(prop_3_start_resurgence >= 0.45 & prop_4_start_resurgence >= 0.8, 2, 
  #                                                    ifelse(prop_3_start_resurgence < 0.45, 3, 4)))]
  tmp4[, prop_start_resurgence_cat := ifelse(prop_3_start_resurgence < 0.45 & prop_4_start_resurgence < 0.825, 1, 
                                             ifelse(prop_3_start_resurgence >= 0.45 & prop_4_start_resurgence >= 0.825, 2, 3))]
  tmp4[, prop_start_resurgence_cat := factor(prop_start_resurgence_cat, levels = c(2,3,1))]
  tmp3 = merge(tmp3, tmp4, by = c('loc_label', 'date', 'age', 'prop_3', 'prop_4'))
  tmp3[, dummy := 1]
  
  # prediction table
  pred_tmp = unique(select(tmp3, loc_label, date_index, prop_3, prop_4, prop_start_resurgence_cat, dummy))
  setnames(pred_tmp, 'loc_label', 'loc_label_pred')
  pred_tmp = merge(pred_tmp, unique(select(tmp3, loc_label, dummy)), allow.cartesian=TRUE, by = 'dummy')
  pred_tmp[, type := 'same intercept - same slope']
  pred_tmp1 = copy(pred_tmp[loc_label_pred == loc_label])
  pred_tmp1[, prop_3 := prop_4]
  pred_tmp1[, prop_start_resurgence_cat := ifelse(min(prop_4) < 0.825, 1, 2)]
  pred_tmp1[, type := 'prop_3 same as prop_4']
  pred_tmp = rbind(pred_tmp1, pred_tmp)
  # pred_tmp1 = merge(unique(select(pred_tmp, -prop_3, -prop_4, -date_index, -loc_label)), unique(select(tmp3, loc_label, prop_3, prop_4, date_index, dummy)), allow.cartesian=TRUE, by = 'dummy')
  # pred_tmp1[, type := 'same intercept']
  # pred_tmp = rbind(pred_tmp1, pred_tmp)
  pred_tmp[, date_index_sigmoid := (1/ (1 + exp(-(date_index - round(max(date_index)) / 2))))]
  pred_tmp[, idx_predict := 1:nrow(pred_tmp)]
  
  # subset(pred_tmp, loc_label_pred %in% c('Florida', 'New York') & type == 'same intercept' & loc_label == 'California')
  tmp3 = subset(tmp3, iterations < 1000)
  tmp3[, date_index_sigmoid := (1/ (1 + exp(-(date_index - round(max(date_index)) / 2))))]
  tmp4 = tmp3[,
     {
    # fit <- lm(log(prop_increase) ~ date_index + date_index_sigmoid + prop_3 + prop_4 + prop_3*date_index + prop_4*date_index + prop_3*prop_4 + prop_start_resurgence_cat)
    fit <- lm(log(prop_increase) ~ date_index + date_index_sigmoid + prop_3 + prop_4 + prop_3*date_index + prop_4*date_index + prop_3*prop_4 )
    
    predict <- exp( predict(fit, newdata = pred_tmp, type = 'response') )
    summary <- summary(fit)$coefficients
    list(predict = predict, idx_predict = pred_tmp$idx_predict, 
         par_prop_3 = summary[rownames(summary) == 'prop_3', 1]/100, 
         par_prop_4 = summary[rownames(summary) == 'prop_4', 1]/100#, 
         # par_prop_3_prop_4 = summary[rownames(summary) == 'prop_3:prop_4', 1]/1000,
         # par_prop_start_resurgence_cat1 = summary[rownames(summary) == 'prop_start_resurgence_cat1', 1]/100, 
         # par_prop_start_resurgence_cat3 = summary[rownames(summary) == 'prop_start_resurgence_cat3', 1]/100, 
         # par_prop_start_resurgence_cat4 = summary[rownames(summary) == 'prop_start_resurgence_cat4', 1]/100, 
         # par_date_index_prop_3 = summary[rownames(summary) == 'date_index:prop_3', 1]/100,
         # par_date_index_prop_4 = summary[rownames(summary) == 'date_index:prop_4', 1]/100
         )
    }, by = c('age', 'iterations')]
  tmp4 = merge(tmp4, pred_tmp, by = 'idx_predict')
  
  
  if(0){ # demonstrate the fit

    df = subset(tmp3, iterations ==1 & age == '65+')
    df[, log_increase := log(prop_increase)]
    df[, date_index_sigmoid := (1/ (1 + exp(-(date_index - round(max(date_index)) / 2))))]
    fit = lm(log_increase ~ date_index + date_index_sigmoid + prop_3 + prop_4, data = df)
    
    # ggplot(subset(df, code == 'TX'), aes(x = date, y = date_index_sigmoid)) + geom_line()
    summary(fit)
    subset(df, code == 'TX')
    df$predict = exp(predict(fit, newdata = df, type = 'response'))
    df[, prop_3 := prop_4]
    df$predict1 = exp(predict(fit, newdata = df, type = 'response'))
    
   ggplot(df, aes(x = date, col = loc_label)) +
      facet_wrap(~age) +
      geom_point(aes(y = prop_increase)) +
      geom_line(aes(y = predict)) +
      geom_line(aes(y = predict1), linetype = 'dashed')
    
    ggplot(df, aes(x = prop_3)) +
      facet_wrap(~age) +
      geom_point(aes(y = prop_increase, col = loc_label)) +
      geom_point(aes(y = predict), col = 'red') + 
      scale_y_log10()
    
    ggplot(df, aes(x = prop_4)) +
      facet_wrap(~age) +
      geom_point(aes(y = prop_increase, col = loc_label)) +
      geom_point(aes(y = predict), col = 'red')

  }
  
  # parameters
  tmp5 = unique(select(tmp4, names(tmp4)[grepl('par', names(tmp4))], iterations, age))
  tmp5 = as.data.table( reshape2::melt(tmp5, id.vars = c('iterations', 'age')))
  
  # predictions
  tmp4 = merge(tmp3, select(tmp4, loc_label, loc_label_pred, date_index, age, predict, type, iterations), by = c('date_index', 'age', 'iterations', 'loc_label'))
  tmp4[, weekly.deaths_predict := weekly_deaths_before_resurgence * predict]
  tmp4[, diff.weekly.deaths := weekly.deaths_predict - value]
  stopifnot(sum(is.na(tmp4$prop_increase_predict)) == 0)
  # subset(tmp4, code =='TX' & type == 'same intercept - same slope' & iterations == 1)
  
  # summarise
  tmp1 = tmp4[, list( 	q= quantile(weekly.deaths_predict, prob=ps, na.rm = T),
                      q_label=paste0(p_labs, '_predict')), 
             by=c('age', 'code', 'loc_label', 'loc_label_pred', 'type', 'date')]	
  tmp1 = dcast(tmp1, code + loc_label + loc_label_pred + date + age + type ~ q_label, value.var = "q")
  tmp3 = tmp4[, list( 	q= quantile(diff.weekly.deaths, prob=ps, na.rm = T),
                       q_label=paste0(p_labs, '_diff_predict')), 
              by=c('age', 'code', 'loc_label', 'loc_label_pred', 'type', 'date')]	
  tmp3 = dcast(tmp3, code + loc_label + loc_label_pred + date + age + type ~ q_label, value.var = "q")
  tmp1 = merge(tmp3, tmp1, by=c('age', 'code', 'loc_label', 'loc_label_pred', 'type', 'date'))
  
  tmp3 = tmp2[, list( 	q= quantile(prop_increase, prob=ps, na.rm = T),
                      q_label=paste0(p_labs)), 
             by=c('age', 'code', 'date', 'loc_label')]	
  tmp3 = dcast(tmp3, loc_label + code + date + age ~ q_label, value.var = "q")
  tmp2 = tmp2[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=paste0(p_labs, '_weekly_deaths')), 
              by=c('age', 'code', 'date', 'loc_label')]	
  tmp2 = dcast(tmp2, loc_label + code + date + age ~ q_label, value.var = "q")
  tmp2 = merge(tmp2, tmp3, by=c('age', 'code', 'date', 'loc_label'))
  
  tmp3 = tmp5[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=paste0(p_labs)), 
              by=c('variable', 'age')]	
  tmp3 = dcast(tmp3, variable + age ~ q_label, value.var = "q")
  
  return(list(tmp1, tmp2, tmp3, tmp))
}

make_forest_plot_table <- function(summary, df_age_vaccination2, df_state, names, math_name, groups, groups_levels){
  tmp <- summary[ grepl(paste(paste0('^',names),collapse = '|'),rownames(summary)) ,]
  
  variables <- rownames(tmp); group = c()
  for(x in 1:length(variables)) group[x] = groups[which(grepl(gsub('(.+)\\[.*', '\\1', variables[x]), names))]
  for(x in 1:length(names)) variables = gsub(names[x], math_name[x], variables)
  for(x in df_age_vaccination2$age_index) variables = gsub(paste0('\\[',x), paste0('\\["', df_age_vaccination2$age[x], '"'), variables) 
  for(x in df_age_vaccination2$age_index) variables[grepl('vac', variables)] = gsub(paste0('\",', x, '\\]'), paste0(', ',df_age_vaccination2$age[x], '"\\]'), variables[grepl('vac', variables)]) 
  for(x in df_state$state_index) variables[!grepl('vac', variables)] = gsub(paste0('\",', x, '\\]'), paste0(', ',df_state$loc_label[x], '"\\]'), variables[!grepl('vac', variables)]) 
  
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

find_vaccine_effects <- function(weekly_deaths, vaccine_data, start_resurgence, pick_resurgence)
{
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  df_age_vaccination = unique(select(weekly_deaths, age_index, age))
  df_age_vaccination[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age_vaccination[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age_vaccination[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age_vaccination[grepl('\\+', age_to), age_to := max(vaccine_data$age)]
  set(df_age_vaccination, NULL, 'age_from', df_age_vaccination[,as.numeric(age_from)])
  set(df_age_vaccination, NULL, 'age_to', df_age_vaccination[,as.numeric(age_to)])
  vaccine_data[, age_index := which(df_age_vaccination$age_from <= age & df_age_vaccination$age_to >= age), by = 'age']
  
  delay = 2*7
  tmp = vaccine_data[, list(prop = unique(prop)), by = c('code', 'date', 'loc_label', 'age_index')]
  tmp[, age_index := paste0('prop_', age_index)]
  tmp = as.data.table( reshape2::dcast(tmp, code + date + loc_label ~ age_index, value.var = 'prop') )
  tmp[, date := date + delay]
  # tmp = subset(tmp, date == start_resurgence-delay)
  
  # population
  tmp1 = vaccine_data[, list(pop = sum(pop)), by = c('code', 'loc_label', 'age_index')]
  tmp2 = tmp1[, list(pop_avg = mean(pop)), by = c('age_index')]
  tmp1 = merge(tmp1, tmp2, by = c('age_index'))
  
  # weekly deaths before resurgence
  tmp2 = subset(weekly_deaths, date %in% c(start_resurgence - 7*1:2))
  tmp2 = tmp2[, list(weekly_deaths_before_resurgence = mean(value)), by = c('code', 'age_index', 'iterations')]
  
  # adjusted weekly deaths
  tmp2 = merge(weekly_deaths, tmp2, by = c('code', 'age_index', 'iterations'))
  # tmp2 = merge(tmp2, tmp1, by = c('code', 'loc_label', 'age_index'))
  tmp2[, prop_increase := value / (weekly_deaths_before_resurgence)]
  tmp2[value == 0 & weekly_deaths_before_resurgence == 0, prop_increase := 1]
  tmp2[value > 0 & weekly_deaths_before_resurgence == 0, prop_increase := value / 0.1]
  stopifnot(sum(is.na(tmp2$prop_increase)) == 0 )
  stopifnot(sum(tmp2$prop_increase == Inf) == 0 )
  tmp2 = subset(tmp2, age %in% c('18-64', '65+'))
  
  # find date index
  tmp3 = merge(tmp2, tmp, by = c('code', 'loc_label', 'date'))
  tmp3 = subset(tmp3, date >= start_resurgence)
  tmp4 = unique(select(tmp3, date))
  tmp4[, min_date := min(date)]
  tmp4[, max_date := max(date)]
  tmp4[, date_index := which(seq.Date(min_date, max_date, 7) == date), by = 'date']
  tmp3 = merge(tmp3, tmp4, by = 'date')
  
  # subset(pred_tmp, loc_label_pred %in% c('Florida', 'New York') & type == 'same intercept' & loc_label == 'California')
  states <- unique(tmp3$loc_label)
  tmp3[, log_increase := log(prop_increase)]
  df = subset(tmp3, iterations == 1)
  df2 = df[, list(min_prop_3 = min(prop_3), min_prop_4 = min(prop_4)), by = c('age', 'loc_label')]
  df = merge(df2, df, by = c('age', 'loc_label'))
  df[, prop_3_adj := prop_3 - min_prop_3]
  df[, prop_4_adj := prop_4 - min_prop_4]
  
  ages = unique(df$age); stan_data = list()
  for(a in 1:length(ages)){
    # a = 1

    prop_3 <- as.matrix(reshape2::dcast(subset(df, age == ages[a]), date_index~ loc_label, value.var = 'prop_3_adj')[,-1])
    prop_4 <- as.matrix(reshape2::dcast(subset(df, age == ages[a]), date_index~ loc_label, value.var = 'prop_4_adj')[,-1])
    
    prop_vac_start_3 = as.numeric(reshape2::dcast(subset(df2, age == ages[a]), .~ loc_label, value.var = 'min_prop_3')[,-1])
    prop_vac_start_4 = as.numeric(reshape2::dcast(subset(df2, age == ages[a]), .~ loc_label, value.var = 'min_prop_4')[,-1])
    
    
    stan_data[[a]] = list(M = length(states),
                          T = max(df$date_index),
                          time = 1:max(df$date_index),
                          C = 2,
                          prop_vac_start = list(prop_vac_start_3, prop_vac_start_4),
                          c_counterfactual = 2,
                          # S = max(as.numeric(df$prop_start_resurgence_cat)),
                          prop_vac = list(prop_3, prop_4))
    
  }


  tmp3 = subset(tmp3, iterations < 1000); tmp4 = list(); tmp5 = list()
  path.to.stan.model = file.path(indir, 'stan-models-union', '211006.stan')
  model = rstan::stan_model(path.to.stan.model)
  
  # fit
  j = 1
  for(a in 1:length(ages)){
    
    for(i in unique(tmp3$iterations)){
      
      # cat('iterations ', i, '\n')
      
      # cat('age ', a, '\n')
      
      stan_data1 <- c(list(y = as.matrix(reshape2::dcast(subset(tmp3, iterations == i & age == ages[a]), date_index ~ loc_label, value.var = 'log_increase')[,-1])), 
                      stan_data[[a]])
      
      if(i == 1){
        fit <- rstan::sampling(model,data=stan_data1,iter=1001,warmup=100,chains=1,seed=JOBID,refresh = 0)
        
      } else{
        fit <- rstan::sampling(model,data=stan_data1,iter=2,warmup=1,chains=1,seed=JOBID,refresh = 0, init=list(stan_init))
      }
      
      samples =  rstan::extract(fit) 
      df = as.data.table( reshape2::melt(samples) )
      df = subset(df, iterations == max(iterations)); stan_init = list()
      for(var in unique(df$L1)){
        stan_init[[var]] = subset(df, L1 == var)$value
      }
      
      variables = c('y_hat', 'y_counterfactual'); df = list(); k = 1
      for(var in variables){
        df[[k]] = as.data.table(reshape2::melt(samples[[var]]))
        df[[k]] = subset(df[[k]], iterations == max(iterations))
        df[[k]][, var := var]
        k = k + 1
      }
      df = do.call('rbind', df)
      df[, iterations := i]
      df[, age := ages[a]]
      tmp4[[j]] = df
      
      variables = c('lambda', 'beta'); df = list(); k = 1
      for(var in variables){
        df[[k]] = as.data.table(reshape2::melt(samples[[var]]))
        df[[k]] = subset(df[[k]], iterations == max(iterations))
        df[[k]][, var := var]
        k = k + 1
      }
      df = do.call('rbind', df)
      df[, iterations := i]
      df[, age := ages[a]]
      tmp5[[j]] = df
      
      j = j + 1
    }
  }
  
  tmp4 = do.call('rbind', tmp4)
  setnames(tmp4, 2:4, c('date_index', 'loc_label_index', 'value_pred'))
  tmp4[, loc_label := states[loc_label_index]]
  tmp4[, value_pred := exp(value_pred)]
  
  tmp5 = do.call('rbind', tmp5)
  setnames(tmp5, 2, c('age_index_param'))
  tmp5[, age_param := ages[age_index_param]]
  
  # predictions
  tmp4 = merge(tmp3, subset(tmp4, var %in% c('y_hat', 'y_counterfactual')), by = c('date_index', 'age', 'iterations', 'loc_label'))
  tmp4[, weekly.deaths_predict := weekly_deaths_before_resurgence * value_pred]
  tmp4[, diff.weekly.deaths := weekly.deaths_predict - value]
  tmp4[, diff.weekly.deaths := cumsum(diff.weekly.deaths), by = c( 'age', 'iterations', 'loc_label')]
  stopifnot(sum(is.na(tmp4$prop_increase_predict)) == 0)

  # summarise
  tmp1 = tmp4[, list( 	q= quantile(weekly.deaths_predict, prob=ps, na.rm = T),
                       q_label=paste0(p_labs, '_predict')), 
              by=c('age', 'code', 'loc_label',  'date', 'var')]	
  tmp1 = dcast(tmp1, code + loc_label  + date + var +age ~ q_label, value.var = "q")
  tmp3 = tmp4[, list( 	q= quantile(diff.weekly.deaths, prob=ps, na.rm = T),
                       q_label=paste0(p_labs, '_diff_predict')), 
              by=c('age', 'code', 'loc_label', 'date', 'var')]	
  tmp3 = dcast(tmp3, code + loc_label  + date + var + age  ~ q_label, value.var = "q")
  tmp1 = merge(tmp3, tmp1, by=c('age', 'code', 'loc_label',  'date', 'var'))
  
  tmp3 = tmp2[, list( 	q= quantile(prop_increase, prob=ps, na.rm = T),
                       q_label=paste0(p_labs)), 
              by=c('age', 'code', 'date', 'loc_label')]	
  tmp3 = dcast(tmp3, loc_label + code + date + age ~ q_label, value.var = "q")
  tmp2 = tmp2[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=paste0(p_labs, '_weekly_deaths')), 
              by=c('age', 'code', 'date', 'loc_label')]	
  tmp2 = dcast(tmp2, loc_label + code + date + age ~ q_label, value.var = "q")
  tmp2 = merge(tmp2, tmp3, by=c('age', 'code', 'date', 'loc_label'))
  
  tmp5 = tmp5[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=paste0(p_labs)), 
              by=c('age', 'var', 'age_param')]	
  tmp5 = dcast(tmp5, var + age  + age_param ~ q_label, value.var = "q")


  return(list(tmp1, tmp2, tmp3, tmp, tmp5))
}

find_vaccine_effects_unscaled <- function(fit, df_week, df_age_continuous, age_groups, var, suffix_var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  # df age
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
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[[paste0(var, suffix_var[1])]]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp2 = as.data.table( reshape2::melt(fit_samples[[paste0(var,  suffix_var[2])]]) )
  setnames(tmp2, c('Var2', 'Var3', 'value'), c('age_index','week_index', 'value_wo_vaccine'))
  tmp1 = merge(tmp1, tmp2, by = c('iterations', 'age_index','week_index'))
  
  # temporary fix 
  if(stan_model == '210923b'){
    df = data.table(vac_start = stan_data$w_vaccination_start, age_index = 1:length(stan_data$w_vaccination_start))
    tmp1 = merge(tmp1, df, by = c('age_index'))
    tmp2 = tmp1[, list(value_before_vac = value[week_index == vac_start - 1]), by = c('age_index', 'iterations')]
    tmp1 = merge(tmp1, tmp2, by = c('age_index', 'iterations'))
    tmp1[vac_start < max(df_week$week_index) & vac_start >= week_index, value := value_before_vac, by = c('age_index', 'iterations', 'week_index')]
  }
  
  # find reduction 
  tmp1[, value1 := (value - value_wo_vaccine) / abs(value_wo_vaccine), by = c('week_index', 'iterations', 'age_index')]
  tmp1[value < 0 & !(value < 0 & value_wo_vaccine < 0), value1 := - (value_wo_vaccine - value) / (value_wo_vaccine), by = c('week_index', 'iterations', 'age_index')]
  tmp1[value_wo_vaccine < 0 & !(value < 0 & value_wo_vaccine < 0), value1 := (value_wo_vaccine - value) / (value_wo_vaccine), by = c('week_index', 'iterations', 'age_index')]
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value1 = median(value1)), by = c('week_index', 'iterations', 'age_state_index')]

  # take quantiles
  tmp1 = tmp1[, list( 	q= quantile(value1, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, week_index + age_state_index ~ q_label, value.var = "q")
  setnames(tmp1, 'age_state_index', 'age_index')
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1 = merge(tmp1, data.table(age_index = 1:length(age_groups), age = age_groups), by = 'age_index')
  
  file =  paste0(outdir, '-', 'VaccineEffects', '_', Code, '.rds')
  saveRDS(tmp1, file)
  
  return(tmp1)
  
}

find_vaccine_effects_unscaled_extrapolated <- function(fit, prop_tab, df_age_continuous, age_groups, var, suffix_var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  # df age
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
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[[paste0(var, suffix_var[1])]]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','prop_index'))
  tmp2 = as.data.table( reshape2::melt(fit_samples[[paste0(var,  suffix_var[2])]]) )
  setnames(tmp2, c('Var2', 'Var3', 'value'), c('age_index','week_index', 'value_wo_vaccine'))
  tmp2 <- subset(tmp2, week_index == max(week_index))
  tmp1 = merge(tmp1, tmp2, by = c('iterations', 'age_index'))
  
  # find reduction 
  tmp1 = tmp1[, value := value / value_wo_vaccine, by = c('prop_index', 'iterations', 'age_index')]
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = mean(value)), by = c('prop_index', 'iterations', 'age_state_index')]
  
  # take quantiles
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('prop_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, prop_index + age_state_index ~ q_label, value.var = "q")
  setnames(tmp1, 'age_state_index', 'age_index')
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, prop_tab, by = 'prop_index')
  tmp1 = merge(tmp1, data.table(age_index = 1:length(age_groups), age = age_groups), by = 'age_index')
  
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

prepare_prop_vac_table <- function(stan_data, df_week2, df_age_vaccination2){
  
  prop_vac = data.table( do.call('rbind', stan_data$prop_vac) ) 
  prop_vac[, age := rep(df_age_vaccination2$age, each = nrow(stan_data$prop_vac[[1]]))]
  prop_vac[, week_index := rep(1:nrow(stan_data$prop_vac[[1]]), length(stan_data$prop_vac))]
  prop_vac = data.table(reshape2::melt(prop_vac, id.vars = c('age', 'week_index')))
  setnames(prop_vac, c('variable', 'value'), c('code', 'prop'))
  prop_vac = merge(prop_vac, df_week2, by = c('week_index', 'code'))
  
  prop_vac_start = data.table( do.call('rbind', stan_data$prop_vac_start) ) 
  colnames(prop_vac_start) = colnames(stan_data$prop_vac[[1]])
  prop_vac_start[, age := df_age_vaccination2$age] 
  prop_vac_start = data.table(reshape2::melt(prop_vac_start, id.vars = c('age')))
  setnames(prop_vac_start, c('variable', 'value'), c('code', 'pre_prop'))
  
  prop_vac = merge(prop_vac, prop_vac_start, by = c('code', 'age'))
  prop_vac = merge(prop_vac, df_age_vaccination2, by = c('age'))
  prop_vac[, prop := prop + pre_prop]
  prop_vac[, cat := paste0('prop_', age_index)]
  prop_vac = data.table(reshape2::dcast(prop_vac, code + date ~ cat, value.var = 'prop'))
  
  return(prop_vac)
}

