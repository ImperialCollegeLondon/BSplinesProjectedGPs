make_predictive_checks_table = function(fit, df_week, df_state_age, data, deaths_predict_var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index', 'week_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_state_age$age[age_index]]
  
  tmp1 = merge(tmp1, data, by = c('date', 'age'))
  tmp1[, age := factor(age, levels = levels(data$age))]
  
  # stat
  tmp1[, inside.CI := weekly.deaths <= CU & weekly.deaths >= CL]
  
  # save
  saveRDS(tmp1, file = paste0(outdir, '-predictive_checks_table_', Code, '.rds'))
  
  return(tmp1)
}

make_convergence_diagnostics_stats = function(fit, outdir)
  {
  
  stopifnot(!is.null(fit))
  
  summary = rstan::summary(fit_cum)$summary
  eff_sample_size_cum = summary[,9][!is.na(summary[,9])]
  Rhat_cum = summary[,10][!is.na(summary[,10])]
  cat("the minimum and maximum effective sample size are ", range(eff_sample_size_cum), "\n")
  cat("the minimum and maximum Rhat are ", range(Rhat_cum), "\n")
  if(min(eff_sample_size_cum) < 500) cat('\nEffective sample size smaller than 500 \n')
  
  tryCatch({
    sampler_params <- get_sampler_params(fit_cum, inc_warmup = FALSE)
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
  re = rstan::extract(fit_cum)
  tryCatch({
    
    if('log_lik' %in% names(re)){
      log_lik <- loo::extract_log_lik(fit_cum)
      log_lik = log_lik[!is.na(log_lik[,1]),]
      .WAIC = loo::waic(log_lik)
      .LOO = loo::loo(log_lik)
      print(.WAIC); print(.LOO)
      WAIC = .WAIC$pointwise
      LOO = .LOO$pointwise
    }} , error = function(e) e)
  
  # time of execution
  time = sum(get_elapsed_time(fit))
  
  # save
  saveRDS(eff_sample_size_cum, file = paste0(outdir, "-eff_sample_size_cum_", Code, ".rds"))
  saveRDS(Rhat_cum, file = paste0(outdir, "-Rhat_cum_", Code, ".rds"))
  saveRDS(.WAIC, file = paste0(outdir, "-WAIC_", Code, ".rds"))
  saveRDS(.LOO, file = paste0(outdir, "-LOO_", Code, ".rds"))
  saveRDS(sampler_diagnostics, file = paste0(outdir, "-sampler_diagnostics_", Code, ".rds"))
  saveRDS(time, file = paste0(outdir, "-time_elapsed_", Code, ".rds"))
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

make_var_by_age_table = function(fit, df_week, df_state_age, var_name, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index', 'week_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  tmp1[, code := Code]

  saveRDS(tmp1, file = paste0(outdir, '-', var_name,  'Table_', Code, '.rds'))
  
  return(tmp1)
}

make_contribution_ref = function(fit, data_10thdeaths, fiveagegroups, data, df_week, df_age_continuous, outdir){
  
  df_week1 = subset(df_week, date >= data_10thdeaths)
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1 = subset(df_week1, month %in% min(df_week1$month):(min(df_week1$month) + 2))
  
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
  
  # empirical estimate
  df_age_reporting[, age_state_index := which(df_age$age_from <= age_from & df_age$age_to >= age_to), by = 'age_index']
  empirical = merge(data, df_age_reporting, by = 'age_from')
  empirical = merge(empirical, df_week1, by = 'date')
  empirical = empirical[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c( 'age_state_index')]
  empirical[, emp_est := weekly.deaths / sum(na.omit(weekly.deaths))]
  setnames(empirical, 'age_state_index', 'age_index')
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[['phi']]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # by selected age groups
  tmp1 = merge(tmp1, df_age_continuous, by = 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('iterations', 'week_index', 'age_state_index')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # take reference
  tmp1 = subset(tmp1, week_index %in% df_week1$week_index)
  tmp1 = tmp1[, list(value = mean(value)), by = c('iterations', 'age_index')]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index')]	
  tmp1 = dcast(tmp1, age_index ~ q_label, value.var = "q")
  
  tmp1[, age := fiveagegroups[age_index]]
  tmp1[, age := factor(age, levels = fiveagegroups)]
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, select(empirical, age_index, emp_est), by = 'age_index')
  
  saveRDS(tmp1, file = paste0(outdir, '-', 'contribution_ref',  'Table_', Code, '.rds'))
  
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
  pop_data1 = subset(pop_data1, code == Code)
  
  # ref week
  df_week1 = subset(df_week, date >= data_10thdeaths)
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1 = subset(df_week1, month %in% min(df_week1$month):(min(df_week1$month) + 2))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[['phi']]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # group by age group
  tmp1 = merge(tmp1, df_age_continuous, by = 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('iterations', 'age_state_index', 'week_index')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # take reference
  tmp1 = subset(tmp1, week_index %in% df_week1$week_index)
  tmp1 = tmp1[, list(value = mean(value)), by = c('iterations', 'age_index')]
  
  # adjust by population proportion
  tmp1 = merge(tmp1, select(pop_data1, age_index, pop_prop, pop_prop_US), by = 'age_index')
  tmp1[, value := value / pop_prop * pop_prop_US]
  tmp1[, value := value / sum(value), by = c('iterations')]
  
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index')]	
  tmp1 = dcast(tmp1, age_index ~ q_label, value.var = "q")
  
  tmp1[, age := df_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_age$age)]
  tmp1[, code := Code]
  
  saveRDS(tmp1, file = paste0(outdir, '-', 'contribution_ref_adj',  'Table_', Code, '.rds'))
  
  return(tmp1)
}


find_contribution_one_age_group = function(fit, df_week, df_age_continuous, df_age_reporting, age_group, date_10thcum, pop_data, data, outdir, with_empirical = F){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
  # find ref week 
  df_week1 = subset(df_week, date >= date_10thcum)
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1 = subset(df_week1, month %in% min(df_week1$month):(min(df_week1$month) + 2))
  
  # df age
  if(grepl('\\+', age_group))
  {
    age_groups = c(paste0('0-', as.numeric(gsub('(.+)\\+', '\\1', age_group))-1), age_group)
  } else if(grepl('0', age_group))
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
  pop_data1 = subset(pop_data1, code == Code)
  stopifnot(nrow(pop_data1) == length(age_groups))
  
  # empirical estimate
  if(with_empirical){
    df_age_reporting[, age_state_index := which(df_age_state$age_from <= age_from & df_age_state$age_to >= age_to), by = 'age_index']
    empirical = merge(data, df_age_reporting, by = 'age_from')
    empirical = empirical[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c('date', 'age_state_index')]
    baseline_emp = merge(empirical, df_week1, by = 'date')
    baseline_emp = baseline_emp[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c( 'age_state_index')]
    baseline_emp[, emp_ref := weekly.deaths / sum(na.omit(weekly.deaths))]
    empirical1 = empirical[, list(weekly.deaths_t = sum(na.omit(weekly.deaths))), by = c( 'date')]
    empirical = merge(empirical, empirical1, by = c( 'date'))
    empirical = merge(empirical, select(baseline_emp, age_state_index, emp_ref), by = c( 'age_state_index'))
    empirical = merge(empirical, pop_data1, by = 'age_state_index')
    empirical[, emp := weekly.deaths / weekly.deaths_t]
    empirical[, emp_rel := emp / emp_ref]
    empirical[, emp_adj := emp * multiplier]
    empirical[, emp_adj := emp_adj / sum(emp_adj), by = c('date')]
    empirical[, emp_ref_adj :=  emp_ref * multiplier]
    empirical[, emp_ref_adj :=  emp_ref_adj / sum(emp_ref_adj), by = c('date')]
    empirical[, emp_rel_adj :=  emp_adj / emp_ref_adj]
    
    empirical = subset(empirical, age_state_index == which(age_groups == age_group))
    setnames(empirical, 'age_state_index', 'age_index')
    
    empirical = select(empirical, date, age_index, emp, emp_rel, emp_adj, emp_rel_adj)
    
  }

  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[['phi']]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'iterations', 'age_state_index')]
  
  # # find reference
  tmp2 = subset(tmp1, week_index %in% df_week1$week_index)
  tmp2 = tmp2[, list(value_ref = mean(value)), by = c('iterations', 'age_state_index')]
  tmp1 = merge(tmp1, tmp2, by = c('iterations', 'age_state_index'))
  tmp1[, value_rel := value / value_ref]
  
  # adjust by population proportion
  tmp1 = merge(tmp1, pop_data1, by = 'age_state_index')
  tmp1[, value_adj := value * multiplier]
  tmp1[, value_adj := value_adj / sum(value_adj), by = c('iterations', 'week_index')]
  tmp1[, value_ref_adj := value_ref * multiplier]
  tmp1[, value_ref_adj := value_ref_adj / sum(value_ref_adj), by = c('iterations', 'week_index')]
  tmp1[, value_rel_adj := value_adj / value_ref_adj ]
  
  # keep only interested age group 
  tmp1 = subset(tmp1, age_state_index == which(age_groups == age_group))
  
  # take quantiles
  tmp2 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index')]	
  tmp2 = dcast(tmp2, week_index ~ q_label, value.var = "q")
  
  tmp3 = tmp1[, list( 	q= quantile(value_rel, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index')]	
  tmp3 = dcast(tmp3, week_index ~ q_label, value.var = "q")
  setnames(tmp3, p_labs, paste0(p_labs, '_rel'))
  tmp2 = merge(tmp2, tmp3, by = 'week_index')
  
  tmp3 = tmp1[, list( 	q= quantile(value_adj, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index')]	
  tmp3 = dcast(tmp3, week_index ~ q_label, value.var = "q")
  setnames(tmp3, p_labs, paste0(p_labs, '_adj'))
  tmp2 = merge(tmp2, tmp3, by = 'week_index')
  
  tmp1 = tmp1[, list( 	q= quantile(value_rel_adj, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index')]	
  tmp1 = dcast(tmp1, week_index ~ q_label, value.var = "q")
  setnames(tmp1, p_labs, paste0(p_labs, '_rel_adj'))
  tmp1 = merge(tmp2, tmp1, by = 'week_index')
  
  tmp1[, age := age_group]
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, after.10thcumdeaths := date >= date_10thcum]
  
  if(with_empirical){
    tmp1= merge(tmp1, empirical, by = c('date'), all.x = T)
    
  }

  saveRDS(tmp1, file = paste0(outdir, '-Contribution_Age_', age_group,'_', Code, '.rds'))

  # return(tmp1)
}


make_death_ratio_table = function(fit, df_week, df_state_age, data, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples$deaths_predict_state_age_strata) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp2 = subset(tmp1, date <= ref_date)
  tmp2 = tmp2[, list(value.ref = mean(value)), by = c('age_index')]
  tmp1 = merge(tmp1, tmp2, by = c('age_index'))
  tmp1[, value := value / value.ref]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_index', 'week_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_state_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age)]
  tmp1[, code := Code]
  
  # find empirical estimate
  tmp = select(data, weekly.deaths, date, age)
  tmp2 = subset(tmp, date <= ref_date)
  tmp2 = tmp2[, list(weekly.deaths.ref = mean(weekly.deaths)), by = c('age')]
  tmp = merge(tmp, tmp2, by = c('age'))
  tmp[, death.ratio := weekly.deaths / weekly.deaths.ref]
  subset(tmp, age %in% unique(data$age))
  
  tmp1 = merge(tmp1, tmp, by = c('age', 'date'), all.x = T)
  
  # save
  saveRDS(tmp1, file = paste0(outdir, '-DeathRatioTable_', Code, '.rds'))
  
  return(tmp1)
}

find_overall_cumulative_deaths = function(fit_cum, df_week, deaths_predict_var){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit_cum)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'iterations')]
  tmp1 = tmp1[, value := cumsum(value), by = c('iterations')]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index')]	
  tmp1 = dcast(tmp1, week_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  
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
  df_missing = as.data.table( reshape2::melt(stan_data$idx_weeks_missing) )[,-1]
  setnames(df_missing, 1:2, c('idx_serie_missing', 'week_index'))
  df_missing = subset(df_missing, week_index != -1)
  
  tmp = data.table(age_index = stan_data$age_missing,
                   min_count_censored = stan_data$min_count_censored,
                   max_count_censored = stan_data$max_count_censored,
                   sum_count_censored = stan_data$sum_count_censored,
                   idx_serie_missing = 1:length(stan_data$age_missing))
  tmp[min_count_censored == -1, min_count_censored := NA]
  tmp[max_count_censored == -1, max_count_censored := NA]
  tmp[sum_count_censored == -1, sum_count_censored := NA]
  df_missing = merge(df_missing, tmp, by = 'idx_serie_missing')
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
  tmp1 = merge(tmp1, df_missing, by.x = c('age_state_index','week_index'), by.y = c('age_index','week_index'))
  
  # take the cum sum
  tmp1 = tmp1[, value := cumsum(value), by = c('iterations', 'age_state_index')]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, week_index + age_state_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1 = merge(tmp1, df_age_state, by.x = 'age_state_index', by.y = 'age_index')
  tmp1 = merge(tmp1, df_missing, by.x = c('age_state_index','week_index'), by.y = c('age_index','week_index'))
  
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
  df_missing = as.data.table( reshape2::melt(stan_data$idx_weeks_missing) )[,-1]
  setnames(df_missing, 1:2, c('idx_serie_missing', 'week_index'))
  df_missing = subset(df_missing, week_index != -1)
  
  tmp = data.table(age_index = stan_data$age_missing,
                   min_count_censored = stan_data$min_count_censored,
                   max_count_censored = stan_data$max_count_censored,
                   sum_count_censored = stan_data$sum_count_censored,
                   idx_serie_missing = 1:length(stan_data$age_missing))
  tmp[min_count_censored == -1, min_count_censored := NA]
  tmp[max_count_censored == -1, max_count_censored := NA]
  tmp[sum_count_censored == -1, sum_count_censored := NA]
  df_missing = merge(df_missing, tmp, by = 'idx_serie_missing')
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
  tmp1 = merge(tmp1, df_missing, by.x = c('age_state_index','week_index'), by.y = c('age_index','week_index'))
  
  # take the cum sum
  tmp1 = tmp1[, list(value = sum(value)), by = c('age_state_index', 'iterations')]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('age_state_index')]	
  tmp1 = dcast(tmp1,age_state_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_age_state, by.x = 'age_state_index', by.y = 'age_index')
  tmp1 = merge(tmp1, unique(select(df_missing, -week_index)),by.x = 'age_state_index', by.y = 'age_index')
  
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # inside CI?
  tmp1[is.na(sum_count_censored), inside.CI := (CL >= min_count_censored & CL <= max_count_censored) | (CU >= min_count_censored & CU <= max_count_censored) | (CL <= min_count_censored & CU >= max_count_censored)]
  tmp1[!is.na(sum_count_censored), inside.CI := sum_count_censored >= CL & sum_count_censored <= CU]
  
  # save
  saveRDS(tmp1, file = paste0(outdir, '-predictive_checks_table2_', Code, '.rds'))
  
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

find_cumulative_deaths_givensum_state_age = function(fit, date_10thcum, df_week, df_age_continuous, data_comp, cum.death.var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  # df age
  df_age_state = data.table(age = unique(subset(data_comp, code == Code)$age))
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
  df_week1 = subset(df_week, date >= date_10thcum)
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1 = subset(df_week1, month %in% min(df_week1$month):(min(df_week1$month) + 2))
  
  # week index
  df_week_10thcum = subset(df_week, date >= date_10thcum)
  tmp2 = as.data.table( subset(data_comp, code == Code))
  tmp2[, date := as.Date(date)]
  tmp2 = tmp2[, list(cum.death = sum(get(cum.death.var))), by = 'date']
  tmp2 = tmp2[order( date)]
  dummy = min(df_week_10thcum$date) %in% tmp2$date 
  last.cum.deaths = max(subset(tmp2, date == ifelse(dummy, min(df_week_10thcum$date) , 
                                                    (min(tmp2$date):(min(tmp2$date)+6))[min(tmp2$date):(min(tmp2$date)+6) %in% df_week_10thcum$date] ) )$cum.death )
  tmp2 = unique(subset(tmp2, date %in% c(max(df_week_10thcum$date) + 7, df_week_10thcum$date)))
  tmp2[, weekly.deaths := c(diff(cum.death), NA)]
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = 'date')
  tmp2 = select(tmp2, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))
  
  # adjust dfweek for cum
  df_week_adj = df_week[, list(date = date + 7), by = 'week_index']
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples['phi']) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # add cumsum from weeks before 10th date
  tmp3 = subset(tmp1, week_index %in% df_week1$week_index)
  tmp3 = tmp3[, list(value_ref = mean(value)), by = c('age_index', 'iterations')]
  tmp3[, value_ref := value_ref * last.cum.deaths]
  
  # multiply by the weekly death
  tmp1 = merge(tmp1, tmp2, by = 'week_index')
  tmp1[, value := value * weekly.deaths]
  
  # find cumulative death
  tmp1 = merge(tmp1, tmp3, by = c('age_index', 'iterations'))
  tmp1[, value := cumsum(value), by = c('iterations', 'age_index')]
  tmp1[, value := value + value_ref]
  
  # take quantiles
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week_adj, by = 'week_index')
  tmp1 = merge(tmp1, df_age_state, by = 'age_index')
  
  saveRDS(list(data_comp, tmp1), file = paste0(outdir, '-CumDeathsComp_', 'ScrapedData', '_', Code, '.rds'))
  
  return(tmp1)
}

find_phi_state_age = function(fit, df_week, df_age_continuous, age_state){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
  # df age
  df_age_state = data.table(age = age_state)
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
  tmp1 = as.data.table( reshape2::melt(fit_samples['phi']) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
  
  # take the cum sum
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, week_index + age_state_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1 = merge(tmp1, df_age_state, by.x = 'age_state_index', by.y = 'age_index')
  
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
  df_week1 = subset(df_week, date >= date_10thcum)
  df_week1[, month := as.numeric(format(date, '%m'))]
  df_week1[grepl('2021', date), month := month + 12]
  df_week1 = subset(df_week1, month %in% min(df_week1$month):(min(df_week1$month) + 2))
  
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
  pop_data1 = subset(pop_data1, code == Code)
  
  # week index
  df_week_10thcum = subset(df_week, date >= date_10thcum)
  data_comp = as.data.table( subset(data_comp, code == Code))
  data_comp[, date := as.Date(date)]
  tmp2 = data_comp[, list(cum.death = sum(get(cum.death.var))), by = 'date']
  tmp2 = tmp2[order( date)]
  dummy = min(df_week_10thcum$date)  %in% tmp2$date
  last.cum.deaths = max(subset(tmp2, date == ifelse(dummy, min(df_week_10thcum$date) , 
                                                    (min(tmp2$date):(min(tmp2$date)+6))[min(tmp2$date):(min(tmp2$date)+6) %in% df_week_10thcum$date] ) )$cum.death )
  tmp2 = unique(subset(tmp2, date %in% c(df_week_10thcum$date, max(df_week_10thcum$date)+7)))
  tmp2[, weekly.deaths := c(diff(cum.death),NA)]
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = 'date')
  tmp2 = select(tmp2, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))
  
  # adjust dfweek for cum
  df_week_adj = df_week[, list(date = date + 7), by = 'week_index']
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples['phi']) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # add cumsum from weeks before 10th date
  tmp3 = subset(tmp1, week_index %in% df_week1$week_index)
  tmp3 = tmp3[, list(value_ref = mean(value)), by = c('age_index', 'iterations')]
  tmp3[, value_ref := value_ref * last.cum.deaths]
  
  # multiply by the weekly death
  tmp1 = merge(tmp1, tmp2, by = 'week_index')
  tmp1[, value := value * weekly.deaths]
  
  # find cumulative death
  tmp1 = merge(tmp1, tmp3, by = c('age_index', 'iterations'))
  tmp1[, value := cumsum(value), by = c('iterations', 'age_index')]
  tmp1[, value := value + value_ref]
  
  # find mortality rate
  tmp1 = merge(tmp1, pop_data1, by = 'age_index')
  tmp1[, value := value / pop]
  
  # quantiles
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week_adj, by = 'week_index')
  tmp1[, age := df_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_age$age)]
  tmp1[, code := Code]
  
  saveRDS(tmp1, file = paste0(outdir, '-', 'MortalityRate', 'Table_', Code, '.rds'))
  
  return(tmp1)
}

make_weekly_death_rate_table = function(fit_cum, fiveagegroups, date_vac, df_week, data_comp, data, df_age_continuous, cum.death.var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit_cum)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
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
  
  # ref week
  data_comp = as.data.table( subset(data_comp, code == Code))
  date_max = subset(data_comp, date >= date_vac)
  date_max = date_max[which.max(daily_deaths), date-7]
  cat('Date max is ', as.character(date_max))
  df_week1  = df_week[ date %in% date_max:(date_max + 6)]
  # df_week1[, month := as.numeric(format(date, '%m'))]
  # df_week1[grepl('2021', date), month := month + 12]
  # df_week1 = subset(df_week1, month %in% min(df_week1$month):(min(df_week1$month) + 1))
  # 
  # week index
  data_comp[, date := as.Date(date)]
  tmp2 = data_comp[, list(cum.death = sum(get(cum.death.var))), by = 'date']
  tmp2 = tmp2[order( date)]
  tmp2 = unique(subset(tmp2, date %in% c(df_week$date, max(df_week$date)+7)))
  tmp2[, weekly.deaths := c(diff(cum.death),NA)]
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = 'date')
  tmp2 = select(tmp2, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))

  # empirical
  df_age_reporting[, age_state_index := which(df_age$age_from <= age_from & df_age$age_to >= age_to), by = 'age_index']
  empirical = merge(data, df_age_reporting, by = 'age_from')
  empirical = merge(empirical, df_week, by = 'date')
  empirical = empirical[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c( 'age_state_index', 'week_index')]
  tmp_emp = empirical[, list(total_deaths = sum(na.omit(weekly.deaths))), by = 'week_index']
  empirical = merge(empirical, tmp_emp, by = c('week_index'))
  empirical[, prop_deaths := weekly.deaths / total_deaths]
  empirical = merge(select(empirical, week_index, age_state_index, prop_deaths), tmp2, by = c('week_index'))
  setnames(empirical, 'weekly.deaths', 'weekly.deaths_t')
  empirical[, weekly.deaths := weekly.deaths_t *  prop_deaths, by = c( 'age_state_index', 'week_index')]
  tmp_emp = subset(empirical, week_index %in% df_week1$week_index)
  tmp_emp = tmp_emp[, list(weekly.deaths_ref = mean(na.omit(weekly.deaths))), by = 'age_state_index']
  empirical= merge(empirical, tmp_emp, by = 'age_state_index')
  empirical[, emp_est := weekly.deaths / weekly.deaths_ref]
  setnames(empirical, 'age_state_index', 'age_index')
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples['phi']) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # multiply by the weekly death
  tmp1 = merge(tmp1, tmp2, by = 'week_index')
  tmp1[, value := value * weekly.deaths]
  
  # find ratio relative to baseline 
  tmp3 = subset(tmp1, week_index %in% df_week1$week_index)
  tmp3 = tmp3[, list(value_ref = mean(value)), by = c('age_index', 'iterations')]
  tmp1 = merge(tmp1, tmp3, by = c('age_index', 'iterations'))
  tmp1[, value := value / value_ref]
  
  # quantiles
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_age$age)]

  tmp1 = merge(tmp1, empirical, by = c('age_index', 'week_index'), all.x = T)
  
  saveRDS(tmp1, file = paste0(outdir, '-', 'DeathRatioWinter', 'Table_', Code, '.rds'))
  
  return(tmp1)
}

make_weekly_death_rate_other_source = function(fit_cum, df_week, data_comp, cum.death.var, var.phi, df_age, outdir, age_groups = NULL, lab = NULL){
  
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
  tmp2 = data_comp[, list(cum.death = sum(get(cum.death.var))), by = 'date']
  tmp2 = tmp2[order( date)]
  tmp2 = unique(subset(tmp2, date %in% c(df_week$date, max(df_week$date)+7)))
  tmp2[, weekly.deaths := c(diff(cum.death),NA)]
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = 'date')
  tmp2 = select(tmp2, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))
  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var.phi]]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  if(!is.null(age_groups)){
    tmp1 = merge(tmp1, df_age_continuous, 'age_index')
    tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
    setnames(tmp1, 'age_state_index', 'age_index')
  }
  
  # multiply by the weekly death
  tmp1 = merge(tmp1, tmp2, by = 'week_index')
  tmp1[, value := value * weekly.deaths]
  
  # quantiles
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), 
              by=c('week_index', 'age_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_age$age[age_index]]
  tmp1[, age := factor(age, levels = df_age$age)]
  
  if(!is.null(lab)){
    file = paste0(outdir, '-', 'DeathByAge', 'Table_', lab, '_', Code, '.rds')
    
  }else{
    file =  paste0(outdir, '-', 'DeathByAge', 'Table_', var.phi, '_', Code, '.rds')
  }

  saveRDS(tmp1, file =file)

  return(tmp1)
}


# 
# add_intervals_missing_data = function(data, stan_data, base_week_idx = 0){
#   
#   # include missing information
#   df_missing = as.data.table( reshape2::melt(stan_data$idx_weeks_missing) )[,-1]
#   setnames(df_missing, 1:2, c('idx_serie_missing', 'week_index'))
#   df_missing = subset(df_missing, week_index != -1)
#   
#   tmp = data.table(age_index = stan_data$age_missing,
#                    min_count_censored = stan_data$min_count_censored,
#                    max_count_censored = stan_data$max_count_censored,
#                    sum_count_censored = stan_data$sum_count_censored,
#                    idx_serie_missing = 1:length(stan_data$age_missing))
#   tmp[min_count_censored == -1, min_count_censored := NA]
#   tmp[max_count_censored == -1, max_count_censored := NA]
#   tmp[sum_count_censored == -1, sum_count_censored := NA]
#   
#   df_missing = merge(df_missing, tmp, by = 'idx_serie_missing')
#   
#   tmp1 = merge(data, df_missing, by = c('week_index', 'age_index'), all.x = T)
#   
#   return(tmp1)
# }

# find_cumulative_deaths_state_age = function(fit, df_week, df_age_continuous, state_age_groups, deaths_predict_var){
#   
#   ps <- c(0.5, 0.025, 0.975)
#   p_labs <- c('M','CL','CU')
#   
#   if(is.null(fit)) return()
#   
#   # extract samples
#   fit_samples = rstan::extract(fit_cum)
#   
#   # df age
#   df_age_state = data.table(age = state_age_groups)
#   df_age_state[, age_index := 1:nrow(df_age_state)]
#   df_age_state[, age_from := gsub('(.+)-.*', '\\1', age)]
#   df_age_state[, age_to := gsub('.*-(.+)', '\\1', age)]
#   df_age_state[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
#   df_age_state[grepl('\\+', age_to), age_to := max(df_age_continuous$age)]
#   df_age_state[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
#   df_age_state[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
#   df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
#   df_age_continuous[, age_state_index := which(df_age_state$age_from_index <= age_index & df_age_state$age_to_index >= age_index), by = 'age_index']
#   
#   
#   # tmp1
#   tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
#   setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
#   
#   # sum by state age group
#   tmp1 = merge(tmp1, df_age_continuous, 'age_index')
#   tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
#   
#   # take the cum sum
#   tmp1 = tmp1[, value := cumsum(value), by = c('iterations', 'age_state_index')]
#   tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
#                        q_label=p_labs), 
#               by=c('week_index', 'age_state_index')]	
#   tmp1 = dcast(tmp1, week_index + age_state_index ~ q_label, value.var = "q")
#   
#   tmp1[, code := Code]
#   
#   tmp1 = merge(tmp1, df_week, by = 'week_index')
#   tmp1 = merge(tmp1, df_age_state, by.x = 'age_state_index', by.y = 'age_index')
#   
#   return(tmp1)
# }
