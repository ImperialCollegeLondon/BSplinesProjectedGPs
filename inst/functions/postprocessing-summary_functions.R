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
  tmp1[, inside.CI := daily.deaths <= CU & daily.deaths >= CL]
  
  # save
  saveRDS(tmp1, file = paste0(outdir, '-predictive_checks_table_', Code, '.rds'))
  
  return(tmp1)
}

add_intervals_missing_data = function(data, stan_data, base_week_idx = 0){
  
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
  
  tmp1 = merge(data, df_missing, by = c('week_index', 'age_index'), all.x = T)

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
  
  # compute WAIC and LOO
  re = rstan::extract(fit_cum)
  if('log_lik' %in% names(re)){
    log_lik <- loo::extract_log_lik(fit_cum)
    .WAIC = loo::waic(log_lik)
    .LOO = loo::loo(log_lik)
    print(.WAIC); print(.LOO)
    WAIC = .WAIC$pointwise
    LOO = .LOO$pointwise
  } 
  
  sampler_params <- get_sampler_params(fit_cum, inc_warmup = FALSE)
  sampler_diagnostics <- data.table()
  for (i in colnames(sampler_params[[1]])) {
    tmp <- data.table(t(sapply(sampler_params, function(x) quantile(x[, i],probs = c(0.025,0.5,0.975)))))
    tmp[, diagnostics:=i ]
    tmp[, chain:= seq_len(length(sampler_params))]
    sampler_diagnostics <- rbind(sampler_diagnostics, tmp)
  }
  print(sampler_diagnostics)
  
  check_all_diagnostics(fit, outdir)
  
  # save
  saveRDS(eff_sample_size_cum, file = paste0(outdir, "-eff_sample_size_cum_", Code, ".rds"))
  saveRDS(Rhat_cum, file = paste0(outdir, "-Rhat_cum_", Code, ".rds"))
  saveRDS(WAIC, file = paste0(outdir, "-WAIC_", Code, ".rds"))
  saveRDS(LOO, file = paste0(outdir, "-LOO_", Code, ".rds"))
  saveRDS(sampler_diagnostics, file = paste0(outdir, "-sampler_diagnostics_", Code, ".rds"))
}

check_all_diagnostics <- function(fit, outdir) {
  check_n_eff(fit)
  check_rhat(fit)
  n_div <- check_div(fit)
  n_treedepth <- check_treedepth(fit,15)
  check_energy(fit)
  
  cat('\nn_div',n_div, '\n' )
  cat('\nn_treedepth',n_treedepth, '\n' )
  saveRDS(list(n_div,n_treedepth), file=paste0(outdir,'-diagnostics_', Code, '.rds'))
}

check_treedepth <- function(fit, max_depth = 10) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  print(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                n, N, max_depth, 100 * n / N))
  if (n > 0)
    print('  Run again with max_depth set to a larger value to avoid saturation')
  return(n)
}

check_energy <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)*2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('E-BFMI indicated no pathological behavior')
  else
    print('  E-BFMI below 0.2 indicates you may need to reparameterize your model')
}

check_n_eff <- function(fit) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
  # remove na neff
  fit_summary <- fit_summary[!is.na(fit_summary[,5]),]
  N <- dim(fit_summary)[[1]]
  
  iter <- dim(rstan::extract(fit)[[1]])[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    ratio <- fit_summary[,5][n] / iter
    if (ratio < 0.001) {
      print(sprintf('n_eff / iter for parameter %s is %s!',
                    rownames(fit_summary)[n], ratio))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('n_eff / iter looks reasonable for all parameters')
  else
    print('  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated')
}

check_div <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  
  print(sprintf('%s of %s iterations ended with a divergence (%s%%)',
                n, N, 100 * n / N))
  if (n > 0)
    print('  Try running with larger adapt_delta to remove the divergences')
  # return iterations ended with a divergence
  return(n)
}

check_rhat <- function(fit) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
  # remove na neff
  fit_summary <- fit_summary[!is.na(fit_summary[,6]),]
  N <- dim(fit_summary)[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[,6][n]
    if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
      print(sprintf('Rhat for parameter %s is %s!',
                    rownames(fit_summary)[n], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('Rhat looks reasonable for all parameters')
  else
    print('  Rhat above 1.1 indicates that the chains very likely have not mixed')
}

make_probability_ratio_table = function(fit, df_week, df_state_age, data, outdir){
  
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
  tmp = select(data, daily.deaths, date, age)
  tmp2 = tmp[, list(total.deaths = sum(na.omit(daily.deaths))), by = 'date']
  tmp = merge(tmp, tmp2, by = 'date')
  tmp[, emp.prob := daily.deaths / total.deaths]
  tmp2 = subset(tmp, date <= ref_date)
  tmp2 = tmp2[, list(emp.prob = mean(emp.prob)), by = c('age')]
  setnames(tmp2, 'emp.prob', 'emp.prob.ref')
  tmp = merge(tmp, tmp2, by = c('age'))
  tmp[, emp.prob.ratio := emp.prob / emp.prob.ref]
  subset(tmp, age %in% unique(data$age))
  
  tmp1 = merge(tmp1, select(tmp, -daily.deaths), by = c('age', 'date'), all.x = T)
  
  # save
  saveRDS(tmp1, file = paste0(outdir, '-ProbabilityRatioTable_', Code, '.rds'))

  return(tmp1)
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
  tmp = select(data, daily.deaths, date, age)
  tmp2 = subset(tmp, date <= ref_date)
  tmp2 = tmp2[, list(daily.deaths.ref = mean(daily.deaths)), by = c('age')]
  tmp = merge(tmp, tmp2, by = c('age'))
  tmp[, death.ratio := daily.deaths / daily.deaths.ref]
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

find_cumulative_deaths_state_age = function(fit, df_week, df_age_continuous, state_age_groups, deaths_predict_var){
  
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
  df_age_continuous[, age_index := 1:nrow(df_age_continuous)]
  df_age_continuous[, age_state_index := which(df_age_state$age_from_index <= age_index & df_age_state$age_to_index >= age_index), by = 'age_index']

  
  # tmp1
  tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'age_state_index', 'iterations')]
  
  # take the cum sum
  tmp1 = tmp1[, value := cumsum(value), by = c('iterations', 'age_state_index')]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, week_index + age_state_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1 = merge(tmp1, df_age_state, by.x = 'age_state_index', by.y = 'age_index')
  
  return(tmp1)
}

find_sum_missing_deaths_state_age = function(fit, df_week, df_age_continuous, state_age_groups, stan_data, deaths_predict_var){
  
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
  tmp1 = tmp1[, value := cumsum(value), by = c('iterations', 'age_state_index')]
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_state_index')]	
  tmp1 = dcast(tmp1, week_index + age_state_index ~ q_label, value.var = "q")
  
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1 = merge(tmp1, df_age_state, by.x = 'age_state_index', by.y = 'age_index')
  tmp1 = merge(tmp1, df_missing, by.x = c('age_state_index','week_index'), by.y = c('age_index','week_index'))
  
  
  return(tmp1)
}


find_sum_bounded_missing_deaths_state_age = function(fit, df_age_continuous, state_age_groups, stan_data, deaths_predict_var){
  
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
  tmp1 = merge(tmp1, df_missing,by.x = 'age_state_index', by.y = 'age_index')
  
  
  return(tmp1)
}
