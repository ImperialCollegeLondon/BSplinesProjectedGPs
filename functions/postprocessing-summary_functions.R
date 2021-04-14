make_predictive_checks_table = function(fit, df_week, df_state_age, data, deaths_predict_var, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[deaths_predict_var]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps),
                       q_label=p_labs), 
              by=c('age_index', 'week_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_state_age$age_cat[age_index]]
  
  tmp1 = merge(tmp1, data, by = c('date', 'age'))
  tmp1[, age := factor(age, levels = levels(data$age))]
  
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
  
  # save
  saveRDS(eff_sample_size_cum, file = paste0(outdir, "-eff_sample_size_cum_", Code, ".rds"))
  saveRDS(Rhat_cum, file = paste0(outdir, "-Rhat_cum_", Code, ".rds"))
  saveRDS(WAIC, file = paste0(outdir, "-WAIC_", Code, ".rds"))
  saveRDS(LOO, file = paste0(outdir, "-LOO_", Code, ".rds"))
}

make_probability_ratio_table = function(fit, df_week, df_state_age, data1, data2, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) stop()
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples$probability_ratio_age_strata) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps),
                       q_label=p_labs), 
              by=c('age_index', 'week_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_state_age$age_cat[age_index]]
  tmp1[, age := factor(age, levels = df_state_age$age_cat)]
  
  tmp1[, code := Code]
  
  # find significant shift
  tmp1[, shift := F]
  tmp1[CL > 1, shift := T]
  tmp1[CU < 1, shift := T]
  tmp2 = tmp1[, list(shift.tot = (sum(shift) > 0)), by = 'age']
  tmp1 = merge(tmp1, tmp2, by = 'age')
  
  # find empirical estimate
  tmp = select(rbind(data1,data2), COVID.19.Deaths, date, age)
  tmp2 = tmp[, list(total.deaths = sum(na.omit(COVID.19.Deaths))), by = 'date']
  tmp = merge(tmp, tmp2, by = 'date')
  tmp[, emp.prob := COVID.19.Deaths / total.deaths]
  tmp2 = subset(tmp, date == min(data2$date))
  setnames(tmp2, 'emp.prob', 'emp.prob.ref')
  tmp = merge(tmp, select(tmp2, -COVID.19.Deaths, -date, -total.deaths), by = c('age'))
  tmp[, emp.prob.ratio := emp.prob / emp.prob.ref]
  subset(tmp, age %in% unique(data2$age))
  
  tmp1 = merge(tmp1, select(tmp, -COVID.19.Deaths), by = c('age', 'date'), all.x = T)
  
  # save
  saveRDS(tmp1, file = paste0(outdir, '-ProbabilityRatioTable_', Code, '.png'))
  
  return(tmp1)
}
