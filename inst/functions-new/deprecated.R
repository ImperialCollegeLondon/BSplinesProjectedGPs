find_contribution_multiple_age_groups = function(fit, df_week, df_age_continuous, age_groups, lab,outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return()
  
  # extract samples
  fit_samples = rstan::extract(fit_cum)
  
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
  tmp1 = as.data.table( reshape2::melt(fit_samples[['phi']]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  
  # sum by state age group
  tmp1 = merge(tmp1, df_age_continuous, 'age_index')
  tmp1 = tmp1[, list(value = sum(value)), by = c('week_index', 'iterations', 'age_state_index')]
  setnames(tmp1, 'age_state_index', 'age_index')
  
  # take quantiles
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('week_index', 'age_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1[, age := age_groups[age_index]]
  tmp1[, code := Code]
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  
  saveRDS(tmp1, file = paste0(outdir, '-Contribution_Age_', lab,'_', Code, '.rds'))
  
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
  date_max = subset(data_comp, date >= date_vac & date <= max(data$date)-7)
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


# 
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

plot_convergence_diagnostics = function(fit, title, suffix, outfile)
{
  
  stopifnot(!is.null(fit))
  
  posterior <- as.array(fit)
  
  pars = list("rho", "sigma", "beta", 'a_age', 'a0', c('gamma', 'tau', 'p'), 'nu', 'lambda', 'tau',
              'alpha_gp', 'delta')
  
  for(j in 1:length(pars))
  {                     
    
    par = pars[[j]]
    
    if(sum(grepl(par, names(fit))) == 0 | sum(grepl(par, names(fit))) > 60) next
    
    p_trace = bayesplot::mcmc_trace(posterior, regex_pars = par) + labs(title = title) 
    ggsave(p_trace, file = paste0(outfile, "-convergence_diagnostics-", "trace_plots_", suffix, '_', Code, "_", par,".png") , w= 10, h = 10)
    
    p_intervals = bayesplot::mcmc_intervals(posterior, regex_pars = par, prob = 0.95,
                                            prob_outer = 0.95) + labs(title = title)
    ggsave(p_intervals, file = paste0(outfile, "-convergence_diagnostics-", "intervals_plots_",  suffix, '_',Code, "_", par,".png") , w=10, h = 10)
    
    if(length(par) > 1){
      p_pairs = gridExtra::arrangeGrob(bayesplot::mcmc_pairs(posterior, regex_pars = par), top = title)
      ggsave(p_pairs, file = paste0(outfile, "-convergence_diagnostics-", "pairs_plots_",  suffix, '_',Code, "_", par,".png"), w= 15, h = 15)
      
    }
    
  }
  
}
plot_ratio_contribution = function(ratio_contribution, outdir)
{
  ggplot(ratio_contribution, aes(x = date, y = age)) + 
    geom_raster(aes(fill = M)) + 
    facet_wrap(~loc_label) + 
    theme_bw() + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90)) + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, trans = 'log', breaks = c(0.2, 1, 20)) + 
    labs(x= '', y = 'Age groups', fill = 'Estimated posterior median')
  ggsave(file = paste0(outdir, '-ProbabilityRatio.png'), w = 9, h = length(unique(ratio_contribution$loc_label)) / 4, limitsize = F)
  
  tmp = subset(ratio_contribution, age %in% c('35-44', '45-54', '55-64', '65-74', '75-84', '85+'))
  # tmp = subset(tmp, date > as.Date('2020-06-01'))
  ggplot(tmp, aes(x = date, y = age)) + 
    geom_raster(aes(fill = M)) + 
    facet_wrap(~loc_label) + 
    theme_bw() + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, trans = 'log', breaks = c(0.2, 1, 5)) + 
    labs(x= '', y = 'Age groups', fill = 'Estimated posterior median')
  ggsave(file = paste0(outdir, '-ProbabilityRatio_ederly.png'), w = 9, h = length(unique(ratio_contribution$loc_label)) / 4, limitsize = F)
  
}
plot_imputed_deaths_by_age = function(tmp1, var_name, data, outdir, discrete = F, wo_ytext = F,  wo_legend = F){
  
  if(!discrete){
    tmp1[, age := as.numeric(age)]
  }
  
  p = ggplot(tmp1, aes(x = date, y = age )) + 
    theme_bw() +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    geom_raster(aes(fill = M))  + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  
  if(discrete){
    limits_wd = range(na.omit(data$weekly.deaths))
    range_wd = sqrt(range(tmp1$M))
    digits_cut= ifelse(range_wd[2]/100 > 1, 2, 1)
    digits_cut= ifelse(range_wd[2]/10 > 1, digits_cut, 0)
    p = p + 
      scale_y_discrete(expand = c(0,0))  + labs(x = '', y = 'Age group', fill = 'Estimated posterior median') + 
      scale_fill_viridis_c(trans = 'sqrt', limits = limits_wd, breaks = round(seq(range_wd[2]/2, range_wd[2], length.out = 3)^2, digits = -digits_cut)) 
  } else {
    p = p + 
      scale_y_continuous(expand = c(0,0))  +labs(x = '', y = 'Age', fill = 'Estimated posterior median') + 
      scale_fill_viridis_c()
  }
  
  if(wo_ytext)
    p = p + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  if(wo_legend)
    p = p + theme(legend.position = 'none')
  
  if(!is.null(outdir)){
    ggsave(p, file = paste0(outdir, "-posterior_", var_name, '_', Code, ".png") , w= 5, h = 5.2, limitsize = FALSE)
    
  }
  
  return(p)
}

compare_CDCestimation_JHU_DoH_plot = function(CDC_data, JHU_data, scraped_data, var.cum.deaths.CDC, outdir)
{
  # prepare JHU data
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  JHUData = subset(JHUData, code == Code)
  JHUData[, date := as.Date(date)]
  JHUData[, CL := NA]
  JHUData[, CU := NA]
  
  # prepare DoH data
  scraped_data = select(as.data.table(scraped_data), code, date, cum.deaths, age)
  scraped_data = scraped_data[, list(cum.deaths = sum(cum.deaths)), by = c('code', 'date')]
  scraped_data = subset(scraped_data, code == Code)
  setnames(scraped_data, 'cum.deaths', 'cumulative_deaths')
  scraped_data[, date := as.Date(date)]
  scraped_data[, CL := NA]
  scraped_data[, CU := NA]
  
  # prepare CDC estimations
  CDC_data = select(as.data.table(CDC_data), code, date, var.cum.deaths.CDC, CL, CU)
  setnames(CDC_data, var.cum.deaths.CDC, 'cumulative_deaths')
  
  # plot
  JHUData[, source := 'JHU']
  scraped_data[, source := 'DoH']
  CDC_data[, source := 'CDC']
  
  tmp2 = rbind(rbind(JHUData, CDC_data), scraped_data)
  tmp2 = subset(tmp2, code %in% unique(CDC_data$code) & date <= max(CDC_data$date))
  
  Limits = range(JHUData$cumulative_deaths)
  
  p = ggplot(tmp2, aes(x = date, y = cumulative_deaths)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = source), alpha = 0.5) +
    geom_line(aes(col = source), size = 1) +
    theme_bw() + 
    scale_y_continuous(limits = Limits) +
    scale_color_viridis_d(option = "B", direction = -1, end = 0.8) + 
    scale_fill_viridis_d(option = "B", direction = -1, end = 0.8) + 
    labs(x = '', y = 'Cumulative COVID-19 \ndeaths', col = '', fill = '') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90)) 
  ggsave(p, file = paste0(outdir, '-comparison_JHU_CDC_DoH_uncertainty_', Code, '.png'), w = 4, h = 4, limitsize = F)
}
plot_covariance_matrix = function(fit_cum, outdir)
{
  
  if(is.null(stan_data$Adj)) return()
  if(!is.null(stan_data$node1)) return()
  if(!is.null(stan_data$num_basis_rows)) return()
  
  
  samples = extract(fit_cum)
  
  D = diag( apply(stan_data$Adj, 2, sum) )
  tau_m = median(samples$tau)
  p_m = median(samples$p)
  cov_m = tau_m * solve(D - p_m * stan_data$Adj)
  
  tmp1 = as.data.table( reshape2::melt( stan_data$Adj ))
  ggplot(tmp1, aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill = as.factor(value)))+
    scale_fill_manual(values = c('beige',"blue")) + 
    labs(x = expression(s[j]), y = expression(s[i]), fill = '') + 
    scale_y_reverse(expand = c(0,0))  +
    scale_x_continuous(expand = c(0,0)) +
    theme(legend.position = 'none')
  ggsave(file = paste0(outdir, '-AdjacencyMatrix_', Code, '.png'), w = 4, h = 4)
  
  tmp1 = as.data.table( reshape2::melt( D ))
  ggplot(tmp1, aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill = value)) +
    labs(x = expression(s[j]), y = expression(s[i]), fill = '') + 
    scale_y_reverse(expand = c(0,0))  +
    scale_x_continuous(expand = c(0,0)) +
    theme(legend.position = 'none') + 
    scale_fill_viridis_c() 
  ggsave(file = paste0(outdir, '-DMatrix_', Code, '.png'), w = 4, h = 4)
  
  tmp1 = as.data.table( reshape2::melt( cov_m ))
  range_value = range(tmp1$value)
  ggplot(tmp1, aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill = value)) +
    labs(x = expression(s[j]), y = expression(s[i]), fill = 'Estimated posterior value') + 
    scale_y_reverse(expand = c(0,0))  +
    scale_x_continuous(expand = c(0,0)) +
    theme(legend.position = 'none') + 
    scale_fill_viridis_c(begin = 0, end = 1, limits = range_value, breaks = seq(0,0.6,0.2)) 
  ggsave(file = paste0(outdir, '-CovarianceMatrix_', Code, '.png'), w = 4, h = 4)
  
  map = as.data.table(expand.grid(idx_week = 1:stan_data$W, idx_basis = 1:stan_data$num_basis ))
  set(map, NULL, 'idx', 1:nrow(map))
  tmp1 = merge(tmp1, map, by.x = 'Var1', by.y = 'idx')
  setnames(tmp1, c('idx_week', 'idx_basis'), c('idx_week_column', 'idx_basis_column'))
  tmp1 = merge(tmp1, map, by.x = 'Var2', by.y = 'idx')
  setnames(tmp1, c('idx_week', 'idx_basis'), c('idx_week_row', 'idx_basis_row'))
  
  tmp2 = subset(tmp1, idx_week_row %in% 1 & idx_week_column %in% 1)
  ggplot(tmp2, aes(x = idx_basis_row, y = idx_basis_column)) + 
    geom_raster(aes(fill = value)) +
    labs(x = 'basis function index', y = 'basis function index', fill = 'Estimated posterior value') + 
    scale_y_reverse(expand = c(0,0), breaks = seq(1, 10, 2))  +
    scale_x_continuous(expand = c(0,0), breaks = seq(1, 10, 2)) +
    theme(legend.position = 'none') + 
    scale_fill_viridis_c(begin = 0, end = 1, limits = range_value, breaks = seq(0,0.6,0.2)) 
  ggsave(file = paste0(outdir, '-CovarianceMatrix_w1_', Code, '.png'), w = 4, h = 4)
  
  tmp2 = subset(tmp1, idx_basis_column %in% 5 & idx_basis_row %in% 5)
  ggplot(tmp2, aes(x = idx_week_row, y = idx_week_column)) + 
    geom_raster(aes(fill = value)) +
    labs(x = 'week index', y = 'week index', fill = 'Estimated posterior value') + 
    scale_y_reverse(expand = c(0,0))  +
    scale_x_continuous(expand = c(0,0)) +
    theme(legend.position = 'none') + 
    scale_fill_viridis_c(begin = 0, end = 1, limits = range_value, breaks = seq(0,0.6,0.2)) 
  ggsave(file = paste0(outdir, '-CovarianceMatrix_k5_', Code, '.png'), w = 4, h = 4)
  
}

plot_probability_ratio = function(probability_ratio_table, df_week, stan_data, outdir)
{
  dates = format( range( subset(df_week, week_index %in% stan_data$w_ref_index)$date ), "%d-%b-%Y")
  tmp = subset(probability_ratio_table, age %in% c('55-64', '65-74', '75-84', '85+'))
  
  # plot
  ggplot(probability_ratio_table, aes(x = date)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line(aes(y = M)) + 
    geom_point(aes(y = emp.prob.ratio), col = 'darkred') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    coord_cartesian(ylim = c(0, max(tmp$CU))) +
    theme_bw() + 
    geom_hline(aes(yintercept = 1), col = 'darkblue') +
    facet_wrap(~age,  ncol =1)+ 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 90)) +
    labs(x = '', y = paste0('Ratio of the age-specific contribution to COVID-19 weekly deaths relative to its mean between ', dates[1], ' and ', dates[2] ))
  ggsave(file = paste0(outdir, '-ProbabilityRatio_', Code, '.png'), w = 4, h = 14)
  
  ggplot(tmp, aes(x = date)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line(aes(y = M)) + 
    geom_point(aes(y = emp.prob.ratio), col = 'darkred') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    theme_bw() + 
    geom_hline(aes(yintercept = 1)) +
    facet_wrap(~age,  ncol =1)+ 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 90)) +
    labs(x = '', y = paste0('Ratio of the age-specific contribution to COVID-19 weekly deaths relative to its mean between ', dates[1], ' and ', dates[2] ))
  ggsave(file = paste0(outdir, '-ProbabilityRatio_elderly_', Code, '.png'), w = 4, h = 8)
  
  tmp = subset(probability_ratio_table, age %in% c('55-64', '65-74', '75-84', '85+'))
  ggplot(tmp, aes(x = date)) + 
    geom_line(aes(y = M, col = age)) + 
    # geom_point(aes(y = emp.prob.ratio, col = age)) +
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.1) +
    theme_bw() + 
    geom_hline(aes(yintercept = 1), col = 'grey30', size = 0.3) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = '', y = paste0('Percentage change in weekly deaths age-specific composition\ncompared to the baseline period'), col = 'Age group', 
         fill = 'Age group') + 
    scale_color_viridis_d(option = 'B', end = 0.9, begin = 0.2) + 
    scale_fill_viridis_d(option = 'B', end = 0.9, begin = 0.2) + 
    theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1), 
          panel.grid.major= element_blank(), legend.box="vertical", 
          legend.margin=margin(), 
          legend.position = 'bottom') + 
    # geom_vline(data = tmp1, aes(xintercept = date, linetype = prop_name)) + 
    guides(color = guide_legend(order=1), fill = guide_legend(order=1)) 
  ggsave(file = paste0(outdir, '-ProbabilityRatio_elderly_median_', Code, '.png'), w = 6, h = 6)
  
}
plot_death_ratio = function(death_ratio_table, outdir)
{
  # plot
  tmp = subset(death_ratio_table, age %in% c( '55-64', '65-74', '75-84', '85+'))
  
  ggplot(tmp, aes(x = date)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line(aes(y = M)) + 
    geom_point(aes(y = death.ratio), col = 'darkred') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    theme_bw() + 
    geom_hline(aes(yintercept = 1), col = 'darkblue') +
    facet_wrap(~age, ncol = 1)+ 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 90))
  ggsave(file = paste0(outdir, '-DeathRatio_', Code, '.png'), w = 4, h = 8)
  
  ggplot(tmp, aes(x = date)) + 
    geom_line(aes(y = M, col = age)) + 
    geom_point(aes(y = death.ratio, col = age)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.1) +
    theme_bw() + 
    geom_hline(aes(yintercept = 1)) 
  ggsave(file = paste0(outdir, '-DeathRatio_elderly_', Code, '.png'), w = 8, h = 8)
}

plot_death_ratio_winter = function(death_ratio_winter, vaccinedata, outdir){
  
  tmp = vaccinedata[prop_vaccinated_1dosep >= 0.5, list(date = min(date), prop = 0.5), by = 'age']
  tmp1 = vaccinedata[prop_vaccinated_1dosep >= 0.75, list(date = min(date), prop = 0.75), by = 'age']
  tmp1 = rbind(tmp, tmp1)
  tmp1 = subset(tmp1, age == '75+')
  tmp1[, prop_name := paste0(prop*100, '%')]
  tmp1[, date := as.Date(date)]
  
  tmp = subset(death_ratio_winter, date >= as.Date('2020-12-01'))
  ggplot(tmp, aes(x = date)) + 
    geom_line(aes(y = M, col = age)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) +
    geom_point(aes(y = emp_est, col = age)) +
    theme_bw() + 
    geom_vline(data = tmp1, aes(xintercept = date, linetype = prop_name)) + 
    scale_y_continuous(labels = scales::percent_format()) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 45,hjust=1,vjust=1), 
          panel.grid.major= element_blank(), legend.box="vertical", 
          legend.margin=margin())+ 
    labs(y = 'Deaths as percentage of winter peak', x = '', col = 'Age group',  fill = 'Age group', 
         linetype = 'Proportion of 75+ vaccinated\nwith at least one dose') + 
    scale_linetype_manual(values = c(2,3)) + 
    scale_color_viridis_d(option = 'B', begin = 0.2, end = 0.8)+ 
    scale_fill_viridis_d(option = 'B', begin = 0.2, end = 0.8) + 
    guides(linetype = guide_legend(order=2), color = guide_legend(order=1), fill = guide_legend(order=1)) 
  ggsave(paste0(outdir, '-DeathsRatioWinter_', Code, '.png'), w = 6, h = 5)
  
  
}

savepdf <- function(fname, width=16, height=10)
{
  pdf(fname, width=width/2.54, height=height/2.54)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(0,0,0,0))
}

plot_posterior_plane = function(fit_cum, df_week, outdir)
{
  fit_samples = extract(fit_cum)
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  tmp1 = as.data.table( reshape2::melt( fit_samples$beta ))
  
  row_name = 'week_index'
  column_name = 'basis_function_index'
  if(max(tmp1$Var3) == stan_data$A) column_name = 'Age'
  
  setnames(tmp1, c('Var2', 'Var3'), c(row_name, column_name))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps),
                       q_label=p_labs), 
              by=c(column_name, row_name)]	
  tmp1 = dcast(tmp1, get(row_name) + get(column_name) ~ q_label, value.var = "q")
  setnames(tmp1, c('row_name', 'column_name'), c(row_name, column_name))
  
  tmp1 = merge(tmp1, df_week, by = row_name)

  ggplot(tmp1, aes(x = date, y = get(column_name))) +
    geom_raster(aes(fill = M))  + 
    labs(x = 'Date', y = column_name, fill = 'Estimated posterior value') + 
    scale_y_reverse(expand = c(0,0), breaks = seq(min(tmp1[,get(column_name)]), max(tmp1[,get(column_name)]), 2))  +
    scale_x_date(expand = c(0,0), breaks = '2 months', date_labels = "%b-%y") + 
    scale_fill_viridis_c(option = "A") + 
    theme(legend.position='bottom')
  ggsave(file = paste0(outdir, '-PlanePosterior_', Code, '.png'), w = 5, h = 5.5)
  
}



plot_posterior_predictive_checks_cumulative = function(tmp, outdir)
{
  
  p = ggplot(tmp, aes(x = date)) + 
    geom_point(aes(y = M)) +
    geom_errorbar(aes(ymin = CL, ymax = CU)) +
    geom_ribbon(aes(ymin = min_count_censored, ymax = max_count_censored), alpha = 0.5, fill = 'red') +
    geom_point(aes(y = sum_count_censored), col = 'blue') +
    facet_wrap(~age, scale = 'free_y') + 
    labs(x = '', y = 'cumulative deaths')
  ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_cum_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
  
}


compare_CDCestimation_JHU_error_plot = function(CDC_data, JHU_data, var.cum.deaths.CDC, outdir)
{
  # prepare JHU data
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  
  # prepare CDC estimations
  CDC_data = select(as.data.table(CDC_data), code, date, var.cum.deaths.CDC)
  setnames(CDC_data, var.cum.deaths.CDC, 'cumulative_deaths')
  
  # plot
  JHUData[, source := 'JHU']
  CDC_data[, source := 'CDC']
  
  tmp2 = rbind(JHUData, CDC_data)
  tmp2 = subset(tmp2, code %in% unique(CDC_data$code) & date <= max(CDC_data$date))
  
  n_code = length(unique(tmp2$code))
  
  p = ggplot(tmp2, aes(x = date, y = cumulative_deaths)) + 
    geom_line(aes(col = source), size = 1) +
    facet_wrap(~code, nrow = length(unique(tmp2$code)), scale = 'free') + 
    theme_bw() + 
    scale_color_viridis_d(option = "B", direction = -1, end = 0.8) + 
    scale_fill_viridis_d(option = "B", direction = -1, end = 0.8)
  ggsave(p, file = paste0(outdir, '-comparison_JHU_CDC_uncertainty_', Code, '.png'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
}

plot_posterior_plane_with_data = function(deathByAge_1, deathByAge_2, outdir)
{
  library(magick)
  
  deathByAge_1 = subset(deathByAge_1, code == Code)
  deathByAge_2 = subset(deathByAge_2, code == Code)
  
  range_d = range(na.omit(c(deathByAge_1$weekly.deaths, deathByAge_2$weekly.deaths)))
  
  p1 = ggplot(deathByAge_1, aes(x = date, y = age)) + 
    geom_raster(aes(fill = weekly.deaths )) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'sqrt', limits = range_d) +
    scale_x_date(expand = c(0,0), date_labels = c("%b %Y"), date_breaks = '1 month') + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 70, vjust  = 0.5), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age groups first specification',
         fill = 'Reported covid-19 deaths')
  
  p2 = ggplot(deathByAge_2, aes(x = date, y = age)) + 
    geom_raster(aes(fill = weekly.deaths )) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'sqrt', limits = range_d) +
    scale_x_date(expand = c(0,0), date_labels = c("%b %Y"), date_breaks = '1 month') + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 70, vjust  = 0.5), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age groups second specification',
         fill = 'Reported covid-19 deaths')
  
  p = ggpubr::ggarrange(p1, p2, nrow = 2,common.legend = T, legend = 'bottom')
  ggsave(p, file = paste0(outdir, '-panel_right_', Code, '.png'), h = 8, w = 6)
  
  panel.left <- magick::image_read(paste0(outdir, '-PlanePosterior_', Code, '.png'))
  panel.right <- magick::image_read(paste0(outdir, '-panel_right_', Code, '.png'))
  
  p = image_append(c( image_scale(panel.left, "2200"), panel.right) )
  
  savepdf(paste0(outdir, '-PlanePosterior_Data_', Code, '.pdf'), w = 14.5*2, h = 10*2)
  plot(p)
  dev.off()
  
}

compare_CDCestimation_DoH_age_plot_compmethod = function(tab_doh, scraped_data, model_name, selected_method)
{
  scraped_data = select(as.data.table(scraped_data), code, date, age, cum.deaths)
  scraped_data = subset(scraped_data, code == unique(tab_doh$code))
  scraped_data[, date := as.Date(date)]
  scraped_data[, CL := NA]
  scraped_data[, CU := NA]
  ages = unique(scraped_data$age)
  scraped_data = scraped_data[order(code, age, date)]
  scraped_data = scraped_data[rep(seq_len(nrow(scraped_data)), length(unique(tab_doh$method))), ]
  scraped_data[, method := rep(unique(tab_doh$method), each = nrow(scraped_data)/length(unique(tab_doh$method)))]
  
  # prepare CDC estimations
  CDC_data = select(as.data.table(tab_doh), code, date, age, M, CL, CU, method)
  setnames(CDC_data, 'M', 'cum.deaths')
  
  # plot
  scraped_data[, source := 'DoH']
  CDC_data[, source := 'Estimated']
  
  tmp2 = rbind(scraped_data, CDC_data)
  tmp2[, age := factor(age, levels = ages)]
  tmp2[, source := factor(source, levels = c('Estimated', 'DoH'))]
  tmp2[, method := factor(method, levels = model_name)]
  
  col = c('#5CC8D7FF', '#00203FFF')
  
  df = CDC_data[method == selected_method, list(maxCL = max(CU)), by = 'age']
  
  
  p = ggplot(tmp2, aes(x = date, y = cum.deaths)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = source), alpha = 0.5) +
    geom_line(aes(col = source), size = 1) +
    theme_bw() + 
    scale_color_manual(values = col) + 
    scale_fill_manual(values = col) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    facet_grid(age~method,  scale = 'free_y') + 
    theme(legend.position = 'bottom',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x = element_text(angle = 45,hjust=1,vjust=1), 
          axis.title.x =  element_blank(), 
          axis.title.y = element_text(size = rel(1.1)),
          strip.text =  element_text(size = rel(1)),
          legend.text = element_text(size = rel(1.1))) + 
    labs(y = 'Cumulative COVID-19 deaths', col = '', fill = '') + 
    coord_cartesian_panels(
      panel_limits = tibble::tribble(
        ~age, ~ymin, ~ymax
        , "0-4" ,    -1,    df$maxCL[1] 
        , "5-14" ,    -1,    df$maxCL[2] 
        , "15-24" ,    -1,    df$maxCL[3] 
        , "25-34" ,    -1,    df$maxCL[4] 
        , "35-44" ,    -1,    df$maxCL[5] 
        , "45-54" ,    -1,    df$maxCL[6] 
        , "55-64" ,    -1,    df$maxCL[7] 
        , "65-74" ,    -1,    df$maxCL[8] 
        , "75-84" ,    -1,    df$maxCL[9] 
        , "85+" ,    -1,    df$maxCL[10] 
      ))
  
  return(p)
}



plot_contribution_comparison_method = function(tab_cd, data, model_name){
  
  tmp = copy(data)
  setnames(tmp, 'prop_deaths', 'M')
  set(tmp, NULL, 'CL', NA_real_)
  set(tmp, NULL, 'CU', NA_real_)
  
  tmp = rbind(select(tab_cd, date, age, M, method), select(tmp, date, age, M, method))
  tmp[, method := factor(method, c('observation', model_name))]
  
  p = ggplot(tmp, aes(x = date, y = age)) + 
    geom_raster(aes(fill = M )) + 
    theme_bw() +
    scale_fill_viridis_c(option = "E") +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    labs(x = '', y = 'Age group',
         fill = 'Contribution to\nweekly\nCOVID-19\ndeaths') + 
    facet_grid(~method) + 
    theme(legend.position = 'left',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white")) 
  
  return(p)
}

plot_death_comparison_method = function(tab_d, data, model_name){
  
  tmp = copy(data)
  setnames(tmp, 'weekly.deaths', 'M')
  set(tmp, NULL, 'CL', NA_real_)
  set(tmp, NULL, 'CU', NA_real_)
  
  tmp = rbind(select(tab_d, date, age, M, method), select(tmp, date, age, M, method))
  tmp[, method := factor(method, c('observation', model_name))]
  
  range_wd = sqrt(range(na.omit(data$weekly.deaths)))
  digits_cut= ifelse(range_wd[2]/100 > 1, 2, 1)
  digits_cut= ifelse(range_wd[2]/10 > 1, digits_cut, 0)
  p =  ggplot(tmp, aes(x = date, y = age)) + 
    geom_raster(aes(fill = M )) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'sqrt', breaks = round(seq(range_wd[2]/2, range_wd[2], length.out = 3)^2, digits = -digits_cut) ) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    labs(x = '', y = 'Age group',
         fill = 'Reported\nCOVID-19\ndeaths') + 
    facet_grid(~method) + 
    theme(legend.position = 'left',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white")) 
  
  return(p)
}


plot_contribution_magnitude_all_states = function(contribution, death2agegroups, selected_states, outdir, nrow = 1, legend.position = 'bottom')
{
  
  tmp = subset(contribution, code %in% selected_states)
  tmp2 = subset(death2agegroups, code %in% selected_states)
  
  tmp1 = unique(select(tmp, code, loc_label))
  tmp1[, code := factor(code, levels = selected_states)]
  tmp1 = tmp1[order(code)]
  
  tmp[, loc_label := factor(loc_label, tmp1$loc_label)]
  tmp2[, loc_label := factor(loc_label, tmp1$loc_label)]
  
  p1 = ggplot(tmp, aes(x= date) ) +
    geom_line(aes(y = M, col = age)) + 
    geom_point(aes(y = emp, col = age), size = 0.5) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
    facet_wrap(~loc_label, nrow = 3) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme_bw() + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 45, hjust =1), 
          panel.grid.major = element_blank(), 
          axis.title.x = element_blank()) +
    labs(y = 'Estimated age-specific contribution\nto COVID-19 weekly deaths') + 
    scale_y_continuous(labels = scales::percent_format()) +
    scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
    scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8)
  
  
  p2 = ggplot(tmp2, aes(x= date) ) +
    geom_line(aes(y = M, col = age)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
    facet_wrap(~loc_label, nrow = 3, scale = 'free_y') +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme_bw() + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 45, hjust =1), 
          panel.grid.major = element_blank(), 
          axis.title.x = element_blank()) +
    labs( y = 'Estimated age-specific COVID-19\nweekly deaths', col = 'Age group', fill = 'Age group') + 
    scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
    scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8)
  
  
  p = ggpubr::ggarrange(p2,p1, common.legend = T, legend = legend.position, nrow = nrow)
  ggsave(p, file = paste0(outdir, '-ContributionMagnitude_panel.png'), w = 7, h = 8)
  
  
  return(p)
}
