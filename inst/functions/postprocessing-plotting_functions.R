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

plot_posterior_predictive_checks = function(data, variable, outdir)
{
  
  data[, PPP := paste0('inside CI: ', round(mean(na.omit(inside.CI))*100, digits = 2), '%'), by = 'date']
  data[, date_ppp := paste0(as.character(date), ' - ', PPP)]
  
  Code = unique(data$code)
  
  # posterior data check
  p = ggplot(data, aes(x = age)) + 
    geom_point(aes(y = M), col = "black", size = 1) +
    geom_errorbar(aes(ymin = CL, ymax = CU), width=0.3, col = "black") +
    geom_point(aes(y = get(variable)), col = "red", size = 1) + 
    theme_bw() +
    labs(y = "Weekly COVID-19 deaths", x = "Age") + 
    facet_wrap(~date_ppp, scale = 'free_y') + 
    theme(axis.text.x = element_text(angle = 90),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) 
  ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_", Code,".png") , w= 15, h = 10, limitsize = FALSE)
  
  # posterior predictive check
  p = ggplot(data, aes(x = date)) + 
    geom_point(aes(y = M), col = "black", size = 1) +
    geom_errorbar(aes(ymin = CL, ymax = CU), width=0.3, col = "black") +
    geom_point(aes(y = get(variable)), col = "red", size = 1) + 
    theme_bw() +
    labs(y = "Weekly COVID-19 deaths", x = "Week") + 
    facet_wrap(~age, ncol = 3, scale = 'free_y') + 
    theme(axis.text.x = element_text(angle = 90),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_2_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
  
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

plot_probability_deaths_age_contribution = function(tmp1, var_name, outdir, discrete = F, limits = NULL, wo_ytext = F, wo_legend = F){
  
  if(!discrete)
    tmp1[, age := as.numeric(age)]
  
  n_row = length(unique(tmp1$date))
  
  p = ggplot(tmp1, aes(x = age)) + 
    theme_bw() +
    labs(y = 'Age-specific contribution to COVID-19 weekly deaths', x = "Age") + 
    facet_wrap(~date) + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  dates = unique(tmp1$date)
  tmp2 = subset(tmp1, date %in% dates[seq(1, length(dates), length.out =3)])
  tmp2[, date_format := format(date, "%d-%b-%y")]
  p1 = ggplot(tmp2, aes(x = age)) + 
    theme_bw() +
    labs(y = 'Age-specific contribution to COVID-19 weekly deaths', x = "Age") + 
    facet_wrap(~date_format, ncol = 1)+ 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  if(discrete){
    p = p + geom_point(aes(y = M)) + geom_errorbar(aes(ymin= CL, ymax = CU)) 
    p1 = p1 + geom_point(aes(y = M)) + geom_errorbar(aes(ymin= CL, ymax = CU)) 
  } else {
    p = p + geom_line(aes(y = M)) + geom_ribbon(aes(ymin= CL, ymax = CU), alpha = 0.5) 
    p1 = p1 + geom_line(aes(y = M)) + geom_ribbon(aes(ymin= CL, ymax = CU), alpha = 0.5) 
  }
  if(is.null(limits)){
    p = p + coord_cartesian(ylim = limits)
  }
    
  if(!is.null(outdir)){
    ggsave(p, file = paste0(outdir, "-continuous_contribution_", var_name, '_', Code, ".png") , w= 10, h = 8, limitsize = FALSE)
    ggsave(p1, file = paste0(outdir, "-continuous_contribution_short_", var_name, '_', Code, ".png") , w= 4, h = 8, limitsize = FALSE)
    
  }

  p = ggplot(tmp1, aes(x = date, y = age)) +
    geom_raster(aes(fill = M))  + 
    theme(legend.position='bottom') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  
  if(is.null(limits)){
    p = p + scale_fill_viridis_c(option = "E") 
  } else {
    p = p + scale_fill_viridis_c(option = "E", limits= limits) 
  }
  
  if(wo_ytext)
    p = p + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  
  if(wo_legend)
    p = p + theme(legend.position = 'none')
  
  if(discrete){
    p = p + scale_y_discrete(expand = c(0,0))  + labs(x = '', y = 'Age group', fill = 'Estimated posterior value') 
  } else {
    p = p + scale_y_continuous(expand = c(0,0))  +labs(x = '', y = 'Age', fill = 'Estimated posterior value') 
  }
  if(!is.null(outdir)){
    ggsave(p, file = paste0(outdir, '-continuous_contribution_allweeks_', var_name, '_', Code, '.png'), w = 5, h = 5.2)
  }
  return(list(p1, p))
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

plot_sum_missing_deaths = function(tmp, outdir)
{
  
  p = ggplot(tmp, aes(x = date)) + 
    geom_point(aes(y = M)) +
    geom_errorbar(aes(ymin = CL, ymax = CU)) +
    geom_ribbon(aes(ymin = min_count_censored, ymax = max_count_censored), alpha = 0.5, fill = 'red') +
    geom_point(aes(y = sum_count_censored), col = 'blue') +
    facet_wrap(~age, scale = 'free_y') + 
    labs(x = '', y = 'cumulative deaths')
  ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_missing_cum_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
  
}

plot_sum_bounded_missing_deaths = function(tmp, outdir)
{
  
  tmp = copy(tmp1)
  tmp[, dummy := 1]
  
  p = ggplot(tmp, aes(x = dummy) ) + 
    geom_point(aes(y = M)) +
    geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.25) +
    geom_errorbar(aes(ymin = min_count_censored, ymax = max_count_censored),  color = 'red', position = 'dodge', width = 0.25) +
    geom_point(aes(y = sum_count_censored), col = 'blue') +
    facet_wrap(~age, scale = 'free_y') + 
    labs(x = '', y = 'cumulative deaths')
  ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_missing_sum_cum_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
  
}

plot_mortality_rate = function(mortality_rate_table, outdir)
{
  
  
  p = ggplot(mortality_rate_table, aes(x = date)) + 
    geom_line(aes(y = M)) +
    geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) +
    facet_wrap(~age, scale = 'free_y') + 
    labs(x = '', y = 'COVID-19 mortality rate', col = '', fill = '') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90))  +
    theme_bw()
  ggsave(p, file = paste0(outdir, '-MortalityRate_', Code, '.png'), w = 7, h =7, limitsize = F)
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

compare_CDCestimation_DoH_age_plot = function(CDC_data, scraped_data, var.cum.deaths.CDC, outdir, overall = F)
{
  scraped_data = select(as.data.table(scraped_data), code, date, age, cum.deaths)
  scraped_data = subset(scraped_data, code == Code)
  scraped_data[, date := as.Date(date)]
  scraped_data[, CL := NA]
  scraped_data[, CU := NA]
  scraped_data = scraped_data[order(code, age, date)]

  # age sorted
  ages = unique(select(CDC_data, age))
  ages[, age_from := gsub('(.+)-.*', '\\1', age)]
  ages[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  ages = ages[order(as.numeric(age_from))]$age
  
  # prepare CDC estimations
  CDC_data = select(as.data.table(CDC_data), code, date, age, var.cum.deaths.CDC )
  setnames(CDC_data, var.cum.deaths.CDC, c('cum.deaths', 'CL', 'CU'))
  CDC_data = subset(CDC_data, date <= max(scraped_data$date))
  
  # plot
  scraped_data[, source := 'DoH']
  CDC_data[, source := 'Estimated']
  
  tmp2 = rbind(scraped_data, CDC_data)
  tmp2[, age := factor(age, levels = ages)]
  tmp2[, source := factor(source, levels = c('Estimated', 'DoH'))]
  
  col = c('#5CC8D7FF', '#00203FFF')
  
  p = ggplot(tmp2, aes(x = date, y = cum.deaths)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = source), alpha = 0.5) +
    geom_line(aes(col = source), size = 1) +
    theme_bw() + 
    scale_color_manual(values = col) + 
    scale_fill_manual(values = col) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    facet_wrap(~age,  scale = 'free_y', ncol = 1) + 
    theme(legend.position = 'bottom',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x = element_text(angle = 45,hjust=1,vjust=1), 
          axis.title.x = element_blank()) + 
    labs(y = 'Cumulative COVID-19 deaths', col = '', fill = '')
  ggsave(p, file = paste0(outdir, '-comparison_DoH_CDC_uncertainty_', Code, '.png'), w = 3, h = 10, limitsize = F)
  
  return(p)
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

plot_posterior_plane = function(fit_cum, df_week, df_age_continuous,stan_data, outdir)
{
  
  # if(is.null(stan_data$num_basis)) return(NULL)
  
  varname = 'beta'
  row_name = 'week_index'
  row_lab = "Date" 
  column_name = 'age_index' 
  column_lab = 'Age'
  if(!is.null(stan_data$num_basis_rows)){
    varname = 'f'
    row_name = 'age_index' 
    column_name = 'week_index' 
  }
  
  # column_name = 'week_index'
  # column_lab = 'Date'
  # row_name = 'basis_function_index_age'
  # row_lab = "Basis function index"
  # n_rows = stan_data$num_basis
  # n_columns = stan_data$W
  # if(is.null(stan_data$num_basis)){
  #   
  # }
  # if(!is.null(stan_data$num_basis_rows)){
  #   column_name = 'basis_function_index_week'
  #   column_lab = 'Week basis function'
  #   row_name = 'basis_function_index_age'
  #   row_lab = "Age basis functions"
  #   n_rows = stan_data$num_basis_rows
  #   n_columns = stan_data$num_basis_columns
  # }
  
  euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))   
  fit_samples = rstan::extract(fit_cum)
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt( fit_samples[[varname]] ))
  setnames(tmp1, c('Var2', 'Var3'), c(row_name, column_name))
  tmp1 = tmp1[, list(q = quantile(value, prob=ps), q_label=p_labs), by=c(column_name, row_name)]
  tmp1 = dcast(tmp1, get(row_name) + get(column_name) ~ q_label, value.var = "q")
  setnames(tmp1, c('row_name', 'column_name'), c(row_name, column_name))
  
  if(!is.null(stan_data$num_basis_rows)){
    row_name_1 = row_name 
    row_name = column_name 
    column_name = row_name_1
  }
  
  tmp1 = merge(tmp1, df_week, by = row_name)

    p = ggplot(tmp1, aes(x = date, y = get(column_name))) + 
      geom_raster(aes(fill = M ), interpolate = TRUE) + 
      theme_bw() +
      theme(legend.position = 'none',
            axis.text.x = element_text(angle = 90), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      labs(y = column_lab, x = row_lab, fill ='') +
      scale_fill_viridis_c(option = 'inferno', begin = 0.55) +
      scale_y_reverse(expand = c(0,0)) + 
      geom_contour(aes(z=M), bins = 12, color = "gray50", 
                   size = 0.5, alpha = 0.5)  + 
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) 
    ggsave(p, file = paste0(outdir, '-PlanePosterior_', Code, '.png'),width = 4, height = 4)
    
    
  return(p)
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



plot_mean_age_death = function(mean_age_death, outdir){
  if(length(unique(mean_age_death$code)) == 1)
  {
    ggplot(mean_age_death, aes(x = date)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) +
      geom_line(aes(y = M)) + 
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90)) +
      labs(x = '', y = 'Mean age death')
    ggsave(file = paste0(outdir, '-MeanAgeDeath_', Code, '.png'), w = 4, h = 4)
  }
  
  else {
    
    mean_age_death[, loc_label := factor(loc_label, levels= sort(unique(mean_age_death$loc_label), decreasing = T))]
    ggplot(mean_age_death, aes(x = date, y = loc_label)) + 
      geom_raster(aes(fill = M)) +
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
      scale_y_discrete(expand = c(0,0)) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90), 
            legend.position = 'bottom') +
      labs(x = '', y = '', fill = 'Mean age death') + 
      scale_fill_viridis_c(option = 'B')
    ggsave(file = paste0(outdir, '-MeanAgeDeath.png'), w = 6, h = length(unique(mean_age_death$loc_label)) / 4, limitsize = F)
  }
}

plot_mortality_rate_all_states = function(mortality_rate, outdir)
{
  for(a in unique(mortality_rate$age_index))
  {
    tmp = subset(mortality_rate, age_index == a)
    age = unique(tmp$age)
    tmp = tmp[order(M)]
    medianM =  tmp[,median(M)]
    tmp[, loc_label := factor(loc_label, tmp$loc_label)]
    ggplot(tmp, aes(x=loc_label, y = M)) + 
      geom_bar(aes(fill = M), stat="identity") +
      geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
      scale_fill_gradient2(low= 'darkturquoise', high = 'darkred', mid = 'beige', midpoint = medianM) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1), 
            legend.position = 'none', 
            panel.grid.major= element_blank()) + 
      scale_y_continuous(expand = c(0,0), labels = scales::percent_format()) +
      coord_cartesian(ylim = c(0,NA)) +
      labs(x ='', y = paste0('Predicted COVID-19 attributable mortality rates\namong individuals ', age, ' as of ', format(unique(tmp$date), '%b %Y')))
    ggsave(paste0(outdir, paste0('-MortalityRate_', age, '.png')), w = 8, h = 5)
  }
}

plot_contribution_all_states = function(contribution, vaccinedata_state, outdir){

  vaccinedata_state = subset(vaccinedata_state, loc_label %in% unique(contribution$loc_label) & date <= max(contribution$date))
  vaccinedata_state[, dummy := '']
  
  tmp = contribution[!is.na(emp)]
  
  ggplot(tmp, aes(x= date) ) +
    geom_line(aes(y = M, col = age)) + 
    geom_point(aes(y = emp, col = age), size = 0.5) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
    facet_wrap(~loc_label, ncol = 6) +
    geom_line(data = vaccinedata_state, aes(y  = prop_vaccinated_1dosep, linetype = dummy), col = 'grey20', size = 0.9) + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme_bw() + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 45, hjust =1), 
          panel.grid.major = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = rel(1.2)), 
          legend.title = element_text(size = rel(1.1)), 
          strip.text = element_text(size = rel(1.1)), 
          legend.position = 'bottom') +
    labs(y = 'Estimated age-specific contribution to COVID-19 weekly deaths', fill = 'Age groups', col = 'Age groups', 
         linetype = 'Proportion of the state population\nvaccinated with at least one dose') + 
    scale_y_continuous(labels = scales::percent_format()) +
    scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
    scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8) + 
    scale_linetype_manual(values = 1)
  ggsave(paste0(outdir, paste0('-Contribution_Vaccination_allStates.png')), w = 9, h = 12)
  
  
}

plot_mortality_all_states = function(death, data, min_vaccine_date, outdir)
{
  codes = unique(death$code)
  ncodes = length(codes)
  mid = round(ncodes/ 2)
  
  df = as.data.table( reshape2::melt(select(death, loc_label, code, date, age, emp), id.vars = c('loc_label', 'code', 'date', 'age')) )
  df[, variable2 := 'CDC data']

  vaccination_start = data.table(date = min_vaccine_date)
  death[, dummy := 'Posterior median prediction\nusing age-aggregated JHU data\nto adjust for reporting delays']
  
  tmp = subset(death, code %in% codes[1:mid])
  df1 =subset(df, code %in% codes[1:mid])
  ggplot(tmp, aes(x= date) ) +
    geom_line(aes(y = M, col = age, linetype = dummy)) + 
    geom_point(data = df1, aes(y = value, shape= variable2, fill = age), size = 1, col = 'grey70', stroke = 0.1) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
    facet_wrap(~loc_label, scale = 'free_y', ncol = 4) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    geom_vline(data = vaccination_start, aes(xintercept = date), linetype = 2, col = 'grey50') +
    theme_bw() + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 45, hjust =1), 
          panel.grid.major = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = rel(1.2)), 
          legend.title = element_text(size = rel(1.1)), 
          legend.text = element_text(size = rel(1.1)), 
          strip.text = element_text(size = rel(1.1)), 
          legend.position = 'bottom') +
    labs( y = 'Age-specific COVID-19 weekly deaths', col = 'Age group', fill = 'Age group', shape = '', linetype = '') + 
    scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
    scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8) +
    scale_shape_manual(values = 21) + 
    guides(shape = guide_legend(override.aes = list(size=1, stroke = 1), order = 1), linetype =  guide_legend(order=2),
           fill = guide_legend(order=3), col = guide_legend(order =3))
  ggsave(paste0(outdir, paste0('-Mortality_allStates_1.png')), w = 9, h = 12)
  
  tmp = subset(death, code %in% codes[(mid+1):ncodes])
  df1 =subset(df, code %in% codes[(mid+1):ncodes])
  ggplot(tmp, aes(x= date) ) +
    geom_line(aes(y = M, col = age, linetype = dummy)) + 
    geom_point(data = df1, aes(y = value, shape= variable2, fill = age), size = 1, col = 'grey70', stroke = 0.1) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
    facet_wrap(~loc_label, scale = 'free_y', ncol = 4) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    geom_vline(data = vaccination_start, aes(xintercept = date), linetype = 2, col = 'grey50') +
    theme_bw() + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 45, hjust =1), 
          panel.grid.major = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = rel(1.2)), 
          legend.title = element_text(size = rel(1.1)), 
          legend.text = element_text(size = rel(1.1)), 
          strip.text = element_text(size = rel(1.1)), 
          legend.position = 'bottom') +
    labs( y = 'Age-specific COVID-19 weekly deaths', col = 'Age group', fill = 'Age group', shape = '', linetype = '') + 
    scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
    scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8)  +
    scale_shape_manual(values = 21) + 
    guides(shape = guide_legend(override.aes = list(size=1, stroke = 1), order = 1), linetype =  guide_legend(order=2),
           fill = guide_legend(order=3), col = guide_legend(order =3))
  ggsave(paste0(outdir, paste0('-Mortality_allStates_2.png')), w = 9, h = 12)
  
}

plot_contribution_continuous_comparison_method = function(tab_cc, tab_d, data, 
                                                          selected_method, model_name, show.method = T, heights= c(0.4,0.6)){
  
  dates = unique(tab_cc$date)
  dates = dates[seq(3, length(dates)-3, length.out =3)]
  
  plot_labels <- format(dates,  "%d-%b-%y")
  
  tmp2 = subset(tab_cc, date %in% dates)
  tmp2[, age_c := as.numeric(age)]
  tmp2[, method := factor(method,  model_name)]
  
  tmp3 = as.data.table(tab_d)
  tmp3[, method := factor(method,  model_name)]
  tmp3[, age_c := as.numeric(age)]
  
  limit_SE = range(subset(tmp2, method == selected_method)$CL, subset(tmp2, method == selected_method)$CU)
  
  tmp= tab_d[, list(sumM = sum(mean)), by = c('date', 'method')]
  
  df = data.frame(value = dates, y = max(tmp$sumM) - max(tmp$sumM)*0.05, 
                  key = as.character(1:length(dates)), 
                  x_prop = 7,
                  y_prop = limit_SE[2] -limit_SE[2]*0.06, date=dates)
  

  p1 = ggplot(tmp2, aes(x = age_c)) + 
    geom_line(aes(y = M)) +
    geom_ribbon(aes(ymin= CL, ymax = CU), alpha = 0.5) +
    theme_bw() +
    labs(y = 'Estimated age-specific contribution\nto COVID-19 weekly deaths', x = "Age") + 
    facet_grid(method~date) +
    coord_cartesian(ylim = limit_SE, xlim = range(tmp2$age_c)) +
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) + 
    theme(panel.border = element_rect(colour = "black", fill = NA), 
          # legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"), 
          panel.grid.major = element_blank(),
          axis.title.y = element_text(size = rel(1.1)), 
          plot.margin = unit(c(5.5,5.5,25,5.5), "pt"), 
          legend.position = 'none')  + 
    geom_point(data = df, aes(x = x_prop, y = y_prop, col = key, shape = key), size =3, stroke = 1.5 ) + 
    scale_shape_manual(name = "Dates", labels = plot_labels, values = c(21, 22, 23)) + 
    scale_colour_manual(name = "Dates", labels = plot_labels, values = gg_color_hue(3)) 
  
  if(show.method==F){
    p1 = p1 + theme(strip.text = element_blank())
  } else{
    p1 = p1 + theme(strip.text.x =  element_blank(),
                    strip.text.y =  element_text(size = rel(1.2)))
  }
  p1 = ggarrange(p1, labels = 'C', font.label = list(size = 20, face = 'bold'), label.x = 0.04, 
                 vjust = 1)


  tmp3 = subset(tmp3, method == selected_method)
  tmp3[, age := factor(age, levels = rev(levels(tmp3$age)))]
  tmp3[, dummy := 'JHU\noverall\nweekly\ndeaths']
  p2 = ggplot(tmp3, aes(x = date)) + 
    theme_bw() +
    geom_step(aes(x = date+ 3.5,y = emp_JHU, linetype = dummy), direction =  "vh") + 
    labs(y = 'Posterior median of the age-specific\nCOVID-19 weekly deaths', x = "", color = '', shape = '', fill = 'Age') + 
    facet_grid(method~.) +
    # coord_cartesian(ylim = limit_SE) + 
    scale_fill_viridis_d(option = 'B',  direction = -1,
    guide = guide_colorsteps(show.limits = TRUE)) +
    # scale_fill_binned(breaks = c(10, 25, 50)) +
    scale_x_date(expand = c(0,0), breaks = '2 months', date_labels = "%b-%y") + 
    geom_segment(data = df, aes(x = value, y = 0, xend = value, yend = y), 
                 linetype = "dashed", colour = "grey30", alpha = 0.75) +
    geom_point(data = df, aes(x = value, y = y, group = key, shape = key, col = key), size = 2, stroke = 1.5 ) +
    scale_shape_manual(name = "", labels = plot_labels, values = c(21, 22, 23)) + 
    scale_colour_manual(name = "", labels = plot_labels, values = gg_color_hue(3)) + 
    geom_bar(aes(y = mean, fill = age), stat = 'identity', width = 7)  +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          strip.background = element_rect(colour="white", fill="white"), 
          panel.grid.major = element_blank(),
          strip.text = element_blank(),
          axis.title.y = element_text(size = rel(1.1)),
          legend.text = element_text(size = rel(1)),
          legend.title = element_text(size = rel(1)),
          axis.text.x = element_text(angle= 40, hjust =1)
    ) +
    guides(col = F, shape = F)
  pl = ggplot(tmp3, aes(x = date, y = M, col = age_c)) + geom_step(aes(x = date+ 3.5,y = emp_JHU, linetype = dummy), direction =  "vh") + 
    geom_point() + scale_color_viridis_c(option = 'B') + labs(color = 'Age', linetype = '')  + theme_bw() + 
    theme(legend.key.height = unit(0.8, "cm"),legend.spacing.y = unit(0.0, 'cm')) + guides(linetype = guide_legend(order=1), color = guide_colorbar(order=2))
  p2 = ggpubr::ggarrange(p2,common.legend = T, legend.grob = get_legend(pl), legend = 'right',
                         labels = 'B', font.label = list(size = 20, face = 'bold'), label.x = 0.04)

  data[, age := factor(age, levels = rev(levels(data$age)))]
  p3 = ggplot(data, aes(x = date)) + 
    theme_bw() +
    labs(y = 'Data age-specific COVID-19\nweekly deaths', x = "", color = '', shape = '', fill = 'Original\nAge\nGroup', linetype = '') + 
    # coord_cartesian(ylim = limit_SE) + 
    scale_fill_viridis_d(option = 'A', end = 0.95, direction = -1) +
    scale_x_date(expand = c(0,0), breaks = '2 months', date_labels = "%b-%y") + 
    geom_bar(aes(y = weekly.deaths, fill = age), stat = 'identity', width = 7)  +
    geom_step(data = tmp3, aes(x = date+ 3.5,y = emp_JHU), direction =  "vh") + 
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          strip.background = element_rect(colour="white", fill="white"), 
          panel.grid.major = element_blank(),
          strip.text = element_blank(),
          axis.title.y = element_text(size = rel(1.1)),
          legend.text = element_text(size = rel(0.8)),
          legend.title = element_text(size = rel(1)),
          # legend.key.height = unit(0.6, "cm"),
          axis.text.x = element_text(angle= 40, hjust =1)) + 
    guides(linetype = guide_legend(order=1), color = guide_legend(order=2))
  p3 = ggarrange(p3, labels = 'A', font.label = list(size = 20, face = 'bold'))
  # p = ggpubr::ggarrange(p2,p1, common.legend = T, legend.grob = get_legend(p3), legend = 'left', widths = c(0.4, 0.6))
              
  p = grid.arrange(p1, p2, p3, layout_matrix = rbind(c(3,2), c(1,1)), heights = heights)
  # ggsave(p, file = '~/Downloads/figure4.png', w = 9, h = 10)
  
  return(p)
  
}



plot_contribution_ref_all_states = function(contribution_ref, contribution_ref_adj, outdir){
  
  tmp = subset(contribution_ref, age == '85+')
  tmp = tmp[order(M)]
  contribution_ref[, loc_label := factor(loc_label, unique(tmp$loc_label))]
  
  ggplot(contribution_ref, aes(x = loc_label, y = M)) + 
    geom_bar(aes(fill = M), stat = 'identity') +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
    facet_grid(age~division,  scales = "free_x", space = 'free_x') + 
    theme_bw() +
    theme(axis.text.x = element_text(angle= 45, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.position = 'none')  +
    scale_fill_viridis(option = 'E', trans = 'sqrt') + 
    scale_y_continuous(labels = scales::percent) +
    labs(x ='', y = paste0('Age-specific contribution to COVID-19 deaths\nduring the baseline period'))
  ggsave(paste0(outdir, '-Contribution_ref.png'), w = 9, h = 6)
  
  ggplot(contribution_ref, aes(x = loc_label, y = M)) + 
    geom_bar(aes(fill = M), stat = 'identity') +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
    geom_point(aes(y = emp_est), col = 'tomato1', size = 0.6) + 
    facet_grid(age~division,  scales = "free_x", space = 'free_x') + 
    theme_bw() +
    theme(axis.text.x = element_text(angle= 45, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.position = 'none')  +
    scale_fill_viridis(option = 'E', trans = 'sqrt') + 
    scale_y_continuous(labels = scales::percent) +
    labs(x ='', y = paste0('Age-specific contribution to COVID-19 deaths\nduring the baseline period'))
  ggsave(paste0(outdir, '-Contribution_ref_empr.png'), w = 9, h = 6)
  
  
  contribution_ref_adj[, loc_label := factor(loc_label, unique(tmp$loc_label))]
  ggplot(contribution_ref_adj, aes(x = loc_label, y = M)) + 
    geom_bar(aes(fill = M), stat = 'identity') +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
    facet_grid(age~division,  scales = "free_x", space = 'free_x') + 
    theme_bw() +
    theme(axis.text.x = element_text(angle= 45, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.position = 'none')  +
    scale_fill_viridis(option = 'E', trans = 'sqrt') + 
    scale_y_continuous(labels = scales::percent) + 
    labs(x ='', y = paste0('Age-specific contribution to COVID-19 deaths in\nage-standardised populations during the baseline period'))
  ggsave(paste0(outdir, '-Contribution_ref_adj.png'), w = 9, h = 6)

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

# savepdf <- function(fname, width=16, height=10)
# {
#   pdf(fname, width=width/2.54, height=height/2.54)
#   par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(0,0,0,0))
# }

# plot_posterior_plane = function(fit_cum, df_week, outdir)
# {
#   fit_samples = extract(fit_cum)
#   
#   ps <- c(0.5, 0.025, 0.975)
#   p_labs <- c('M','CL','CU')
#   tmp1 = as.data.table( reshape2::melt( fit_samples$beta ))
#   
#   row_name = 'week_index'
#   column_name = 'basis_function_index'
#   if(max(tmp1$Var3) == stan_data$A) column_name = 'Age'
#   
#   setnames(tmp1, c('Var2', 'Var3'), c(row_name, column_name))
#   tmp1 = tmp1[, list( 	q= quantile(value, prob=ps),
#                        q_label=p_labs), 
#               by=c(column_name, row_name)]	
#   tmp1 = dcast(tmp1, get(row_name) + get(column_name) ~ q_label, value.var = "q")
#   setnames(tmp1, c('row_name', 'column_name'), c(row_name, column_name))
#   
#   tmp1 = merge(tmp1, df_week, by = row_name)
# 
#   ggplot(tmp1, aes(x = date, y = get(column_name))) +
#     geom_raster(aes(fill = M))  + 
#     labs(x = 'Date', y = column_name, fill = 'Estimated posterior value') + 
#     scale_y_reverse(expand = c(0,0), breaks = seq(min(tmp1[,get(column_name)]), max(tmp1[,get(column_name)]), 2))  +
#     scale_x_date(expand = c(0,0), breaks = '2 months', date_labels = "%b-%y") + 
#     scale_fill_viridis_c(option = "A") + 
#     theme(legend.position='bottom')
#   ggsave(file = paste0(outdir, '-PlanePosterior_', Code, '.png'), w = 5, h = 5.5)
#   
# }


# 
# plot_posterior_predictive_checks_cumulative = function(tmp, outdir)
# {
#   
#   p = ggplot(tmp, aes(x = date)) + 
#     geom_point(aes(y = M)) +
#     geom_errorbar(aes(ymin = CL, ymax = CU)) +
#     geom_ribbon(aes(ymin = min_count_censored, ymax = max_count_censored), alpha = 0.5, fill = 'red') +
#     geom_point(aes(y = sum_count_censored), col = 'blue') +
#     facet_wrap(~age, scale = 'free_y') + 
#     labs(x = '', y = 'cumulative deaths')
#   ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_cum_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
#   
# }

# 
# compare_CDCestimation_JHU_error_plot = function(CDC_data, JHU_data, var.cum.deaths.CDC, outdir)
# {
#   # prepare JHU data
#   JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
#   
#   # prepare CDC estimations
#   CDC_data = select(as.data.table(CDC_data), code, date, var.cum.deaths.CDC)
#   setnames(CDC_data, var.cum.deaths.CDC, 'cumulative_deaths')
#   
#   # plot
#   JHUData[, source := 'JHU']
#   CDC_data[, source := 'CDC']
#   
#   tmp2 = rbind(JHUData, CDC_data)
#   tmp2 = subset(tmp2, code %in% unique(CDC_data$code) & date <= max(CDC_data$date))
#   
#   n_code = length(unique(tmp2$code))
#   
#   p = ggplot(tmp2, aes(x = date, y = cumulative_deaths)) + 
#     geom_line(aes(col = source), size = 1) +
#     facet_wrap(~code, nrow = length(unique(tmp2$code)), scale = 'free') + 
#     theme_bw() + 
#     scale_color_viridis_d(option = "B", direction = -1, end = 0.8) + 
#     scale_fill_viridis_d(option = "B", direction = -1, end = 0.8)
#   ggsave(p, file = paste0(outdir, '-comparison_JHU_CDC_uncertainty_', Code, '.png'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
# }

# plot_posterior_plane_with_data = function(deathByAge_1, deathByAge_2, outdir)
# {
#   library(magick)
#   
#   deathByAge_1 = subset(deathByAge_1, code == Code)
#   deathByAge_2 = subset(deathByAge_2, code == Code)
#   
#   range_d = range(na.omit(c(deathByAge_1$weekly.deaths, deathByAge_2$weekly.deaths)))
#   
#   p1 = ggplot(deathByAge_1, aes(x = date, y = age)) + 
#     geom_raster(aes(fill = weekly.deaths )) + 
#     theme_bw() +
#     scale_fill_viridis_c(trans = 'sqrt', limits = range_d) +
#     scale_x_date(expand = c(0,0), date_labels = c("%b %Y"), date_breaks = '1 month') + 
#     scale_y_discrete(expand = c(0,0)) + 
#     theme(legend.position = 'bottom',
#           axis.text.x = element_text(angle = 70, vjust  = 0.5), 
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank()) +
#     labs(x = '', y = 'Age groups first specification',
#          fill = 'Reported covid-19 deaths')
#   
#   p2 = ggplot(deathByAge_2, aes(x = date, y = age)) + 
#     geom_raster(aes(fill = weekly.deaths )) + 
#     theme_bw() +
#     scale_fill_viridis_c(trans = 'sqrt', limits = range_d) +
#     scale_x_date(expand = c(0,0), date_labels = c("%b %Y"), date_breaks = '1 month') + 
#     scale_y_discrete(expand = c(0,0)) + 
#     theme(legend.position = 'bottom',
#           axis.text.x = element_text(angle = 70, vjust  = 0.5), 
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank()) +
#     labs(x = '', y = 'Age groups second specification',
#          fill = 'Reported covid-19 deaths')
#   
#   p = ggpubr::ggarrange(p1, p2, nrow = 2,common.legend = T, legend = 'bottom')
#   ggsave(p, file = paste0(outdir, '-panel_right_', Code, '.png'), h = 8, w = 6)
#   
#   panel.left <- magick::image_read(paste0(outdir, '-PlanePosterior_', Code, '.png'))
#   panel.right <- magick::image_read(paste0(outdir, '-panel_right_', Code, '.png'))
#   
#   p = image_append(c( image_scale(panel.left, "2200"), panel.right) )
#   
#   savepdf(paste0(outdir, '-PlanePosterior_Data_', Code, '.pdf'), w = 14.5*2, h = 10*2)
#   plot(p)
#   dev.off()
#   
# }


compare_CDCestimation_DoH_age_prop_plot = function(tmp, outdir)
{
  tmp1 = select(tmp, date, age, prop.weekly.deaths)
  tmp1[, CL_prop := NA]
  tmp1[, CU_prop := NA]
  tmp1[, source := 'DoH']
  
  tmp2 = select(tmp, date, age, CL_prop,  CU_prop,  M_prop)
  setnames(tmp2, 'M_prop', 'prop.weekly.deaths')
  tmp2[, source := 'estimated']

  tmp1 = rbind(tmp1, tmp2)

  col = c('#5CC8D7FF', '#00203FFF')
  
  p = ggplot(tmp1, aes(x = date, y = prop.weekly.deaths)) +
    geom_ribbon(aes(ymin = CL_prop, ymax = CU_prop, fill = source), alpha = 0.5) +
    geom_line(aes(col = source), size = 0.5) +
    theme_bw() +
    scale_color_manual(values = col) +
    scale_fill_manual(values = col) +
    facet_wrap(~age,  scale = 'free', ncol = 1) +
    theme(legend.position = 'bottom',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x = element_text(angle = 90)) +
    labs(y = 'Age-specific contribution to COVID-19 weekly deaths', col = '', fill = '', x = '')
  ggsave(p, file = paste0(outdir, '-comparison_DoH_CDC_prop_', Code, '.png'), w = 5, h = 12, limitsize = F)

}

compare_CDCestimation_DoH_age_weekly_plot = function(tmp, outdir)
{
  tmp1 = select(tmp, date, age, weekly.deaths)
  tmp1[, CL_abs_weekly := NA]
  tmp1[, CU_abs_weekly := NA]
  tmp1[, source := 'DoH']
  
  tmp2 = select(tmp, date, age, CL_abs_weekly,  CU_abs_weekly,  M_abs_weekly)
  setnames(tmp2, 'M_abs_weekly', 'weekly.deaths')
  tmp2[, source := 'estimated']
  
  tmp1 = rbind(tmp1, tmp2)
  
  col = c('#5CC8D7FF', '#00203FFF')
  
  p = ggplot(tmp1, aes(x = date, y = weekly.deaths)) +
    geom_ribbon(aes(ymin = CL_abs_weekly, ymax = CU_abs_weekly, fill = source), alpha = 0.5) +
    geom_line(aes(col = source), size = 0.5) +
    theme_bw() +
    scale_color_manual(values = col) +
    scale_fill_manual(values = col) +
    facet_wrap(~age,  scale = 'free', ncol = 1) +
    theme(legend.position = 'bottom',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x = element_text(angle = 90)) +
    labs(y = 'Age-specific COVID-19 weekly deaths', col = '', fill = '', x = '')
  ggsave(p, file = paste0(outdir, '-comparison_DoH_CDC_weekly_', Code, '.png'), w = 5, h = 12, limitsize = F)
  
}

# 
# plot_contribution_comparison_method = function(tab_cd, data, model_name){
#   
#   tmp = copy(data)
#   setnames(tmp, 'prop_deaths', 'M')
#   set(tmp, NULL, 'CL', NA_real_)
#   set(tmp, NULL, 'CU', NA_real_)
#   
#   tmp = rbind(select(tab_cd, date, age, M, method), select(tmp, date, age, M, method))
#   tmp[, method := factor(method, c('observation', model_name))]
#   
#   p = ggplot(tmp, aes(x = date, y = age)) + 
#     geom_raster(aes(fill = M )) + 
#     theme_bw() +
#     scale_fill_viridis_c(option = "E") +
#     scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
#     scale_y_discrete(expand = c(0,0)) + 
#     labs(x = '', y = 'Age group',
#          fill = 'Contribution to\nweekly\nCOVID-19\ndeaths') + 
#     facet_grid(~method) + 
#     theme(legend.position = 'left',
#           axis.text.x = element_text(angle = 90), 
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(), 
#           legend.key = element_blank(), 
#           strip.background = element_rect(colour="white", fill="white")) 
#   
#   return(p)
# }

# plot_death_comparison_method = function(tab_d, data, model_name){
#   
#   tmp = copy(data)
#   setnames(tmp, 'weekly.deaths', 'M')
#   set(tmp, NULL, 'CL', NA_real_)
#   set(tmp, NULL, 'CU', NA_real_)
#   
#   tmp = rbind(select(tab_d, date, age, M, method), select(tmp, date, age, M, method))
#   tmp[, method := factor(method, c('observation', model_name))]
#   
#   range_wd = sqrt(range(na.omit(data$weekly.deaths)))
#   digits_cut= ifelse(range_wd[2]/100 > 1, 2, 1)
#   digits_cut= ifelse(range_wd[2]/10 > 1, digits_cut, 0)
#   p =  ggplot(tmp, aes(x = date, y = age)) + 
#     geom_raster(aes(fill = M )) + 
#     theme_bw() +
#     scale_fill_viridis_c(trans = 'sqrt', breaks = round(seq(range_wd[2]/2, range_wd[2], length.out = 3)^2, digits = -digits_cut) ) +
#     scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
#     scale_y_discrete(expand = c(0,0)) + 
#     labs(x = '', y = 'Age group',
#          fill = 'Reported\nCOVID-19\ndeaths') + 
#     facet_grid(~method) + 
#     theme(legend.position = 'left',
#           axis.text.x = element_text(angle = 90), 
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(), 
#           legend.key = element_blank(), 
#           strip.background = element_rect(colour="white", fill="white")) 
#   
#   return(p)
# }

# 
# plot_contribution_magnitude_all_states = function(contribution, death2agegroups, selected_states, outdir, nrow = 1, legend.position = 'bottom')
# {
#   
#   tmp = subset(contribution, code %in% selected_states)
#   tmp2 = subset(death2agegroups, code %in% selected_states)
#   
#   tmp1 = unique(select(tmp, code, loc_label))
#   tmp1[, code := factor(code, levels = selected_states)]
#   tmp1 = tmp1[order(code)]
#   
#   tmp[, loc_label := factor(loc_label, tmp1$loc_label)]
#   tmp2[, loc_label := factor(loc_label, tmp1$loc_label)]
#   
#   p1 = ggplot(tmp, aes(x= date) ) +
#     geom_line(aes(y = M, col = age)) + 
#     geom_point(aes(y = emp, col = age), size = 0.5) + 
#     geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
#     facet_wrap(~loc_label, nrow = 3) +
#     scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
#     theme_bw() + 
#     theme(strip.background = element_blank(),
#           panel.border = element_rect(colour = "black", fill = NA), 
#           axis.text.x = element_text(angle = 45, hjust =1), 
#           panel.grid.major = element_blank(), 
#           axis.title.x = element_blank()) +
#     labs(y = 'Estimated age-specific contribution\nto COVID-19 weekly deaths') + 
#     scale_y_continuous(labels = scales::percent_format()) +
#     scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
#     scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8)
#   
#   
#   p2 = ggplot(tmp2, aes(x= date) ) +
#     geom_line(aes(y = M, col = age)) + 
#     geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
#     facet_wrap(~loc_label, nrow = 3, scale = 'free_y') +
#     scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
#     theme_bw() + 
#     theme(strip.background = element_blank(),
#           panel.border = element_rect(colour = "black", fill = NA), 
#           axis.text.x = element_text(angle = 45, hjust =1), 
#           panel.grid.major = element_blank(), 
#           axis.title.x = element_blank()) +
#     labs( y = 'Estimated age-specific COVID-19\nweekly deaths', col = 'Age group', fill = 'Age group') + 
#     scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
#     scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8)
#   
#   
#   p = ggpubr::ggarrange(p2,p1, common.legend = T, legend = legend.position, nrow = nrow)
#   ggsave(p, file = paste0(outdir, '-ContributionMagnitude_panel.png'), w = 7, h = 8)
#   
#   
#   return(p)
# }
