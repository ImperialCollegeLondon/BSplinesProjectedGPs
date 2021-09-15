
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
    geom_line(aes(y = M, col = age)) + 
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
    labs( y = 'Predicted age-specific COVID-19 attributable weekly deaths', col = 'Age group', fill = 'Age group', shape = '') + 
    scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
    scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8) +
    scale_shape_manual(values = 21) + 
    guides(shape = guide_legend(override.aes = list(size=1, stroke = 1), order = 1), 
           fill = guide_legend(order=2), col = guide_legend(order =2))
  ggsave(paste0(outdir, paste0('-Mortality_allStates_1.png')), w = 9, h = 12)
  
  tmp = subset(death, code %in% codes[(mid+1):ncodes])
  df1 =subset(df, code %in% codes[(mid+1):ncodes])
  ggplot(tmp, aes(x= date) ) +
    geom_line(aes(y = M, col = age)) + 
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
    labs( y = 'Predicted age-specific COVID-19 attributable weekly deaths', col = 'Age group', fill = 'Age group', shape = '') + 
    scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
    scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8)  +
    scale_shape_manual(values = 21) + 
    guides(shape = guide_legend(override.aes = list(size=1, stroke = 1), order = 1), 
           fill = guide_legend(order=2), col = guide_legend(order =2))
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


  tmp3[, age := factor(age, levels = rev(levels(tmp3$age)))]
  tmp3[, dummy := 'JHU\noverall\nweekly\ndeaths']
  p2 = ggplot(tmp3, aes(x = date)) + 
    theme_bw() +
    geom_step(aes(x = date+ 3.5,y = emp_JHU, linetype = dummy), direction =  "vh") + 
    labs(y = 'Posterior mean of the age-specific\nCOVID-19 weekly deaths', x = "", color = '', shape = '', fill = 'Age') + 
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
  contribution_ref_adj[, loc_label := factor(loc_label, unique(tmp$loc_label))]
  
  uscensus = sort(unique(contribution_ref$division))
  contribution_ref[, dummy := ifelse(division %in% uscensus[c(1,8,2,6)], 1, 2)]
  contribution_ref_adj[, dummy := ifelse(division %in% uscensus[c(1,8,2,6)], 1, 2)]
  
  limits_ref = range(c(0, contribution_ref$CU))
  limits_ref_adj = range(c(0, contribution_ref$CU))
  
  for(a in unique(contribution_ref$age)){
    
    tmp = subset(contribution_ref, age == a)
    tmp = tmp[order(M)]
    tmp[, loc_label := factor(loc_label, levels = tmp$loc_label)]
    limits = c(0, max(tmp$CU) + max(tmp$CU) * 0.01)
    p1 = ggplot(subset(tmp, dummy == 1), aes(x = loc_label, y = M)) + 
      geom_bar(aes(fill = M), stat = 'identity') +
      geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
      facet_wrap(~division, scale = 'free_x', nrow = 1) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle= 45, hjust = 1), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"), 
            legend.position = 'none', axis.title = element_blank())  +
      scale_fill_viridis(option = 'E', trans = 'sqrt', limits = limits_ref) + 
      scale_y_continuous(labels = scales::percent, limits = limits, expand = c(0,0)) 
    # p1 = adjust_facet_size(p1)
    
    p2 = ggplot(subset(tmp, dummy == 2), aes(x = loc_label, y = M)) + 
      geom_bar(aes(fill = M), stat = 'identity') +
      geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
      facet_wrap(~division, scale = 'free_x', nrow = 1) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle= 45, hjust = 1), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"), 
            legend.position = 'none', axis.title = element_blank())  +
      scale_fill_viridis(option = 'E', trans = 'sqrt', limits = limits_ref) + 
      scale_y_continuous(labels = scales::percent, limits = limits, expand = c(0,0)) 
    # p2 = adjust_facet_size(p2)
    
    p = gridExtra::grid.arrange(p1, p2, nrow = 2, 
                                left=   textGrob(paste0('Age-specific contribution from individuals aged ', a, '\nto COVID-19 deaths during the baseline period'), gp = gpar(fontsize = 10), rot = 90))
    ggsave(p, file = paste0(outdir, '-Contribution_ref_',a,'.png'), w = 6, h = 5)
    
    
    tmp = subset(contribution_ref, age == a)
    tmp = tmp[order(M)]
    tmp[, loc_label := factor(loc_label, levels = tmp$loc_label)]
    limits = c(0, max(tmp$CU) + max(tmp$CU) * 0.01)
    p1 = ggplot(subset(tmp, dummy == 1), aes(x = loc_label, y = M)) + 
      geom_bar(aes(fill = M), stat = 'identity') +
      geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
      geom_point(aes(y = emp_est), col = 'tomato1', size = 0.6) + 
      facet_wrap(~division, scale = 'free_x', nrow = 1) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle= 45, hjust = 1), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"), 
            legend.position = 'none', axis.title = element_blank())  +
      scale_fill_viridis(option = 'E', trans = 'sqrt', limits = limits_ref) + 
      scale_y_continuous(labels = scales::percent, limits = limits, expand = c(0,0)) 
    # p1 = adjust_facet_size(p1)
    
    p2 = ggplot(subset(tmp, dummy == 2), aes(x = loc_label, y = M)) + 
      geom_bar(aes(fill = M), stat = 'identity') +
      geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
      geom_point(aes(y = emp_est), col = 'tomato1', size = 0.6) + 
      facet_wrap(~division, scale = 'free_x', nrow = 1) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle= 45, hjust = 1), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"), 
            legend.position = 'none', axis.title = element_blank())  +
      scale_fill_viridis(option = 'E', trans = 'sqrt', limits = limits_ref) + 
      scale_y_continuous(labels = scales::percent, limits = limits, expand = c(0,0)) 
    # p2 = adjust_facet_size(p2)
    
    p = gridExtra::grid.arrange(p1, p2, nrow = 2, 
                                left= textGrob(paste0('Age-specific contribution from individuals aged ', a, '\nto COVID-19 deaths during the baseline period'), gp = gpar(fontsize = 10), rot = 90))
    ggsave(p, file = paste0(outdir, '-Contribution_ref_empr_',a,'.png'), w = 6, h = 5)
    
    
    tmp = subset(contribution_ref_adj, age == a)
    tmp = tmp[order(M)]
    tmp[, loc_label := factor(loc_label, levels = tmp$loc_label)]
    limits = c(0, max(tmp$CU) + max(tmp$CU) * 0.01)
    p1 = ggplot(subset(tmp, dummy == 1), aes(x = loc_label, y = M)) + 
      geom_bar(aes(fill = M), stat = 'identity') +
      geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
      facet_wrap(~division, scale = 'free_x', nrow = 1) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle= 45, hjust = 1), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"), 
            legend.position = 'none', axis.title = element_blank())  +
      scale_fill_viridis(option = 'E', trans = 'sqrt', limits = limits_ref_adj) + 
      scale_y_continuous(labels = scales::percent, limits = limits, expand = c(0,0)) 
    # p1 = adjust_facet_size(p1)
    
    p2 = ggplot(subset(tmp, dummy == 2), aes(x = loc_label, y = M)) + 
      geom_bar(aes(fill = M), stat = 'identity') +
      geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
      facet_wrap(~division, scale = 'free_x', nrow = 1) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle= 45, hjust = 1), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"), 
            legend.position = 'none', axis.title = element_blank())  +
      scale_fill_viridis(option = 'E', trans = 'sqrt', limits = limits_ref_adj) + 
      scale_y_continuous(labels = scales::percent, limits = limits, expand = c(0,0)) 
    # p2 = adjust_facet_size(p2)
    
    p = gridExtra::grid.arrange(p1, p2, nrow = 2, 
                                left= textGrob(paste0('Contribution from individuals aged ', a, ' to COVID-19 deaths in\nage-standardised populations during the baseline period'), gp = gpar(fontsize = 10), rot = 90))
    ggsave(p, file = paste0(outdir, '-Contribution_ref_adj_',a,'.png'), w = 6, h = 5)
    
  }
  
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

plot_vaccine_effects <- function(vaccine_data, weeklydv, weeklyf, weeklyphi, outdir){
  
  vaccine_data[, age_index := which(df_age_vaccination$age_from <= age & df_age_vaccination$age_to >= age), by = 'age']
  tmp <- unique(select(vaccine_data, code, date, prop, age_index))
  
  delay = 7*2
  weeklydv[, date := date + delay]
  weeklyf[, date := date + delay]
  weeklyphi[, date := date + delay]
  
  tmp1 <- merge(tmp, weeklydv, by = c('code', 'date', 'age_index'))
  p <- ggplot(tmp1, aes(x = prop)) + 
    geom_point(aes(y = M)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU)) + 
    geom_point(aes(y = emp), col = 'red') + 
    facet_wrap(~age, nrow = nrow(df_age_vaccination)) + 
    labs(x = 'Proportion of vaccinated', y = 'Weekly deaths 2 weeks later') + 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white")) 
  ggsave(p, file = paste0(outdir, "-vaccination_effects_weekly_deaths_", Code,".png") , w= 7, h = 8, limitsize = FALSE)
  
  
  tmp1 <- merge(tmp, weeklyf, by = c('code', 'date', 'age_index'))
  p <- ggplot(tmp1, aes(x = date, col = prop)) + 
    geom_point(aes(y = M, col = prop)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = prop)) + 
    facet_wrap(~age, nrow = nrow(df_age_vaccination)) + 
    labs(col = 'Proportion of\nvaccinated', y = 'f 2 weeks later', x = '')+ 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white")) 
  ggsave(p, file = paste0(outdir, "-vaccination_effects_f_", Code,".png") , w= 7, h = 8, limitsize = FALSE)
  
  tmp1 <- merge(tmp, weeklyphi, by = c('code', 'date', 'age_index'))
  p <- ggplot(tmp1, aes(x = date, col = prop)) + 
    geom_point(aes(y = M, col = prop)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = prop)) + 
    geom_point(aes(y = prop_deaths), col = 'red') + 
    facet_wrap(~age, nrow = nrow(df_age_vaccination)) + 
    labs(col = 'Proportion of\nvaccinated', y = 'Contribution to weekly deaths 2 weeks later', x = '')+ 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white")) 
  ggsave(p, file = paste0(outdir, "-vaccination_effects_contribution_", Code,".png") , w= 7, h = 8, limitsize = FALSE)
  
  p <- ggplot(tmp1, aes(x = prop)) + 
    geom_point(aes(y = M)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU)) + 
    geom_point(aes(y = prop_deaths), col = 'red') + 
    facet_wrap(~age, nrow = nrow(df_age_vaccination)) + 
    labs(x = 'Proportion of vaccinated', y = 'Contribution to weekly deaths 2 weeks later')+ 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white")) 
  ggsave(p, file = paste0(outdir, "-vaccination_effects_contribution2_", Code,".png") , w= 7, h = 8, limitsize = FALSE)
  
}

