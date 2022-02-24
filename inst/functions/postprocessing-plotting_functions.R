plot_posterior_predictive_checks = function(data, variable, outdir)
{
  
  data[, PPP := paste0('inside CI: ', round(mean(na.omit(inside.CI))*100, digits = 2), '%'), by = 'date']
  data[, date_ppp := paste0(as.character(date), ' - ', PPP)]
  
  for(Code in unique(data$code)){
    # posterior data check
    p = ggplot(subset(data, code == Code), aes(x = age)) + 
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
    p = ggplot(subset(data, code == Code), aes(x = date)) + 
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
  
}

plot_probability_deaths_age_contribution = function(tmp1, var_name, outdir, discrete = F, limits = NULL, wo_ytext = F, wo_legend = F){
  
  if(!discrete)
    tmp1[, age := as.numeric(age)]
  
  n_row = length(unique(tmp1$date))
  
  df = copy(tmp1)
  
  for(Code in unique(df$code)){
    
    tmp1 = subset(df, code == Code)
    
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
  }
  
  return(list(p1, p))
}

plot_sum_missing_deaths = function(tmp, outdir)
{
  
  for(Code in unique(tmp$code)){
    p = ggplot(subset(tmp, code == Code), aes(x = date)) + 
      geom_point(aes(y = M)) +
      geom_errorbar(aes(ymin = CL, ymax = CU)) +
      geom_ribbon(aes(ymin = min_count_censored, ymax = max_count_censored), alpha = 0.5, fill = 'red') +
      geom_point(aes(y = sum_count_censored), col = 'blue') +
      facet_wrap(~age, scale = 'free_y') + 
      labs(x = '', y = 'cumulative deaths')
    ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_missing_cum_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
    
  }

}

plot_sum_bounded_missing_deaths = function(tmp, outdir)
{
  
  tmp[, dummy := 1]
  
  for(Code in unique(tmp$code)){
    p = ggplot(subset(tmp, code == Code), aes(x = dummy) ) + 
      geom_point(aes(y = M)) +
      geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.25) +
      geom_errorbar(aes(ymin = min_count_censored, ymax = max_count_censored),  color = 'red', position = 'dodge', width = 0.25) +
      geom_point(aes(y = sum_count_censored), col = 'blue') +
      facet_wrap(~age, scale = 'free_y') + 
      labs(x = '', y = 'cumulative deaths')
    ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_missing_sum_cum_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
    
  }

}

plot_mortality_rate = function(mortality_rate_table, mortality_rate_table_continuous, outdir)
{
  mortality_rate_table_continuous[, age_cat := as.character(age)]
  mortality_rate_table_continuous[age_cat == '85', age_cat := '85+']
  
  for(Code in unique(mortality_rate_table$code)){
    p = ggplot(subset(mortality_rate_table, code == Code), aes(x = date)) + 
      geom_line(aes(y = M)) +
      geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) +
      facet_wrap(~age, scale = 'free_y') + 
      labs(x = '', y = 'COVID-19 mortality rate', col = '', fill = '') + 
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
      theme(legend.position = 'bottom',
            axis.text.x = element_text(angle = 90))  +
      theme_bw()
    ggsave(p, file = paste0(outdir, '-MortalityRate_', Code, '.png'), w = 7, h =7, limitsize = F)
    
    tmp <- subset(mortality_rate_table_continuous, code == Code & date == max(date))
    p = ggplot(subset(tmp, age != '85'), aes(x = age)) + 
      geom_line(aes(y = M)) +
      geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) +
      labs(y = 'COVID-19 mortality rate', col = '', fill = '') +
      theme_bw() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), limits = range(c(tmp$CL, tmp$CU + 0.001)), breaks = seq(0, max(c( tmp$CU)), 0.01), labels = scales::percent) +
      theme(axis.title.x = element_blank(),
            plot.margin=unit(c(5.5, 0, 5.5, 5.5), "pt")) 
    
    
    p1 <- ggplot(subset(tmp, age == '85'), aes(x = age_cat)) + 
      theme_bw() + 
      geom_point(aes(y = M)) +
      geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5, width = 0) +
      labs(x = 'age', y = 'COVID-19 mortality rate', col = '', fill = '') + 
      scale_x_discrete(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), limits = range(c(tmp$CL, tmp$CU + 0.001)), breaks = seq(0, max(c( tmp$CU)), 0.01), labels = scales::percent) +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(),
            plot.margin=unit(c(5.5, 5.5, 5.5, 0), "pt")) 
    
    p <- grid.arrange(p, p1, nrow = 1, widths = c(1, 0.1), bottom = 'Age')

    ggsave(p, file = paste0(outdir, '-MortalityRateContinuous_', Code, '.png'), w = 7, h =7, limitsize = F)
  }

}

compare_CDCestimation_DoH_age_plot = function(CDC_data, scraped_data, var.cum.deaths.CDC, outdir, overall = F)
{
  
  if('GA' %in% Code){
    scraped_data_ga <- reduce_agebands_scrapedData_GA(subset(scraped_data, code == 'GA'))
    scraped_data = rbind(subset(scraped_data, code != 'GA'), scraped_data_ga)
  }
  
  
  scraped_data = select(as.data.table(scraped_data), code, date, age, cum.deaths)
  scraped_data = subset(scraped_data, code %in% Code)
  scraped_data[, date := as.Date(date)]
  scraped_data[, CL := NA]
  scraped_data[, CU := NA]
  scraped_data = merge(scraped_data, unique(select(CDC_data, code, loc_label)), by = 'code')
  scraped_data = scraped_data[order(loc_label, code, age, date)]


  # age sorted
  ages = unique(select(CDC_data, age))
  ages[, age_from := gsub('(.+)-.*', '\\1', age)]
  ages[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  ages = ages[order(as.numeric(age_from))]$age
  
  # prepare CDC estimations
  CDC_data = select(as.data.table(CDC_data), code, loc_label, date, age, var.cum.deaths.CDC )
  setnames(CDC_data, var.cum.deaths.CDC, c('cum.deaths', 'CL', 'CU'))
  CDC_data = subset(CDC_data, date <= max(scraped_data$date))
  
  # plot
  scraped_data[, source := 'DoH']
  CDC_data[, source := 'Estimated']
  
  tmp2 = rbind(scraped_data, CDC_data)
  tmp2[, age := factor(age, levels = ages)]
  tmp2[, source := factor(source, levels = c('Estimated', 'DoH'))]
  
  col = c('#5CC8D7FF', '#00203FFF')
  
  for(Code in unique(tmp2$code)){
    p = ggplot(subset(tmp2, code == Code), aes(x = date, y = cum.deaths)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = source), alpha = 0.5) +
      geom_line(aes(col = source), size = 1) +
      theme_bw() + 
      scale_color_manual(values = col) + 
      scale_fill_manual(values = col) +
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
      facet_grid(age~loc_label,  scale = 'free_y') + 
      theme(legend.position = 'bottom',
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.text.x = element_text(angle = 45,hjust=1,vjust=1), 
            axis.title.x = element_blank()) + 
      labs(y = 'Cumulative COVID-19 deaths', col = '', fill = '')
    ggsave(p, file = paste0(outdir, '-comparison_DoH_CDC_uncertainty_', Code, '.png'), w = 3, h = 10, limitsize = F)
    
  }

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
  setnames(tmp1, 2:4, c('state_index', row_name, column_name))
  tmp1 = tmp1[, list(q = quantile(value, prob=ps), q_label=p_labs), by=c('state_index', column_name, row_name)]
  tmp1 = dcast(tmp1, state_index + get(row_name) + get(column_name) ~ q_label, value.var = "q")
  setnames(tmp1, c('row_name', 'column_name'), c(row_name, column_name))
  
  if(!is.null(stan_data$num_basis_rows)){
    row_name_1 = row_name 
    row_name = column_name 
    column_name = row_name_1
  }
  
  tmp1 = merge(tmp1, df_week, by = row_name)
  tmp1 = merge(tmp1, df_state, by = 'state_index')
  
  for(Code in unique(tmp1$code)){
    p = ggplot(subset(tmp1, code == Code), aes(x = date, y = get(column_name))) + 
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
    
  }
  
  
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
  
  tmp =   subset(mortality_rate, age == '85+')
  age = unique(tmp$age)
  tmp = tmp[order(M)]
  medianM =  tmp[,median(M)]
  mortality_rate[, loc_label := factor(loc_label, tmp$loc_label)]
  
  mortality_rate <- mortality_rate[age != '0-24']
  mortality_rate[, `Age group` := age]

  p <- ggplot(mortality_rate, aes(x=loc_label, y = M)) + 
    geom_bar(aes(fill = M), stat="identity") +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
    scale_fill_gradient2(low= 'darkturquoise', high = 'darkred', mid = 'beige', midpoint = medianM, labels = scales::percent_format()) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 70,hjust=1,vjust=1), 
          legend.position = 'bottom', 
          panel.grid.minor= element_blank(), 
          panel.grid.major.x= element_blank(), 
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + 
    labs(x = '', y = paste0('Predicted COVID-19 attributable mortality rates as of ', format(unique(mortality_rate$date), '%B %Y')), 
         fill = '') + 
    scale_y_continuous(expand =expansion(mult = c(0, .05)), labels = scales::percent_format())+ 
    guides(fill = guide_colourbar(barwidth = 10,  barheight = 1))
  
  
  if(length(unique(mortality_rate$code)) <= 10){
    p <- p +  facet_grid(.~`Age group`, label = 'label_both')
    ggsave(p, file = paste0(outdir, paste0('-MortalityRate_allages.png')), w = 6.5, h = 3.5)
    
  } else{
    p <- p +  facet_wrap(~`Age group`, nrow =4, label = 'label_both', scales = 'free_y') 
     
    ggsave(p, file = paste0(outdir, paste0('-MortalityRate_allages_large.png')), w = 7, h = 9)
    
  }
  # +
  #   labs(y = '')

}

plot_mortality_rate_continuous_all_states = function(mortality_rate, outdir)
{
  
  mortality_rate[, age := as.numeric(age)]
  mortality_rate[, age_cat := as.character(age)]
  mortality_rate[age_cat == '85', age_cat := '85+']
  
  tmp <- subset(mortality_rate, code %in% selected_codes)
  
  p <- ggplot(subset(tmp, age != '85'), aes(x=age)) + 
    geom_line(aes(y = M, col = loc_label)) +
    geom_ribbon(aes(ymin=CL, ymax=CU, fill = loc_label), alpha = 0.4) + 
    theme_bw() +
    theme(legend.position = 'bottom', 
          panel.grid.minor= element_blank(), 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          plot.margin = unit(c(5.5,0,5.5,5.5), "pt")) + 
    scale_color_jcolors('pal8') + 
    scale_fill_jcolors('pal8') + 
    scale_y_continuous(expand = c(0,0), labels = scales::percent_format(), limits = range(c(tmp$CL, tmp$CU + 0.001)), 
                       breaks = seq(0, max(tmp$CU), 0.01)) +
    scale_x_continuous(expand = c(0,0)) +
    labs(y = paste0('Predicted COVID-19 attributable mortality\nrates as of ', format(unique(mortality_rate$date), '%b %Y')),
         col = '', fill = '', x = 'Age')
  
  p1 <- ggplot(subset(tmp, age == '85'), aes(x=age_cat)) + 
    geom_errorbar(aes(ymin=CL, ymax=CU), width = 0) + 
    geom_point(aes(y = M, col = loc_label)) +
    theme_bw() +
    theme(
      panel.grid.minor= element_blank(), 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.title.x = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          plot.margin = unit(c(5.5,5.5,18,0), "pt"),
          legend.position = 'none',
          ) + 
    scale_color_jcolors('pal8') + 
    scale_fill_jcolors('pal8') + 
    scale_y_continuous(expand = c(0,0), labels = scales::percent_format(), limits = range(c(tmp$CL, tmp$CU + 0.001)), 
                       breaks = seq(0, max(tmp$CU), 0.01)) +
    scale_x_discrete(expand = c(0,0)) +
    labs(y = paste0('Predicted COVID-19 attributable mortality\nrates as of ', format(unique(mortality_rate$date), '%b %Y')),
         col = '')
  
  p <- ggarrange(p, p1, nrow = 1, common.legend = T, legend = 'bottom', widths = c(1, 0.1))
  ggsave(p, file = paste0(outdir, paste0('-MortalityRateContinuous_allages_selectedstates.png')), w = 6, h = 4)
  
  
  ###
  
  tmp <- subset(mortality_rate, code == 'NY')
  
  if(nrow(tmp) == 0){
    tmp <- subset(mortality_rate, code == unique(mortality_rate)[1])
  }
  
  tmp[, State := loc_label]
  p <- ggplot(subset(tmp, age != '85'), aes(x=age)) + 
    geom_line(aes(y = M)) +
    geom_ribbon(aes(ymin=CL, ymax=CU), alpha = 0.4) + 
    theme_bw() + 
    facet_wrap(~State, label = 'label_both') +
    theme(legend.position = 'bottom', 
          # axis.title = element_text(size = rel(1.2)),
          # axis.text = element_text(size = rel(1.1)),
          axis.title.y =element_text(size = rel(1)),
          strip.text =element_text(size = rel(1)),
          axis.title.x =element_text(size = rel(1)),
          panel.grid.minor= element_blank(), 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + 
    scale_y_continuous(expand = c(0,0), labels = scales::percent_format()) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = 'Age', y = paste0('Predicted COVID-19 attributable mortality\nrates as of ', format(unique(mortality_rate$date), '%B %Y'))) 
   ggsave(p, file = paste0(outdir, paste0('-MortalityRateContinuous_allages_NY.png')), w = 6.5, h = 3)
  
}

plot_mortality_all_states = function(death, lab = 'allStates', outdir)
{
  
  df = as.data.table( reshape2::melt(select(death, loc_label, code, date, age, emp), id.vars = c('loc_label', 'code', 'date', 'age')) )
  df[, variable2 := 'CDC data']
  df[, `Age group` := age]
  
  death[, dummy := 'Posterior median prediction\nusing age-aggregated JHU data\nto adjust for reporting delays']
  # death[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  

  death[, `Age group` := age]
  
  p <- ggplot(subset(death) ) +
    facet_grid(loc_label~`Age group`, scale = 'free') +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
    theme_bw() + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 45, hjust =1), 
          panel.grid.major = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = rel(1.1)), 
          legend.text = element_text(size = rel(1.1)), 
          strip.text = element_text(size = rel(1)), 
          legend.position = 'bottom') +
    labs( y = 'Predicted age-specific COVID-19 attributable weekly deaths', shape = '',
          linetype = '', col = '', fill = '') + 
    scale_shape_manual(values = 16) +
    scale_linetype_manual(values = 2) + 
    guides(shape = guide_legend(override.aes = list(size=1, stroke =1.5), order = 2),
           linetype = guide_legend(order = 3), 
           col = guide_legend(order = 1),
           fill = guide_legend(order = 1))
  
  if(all(unique(death$code) %in% selected_codes) | all(unique(death$code) %in% selected_10_codes) ){
    
    if(all(unique(death$code) %in% selected_codes))
      colfunc <- jcolors("pal8")[1:length(unique(death$code))]

    if(all(unique(death$code) %in% selected_10_codes))
      colfunc <- jcolors("pal8")[5:10][1:length(unique(death$code))]
    
    p <- p + 
      geom_point(data = subset(df), aes(x= date, y = value, shape= variable2), col = 'black', fill = 'black', size = 0.9, alpha = 0.7) + 
      geom_line(aes(x= date, y = M, col = loc_label), show.legend = F) +
      geom_ribbon(aes(x= date, ymin = CL, ymax = CU, fill = loc_label), alpha = 0.5, show.legend = F) + 
      scale_color_manual(values = as.character(colfunc)) + 
      scale_fill_manual(values = as.character(colfunc)) 
    ggsave(p, file = paste0(outdir, paste0('-Mortality_', lab, '.png')), w = 7.5, h = 8)
    
  } else {
    p <- p + 
      geom_point(data = subset(df), aes(x= date, y = value, shape= variable2), col = 'darkred', fill = 'darkred', size = 0.9, alpha = 0.6)+ 
      geom_line(aes(x= date, y = M, col = loc_label), show.legend = F) +
      geom_ribbon(aes(x= date, ymin = CL, ymax = CU, fill = loc_label), alpha = 0.4, show.legend = F) +
    scale_color_manual(values = rep('black', length(unique(death$code)))) + 
      scale_fill_manual(values = rep('black', length(unique(death$code)))) 
    ggsave(p, file = paste0(outdir, paste0('-Mortality_', lab, '.png')), w = 7.5, h = 4 + 0.72 * length(unique(death$code)))
    
  }


}


base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

plot_contribution_continuous_comparison_method_with_data = function(tab_cc, tab_d, data, 
                                                          selected_method, model_name, show.method = T, heights= c(0.4,0.6), outdir = NULL){
  
  dates = unique(tab_cc$date)
  dates = dates[seq(3, length(dates)-3, length.out =3)]
  
  plot_labels <- format(dates,  "%d-%b-%y")
  
  tmp2 = subset(tab_cc, date %in% dates)
  tmp2[, age_c := as.numeric(age)]
  tmp2[, method := factor(method,  model_name)]
  tmp22 = copy(tmp2)
  
  tmp3 = as.data.table(tab_d)
  tmp3[, age_c := as.numeric(age)]
  tmp33 = copy(tmp3)
  
  tmp= tab_d[, list(sumM = sum(mean)), by = c('date', 'method', 'code')]
  tmpp = copy(tmp)
  
  dataa = copy(data)
  
  for(Code in unique(tmp22$code)){
    
    tmp = subset(tmpp, code == Code)
    tmp2 = subset(tmp22, code == Code)
    tmp3 = subset(tmp33, code == Code)
    data = subset(dataa, code == Code)
    
    limit_SE = range(subset(tmp2, method == selected_method)$CL, subset(tmp2, method == selected_method)$CU)
    
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
      # coord_cartesian(ylim = limit_SE, xlim = range(tmp2$age_c)) +
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
    tmp3[, dummy_new := 'JHU\noverall\nweekly\ndeaths']
    p2 = ggplot(tmp3, aes(x = date)) + 
      theme_bw() +
      geom_step(aes(x = date+ 3.5,y = emp_JHU, linetype = dummy_new), direction =  "vh") + 
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
    pl = ggplot(tmp3, aes(x = date, y = M, col = age_c)) + geom_step(aes(x = date+ 3.5,y = emp_JHU, linetype = dummy_new), direction =  "vh") + 
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
    
    if(!is.null(outdir))
      ggsave(p, file = paste0(outdir, '-panel_plot_1_', Code, '.png'), w = 9, h = 7)
  }
  
  return(p)
  
}


plot_contribution_continuous_comparison_method = function(tab_cc, selected_method, model_name){
  
  dates = unique(tab_cc$date)
  dates = dates[seq(3, length(dates)-3, length.out =3)]
  
  plot_labels <- format(dates,  "%d-%b-%y")
  
  tmp2 = subset(tab_cc, date %in% dates)
  tmp2[, age_c := as.numeric(age)]
  tmp2[, method := factor(method,  model_name)]
  tmp2[, date_name := format(date, c("%B %d, %Y"))]
  
  limit_SE = range(subset(tmp2, method == selected_method)$CL, subset(tmp2, method == selected_method)$CU)
  
  p1 = ggplot(tmp2, aes(x = age_c)) + 
    geom_line(aes(y = M)) +
    geom_ribbon(aes(ymin= CL, ymax = CU), alpha = 0.5) +
    theme_bw() +
    labs(y = 'Estimated age-specific contribution to COVID-19 weekly deaths', x = "Age") + 
    facet_grid(date_name~method) +
    # coord_cartesian(ylim = limit_SE, xlim = range(tmp2$age_c)) +
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) + 
    theme(panel.border = element_rect(colour = "black", fill = NA), 
          # legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"), 
          panel.grid.major = element_blank(),
          axis.title.y = element_text(size = rel(1.1)), 
          plot.margin = unit(c(5.5,5.5,25,5.5), "pt"), 
          legend.position = 'none')
  
  return(p1)
  
}

plot_contribution_ref_all_states = function(contribution_ref, contribution_ref_adj, outdir){
  
  # tmp = subset(contribution_ref, age == '85+')
  # tmp = tmp[order(M)]
  
  # contribution_ref[, loc_label := factor(loc_label, unique(tmp$loc_label))]
  # contribution_ref_adj[, loc_label := factor(loc_label, unique(tmp$loc_label))]
  
  limits = c(0, max(contribution_ref$CU) + max(contribution_ref$CU) * 0.1)
  p1 = ggplot(contribution_ref, aes(x = loc_label, y = M)) + 
    geom_bar(aes(fill = M), stat = 'identity') +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
    facet_grid(.~age, scale = 'free_x') + 
    theme_bw() +
    theme(axis.text.x = element_text(angle= 45, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.position = 'none', axis.title.x = element_blank())  +
    scale_fill_viridis(option = 'E', trans = 'sqrt') + 
    scale_y_continuous(labels = scales::percent, limits = limits, expand = c(0,0)) + 
    labs(y = paste0('Age-specific contribution to COVID-19 deaths\nduring the baseline period'))
  ggsave(p1, file = paste0(outdir, '-Contribution_ref.png'), w = 8, h = 3.5)
  
  p1 = ggplot(contribution_ref, aes(x = loc_label, y = M)) + 
    geom_bar(aes(fill = M), stat = 'identity') +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
    facet_grid(.~age, scale = 'free_x') + 
    geom_point(aes(y = emp_est), col = 'tomato1', size = 1) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle= 45, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.position = 'none', axis.title.x = element_blank())  +
    scale_fill_viridis(option = 'E', trans = 'sqrt') + 
    scale_y_continuous(labels = scales::percent, limits = limits, expand = c(0,0)) + 
    labs(y = paste0('Age-specific contribution to COVID-19 deaths\nduring the baseline period'))
  ggsave(p1, file = paste0(outdir, '-Contribution_ref_empr.png'), w = 8, h = 3.5)
  
  limits = c(0, max(contribution_ref_adj$CU) + max(contribution_ref_adj$CU) * 0.1)
  
  p1 = ggplot(contribution_ref_adj, aes(x = loc_label, y = M)) + 
    geom_bar(aes(fill = M), stat = 'identity') +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
    facet_grid(.~age, scale = 'free_x') + 
    theme_bw() +
    theme(axis.text.x = element_text(angle= 45, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.position = 'none', axis.title.x = element_blank())  +
    scale_fill_viridis(option = 'E', trans = 'sqrt') + 
    scale_y_continuous(labels = scales::percent, limits = limits, expand = c(0,0)) + 
    labs(y = paste0('Age-specific contribution to COVID-19 deaths\nin age-standardised populations\nduring the baseline period'))
  ggsave(p1, file = paste0(outdir, '-Contribution_ref_adj.png'), w = 8, h = 3.5)
  
  
}


compare_CDCestimation_DoH_age_prop_plot = function(tmp, outdir)
{
  tmp1 = select(tmp, code, loc_label, date, age, prop.weekly.deaths)
  tmp1[, CL_prop := NA]
  tmp1[, CU_prop := NA]
  tmp1[, source := 'DoH']
  
  tmp2 = select(tmp, code, loc_label, date, age, CL_prop,  CU_prop,  M_prop)
  setnames(tmp2, 'M_prop', 'prop.weekly.deaths')
  tmp2[, source := 'estimated']

  tmp1 = rbind(tmp1, tmp2)

  col = c('#5CC8D7FF', '#00203FFF')
  
  for(Code in unique(tmp1$code)){
    
    p = ggplot(subset(tmp1, code == Code), aes(x = date, y = prop.weekly.deaths)) +
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

}

compare_CDCestimation_DoH_age_weekly_plot = function(tmp, outdir)
{
  tmp1 = select(tmp, code, loc_label, date, age, weekly.deaths)
  tmp1[, CL_abs_weekly := NA]
  tmp1[, CU_abs_weekly := NA]
  tmp1[, source := 'DoH']
  
  tmp2 = select(tmp,code, loc_label,  date, age, CL_abs_weekly,  CU_abs_weekly,  M_abs_weekly)
  setnames(tmp2, 'M_abs_weekly', 'weekly.deaths')
  tmp2[, source := 'estimated']
  
  tmp1 = rbind(tmp1, tmp2)
  
  col = c('#5CC8D7FF', '#00203FFF')
  
  for(Code in unique(tmp1$code)){
    
    p = ggplot(subset(tmp1, code == Code), aes(x = date, y = weekly.deaths)) +
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

}

plot_lambda_table <- function(lambda_table, outdir){
  
  state_indices = unique(lambda_table$state_index)
  
  if(length(state_indices) < 10){
    
    p <- ggplot(lambda_table, aes(x = date, col = type)) + 
      geom_point(aes(y = M), position = position_dodge(5), size = 0.75) + 
      geom_errorbar(aes(ymin = CL, ymax = CU), position = position_dodge(5), width = 0) + 
      facet_grid(loc_label~., scales = 'free_y') + 
      theme_bw() +
      labs(y = expression(lambda), col = '') +
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
      theme(legend.position = 'bottom',
            axis.text.x = element_text(angle = 70, hjust =1), 
            strip.background = element_blank(), 
            axis.title.x = element_blank()) 
    ggsave(p, file = paste0(outdir, '-lambda_prior_posterior.png'), w = 9, h = 8)
    
  } else{
    
    mid_point = length(state_indices) / 2
    state_indices_list = list(state_indices[1:mid_point], 
                              state_indices[(mid_point+1):length(state_indices)])
    
    for(i in 1:2){
      tmp <- subset(lambda_table, state_index %in% state_indices_list[[i]])
      # tmp <- subset(tmp, date < as.Date('2020-12-01'))
      p <- ggplot(tmp, aes(x = date, col = type)) + 
        geom_point(aes(y = M), position = position_dodge(5), size = 0.75) + 
        geom_errorbar(aes(ymin = CL, ymax = CU), position = position_dodge(5), width = 0) + 
        facet_grid(loc_label~., scales = 'free_y') + 
        theme_bw() +
        labs(y = expression(lambda), col = '') +
        scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
        theme(legend.position = 'bottom',
              axis.text.x = element_text(angle = 70, hjust =1), 
              strip.background = element_blank(), 
              axis.title.x = element_blank()) 
      ggsave(p, file = paste0(outdir, '-lambda_prior_posterior_part', i, '.png'), w = 9, h = 8)
    }
    
  }
  
  
}

plot_var_base_model_table <- function(loc_label, outdir){
  
  scales_y <- list(
    `zeta` = scale_y_log10(),
    `gamma[1]` = scale_y_continuous(),
    `gamma[2]` = scale_y_continuous(),
    `nu` = scale_y_continuous()
  )  
  
  p <- ggplot(var_base_model_table, aes(x = loc_label, col = type)) + 
    geom_point(aes(y = M), position = position_dodge(0.5), size = 0.75) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), position = position_dodge(0.5), width = 0.1) + 
    facet_grid_sc(rows = vars(math_name), cols = NULL, labeller = label_parsed,scales = list(y = scales_y)) +
    theme_bw() +
    labs(y = '', col = '') +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 70, hjust =1), 
          strip.background = element_blank(), 
          axis.title.x = element_blank(), axis.title.y = element_blank()) 
  scale_y_discrete(labels = label_parse(), limits=rev) 
  ggsave(p, file = paste0(outdir, '-var_base_model_prior_posterior.png'), w = 6, h = 6)
  
}

plot_lambda_table_sensitivity <- function(lambda_table, outdir){
  p <- ggplot(lambda_table, aes(x = date, col = model)) + 
    geom_point(aes(y = M), position = position_dodge(7), size = 0.75) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), position = position_dodge(7), width = 0) + 
    facet_wrap(~type, nrow =2) + 
    theme_bw() +
    labs(y = expression(lambda), col = '') +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 70, hjust =1, size = rel(1.2)), 
          strip.background = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.y = element_text(size = rel(1.2)),
          axis.title.y = element_text(size = rel(1.3)),
          legend.text =  element_text(size = rel(1.3)),
          strip.text =  element_text(size = rel(1.3))) + 
    scale_color_manual(values = colors) 
  ggsave(p, file = paste0(outdir, '-lambda_prior_posterior_sensitivity.png'), w = 9, h = 6)
}


plot_var_base_model_table_sensitivity <- function(table, lab, outdir){
  
  scales_y <- list(
    `zeta` = scale_y_log10(),
    `gamma[1]` = scale_y_continuous(),
    `gamma[2]` = scale_y_continuous(),
    `nu` = scale_y_continuous(trans = 'log1p')
  )
  
  p <- ggplot(table, aes(x = 1, col = model)) + 
    geom_point(aes(y = M), position = position_dodge(0.5), size = 0.75) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), position = position_dodge(0.5), width = 0.1) + 
    facet_grid_sc(rows = vars(math_name), cols = vars(type), labeller = label_parsed,scales = list(y = scales_y)) +
    theme_bw() +
    labs(y = '', col = '') +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(size = rel(1.2)), 
          strip.background = element_blank(), 
          axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    scale_x_discrete(labels = function(l) parse(text=l))+ 
    scale_color_manual(values = colors) + 
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  
  ggsave(p, file = paste0(outdir, paste0('-var_base_model_prior_posterior_sensitivity_', lab, '.png')), w = 5.5, h = 2 + 1 * length(unique(table$math_name)))
  
  return(p)
}

