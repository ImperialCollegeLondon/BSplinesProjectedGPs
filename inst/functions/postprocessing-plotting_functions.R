
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

plot_estimate_vaccine <- function(data_res5, label, lab_short, outdir){
  # data_res5[, var := paste0(var,'_', age_param)]

  ages = unique(data_res5$age_index_recipient)
  data_res5[, var2:= paste0('Effect of the ', label,' of\nfully vaccinated individuals aged ', age_index_source) ]
  data_res5[, var3:= factor(paste0('among individuals\naged ', age_index_recipient), 
                            levels = c(paste0('among individuals\naged ', rev(ages))))]
  
  p = ggplot(data_res5, aes(y = var3)) + 
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'grey50') + 
    geom_errorbarh(aes(xmin = CL, xmax = CU),  height = 0.2) + 
    geom_point(aes(x = M, col = age_index_recipient, shape = age_index_source), size = 2, stroke = 0.1) + 
    theme_bw() + 
    facet_wrap(~var2, nrow = 2) +
    theme(axis.title.x = element_blank(), 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'none') + 
    labs(y = 'On relative COVID-19\nweekly deaths during the\nSummer 2021 resurgence period')
  ggsave(p, file = paste0(outdir, '-vaccine_effect_estimate_', lab_short, '.png'), w = 6, h = 2.5)
  
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
            panel.grid.major= element_blank(), 
            axis.title.x = element_blank()) + 
      scale_y_continuous(expand = c(0,0), labels = scales::percent_format()) +
      coord_cartesian(ylim = c(0,max(tmp$CU)+0.001)) +
      labs(y = paste0('Predicted COVID-19 attributable mortality rates\namong individuals ', age, ' as of ', format(unique(tmp$date), '%b %Y')))
    ggsave(paste0(outdir, paste0('-MortalityRate_', age, '.png')), w = 8, h = 5)
  }
  
  tmp =   subset(mortality_rate, age == '85+')
  age = unique(tmp$age)
  tmp = tmp[order(M)]
  medianM =  tmp[,median(M)]
  mortality_rate[, loc_label := factor(loc_label, tmp$loc_label)]
  ggplot(mortality_rate, aes(x=loc_label, y = M)) + 
    geom_bar(aes(fill = M), stat="identity") +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
    scale_fill_gradient2(low= 'darkturquoise', high = 'darkred', mid = 'beige', midpoint = medianM) + 
    theme_bw() +
    facet_grid(.~age)+
    theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1), 
          legend.position = 'none', 
          panel.grid.major= element_blank(), 
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + 
    scale_y_continuous(expand = c(0,0), labels = scales::percent_format()) +
    coord_cartesian(ylim = c(0,max(tmp$CU)+0.001)) +
    labs(y = paste0('Predicted COVID-19 attributable mortality\nrates among individuals as of ', format(unique(tmp$date), '%b %Y')))
  ggsave(paste0(outdir, paste0('-MortalityRate_allages.png')), w = 8, h = 3.5)
  
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
          panel.grid.major= element_blank(), 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          plot.margin = unit(c(5.5,0,5.5,5.5), "pt")) + 
    scale_color_jcolors('pal8') + 
    scale_fill_jcolors('pal8') + 
    scale_y_continuous(expand = c(0,0), labels = scales::percent_format(), limits = range(c(tmp$CL, tmp$CU + 0.001)), 
                       breaks = seq(0, max(tmp$CU), 0.01)) +
    scale_x_continuous(expand = c(0,0)) +
    labs(y = paste0('Predicted COVID-19 attributable mortality\nrates among individuals as of ', format(unique(mortality_rate$date), '%b %Y')),
         col = '', fill = '', x = 'Age')
  
  p1 <- ggplot(subset(tmp, age == '85'), aes(x=age_cat)) + 
    geom_errorbar(aes(ymin=CL, ymax=CU), width = 0) + 
    geom_point(aes(y = M, col = loc_label)) +
    theme_bw() +
    theme(
          panel.grid.major= element_blank(), 
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
    labs(y = paste0('Predicted COVID-19 attributable mortality\nrates among individuals as of ', format(unique(mortality_rate$date), '%b %Y')),
         col = '')
  
  ggarrange(p, p1, nrow = 1, common.legend = T, legend = 'bottom', widths = c(1, 0.1))
  ggsave(paste0(outdir, paste0('-MortalityRateContinuous_allages_selectedstates.png')), w = 7, h = 5)
  
  
  
  ########
  
  p <- ggplot(subset(mortality_rate, age != 85), aes(x=age)) + 
    geom_line(aes(y = M, col = loc_label)) +
    geom_ribbon(aes(ymin=CL, ymax=CU, fill = loc_label), alpha = 0.4) + 
    theme_bw() +
    theme(legend.position = 'bottom', 
          panel.grid.major= element_blank(), 
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + 
    scale_color_jcolors('pal8') + 
    scale_fill_jcolors('pal8') + 
    facet_grid(loc_label~.) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent_format(), limits = range(c(mortality_rate$CL, mortality_rate$CU + 0.001)), 
                       breaks = seq(0, max(mortality_rate$CU), 0.01)) +
    scale_x_continuous(expand = c(0,0)) +
    labs(y = paste0('Predicted COVID-19 attributable mortality\nrates among individuals as of ', format(unique(mortality_rate$date), '%b %Y')),
         col = '', fill = '', x = 'Age')
  
  p1 <- ggplot(subset(mortality_rate, age == 85), aes(x=age_cat)) + 
    geom_errorbar(aes(ymin=CL, ymax=CU), width = 0) + 
    geom_point(aes(y = M, col = loc_label)) +
    theme_bw() +
    theme(
      panel.grid.major= element_blank(), 
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA), 
      axis.title.x = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_blank(), 
      plot.margin = unit(c(5.5,5.5,18,0), "pt"),
      legend.position = 'none',
    ) + 
    facet_grid(loc_label~.) +
    scale_color_jcolors('pal8') + 
    scale_fill_jcolors('pal8') + 
    scale_y_continuous(expand = c(0,0), labels = scales::percent_format(), limits = range(c(mortality_rate$CL, mortality_rate$CU + 0.001)), 
                       breaks = seq(0, max(mortality_rate$CU), 0.01)) +
    scale_x_discrete(expand = c(0,0)) +
    labs(y = paste0('Predicted COVID-19 attributable mortality\nrates among individuals as of ', format(unique(mortality_rate$date), '%b %Y')),
         col = '')
  
  ggarrange(p, p1, nrow = 1, common.legend = T, legend = 'bottom', widths = c(1, 0.1))
  ggsave(paste0(outdir, paste0('-MortalityRateContinuous_allages.png')), w = 7, h = 5 + 2*(length(unique(mortality_rate$code))/4))
  
}

plot_contribution_all_states <- function(contribution, vaccinedata, outdir){
  
  tmp = subset(vaccinedata, date %in% unique(vaccinationeffect$date))
  tmp1 = subset(vaccinationeffect, date == max(tmp$date))
  
  ggplot(tmp1, aes(y = loc_label)) + 
    geom_vline(xintercept = 0, linetype = 'dashed', col = 'grey50') +
    geom_errorbarh(aes(xmin = CL, xmax = CU, group =  as.factor(age)), 
                   position = position_dodge(width = 0.5), height = 0.2) +
    geom_point(aes(x = M, col = as.factor(age)), position = position_dodge(width = 0.5)) 
      
  
  
}


# plot_contribution_all_states = function(contribution, vaccinedata_state, outdir){
# 
#   vaccinedata_state = subset(vaccinedata_state, loc_label %in% unique(contribution$loc_label) & date <= max(contribution$date))
#   vaccinedata_state[, dummy := '']
#   
#   tmp = contribution[!is.na(emp)]
#   
#   ggplot(tmp, aes(x= date) ) +
#     geom_line(aes(y = M, col = age)) + 
#     geom_point(aes(y = emp, col = age), size = 0.5) + 
#     geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
#     facet_wrap(~loc_label, ncol = 6) +
#     geom_line(data = vaccinedata_state, aes(y  = prop_vaccinated_1dosep, linetype = dummy), col = 'grey20', size = 0.9) + 
#     scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
#     theme_bw() + 
#     theme(strip.background = element_blank(),
#           panel.border = element_rect(colour = "black", fill = NA), 
#           axis.text.x = element_text(angle = 45, hjust =1), 
#           panel.grid.major = element_blank(), 
#           axis.title.x = element_blank(), 
#           axis.title.y = element_text(size = rel(1.2)), 
#           legend.title = element_text(size = rel(1.1)), 
#           strip.text = element_text(size = rel(1.1)), 
#           legend.position = 'bottom') +
#     labs(y = 'Estimated age-specific contribution to COVID-19 weekly deaths', fill = 'Age groups', col = 'Age groups', 
#          linetype = 'Proportion of the state population\nvaccinated with at least one dose') + 
#     scale_y_continuous(labels = scales::percent_format()) +
#     scale_color_viridis_d(option = 'B', begin = 0.4, end = 0.8)+
#     scale_fill_viridis_d(option = 'B', begin = 0.4, end = 0.8) + 
#     scale_linetype_manual(values = 1)
#   ggsave(paste0(outdir, paste0('-Contribution_Vaccination_allStates.png')), w = 9, h = 12)
#   
#   
# }

plot_mortality_all_states = function(death, resurgence_dates, lab = 'allStates', outdir)
{
  
  df = as.data.table( reshape2::melt(select(death, loc_label, code, date, age, emp), id.vars = c('loc_label', 'code', 'date', 'age')) )
  df[, variable2 := 'CDC data']
  df[, `Age group` := age]
  
  death[, dummy := 'Posterior median prediction\nusing age-aggregated JHU data\nto adjust for reporting delays']
  # death[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  colfunc <- jcolors("pal8")[1:length(locs)]
  colfunc <- colfunc[which(locs %in% unique(death$code))]
  
  death[, `Age group` := age]
  
  dummy.dt = merge(resurgence_dates, unique(select(death, code, loc_label)), by = 'code')
  dummy.dt[, text := 'Beginning of Summer 2021 resurgence period']
  
  p <- ggplot(subset(death), aes(x= date) ) +
    geom_point(data = subset(df), aes(y = value, shape= variable2), col = 'black', fill = 'black', size = 0.9, alpha = 0.7) + 
    geom_line(aes(y = M, fill = loc_label), show.legend = F) +
    geom_vline(data = dummy.dt, aes(xintercept = start_resurgence, linetype = text), col = 'grey50') +
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = loc_label), alpha = 0.5, show.legend = F) + 
    facet_grid(loc_label~`Age group`, scale = 'free') +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    scale_y_continuous(expand = c(0,0)) + 
    theme_bw() + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 45, hjust =1), 
          panel.grid.major = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          legend.title = element_text(size = rel(1.1)), 
          legend.text = element_text(size = rel(1.1)), 
          strip.text = element_text(size = rel(1)), 
          legend.position = 'bottom') +
    labs( y = 'Predicted age-specific COVID-19 attributable weekly deaths', shape = '',
          linetype = '', col = '', fill = '') + 
    scale_shape_manual(values = 16) +
    scale_linetype_manual(values = 2) + 
    scale_color_manual(values = as.character(colfunc)) + 
    scale_fill_manual(values = as.character(colfunc)) + 
    guides(shape = guide_legend(override.aes = list(size=1, stroke =1.5), order = 2),
           linetype = guide_legend(order = 3), 
           col = guide_legend(order = 1),
           fill = guide_legend(order = 1))
  ggsave(p, file = paste0(outdir, paste0('-Mortality_', lab, '.png')), w = 7.5, h = 8)
  
}

plot_vaccine_effects_counterfactual <- function(data_res1, data_res2, resurgence_dates, var, lab, outdir){
  
  label_fit <- 'Fit to observed data'
  data_res2[, label_counterfactual := label_fit]
  
  data_res1 = merge(data_res1, select(resurgence_dates, code, start_resurgence), by = 'code')
  data_res1 = data_res1[date >= start_resurgence ]
  data_res1 = select(data_res1, -start_resurgence, -counterfactual_index)
  
  tmp = rbind(data_res1, data_res2)
  
  tmp = subset(tmp, date >= as.Date('2021-01-01'))
  tmp = merge(tmp, select(resurgence_dates, code, stop_resurgence), by = 'code')
  tmp = tmp[date <= stop_resurgence]
  tmp[, label_counterfactual := factor(label_counterfactual, levels = c(label_fit, df_counterfactual$label_counterfactual))]

  dummy.dt = merge(resurgence_dates, df_state, by = 'code')
  dummy.dt[, text := 'Beginning of Summer 2021 resurgence period']
  
  # values_col = c('grey50', 'darkorchid4', )
  
  cols <- viridisLite::viridis(length(unique(tmp$label_counterfactual)), direction = -1, begin = 0.1)
  ages = unique(tmp$age)
  
  for(j in df_counterfactual$counterfactual_index){
    tmp1 <- subset(tmp, label_counterfactual %in% c(label_fit, df_counterfactual[counterfactual_index == j, label_counterfactual]))
    
    p = list()
    for(i in 1:length(ages)){
      # i = 1
      Age = ages[i]
      
      p[[i]] = ggplot(subset(tmp1, age == Age), aes(x = date)) + 
        geom_line(aes(y = M, col = label_counterfactual)) + 
        geom_ribbon(aes(ymin = CL, ymax = CU, fill = label_counterfactual), alpha = 0.5) + 
        facet_grid(loc_label~age) + 
        scale_x_date(expand = expansion(mult = c(0.05,0)), date_labels = c("%b-%y"), breaks = '1 month') + 
        theme_bw() + 
        geom_vline(data = dummy.dt, aes(xintercept = start_resurgence, linetype = text), col = 'grey50') +
        theme(strip.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA), 
              legend.position = 'bottom', 
              axis.text.x = element_text(angle = 70, hjust = 1),
              axis.title.x = element_blank()) + 
        # scale_color_manual(values = values_col) +
        # scale_fill_manual(values = values_col) + 
        scale_color_manual(values = c(cols[1], cols[j + 1])) + 
        scale_fill_manual(values = c(cols[1], cols[j + 1])) + 
        labs(col = '', y = paste0('Predicted age-specific ',var,' COVID-19 attributable weekly deaths'),
             fill = '', linetype = '') +
        scale_linetype_manual(values = 2) +
        guides(fill=guide_legend(nrow=2,byrow=TRUE, order =1), col=guide_legend(nrow=2,byrow=TRUE, order =1), 
               linetype = guide_legend(order=2))
      
      if(i != length(ages)){
        p[[i]] = p[[i]] + theme(strip.text.y = element_blank())
      }
      
      if(i != 1){
        p[[i]] = p[[i]] + theme(axis.title.y = element_blank())
      }
    }
    
    p = ggarrange(plotlist = p, common.legend = T, legend = 'bottom', nrow = 1, widths = c(1.1, 1.1))
    ggsave(p, file = paste0(outdir, '-predicted_weekly_deaths_vaccine_coverage_counterfactual', j, '_', lab, '.png'), w = 6 + length(unique(data_res1$code))/6, h = 5 + 2*(length(unique(data_res1$code))/4))
    
  }
  
  tmp1 <- tmp[, list(max_date = max(date)), by = 'code']
  tmp1 <- merge(tmp, tmp1, by = 'code')
  tmp1 <- tmp1[date == max_date]
  
  tmp2 <- subset(tmp1, label_counterfactual == label_fit)
  setnames(tmp2, c('M', 'CL', 'CU'), c('M_fit', 'CL_fit', "CU_fit"))
  tmp1 <- merge(tmp1, tmp2[, .(M_fit, CL_fit, CU_fit, code, age_index)], by = c('code', 'age_index'))

# tmp1[, CU := CU + M*2]
# tmp1[, CL := CL - M*2]
  p = list()
  for(i in 1:length(ages)){
    # i = 1
    Age = ages[i]
    
    
    p[[i]] =  ggplot(subset(tmp1, age == Age), aes(x = label_counterfactual)) + 
      geom_hline(aes(yintercept=M_fit), col = cols[1], linetype = 'dashed') +
      geom_rect( aes(ymin = CL_fit, ymax = CU_fit), xmin = -Inf, xmax = Inf, fill = cols[1], alpha = 0.05) +
      geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.9, width = 0, col = 'grey40') + 
      geom_point(aes(y = M, col = label_counterfactual)) + 
      facet_grid(loc_label~age) +
      scale_color_manual(values = cols) + 
      scale_fill_manual(values = cols) + 
      theme_bw() +
      theme(strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA), 
            legend.position = 'bottom', 
            axis.text.x = element_blank(),
            axis.title.x = element_blank()) + 
      labs(col = '', y = paste0('Predicted age-specific ',var,' COVID-19 attributable weekly deaths\nat the end of the resurgence period'),
           fill = '', linetype = '') +
      guides(fill=guide_legend(nrow=1 + stan_data$N_COUNTERFACTUAL,byrow=TRUE, order =1), 
             col=guide_legend(nrow=1 + stan_data$N_COUNTERFACTUAL,byrow=TRUE, order =1), 
             linetype = guide_legend(order=2)) 
    
    if(i != length(ages)){
      p[[i]] = p[[i]] + theme(strip.text.y = element_blank())
    }
    
    if(i != 1){
      p[[i]] = p[[i]] + theme(axis.title.y = element_blank())
    }

  }
  
  p = ggarrange(plotlist = p, common.legend = T, legend = 'bottom', nrow = 1, widths = c(1.1, 1.05))
  ggsave(p, file = paste0(outdir, '-predicted_weekly_deaths_vaccine_coverage_', lab, '.png'), w = 5, h = 5 + 2*(length(unique(data_res1$code))/4))
  
}
  
plot_forest_plot <- function(tmp, outdir){
  
  tmp1 = subset(tmp, grepl('\\["18-64', variable) | (grepl('18-64"\\]', variable) & !grepl('\\["65', variable)))
  p1 <- ggplot(tmp1, aes(y = variable)) + 
    geom_point(aes(x = M)) + 
    geom_errorbarh(aes(xmin = CL, xmax = CU), height = 0.2) + 
    theme_bw() + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    facet_grid(group~., scales= "free", space="free") +
    scale_y_discrete(labels = label_parse(), limits=rev) + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          title = element_text(size = rel(0.8)),
          plot.title = element_text(hjust = 0.5)) +
  ggtitle("Effects on the relative rates\namong 18-64")
  
  tmp1 = subset(tmp, grepl('\\["65', variable) | (grepl('65\\+"\\]', variable) & !grepl('\\["18-64', variable)))
  p2 <- ggplot(tmp1, aes(y = variable)) + 
    geom_point(aes(x = M)) + 
    geom_errorbarh(aes(xmin = CL, xmax = CU), height = 0.2) + 
    theme_bw() + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    facet_grid(group~., scales= "free", space="free") +
    scale_y_discrete(labels = label_parse(), limits=rev) + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          title = element_text(size = rel(0.8)),
          plot.title = element_text(hjust = 0.5))+
    ggtitle("Effects on the relative rates\namong 65+")
  
  p <- ggarrange(p1, p2, ncol = 2)
  ggsave(p, file = paste0(outdir, '-forest_plot.png'), w = 8, h = 7)
  
  return(p)
}

plot_forest_plot_with_common_effect <- function(forest_plot_without_common_effect, tmp, outdir){
  
  p1 <- ggplot(tmp, aes(y = variable)) + 
    geom_point(aes(x = M)) + 
    geom_errorbarh(aes(xmin = CL, xmax = CU), height = 0.2) + 
    theme_bw() + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    scale_y_discrete(labels = label_parse(), limits=rev) + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          title = element_text(size = rel(0.8)),
          plot.title = element_text(hjust = 0.5))+
    ggtitle("Common vaccine effects on the relative\nrates across 18-64 and 65+")
  
  p2 <- grid.arrange(forest_plot_without_common_effect, p1, 
                     layout_matrix = rbind(c(1, 1, 1),
                                           c(NA, 2, NA)),
                     heights = c(1, 0.2))
  ggsave(p2, file = paste0(outdir, '-forest_plot.png'), w = 8, h = 8)
  
  
}

plot_relative_resurgence_vaccine <- function(data_res1, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, outdir){
  
  data_res = merge(data_res1, prop_vac, by = c('code', 'date'))
  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  data_res = merge(data_res, resurgence_dates, by = 'code')
  
  p1 <- ggplot(data_res, aes(x = prop_1)) + 
    geom_line(aes(y = M, col = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = loc_label), alpha = 0.5) +
    facet_grid(`Age group`~., label = 'label_both') +
    geom_point(data = subset(data_res, date == start_resurgence), aes(shape = '', y = M), size = 2, stroke = 1, col = 'grey50')  + 
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = 'Proportion of individuals aged 18-64\nfully vaccinated two weeks before', shape = 'Beginning of Summer 2021 resurgence period', col = '', fill = '') + 
    theme_bw() +
    scale_x_continuous(labels = scales::percent) + 
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1)),
          legend.spacing.y = unit(-0, "cm")) +
    scale_shape_manual(values = 4) + 
      scale_color_jcolors(palette = "pal8") + 
      scale_fill_jcolors(palette = "pal8") +
    guides(color = guide_legend(order=1), fill = guide_legend(order=1), shape = guide_legend(order=2)) 

  p2 = ggplot(data_res, aes(x = prop_2)) + 
    geom_line(aes(y = M, col = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = loc_label), alpha = 0.5) +
    facet_grid(`Age group`~., label = 'label_both') +
    geom_point(data = subset(data_res, date == start_resurgence), aes(shape = '', y = M), size = 2, stroke = 1, col = 'grey50')     + 
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = 'Proportion of individuals aged 65+\nfully vaccinated two weeks before', shape = 'Beginning of Summer 2021 resurgence period', col = '', fill = '') + 
    theme_bw() +
    scale_x_continuous(labels = scales::percent) + 
    theme(strip.background = element_blank(),
          panel.border =  element_rect(colour = "black", fill = NA), 
          strip.text = element_text(size = rel(0.85)),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = rel(0.9)),
          legend.box="vertical") + 
    scale_shape_manual(values = 4)+ 
      scale_color_jcolors(palette = "pal8") + 
      scale_fill_jcolors(palette = "pal8") +
    guides(color = guide_legend(order=1), fill = guide_legend(order=1), shape = guide_legend(order=2)) 
  
  p = ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom', widths = c(1.1,1)) # legend.grob = get_legend(p0)
  ggsave(p, file = paste0(outdir, '-relative_deaths_vaccine_coverage.png'), w = 6.5, h = 5)
  
  
}



plot_relative_resurgence_vaccine_indicator <- function(data_res1, prop_vac_indicator, df_age_vaccination2, df_week2, outdir){
  
  data_res = merge(data_res1, select(prop_vac_indicator, -date), by = c('code'))
  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  prop_vac_indicator[indic1 & indic2, group := '1']
  prop_vac_indicator[indic1 & !indic2, group := '2']
  prop_vac_indicator[!indic1 & indic2, group := '3']
  prop_vac_indicator[!indic1 & !indic2, group := '4']
  
  lab = function(Age, prop) paste0('Proportion of individuals aged ', Age, ' fully vaccinated\ntwo weeks before the beginning of Summer 2021\nresurgence period greater than ', prop, '%')
  
  p1 <- ggplot(data_res, aes(x = week_index)) + 
    geom_line(aes(y = M, col = indic1, linetype = indic2, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = indic1, group = loc_label), alpha = 0.5, width = 0.1) +
    facet_grid(`Age group`~., label = 'label_both') +
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = '', shape = 'Beginning of Summer 2021 resurgences', 
         col = lab('18-64', cutoff_1864*100), fill = lab('18-64', cutoff_1864*100), 
         linetype = lab('64+', cutoff_65p*100)) + 
    theme_bw() +
    geom_text(data = subset(data_res, week_index == max(week_index)), aes(label = loc_label, x = week_index, y = M), hjust = -.1, size = 3) +
    # scale_x_date(breaks = '2 weeks', expand=  expansion(mult = c(0,0.15)), date_labels = "%b-%y") + 
    scale_x_continuous(expand=  expansion(mult = c(0,0.15))) + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1)),
          strip.text = element_text(size = rel(0.9)),
          legend.spacing.y = unit(-0, "cm"), 
          legend.position = 'bottom') +
    scale_color_jcolors(palette = "pal8") + 
    scale_fill_jcolors(palette = "pal8") +
    guides(color = guide_legend(order=1), fill = guide_legend(order=1), linetype = guide_legend(order=2)) 

  ggsave(p1, file = paste0(outdir, '-relative_deaths_vaccine_coverage_indicator.png'), w = 6.5, h = 7)
  
}

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

plot_slope_by_prop <- function(log_r_pdeaths, prop_vac, outdir){
  tmp <- log_r_pdeaths[, {
    fit = lm(M ~ week_index)
    list(beta = fit$coefficient[2])
  }, by = c('code', 'age')]
  
  tmp1 = prop_vac[, list(prop_1 = prop_1[date == min(date)], prop_2 = prop_2[date == min(date)]), by = 'code']
  tmp = merge(tmp,tmp1, by = 'code')
  
  p1 <- ggplot(tmp, aes(x = prop_1, y = beta)) + 
    geom_point() + 
    facet_wrap(~age, nrow = 2) + 
    geom_smooth(method = 'lm')
  
  p2 <- ggplot(tmp, aes(x = prop_2, y = beta)) + 
    geom_point() + 
    facet_wrap(~age, nrow = 2) + 
    geom_smooth(method = 'lm')
  
  p = ggpubr::ggarrange(p1, p2, ncol = 2)
  ggsave(p, file = paste0(outdir, '-relative_deaths_slope_to_prop.png'), w = 8, h = 8)
  
}


plot_relative_resurgence_vaccine2 <- function(data_res1, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, log_transform, outdir, lab.fig = ''){

  prop_vac_init = prop_vac[, list(prop_1_init = prop_1[date == min(date)], prop_2_init = prop_2[date == min(date)]), by = 'code']
  data_res = merge(data_res1, prop_vac_init, by = 'code')

  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]

  data_res = as.data.table( merge(data_res, resurgence_dates, by = 'code') )
  
  data_res[, max_week_index := max(week_index), by = 'code']
  data_text = data_res[ week_index == max_week_index]
  data_text[, week_index := max(week_index)]
  data_text[code == 'CA' & age == '18-64', M := M + 0.1]
  data_text[code == 'CA' & age == '65+', M := M + 0.05]
  
  lab = function(Age) paste0('Proportion of individuals aged ', Age, '\nfully vaccinated two weeks before\nthe beginning of Summer 2021\nresurgence period')
  
  p1 <- ggplot(data_res, aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_1_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_1_init, group = loc_label), alpha = 0.5) +
    facet_grid(`Age group`~., label = 'label_both') +
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = 'Week index of Summer 2021 resurgences period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('18-64'), fill = lab('18-64')) + 
    theme_bw() +
    geom_text(data = data_text, aes(label = loc_label, x = week_index, y = M, col = prop_1_init), hjust = -.1, size = 3) +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    scale_x_continuous(expand=  expansion(mult = c(0,0.25)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1)),
          strip.text = element_blank(),
          legend.spacing.x = unit(0.3, "cm"), 
          legend.position = 'bottom') +
    scale_color_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init))) + 
    scale_fill_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)))
  
  p2 <- ggplot(data_res, aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_2_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_2_init, group = loc_label), alpha = 0.5, width = 0.1) +
    facet_grid(`Age group`~., label = 'label_both') +
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = 'Week index of Summer 2021 resurgence period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('65+'), fill = lab('65+')) + 
    theme_bw() +
    geom_text(data = data_text, aes(label = loc_label, x = week_index, y = M, col = prop_2_init), hjust = -.1, size = 3) +
    # scale_x_date(breaks = '2 weeks', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") +
    scale_x_continuous(expand=  expansion(mult = c(0,0.25)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = rel(0.9)),
          legend.spacing.x = unit(0.3, "cm"), 
          legend.position = 'bottom') +
    scale_color_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init))) + 
    scale_fill_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)))
  
  if(log_transform){
    p1 = p1 + scale_y_continuous(trans = 'log', breaks = base_breaks()) + labs(y = 'log relative COVID-19 attributable weekly deaths')
    p2 = p2 + scale_y_continuous(trans = 'log', breaks = base_breaks()) 
  }
  
  p = grid.arrange(p1, p2, ncol = 2)
  
  if(log_transform){
    file =  paste0(outdir, '-log_relative_deaths_vaccine_coverage2', lab.fig, '.png')
  } else{
    file =  paste0(outdir, '-relative_deaths_vaccine_coverage2', lab.fig, '.png')
  }
  
  ggsave(p, file =file, w = 8, h = 5)
  
}

plot_relative_resurgence_vaccine2_long <- function(data_res1, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, outdir){
  
  prop_vac_init = prop_vac[, list(prop_1_init = prop_1[date == min(date)], prop_2_init = prop_2[date == min(date)]), by = 'code']
  data_res = merge(data_res1, prop_vac_init, by = 'code')
  
  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  data_res = merge(data_res, resurgence_dates, by = 'code')
  
  data_text = subset(data_res, week_index == max(week_index))
  data_text[code == 'CA' & age == '18-64', M := M + 0.1]
  data_text[code == 'CA' & age == '65+', M := M + 0.05]
  
  lab = function(Age) paste0('Proportion of individuals aged ', Age, ' fully vaccinated two weeks\nbefore the beginning of Summer 2021 resurgence period')
  
  ## 18-64
  mid_code = round(length(Code) / 2)
  Code_ordered = prop_vac_init[order(prop_1_init)]$code
  tmp <- unique(select(data_res, code, loc_label))
  tmp[, code := factor(code, levels = Code_ordered)]
  
  data_res = data_res[, loc_label := factor(loc_label, levels = tmp[order(code)]$loc_label)]
  p1 <- ggplot(subset(data_res, code %in% Code_ordered[1:mid_code]), aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_1_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_1_init, group = loc_label), alpha = 0.5) +
    facet_grid(`Age group`~`loc_label`,  scale = 'free_y') +
    labs(y = '', 
         x = 'Week index of Summer 2021 resurgences period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('18-64'), fill = lab('18-64')) + 
    theme_bw() +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    # scale_x_continuous(expand=  expansion(mult = c(0,0.25)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_blank(),
          legend.spacing.x = unit(0.5, "cm"), 
          legend.position = 'bottom',
          legend.key.width = unit(1 , "cm")) +
    scale_y_continuous( limits =range(c(data_res$CL, data_res$CU))) +
    scale_x_continuous(breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    scale_color_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)),
                          limits = range(prop_vac_init$prop_1_init)) + 
    scale_fill_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)),
                         limits = range(prop_vac_init$prop_1_init))
  
  p2 <- ggplot(subset(data_res, code %in% Code_ordered[(mid_code +1):(mid_code*2)]), aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_1_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_1_init, group = loc_label), alpha = 0.5) +
    facet_grid(`Age group`~`loc_label`,  scale = 'free_y') +
    labs(y = '', 
         x = 'Week index of Summer 2021 resurgences period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('18-64'), fill = lab('18-64')) + 
    theme_bw() +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    # scale_x_continuous(expand=  expansion(mult = c(0,0.25)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    scale_x_continuous( breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    scale_y_continuous( limits =range(c(data_res$CL, data_res$CU))) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1)),
          legend.spacing.x = unit(0.5, "cm"), 
          legend.position = 'bottom') +
    scale_color_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)), 
                          limits = range(prop_vac_init$prop_1_init)) + 
    scale_fill_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)), 
                         limits = range(prop_vac_init$prop_1_init))
  
  p <- ggarrange(p1, p2, ncol = 1, common.legend = T, legend = 'bottom', heights = c(0.9, 1))
  p <- grid.arrange(p, left = 'Relative COVID-19 attributable weekly deaths')
  file =  paste0(outdir, '-relative_deaths_vaccine_coverage2', '_extended_part1', '.png')
  ggsave(p, file =file, w = 8, h = 6)
  
  ## 65+
  mid_code = round(length(Code) / 2)
  Code_ordered = prop_vac_init[order(prop_2_init)]$code
  tmp <- unique(select(data_res, code, loc_label))
  tmp[, code := factor(code, levels = Code_ordered)]
  
  data_res = data_res[, loc_label := factor(loc_label, levels = tmp[order(code)]$loc_label)]
  p1 <- ggplot(subset(data_res, code %in% Code_ordered[1:mid_code]), aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_2_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_2_init, group = loc_label), alpha = 0.5) +
    facet_grid(`Age group`~`loc_label`,  scale = 'free_y') +
    labs(y = '', 
         x = 'Week index of Summer 2021 resurgences period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('65+'), fill = lab('65+')) + 
    theme_bw() +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    # scale_x_continuous(expand=  expansion(mult = c(0,0.25)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(1)),
          legend.spacing.x = unit(0.5, "cm"), 
          legend.position = 'bottom',
          legend.key.width = unit(1 , "cm")) +
    scale_y_continuous( limits =range(c(data_res$CL, data_res$CU))) +
    scale_x_continuous(breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    scale_color_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)),
                          limits = range(prop_vac_init$prop_2_init)) + 
    scale_fill_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)),
                         limits = range(prop_vac_init$prop_2_init))
  
  p2 <- ggplot(subset(data_res, code %in% Code_ordered[(mid_code +1):(mid_code*2)]), aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_2_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_2_init, group = loc_label), alpha = 0.5) +
    facet_grid(`Age group`~`loc_label`,  scale = 'free_y') +
    labs(y = '', 
         x = 'Week index of Summer 2021 resurgences period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('64+'), fill = lab('65+')) + 
    theme_bw() +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    # scale_x_continuous(expand=  expansion(mult = c(0,0.25)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    scale_x_continuous( breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    scale_y_continuous( limits =range(c(data_res$CL, data_res$CU))) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1)),
          legend.spacing.x = unit(0.5, "cm"), 
          legend.position = 'bottom') +
    scale_color_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)), 
                          limits = range(prop_vac_init$prop_2_init)) + 
    scale_fill_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)), 
                         limits = range(prop_vac_init$prop_2_init))
  
  p <- ggarrange(p1, p2, ncol = 1, common.legend = T, legend = 'bottom', heights = c(0.9, 1))
  p <- grid.arrange(p, left = 'Relative COVID-19 attributable weekly deaths')
  file =  paste0(outdir, '-relative_deaths_vaccine_coverage2', '_extended_part2', '.png')
  ggsave(p, file =file, w = 8, h = 6)
  
}


plot_relative_resurgence_vaccine_no_time <- function(data_res1, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, log_transform, outdir){
  
  prop_vac_init = prop_vac[, list(prop_1_init = prop_1[date == min(date)], prop_2_init = prop_2[date == min(date)]), by = 'code']
  data_res = merge(data_res1, prop_vac_init, by = 'code')
  
  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  data_res = merge(data_res, resurgence_dates, by = 'code')
  
  lab = function(Age) paste0('Proportion of individuals aged ', Age, '\nfully vaccinated two weeks before\nthe beginning of Summer 202\nresurgences')
  
  p1 <- ggplot(data_res, aes(x = prop_1_init)) + 
    geom_point(aes(y = M, col = week_index, shape = loc_label, group = loc_label)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = loc_label), alpha = 0.5) +
    facet_grid(`Age group`~., label = 'label_both') +
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = '', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('18-64'), fill = lab('18-64')) + 
    theme_bw() +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    scale_x_continuous(expand=  expansion(mult = c(0,0.25))) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1)),
          strip.text = element_blank(),
          legend.spacing.y = unit(-0, "cm"), 
          legend.position = 'bottom') 
  
  p2 <- ggplot(data_res, aes(x = prop_2_init)) + 
    geom_point(aes(y = M, col = week_index, shape = loc_label, group = loc_label)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = loc_label), alpha = 0.5) +
    facet_grid(`Age group`~., label = 'label_both') +
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = '', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('65+'), fill = lab('65+')) + 
    theme_bw() +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    scale_x_continuous(expand=  expansion(mult = c(0,0.25))) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1)),
          strip.text = element_blank(),
          legend.spacing.y = unit(-0, "cm"), 
          legend.position = 'bottom') 
  
  if(log_transform){
    p1 = p1 + scale_y_continuous(trans = 'log', breaks = base_breaks()) + labs(y = 'log relative COVID-19 attributable weekly deaths')
    p2 = p2 + scale_y_continuous(trans = 'log', breaks = base_breaks()) 
  }
  
  p = grid.arrange(p1, p2, ncol = 2)
  
  if(log_transform){
    file =  paste0(outdir, '-log_relative_deaths_vaccine_coverage_no_time.png')
  } else{
    file =  paste0(outdir, '-relative_deaths_vaccine_coverage_no_time.png')
  }
  
  ggsave(p, file =file, w = 8, h = 5)
  
}


plot_vaccination_effect_prediction <- function(prediction, var, outdir){
  
  p = list()
  for(i in 1:length(unique(prediction$age_recipient))){
    # Age="18-64"
    Age = unique(prediction$age_recipient)[i]
    tmp <- subset(prediction, age_source == Age)
    p[[i]] = ggplot(tmp, aes(x = prop)) + 
      geom_line(aes(y = M, col = age_recipient)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = age_recipient), alpha = 0.3) +
      labs(y = 'Mitigation on the magnitude of\nrelative COVID-19 attributable weekly deaths', 
           x = paste0('Proportion of individuals aged ', Age, '\nfully vaccinated two weeks before the\nbeginning of Summer 2021 resurgence period'), 
           col = 'Age group', fill = 'Age group') + 
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
            legend.title = element_text(size = rel(0.85)),
            axis.title.x = element_text(size = rel(0.9)),
            axis.title.y = element_text(size = rel(1)),
            legend.spacing.y = unit(-0, "cm")) +
      scale_color_jcolors(palette = "pal6") + 
      scale_fill_jcolors(palette = "pal6") 
    
    if(i > 1){
      p[[i]] <- p[[i]] + theme(axis.title.y = element_blank())
    }
  }
  
  
  p = ggarrange(plotlist = p, ncol = 2, common.legend = T, legend = 'bottom', widths = c(1.1,1)) # legend.grob = get_legend(p0)
  ggsave(p, file = paste0(outdir, '-vaccine_effect_prediction_', var, '.png', collapse= '_'), w = 6.5, h = 5)
  
  return(p)
}


plot_PPC_relative_resurgence <- function(data_res1, data_res2, lab, outdir){
  
  data_res1[, type := 'Fit to observed data']
  data_res2[, type := 'Predicted with vaccination rates']
  data_res = rbind(data_res1, data_res2)
  
  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  p1 <- ggplot(data_res, aes(x = date)) + 
    geom_line(aes(y = M, col = type)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = type), alpha = 0.5) +
    facet_grid(`Age group`~loc_label) +
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         shape = 'Beginning of Summer 2021 resurgence period', col = '', fill = '') + 
    theme_bw() +
    scale_x_date(expand = c(0,0), breaks = '1 month', date_labels = "%b-%y") + 
    theme(strip.background = element_blank(),
          strip.text = element_text(size = rel(0.85)),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.y = element_text(size = rel(0.9)),
          axis.title.x = element_blank(),
          legend.spacing.y = unit(-0, "cm"), 
          legend.position = 'bottom',
          axis.text.x = element_text(angle = 70, hjust = 1)) +
    guides(color = guide_legend(order=1), fill = guide_legend(order=1)) 
  
  ggsave(p1, file = paste0(outdir, '-relative_deaths_vaccine_coverage_PPC', lab, '.png'), w = 8, h = 5)
  
}

plot_contribution_vaccine <- function(contribution, vaccine_data, resurgence_dates, lab, outdir){
  
  delay = 7*2
  df_age_vaccination = unique(select(contribution, age_index, age))
  df_age_vaccination[, age_from := gsub('(.+)-.*', '\\1', age)]
  df_age_vaccination[, age_to := gsub('.*-(.+)', '\\1', age)]
  df_age_vaccination[grepl('\\+', age_from), age_from := gsub('(.+)\\+', '\\1', age)]
  df_age_vaccination[grepl('\\+', age_to), age_to := max(vaccine_data$age)]
  set(df_age_vaccination, NULL, 'age_from', df_age_vaccination[,as.numeric(age_from)])
  set(df_age_vaccination, NULL, 'age_to', df_age_vaccination[,as.numeric(age_to)])
  vaccine_data[, age_index := which(df_age_vaccination$age_from <= age & df_age_vaccination$age_to >= age), by = 'age']
  vaccine_data = vaccine_data[!is.na(age_index), list(prop = unique(prop)), by = c('code', 'date', 'loc_label', 'age_index')]
  vaccine_data[, date := date + delay]
  
  tmp <- merge(contribution, vaccine_data, by = c('code', 'date', 'loc_label', 'age_index'), all.x = T)
  tmp <- merge(tmp, resurgence_dates, by = 'code')
    
  p1 = ggplot(tmp, aes(x = date)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
    geom_line(aes(y = M, col = as.factor(age))) + 
    geom_point(data = subset(tmp, date == start_resurgence), aes(shape = '', y = M), size = 2, stroke = 1.1) + 
    facet_wrap(~loc_label, nrow = length(unique(tmp$loc_label))) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    labs(y = paste0("Estimated contribution to COVID-19 weekly deaths"), 
         col = "Age group", fill = "Age group", x = 'Date',
         shape = 'Beginning of Summer 2021 resurgence period') + 
    scale_shape_manual(values = 4)+ 
    scale_x_date(expand = c(0,0), breaks = '3 months', date_labels = "%b-%y") 
  

  p2 = ggplot(tmp, aes(x = prop)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
    geom_line(aes(y = M, col = as.factor(age))) + 
    geom_point(data = subset(tmp, date == start_resurgence), aes(shape = '', y = M), size = 2, stroke = 1.1) + 
    facet_wrap(~loc_label, nrow = length(unique(tmp$loc_label))) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
    labs(y = paste0(""), 
         x = 'Proportion of fully vaccinated individuals',
         col = "Age group", fill = "Age group",
         shape = 'Beginning of Summer 2021 resurgence period') + 
    scale_shape_manual(values = 4) 
  
  p = ggarrange(p1, p2,common.legend = T, legend = 'bottom')
  ggsave(p, file = paste0(outdir, '-contribution_vaccine_coverage_', lab, '.png'), w = 7, h = 8)
  
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

plot_vaccine_effects <- function(vaccine_data, weeklydv, weeklyf, weeklyphi, outdir){
  
  vaccine_data[, age_index := which(df_age_vaccination$age_from <= age & df_age_vaccination$age_to >= age), by = 'age']
  tmp <- unique(select(vaccine_data, code, date, prop, age_index))
  
  delay = 7*2
  weeklydv[, date := date - delay]
  weeklyphi[, date := date - delay]
  
  tmp1 <- merge(tmp, weeklydv, by = c('code', 'date', 'age_index'))
  p <- ggplot(tmp1, aes(x = prop)) + 
    geom_point(aes(y = M)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU)) + 
    facet_wrap(~age, nrow = nrow(df_age_vaccination)) + 
    labs(x = 'Proportion of vaccinated', y = 'Weekly deaths 2 weeks later') + 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white")) 
  ggsave(p, file = paste0(outdir, "-vaccination_effects_weekly_deaths_", Code,".png") , w= 7, h = 8, limitsize = FALSE)
  
  p <- ggplot(weeklyf, aes(x = date)) + 
    geom_point(aes(y = M)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU)) + 
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

plot_vaccine_effects2 <- function(vaccine_effects, phi, phi_wo_vaccine, lab, outdir){
  
  isunscaled = any(phi$M < 0 | phi$M > 1)
  isunscaled = ifelse(isunscaled, 'unscaled ', '')
    
  # tmp1 <- subset(vaccine_effects, age_index %in% rev(unique(sort(vaccine_effects$age_index)))[1:3])
  tmp1 <- subset(vaccine_effects, !age_index %in% unique(sort(vaccine_effects$age_index))[1])
  tmp1 <- subset(tmp1, date >= as.Date('2021-01-01'))
  p <- ggplot(tmp1, aes(x = date)) + 
    geom_line(aes(y = M , col = age)) + 
    geom_ribbon(aes(ymin = CL , ymax = CU , fill = age), alpha = 0.25) +
    theme_bw() +
    labs(y = paste0("Vaccine effect \n (% change ", isunscaled, "contribution to weekly deaths)"), 
         x = "", col = 'Age groups', fill = 'Age groups') + 
    scale_color_viridis_d(option = 'C', end = 0.9)+ 
    scale_fill_viridis_d(option = 'C', end = 0.9) + 
    scale_y_continuous(labels = scales::percent) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) 
  ggsave(p, file = paste0(outdir, "-vaccine_effects_perc_contribution_", lab, '_', Code,".png") , w= 6, h = 5, limitsize = FALSE)
  
  tmp <- reshape2::melt(phi, id.vars = c('age_index', 'week_index', 'code', 'date', 'age'))
  setnames(tmp, c('value', 'variable'), c('with vaccination', 'stat'))
  tmp1 <- reshape2::melt(phi_wo_vaccine, id.vars = c('age_index', 'week_index', 'code', 'date', 'age'))
  setnames(tmp1, c('value', 'variable'), c('without vaccination', 'stat'))
  tmp <- merge(tmp, tmp1, by = c('age_index', 'week_index', 'code', 'date', 'age', 'stat'))
  tmp <- reshape2::melt(tmp, id.vars = c('age_index', 'week_index', 'code', 'date', 'age', 'stat'))
  tmp <- subset(tmp, !age_index %in% unique(sort(tmp$age_index))[1])
  tmp <- reshape2::dcast(tmp, date + age + variable ~ stat, value.var = 'value')
  
  p <- ggplot(tmp, aes(x = date)) + 
    geom_line(aes(y = M, col = variable)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = variable), alpha = 0.5) +
    facet_wrap(~age, nrow = length(unique(tmp$age))) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    labs(y = paste0("Estimated ", isunscaled, "contribution to COVID-19 weekly deaths"), x = "", col = '', fill = '') 
  ggsave(p, file = paste0(outdir, "-vaccine_effects_compcont_", lab, '_', Code,".png") , w= 5, h = 7, limitsize = FALSE)
  
}

plot_lambda_table <- function(lambda_table, outdir){
  
  state_indices = unique(lambda_table$state_index)
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

plot_var_base_model_table <- function(loc_label, outdir){
  
  p <- ggplot(var_base_model_table, aes(x = loc_label, col = type)) + 
    geom_point(aes(y = M), position = position_dodge(0.5), size = 0.75) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), position = position_dodge(0.5), width = 0.1) + 
    facet_grid(math_name~., labeller = label_parsed, scales = 'free_y') + 
    theme_bw() +
    labs(y = '', col = '') +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 70, hjust =1), 
          strip.background = element_blank(), 
          axis.title.x = element_blank(), axis.title.y = element_blank()) 
  scale_y_discrete(labels = label_parse(), limits=rev) 
  ggsave(p, file = paste0(outdir, '-var_base_model_prior_posterior.png'), w = 6, h = 6)
  
}


plot_vaccine_effects_counterfactual_perc <- function(data_res, prop_vac_counterfactual, lab, outdir){
  
  prop_vac_counterfactual_df <- copy(prop_vac_counterfactual)
  setnames(prop_vac_counterfactual_df, 'age_index', 'age_index_counterfactual')
  df <- copy(df_age_vaccination2[, .(age, age_index)])
  setnames(df, 1:2, c('age_counterfactual', 'age_index_counterfactual'))
  prop_vac_counterfactual_df <- merge(prop_vac_counterfactual_df, df, by = 'age_index_counterfactual')
  
  tmp1 <- data_res[, list(max_date = max(date)), by = 'code']
  tmp1 <- merge(data_res, tmp1, by = 'code')
  tmp1 <- tmp1[date == max_date]
  
  tmp1 <- merge(tmp1, prop_vac_counterfactual_df, by = c('state_index', 'counterfactual_index'))
  tmp1[, age_counterfactual2 := gsub('.* aged (.+)', '\\1', label_counterfactual)]
  tmp1[, label_age_counterfactual := paste0('Counterfactual analysis with a change in the\nvaccine coverage among individuals aged ', age_counterfactual)]
  
  tmp1 <- tmp1[age_counterfactual2 != '18-64 and 65+']
  tmp1 <- tmp1[age_counterfactual == age_counterfactual2]
  
  cols <- viridisLite::viridis(length(unique(tmp1$label_counterfactual)) , direction = -1, begin = 0.1)
  
  p <- ggplot(tmp1, aes(x = diff_value)) + 
    geom_hline(aes(yintercept=0), linetype = 'dashed') +
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.9, width = 0, col = 'grey40') + 
    geom_point(aes(y = M, col = label_counterfactual)) + 
    facet_grid(loc_label~age) +
    scale_color_manual(values = cols) + 
    scale_fill_manual(values = cols) + 
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom') + 
    labs(col = '', y = paste0('Change in age-specific COVID-19 attributable weekly deaths\nat the end of the resurgence period'),
         fill = '', linetype = '', 
         x = 'Change in vaccination coverage') +
    guides(fill=guide_legend(nrow=2,byrow=TRUE, order =1), 
           col=guide_legend(nrow=2,byrow=TRUE, order =1), 
           linetype = guide_legend(order=2)) 
  ggsave(p, file = paste0(outdir, '-predicted_change_weekly_deaths_vaccine_coverage_', lab, '.png'), w = 7.5, h = 5 + 2*(length(unique(data_res$code))/4))
  
  p <- ggplot(tmp1, aes(x = diff_value)) + 
    geom_hline(aes(yintercept=0), linetype = 'dashed') +
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.9, width = 0, col = 'grey40') + 
    geom_point(aes(y = M, col = label_age_counterfactual)) + 
    facet_grid(loc_label~age) +
    scale_color_manual(values = cols) + 
    scale_fill_manual(values = cols) + 
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom') + 
    labs(col = '', y = paste0('Change in age-specific COVID-19 attributable weekly deaths\nat the end of the resurgence period'),
         fill = '', linetype = '', 
         x = 'Change in vaccination coverage') +
    guides(fill=guide_legend(nrow=1,byrow=TRUE, order =1), 
           col=guide_legend(nrow=1,byrow=TRUE, order =1), 
           linetype = guide_legend(order=2)) 
  ggsave(p, file = paste0(outdir, '-predicted_change_weekly_deaths_vaccine_coverage_', lab, '2.png'), w = 7, h = 5 + 2*(length(unique(data_res$code))/4))
  
  p <- ggplot(tmp1, aes(x = diff_value)) + 
    geom_hline(aes(yintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_vline(aes(xintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.9, width = 0, col = 'grey40') + 
    geom_point(aes(y = M, col = label_age_counterfactual)) + 
    facet_grid(loc_label~age) +
    scale_color_manual(values = cols) + 
    scale_fill_manual(values = cols) + 
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom') + 
    labs(col = '', y = paste0('Change in age-specific COVID-19 attributable weekly deaths\nat the end of the resurgence period'),
         fill = '', linetype = '', 
         x = 'Change in vaccination coverage') +
    guides(fill=guide_legend(nrow=1,byrow=TRUE, order =1), 
           col=guide_legend(nrow=1,byrow=TRUE, order =1), 
           linetype = guide_legend(order=2)) 
  ggsave(p, file = paste0(outdir, '-predicted_change_weekly_deaths_vaccine_coverage_', lab, '3.png'), w = 7, h = 5 + 2*(length(unique(data_res$code))/4))
  
}


