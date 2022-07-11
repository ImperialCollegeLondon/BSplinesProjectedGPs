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

plot_mortality_rate_all_states = function(mortality_rate, crude_mortality_rate, outdir)
{
  
  tmp =   subset(mortality_rate, age == '85+')
  age = unique(tmp$age)
  tmp = tmp[order(M)]
  medianM =  tmp[,median(M)]
  mortality_rate[, loc_label := factor(loc_label, tmp$loc_label)]
  
  mortality_rate <- mortality_rate[age != '0-24']
  mortality_rate[, `Age group` := age]
  
  crude_mortality_rate <- crude_mortality_rate[age != '0-24']
  crude_mortality_rate[, `Age group` := age]
  
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
  
  p <- p + 
    geom_point(data = crude_mortality_rate, aes(x = loc_label, y = crude_mortality_rate, col = 'CDC crude estimate')) + 
    labs(col = '')
  ggsave(p, file = paste0(outdir, paste0('-MortalityRate_allages_withcrude.png')), w = 7, h = 9)
  
  
  # +
  #   labs(y = '')
  
}

plot_mortality_rate_all_states2 = function(mortality_rate, outdir)
{
  
  tmp =   subset(mortality_rate, age == '85+')
  age = unique(tmp$age)
  tmp = tmp[order(M)]
  mortality_rate[, loc_label := factor(loc_label, tmp$loc_label)]
  
  mortality_rate <- mortality_rate[age != '0-24']
  mortality_rate[, `Age group` := age]
  
  mortality_rate[, medianM := median(M), by = 'age']
  
  crude_mortality_rate <- crude_mortality_rate[age != '0-24']
  crude_mortality_rate[, `Age group` := age]

  my_plot <- function(tmp, with_axis){
    
    pp <- ggplot(tmp, aes(x=loc_label, y = M)) + 
      geom_hline(aes(yintercept = tmp[, unique(medianM)], linetype = ''), alpha = 0.75, col = 'darkblue') + 
      geom_bar(aes(fill = M), stat="identity") +
      geom_errorbar(aes(ymin=CL, ymax=CU), width=.2, position=position_dodge(.9), color = 'grey30') + 
      scale_fill_gradient2(low= 'darkturquoise', high = 'darkred', mid = 'beige',
                           midpoint = tmp[, unique(medianM)], labels = scales::percent_format()) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 70,hjust=1,vjust=1), 
            legend.position = 'bottom', 
            panel.grid.minor= element_blank(), 
            panel.grid.major.x= element_blank(), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.text = element_text(colour = 'white'),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA)) + 
      labs(x = '', y = '', fill = '', linetype = 'National median') + 
      scale_linetype_manual(values = 'dashed') +
      scale_y_continuous(expand =expansion(mult = c(0, .05)), labels = scales::percent_format())+ 
      guides(fill = guide_colourbar(barwidth = 10,  barheight = 1, order = 1)) +  
      facet_grid(`Age group`~., label = 'label_both') 
    
    if(!with_axis){
      pp <- pp + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
    
    }
    return(pp)
  }
  
  p1 <- my_plot(mortality_rate[age == sort(unique(mortality_rate$age))[1]], F)
  p2 <- my_plot(mortality_rate[age == sort(unique(mortality_rate$age))[2]], F)
  p3 <- my_plot(mortality_rate[age == sort(unique(mortality_rate$age))[3]], F)
  p4 <- my_plot(mortality_rate[age == sort(unique(mortality_rate$age))[4]], T)

  p <- ggarrange(p1, p2, p3, p4,common.legend = T, legend = 'bottom', ncol= 1, heights = c(0.2, 0.2, 0.2, 0.3))
  p <- grid.arrange(p, left = paste0('                             Predicted COVID-19 attributable mortality rates as of ', format(unique(mortality_rate$date), '%B %Y')),
                    bottom = text_grob('Lower <-> Higher   \nthan national median', hjust = 1.15, vjust = -0.1, size = 10))
  ggsave(p, file = paste0(outdir, paste0('-MortalityRate_allages2.png')), w = 7, h = 9)
  
}

plot_mortality_rate_all_states_map <- function(mortality_rate, mortality_rateJ21, outdir){
  
  mortality_rate <- mortality_rate[age != '0-24']
  
  mortality_rate[, state := code]
  mortality_rateJ21[, state := code]

  mortality_rate <- mortality_rate[, .(M, state, age)]
  mortality_rateJ21 <- mortality_rateJ21[, .(state, M_rel, age)]
  
  to_include <- mortality_rate[, unique(state)]
  
  mortality_rate[, `Age group`:=age]
  mortality_rateJ21[, `Age group`:=age]
  
  tmp <- mortality_rate[age == '85+']
  MedianM <- tmp[, unique(median(M))]
  p1 <- plot_usmap(data =tmp, values = 'M', include = to_include, color='white') +
    # facet_wrap(~`Age group`, label = 'label_both') + 
    theme(legend.position = "right",
          legend.text = element_text(size = rel(1)), 
          strip.text = element_text(size = rel(1)),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "white", fill = NA)) + 
    labs(fill = 'Predicted COVID-19\nattributable mortality\nrate in 85+') + 
    scale_fill_gradient(low= '#FBF4F3', high = 'firebrick4', labels = scales::percent_format()) + 
    guides(fill = guide_colorbar(barwidth = 0.4))
  
  my_plots <- function(tmp){
    plot_usmap(data = tmp, values = 'M_rel', include = to_include, color='white') +
      # facet_wrap(~`Age group`, label = 'label_both') + 
      theme(legend.position = "right",
            legend.text = element_text(size = rel(1)), 
            strip.text = element_text(size = rel(1)),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "white", fill = NA))  + 
      scale_fill_gradient(high= '#013D47', low = '#E4F3F4') + 
      guides(fill = guide_colorbar(barwidth = 0.4)) + 
      labs(fill = 'Predicted COVID-19\nattributable mortality\nrate ratio in 85+\nrelative to 55-84')
  }

  # p2 <- my_plots(copy(mortality_rate[age == '25-54']))
  p3 <- my_plots(copy(mortality_rateJ21[age == '55-84']))
  # p4 <- my_plots(copy(mortality_rate[age == '75-84']))
  
  # p <- grid.arrange(p1, p3, layout_matrix = rbind(c(1, 1, 2), c(NA, 3, 4)), widths = c(0.04, 0.48, 0.48))
  p <- grid.arrange(p1, p3, nrow = 2)
  ggsave(p, file = paste0(outdir, paste0('-MortalityRate_map.png')), w = 9, h = 7)
  
  # p <- grid.arrange(p1, p3, p4, layout_matrix = rbind(c(1, 3, 4)), widths = c(0.35, 0.3, 0.3))
  # ggsave(p, file = paste0(outdir, paste0('-MortalityRate_map2.png')), w = 12, h = 4)
  # p <- grid.arrange(p1, p3, p4, ncol = 1)
  # ggsave(p, file = paste0(outdir, paste0('-MortalityRate_map3.png')), w = 4, h = 7)
  
  return(p)
}

plot_contributiondiff_map <- function(contributiondiff, var, outdir, lab = NULL){
  
  if(var == 'diff1'){
    label = 'Estimated change in the contribution of 65+\nto COVID-19 weekly deaths between\nMay 2020 to when vaccination started'
  }
  
  if(var == 'diff2'){
    label = 'Estimated change in the contribution of 65+\nto COVID-19 weekly deaths two months\nafter vacccination started'
    
  }
  contributiondiff <- contributiondiff[age == '65+' & variable==var]
  
  contributiondiff[, state := code]
  contributiondiff[, change := 'no significant change']
  contributiondiff[CL < 0 & CU < 0, change := 'significant decrease']
  contributiondiff[CL > 0 & CU > 0, change := 'significant increase']
  contributiondiff[, change := factor(change, levels = c('significant increase', 'no significant change', 'significant decrease'))]
  contributiondiff <- contributiondiff[, .(change, state)]
  
  to_include <- contributiondiff[, unique(state)]
  
  cols <- ggsci::pal_npg()(10)[c(10,5,6)]
  
  p <- plot_usmap(data =contributiondiff, values = 'change', include = to_include, color='white') +
    # facet_wrap(~`Age group`, label = 'label_both') + 
    theme(legend.position = "right",
          legend.text = element_text(size = rel(1)), 
          strip.text = element_text(size = rel(1)),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "white", fill = NA)) + 
    labs(fill = label) + 
    scale_fill_manual(values = cols, drop = FALSE)  +
    bgcolor('white')
  

  ggsave(p, file = paste0(outdir, paste0('-contribution', var, lab, '_map.png')), w = 9, h = 7)
  ggsave(p, file = paste0(outdir, paste0('-contribution', var, lab, '_map.pdf')), w = 9, h = 7)
  
  return(p)
}


plot_mortality_rate_continuous_all_states = function(mortality_rate, selected_codes, outdir)
{
  
  mortality_rate[, age := as.numeric(age)]
  mortality_rate[, age_cat := as.character(age)]
  mortality_rate[age_cat == '85', age_cat := '85+']
  
  tmp <- mortality_rate[code %in% selected_codes]
  
  p <- ggplot(subset(tmp, age != '85'), aes(x=age)) + 
    geom_line(aes(y = M, col = loc_label)) +
    geom_ribbon(aes(ymin=CL, ymax=CU, fill = loc_label), alpha = 0.4) + 
    theme_bw() +
    theme(legend.position = 'bottom', 
          panel.grid.minor= element_blank(), 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + 
    scale_color_jcolors('pal8') + 
    scale_fill_jcolors('pal8') + 
    scale_y_continuous(expand = c(0,0), labels = scales::percent_format(), 
                       breaks = seq(0, max(tmp$CU), 0.01)) +
    scale_x_continuous(expand = c(0,0)) +
    labs(y = paste0('Predicted COVID-19 attributable\nmortality rates as of ', format(unique(mortality_rate$date), '%B %Y')),
         col = '', fill = '', x = 'Age')
  
  ggsave(p, file = paste0(outdir, paste0('-MortalityRateContinuous_allages_selectedstates_wo85p.png')), w = 6.5, h = 3.3)
  
  p <- p + theme(plot.margin = unit(c(5.5,0,5.5,5.5), "pt"))
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
    labs(y = paste0('Predicted COVID-19 attributable\nmortality rates as of ', format(unique(mortality_rate$date), '%b %Y')),
         col = '')
  
  p <- ggarrange(p, p1, nrow = 1, common.legend = T, legend = 'bottom', widths = c(1, 0.1))
  ggsave(p, file = paste0(outdir, paste0('-MortalityRateContinuous_allages_selectedstates.png')), w = 6, h = 4)
  
  
  ###
  
  if('NY' %in% mortality_rate[, unique(code)]){
    tmp <- subset(mortality_rate, code == 'NY')
    
    if(nrow(tmp) == 0){
      tmp <- subset(mortality_rate, code == unique(mortality_rate$code)[1])
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
      labs(x = 'Age', y = paste0('Predicted COVID-19 attributable\nmortality rates as of ', format(unique(mortality_rate$date), '%B %Y'))) 
    ggsave(p, file = paste0(outdir, paste0('-MortalityRateContinuous_allages_NY.png')), w = 6.5, h = 3)
    
  }

}

plot_mortality_all_states = function(death, resurgence_dates, lab = 'allStates', outdir)
{
  
  df = as.data.table( reshape2::melt(select(death, loc_label, code, date, age, emp), id.vars = c('loc_label', 'code', 'date', 'age')) )
  df[, variable2 := 'CDC data']
  df[, `Age group` := age]
  
  death[, dummy := 'Posterior median prediction\nusing age-aggregated JHU data\nto adjust for reporting delays']
  # death[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  colfunc <- jcolors("pal8")[1:length(death[, unique(code)])]

  if(length(death[, unique(code)]) > 4){
    colfunc <- jcolors("pal8")[1:length(locs)]
    colfunc <- colfunc[-(1:4)]
  }
  
  death[, `Age group` := age]
  
  p <- ggplot(subset(death), aes(x= date) ) +
    geom_point(data = subset(df), aes(y = value, shape= variable2), col = 'black', fill = 'black', size = 0.9, alpha = 0.7) + 
    geom_line(aes(y = M, fill = loc_label), show.legend = F) +
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = loc_label), alpha = 0.5, show.legend = F) + 
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
          legend.title = element_text(size = rel(1.1)), 
          legend.text = element_text(size = rel(1.1)), 
          strip.text = element_text(size = rel(1)), 
          legend.position = 'bottom') +
    labs( y = 'Predicted COVID-19 attributable weekly deaths', shape = '',
          linetype = '', col = '', fill = '', alpha = '') + 
    scale_shape_manual(values = 16) +
    scale_linetype_manual(values = 2) + 
    scale_color_manual(values = as.character(colfunc)) + 
    scale_fill_manual(values = as.character(colfunc)) + 
    guides(shape = guide_legend(override.aes = list(size=1, stroke =1.5), order = 2),
           linetype = guide_legend(order = 3), 
           col = guide_legend(order = 1),
           fill = guide_legend(order = 1))

  if(!is.null(resurgence_dates)){
    dummy.dt <- resurgence_dates[, list(date = seq.Date(start_resurgence, stop_resurgence, by = 'week')), by = 'code']
    dummy.dt = merge(dummy.dt, unique(select(death, code, loc_label)), by = 'code')
    dummy.dt[, text := 'Summer 2021 resurgence period']
    p = p +     geom_ribbon(data = dummy.dt, aes(ymin = -Inf, ymax = Inf, alpha = text)) 
  }
  
  ggsave(p, file = paste0(outdir, paste0('-Mortality_', lab, '.png')), w = 7.5, h = 8)
  
  
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

plot_relative_resurgence_vaccine2 <- function(data_res1, log_transform, outdir, lab.fig = '', 
                                                  selected_codes = NULL, withlimits = T){
  
  prop_vac_init = unique(data_res1[, list(prop_1_init = prop_1[date == min(date)], prop_2_init = prop_2[date == min(date)]), by = 'code'])
  
  data_res = merge(data_res1, prop_vac_init, by = 'code')
  
  if(!is.null(selected_codes)){
    data_res <- subset(data_res, code %in% selected_codes)
  }
  
  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  # data_res = as.data.table( merge(data_res, resurgence_dates, by = 'code') )
  
  data_res[, max_week_index := max(week_index), by = 'code']
  data_text = data_res[ week_index == max_week_index]
  data_text[, week_index := max(week_index)]
  data_text[code == 'CA' & age == '18-64', M := M + 0.1]
  data_text[code == 'TX' & age == '65+', M := M + 0.07]
  data_text[code == 'PA' & age == '65+', M := M + 0.05]
  
  if(withlimits){
    data_text[code == 'CA' & age == '65+', M := M + 0.13]
  }
  
  data_res2 <- data_res[week_index == max_week_index]
  
  lims = data_res1[, range(c(CL, CU))]
  lim1864 = prop_vac_init[, range(prop_1_init)]
  lim65p =  prop_vac_init[, range(prop_2_init)]
  lab = function(Age) paste0('Pre-resurgence\nvaccination rate in ', Age)
  lab_right = function(Age) paste0('Pre-resurgence\nvaccination rate\nin ', Age)
  
  space_mult <- 0.05
  if(any(nchar(data_res[, unique(loc_label)]) > 10)){
    space_mult <- 0.3
  }
  
  p1 <- ggplot(data_res[age == '18-64'], aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_1_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_1_init, group = loc_label), alpha = 0.5) +
    facet_grid(.~`Age group`, label = 'label_both') +
    labs(y = 'Relative COVID-19\nweekly deaths', 
         x = 'Week index of Summer 2021 resurgences period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('18-64'), fill = lab('18-64')) + 
    theme_bw() +
    geom_text(data = data_text[age == '18-64'], aes(label = code, x = week_index, y = M), hjust = -.1, size = 2.5) +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    scale_x_continuous(expand=  expansion(mult = c(0,space_mult)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    scale_y_continuous(limits = lims) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1)),
          legend.spacing.x = unit(0.3, "cm"), 
          strip.text = element_text(size = rel(0.9))) +
    scale_color_gradient(high = '#33ccff', low = '#ff6600', 
                          labels = scales::percent_format(accuracy = 1), limits = lim1864) + 
    scale_fill_gradient(high = '#33ccff', low = '#ff6600', 
                         labels = scales::percent_format(accuracy = 1), limits = lim1864) + 
    # scale_color_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)), 
    #                       labels = scales::percent_format(accuracy = 1), limits = lim1864) + 
    # scale_fill_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)), 
    #                      labels = scales::percent_format(accuracy = 1), limits = lim1864) + 
    guides(shape = 'none')
  
  p2 <- ggplot(data_res[age == '65+'], aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_2_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_2_init, group = loc_label), alpha = 0.5, width = 0.1) +
    facet_grid(.~`Age group`, label = 'label_both') +
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = 'Week index of Summer 2021 resurgence period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('65+'), fill = lab('65+')) + 
    theme_bw() +
    geom_text(data = data_text[age == '65+'], aes(label = code, x = week_index, y = M), hjust = -.1, size = 2.5) +
    # scale_x_date(breaks = '2 weeks', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") +
    scale_x_continuous(expand=  expansion(mult = c(0,space_mult)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_blank(),
          strip.text = element_text(size = rel(0.9)),
          legend.spacing.x = unit(0.3, "cm")) +
    scale_color_gradient(low = '#006666', high = '#ffccff',
                          labels = scales::percent_format(accuracy = 1), limits = lim65p) + 
    scale_fill_gradient(low = '#006666', high = '#ffccff', 
                         labels = scales::percent_format(accuracy = 1), limits = lim65p) + 
    # scale_color_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)), 
    #                       labels = scales::percent_format(accuracy = 1), limits = lim65p) + 
    # scale_fill_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)),
    #                      labels = scales::percent_format(accuracy = 1), limits = lim65p) + 
    guides(shape = 'none')
  
  if(withlimits){
    p1 <- p1 + theme(legend.position = 'bottom')
    p2 <- p2 + scale_y_continuous(limits = lims) + 
      theme(          axis.text.y = element_text(colour = 'white'),
                      axis.ticks.y = element_line(colour = 'white'), 
                      legend.position = 'bottom')
  }else{
    p1 <- p1 + guides(fill = guide_colorbar(barwidth = 0.4), col = guide_colorbar(barwidth = 0.4)) + 
      labs(col = lab_right('18-64'), fill = lab_right('18-64'), 
           y = 'Relative COVID-19 attributable\nweekly deaths') + theme(legend.position = 'right')
    p2 <- p2 + guides(fill = guide_colorbar(barwidth = 0.4), col = guide_colorbar(barwidth = 0.4))+ 
      labs(col = lab_right('65+'), fill = lab_right('65+')) + theme(legend.position = 'right')
  }
  if(log_transform){
    p1 = p1 + scale_y_continuous(trans = 'log', breaks = base_breaks()) + labs(y = 'log relative COVID-19 attributable weekly deaths')
    p2 = p2 + scale_y_continuous(trans = 'log', breaks = base_breaks()) 
  }
  
  p = ggarrange(p1, p2, ncol = 2, widths = c(0.51, 0.49))
  
  if(log_transform){
    file =  paste0(outdir, '-log_relative_deaths_vaccine_coverage2', lab.fig, '.png')
  } else{
    file =  paste0(outdir, '-relative_deaths_vaccine_coverage2', lab.fig, '.png')
  }
  
  ggsave(p, file =file, w = 10.5, h = 3)
  
  p1 <- p1 + theme(axis.title.x = element_text(hjust = 1.04)) + 
    labs(x ='', y  ="\nRelative COVID-19 attributable\nweekly deaths") + theme(legend.position = 'none')
  p2 <- p2 + theme(axis.title.y = element_blank(), axis.title.x = element_text(hjust = -0.05)) + 
    labs(x='', y  ="Over time") + theme(legend.position = 'none')

  return(list(p1, p2))
}

plot_relative_resurgence_vaccine_panel <- function(p4, p_all, lab, outdir){
  p4[[1]] <- ggarrange(p4[[1]], labels = 'A')
  p_all[[1]] <-  ggarrange(p_all[[1]], labels = 'B', label.y = 1.07, label.x = 0.03)
  p <- grid.arrange(grobs = c(p4, p_all), layout_matrix = rbind(c(NA, 1, NA, 2), c(3, 3, 4, 4)), 
                    heights = c(0.46, 0.54), widths = c(0.0155034, 0.5121270, 0.01664, 0.4510599), 
                    bottom = text_grob('Week index of the summer 2021 resurgence period', vjust = -36.5, size =10))
  ggsave(p, file =paste0(outdir, '-relative_deaths_vaccine_coverage_panel_', lab, '.png'), w = 8, h = 6.5)
}


plot_relative_resurgence_vaccine2_old <- function(data_res1, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, log_transform, outdir, lab.fig = '', 
                                              selected_codes = NULL, withlimits = T){
  
  prop_vac_init = prop_vac[, list(prop_1_init = prop_1[date == min(date)], prop_2_init = prop_2[date == min(date)]), by = 'code']

  data_res = merge(data_res1, prop_vac_init, by = 'code')
  
  if(!is.null(selected_codes)){
    data_res <- subset(data_res, code %in% selected_codes)
  }
  
  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  data_res = as.data.table( merge(data_res, resurgence_dates, by = 'code') )
  
  data_res[, max_week_index := max(week_index), by = 'code']
  data_text = data_res[ week_index == max_week_index]
  data_text[, week_index := max(week_index)]
  data_text[code == 'CA' & age == '18-64', M := M + 0.1]
  data_text[code == 'TX' & age == '65+', M := M + 0.07]
  data_text[code == 'PA' & age == '65+', M := M + 0.05]
  
  if(withlimits){
    data_text[code == 'CA' & age == '65+', M := M + 0.13]
  }
  
  data_res2 <- data_res[week_index == max_week_index]

  lims = data_res1[, range(c(CL, CU))]
  lim1864 = prop_vac_init[, range(prop_1_init)]
  lim65p =  prop_vac_init[, range(prop_2_init)]
  lab = function(Age) paste0('Pre-resurgence\nvaccination rate in ', Age)
  lab_right = function(Age) paste0('Pre-resurgence\nvaccination rate\nin ', Age)
  
  space_mult <- 0.25
  if(any(nchar(data_res[, unique(loc_label)]) > 10)){
    space_mult <- 0.3
  }
  
  p1 <- ggplot(data_res[age == '18-64'], aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_1_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_1_init, group = loc_label), alpha = 0.5) +
    geom_point(data = data_res2[age == '18-64'], aes(y = M, col = prop_1_init, shape = loc_label)) + 
    facet_grid(.~`Age group`, label = 'label_both') +
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = 'Week index of Summer 2021 resurgences period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('18-64'), fill = lab('18-64')) + 
    theme_bw() +
    geom_text(data = data_text[age == '18-64'], aes(label = loc_label, x = week_index, y = M, col = prop_1_init), hjust = -.1, size = 3) +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    scale_x_continuous(expand=  expansion(mult = c(0,space_mult)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    scale_y_continuous(limits = lims) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1)),
          legend.spacing.x = unit(0.3, "cm"), 
          strip.text = element_text(size = rel(0.9))) +
    scale_color_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)), 
                          labels = scales::percent_format(accuracy = 1), limits = lim1864) + 
    scale_fill_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)), 
                         labels = scales::percent_format(accuracy = 1), limits = lim1864) + 
    scale_shape_manual(values = c(15, 17, 20, 4, 3, 10, 11, 12, 13, 14)[1:length(data_res$code)]) + 
    guides(shape = 'none')
  
  p2 <- ggplot(data_res[age == '65+'], aes(x = week_index)) + 
    geom_line(aes(y = M, col = prop_2_init, group = loc_label)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = prop_2_init, group = loc_label), alpha = 0.5, width = 0.1) +
    geom_point(data = data_res2[age == '65+'], aes(y = M, col = prop_2_init, shape = loc_label)) + 
    facet_grid(.~`Age group`, label = 'label_both') +
    labs(y = 'Relative COVID-19 attributable weekly deaths', 
         x = 'Week index of Summer 2021 resurgence period', shape = 'Beginning of Summer 2021 resurgence period', 
         col = lab('65+'), fill = lab('65+')) + 
    theme_bw() +
    geom_text(data = data_text[age == '65+'], aes(label = loc_label, x = week_index, y = M, col = prop_2_init), hjust = -.1, size = 3) +
    # scale_x_date(breaks = '2 weeks', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") +
    scale_x_continuous(expand=  expansion(mult = c(0,space_mult)), breaks = min(data_res$week_index) + 2*c(0:floor((max(data_res$week_index)-1)/2) )) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_blank(),
          strip.text = element_text(size = rel(0.9)),
          legend.spacing.x = unit(0.3, "cm")) +
    scale_color_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)), 
                          labels = scales::percent_format(accuracy = 1), limits = lim65p) + 
    scale_fill_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)),
                         labels = scales::percent_format(accuracy = 1), limits = lim65p) + 
    scale_shape_manual(values = c(15, 17, 20, 4, 3, 10, 11, 12, 13, 14)[1:length(data_res$code)])+ 
    guides(shape = 'none')
  
  if(withlimits){
    p1 <- p1 + theme(legend.position = 'bottom')
    p2 <- p2 + scale_y_continuous(limits = lims) + 
      theme(          axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(), 
                      legend.position = 'bottom')
  }else{
    p1 <- p1 + guides(fill = guide_colorbar(barwidth = 0.4), col = guide_colorbar(barwidth = 0.4)) + 
      labs(col = lab_right('18-64'), fill = lab_right('18-64'), 
           y = 'Relative COVID-19 attributable\nweekly deaths') + theme(legend.position = 'right')
    p2 <- p2 + guides(fill = guide_colorbar(barwidth = 0.4), col = guide_colorbar(barwidth = 0.4))+ 
      labs(col = lab_right('65+'), fill = lab_right('65+')) + theme(legend.position = 'right')
  }
  if(log_transform){
    p1 = p1 + scale_y_continuous(trans = 'log', breaks = base_breaks()) + labs(y = 'log relative COVID-19 attributable weekly deaths')
    p2 = p2 + scale_y_continuous(trans = 'log', breaks = base_breaks()) 
  }
  
  p = ggarrange(p1, p2, ncol = 2, widths = c(0.51, 0.49))
  
  if(log_transform){
    file =  paste0(outdir, '-log_relative_deaths_vaccine_coverage2', lab.fig, '.png')
  } else{
    file =  paste0(outdir, '-relative_deaths_vaccine_coverage2', lab.fig, '.png')
  }

  ggsave(p, file =file, w = 10.5, h = 3)
  
  p1 <- p1 +  guides(colour = 'none', fill ='none') + theme(axis.title.x = element_text(hjust = 1.05)) + 
    labs(x ='Week index of the summer', y  ="\nRelative COVID-19 attributable\nweekly deaths")
  p2 <- p2 +  guides(colour = 'none', fill ='none') + theme(axis.title.y = element_blank(), axis.title.x = element_text(hjust = -0.05)) + 
    labs(x='2021 resurgence period', y  ="Over time")
  p = ggarrange(p1, p2, ncol = 2, widths = c(0.53, 0.47))
  
  return(p)
}

plot_relative_resurgence_vaccine_end_2 <- function(data_res1, prop_vac, df_age_vaccination2, df_week2, resurgence_dates, log_transform, outdir, lab.fig = ''){
  
  prop_vac_init = prop_vac[, list(prop_1_init = prop_1[date == min(date)], prop_2_init = prop_2[date == min(date)]), by = 'code']
  data_res = merge(data_res1, prop_vac_init, by = 'code')
  
  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  data_res = as.data.table( merge(data_res, resurgence_dates, by = 'code') )
  
  data_res[, max_week_index := max(week_index), by = 'code']
  
  # last week
  data_res <- data_res[week_index == max_week_index]

  lab = function(Age) paste0('Pre-resurgence\nvaccination rate in ', Age)
  lab_long = function(Age) paste0('Pre-resurgence vaccination rate in ', Age)
  
  p1 <- ggplot(data_res[age == '18-64'], aes(x = prop_1_init)) + 
    geom_point(aes(y = M, col = prop_1_init, shape = loc_label)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5, width = 0) +
    labs(y = 'Relative COVID-19 attributable\nweekly deaths at the peak of the\nsummer 2021 resurgence wave', 
         x = lab_long('18-64'), shape = '', 
         col = lab('18-64')) +
    theme_bw() +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    scale_x_continuous(labels = scales::percent_format()) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(),
          strip.text = element_blank(),
          legend.spacing.x = unit(0.3, "cm"), 
          legend.position = 'bottom') +
    scale_color_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)), 
                          labels = scales::percent_format(accuracy = 1)) + 
    scale_shape_manual(values = c(15, 17, 3, 10, 11, 20, 12, 13, 14, 4)[1:length(data_res$code)]) + 
    guides(col = 'none')

  p2 <- ggplot(data_res[age == '65+'], aes(x = prop_2_init)) + 
    geom_point(aes(y = M, col = prop_2_init, shape = loc_label, fill = prop_1_init)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5, width = 0) +
    labs(x = lab_long('65+'), shape = '', 
         col = lab('65+'), fill = lab('18-64')) + 
    theme_bw() +
    # scale_x_date(breaks = '2 weeks', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") +
    scale_x_continuous(labels = scales::percent_format()) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          # legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_blank(),
          # axis.text.y = element_blank(),
          strip.text = element_blank(),
          # strip.text = element_text(size = rel(0.9)),
          # legend.spacing.x = unit(0.3, "cm"), 
          legend.position = 'bottom') +
    scale_color_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)), 
                          labels = scales::percent_format(accuracy = 1)) + 
    scale_fill_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)), 
                         labels = scales::percent_format(accuracy = 1)) +
    scale_shape_manual(values = c(15, 17, 3, 10, 11, 20, 12, 13, 14, 4)[1:length(data_res$code)]) + 
    guides(fill = guide_colorbar(order = 1, row =1), 
           col = guide_colorbar(order = 2, row = 1), 
           shape = 'none')
  
  if(log_transform){
    p1 = p1 + scale_y_continuous(trans = 'log', breaks = base_breaks()) + labs(y = 'log relative COVID-19 attributable weekly deaths')
    p2 = p2 + scale_y_continuous(trans = 'log', breaks = base_breaks()) 
  }
  
  p = ggarrange(p1, p2, nrow = 1, widths = c(0.51, 0.49), legend.grob = get_legend(p2), legend = 'bottom')
  p = ggarrange(p, legend.grob = get_legend(p1), legend = 'bottom')
  
  if(log_transform){
    file =  paste0(outdir, '-log_relative_deaths_vaccine_coverage2_all', lab.fig, '.png')
  } else{
    file =  paste0(outdir, '-relative_deaths_vaccine_coverage2_all', lab.fig, '.png')
  }
  
  ggsave(p, file =file, w = 8, h = 5)
  
  return(p)
}

plot_relative_resurgence_vaccine_end_3 <- function(data_res1,log_transform, outdir, lab.fig = '', selected_codes= NULL){
  
  prop_vac_init = unique(data_res1[, list(prop_1_init = prop_1[date == min(date)], prop_2_init = prop_2[date == min(date)]), by = 'code'])
  data_res = merge(data_res1, prop_vac_init, by = 'code')
  
  data_res[, `Age group` := age]
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  data_res[, max_week_index := max(week_index), by = 'code']
  
  # last week
  data_res <- data_res[week_index == max_week_index]
  
  # code 
  if(!is.null(selected_codes)){
    data_res <- data_res[code %in% selected_codes]
  }
  
  # data_text <- data_res[(age == '18-64' & (prop_1_init < 0.275 | prop_1_init > 0.5 | M > 2 | M < 0.2)) | (age == '65+' & (prop_2_init < 0.7 | prop_2_init > 0.875 | M > 0.7 | M < 0.08))]
  data_text <- copy(data_res)
  
  lab = function(Age) paste0('Pre-resurgence\nvaccination rate in ', Age)
  lab_long = function(Age) paste0('Pre-resurgence vaccination rate in ', Age)
  
  p1 <- ggplot(data_res[age == '18-64'], aes(x = prop_1_init)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5, width = 0) +
    geom_point(aes(y = M, col = prop_1_init)) + 
    labs(y = 'Relative COVID-19 attributable\nweekly deaths at the peak of the\nsummer 2021 resurgence wave', 
         x = lab_long('18-64'), 
         col = lab('18-64'), fill = lab('18-64')) +
    theme_bw() +
    # scale_x_date(breaks = '1 month', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") + 
    scale_x_continuous(labels = scales::percent_format()) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(),
          strip.text = element_blank(),
          legend.spacing.x = unit(0.3, "cm"), 
          legend.position = 'bottom') +
    scale_color_gradient(high = '#33ccff', low = '#ff6600', 
                         labels = scales::percent_format(accuracy = 1)) + 
    scale_fill_gradient(high = '#33ccff', low = '#ff6600', 
                        labels = scales::percent_format(accuracy = 1)) + 
    geom_text_repel(data = data_text[age == '18-64'], aes(x = prop_1_init, y = M, label = code), size = 2)
  
  p2 <- ggplot(data_res[age == '65+'], aes(x = prop_2_init)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5, width = 0) +
    geom_point(aes(y = M, col = prop_2_init)) + 
    labs(x = lab_long('65+'), shape = '', 
         col = lab('65+'), fill = lab('65+')) + 
    theme_bw() +
    # scale_x_date(breaks = '2 weeks', expand=  expansion(mult = c(0,0.25)), date_labels = "%b-%y") +
    scale_x_continuous(labels = scales::percent_format()) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          # legend.box="vertical", 
          legend.title = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_blank(),
          # axis.text.y = element_blank(),
          strip.text = element_blank(),
          # strip.text = element_text(size = rel(0.9)),
          # legend.spacing.x = unit(0.3, "cm"), 
          legend.position = 'bottom') +
    geom_text_repel(data = data_text[age == '65+'], aes(x = prop_2_init, y = M, label = code), size = 2)+
    scale_color_gradient(low = '#006666', high = '#ffccff',
                         labels = scales::percent_format(accuracy = 1)) + 
    scale_fill_gradient(low = '#006666', high = '#ffccff', 
                        labels = scales::percent_format(accuracy = 1)) 
  
  if(log_transform){
    p1 = p1 + scale_y_continuous(trans = 'log', breaks = base_breaks()) #+ labs(y = 'log relative COVID-19 attributable weekly deaths')
    p2 = p2 + scale_y_continuous(trans = 'log', breaks = base_breaks()) 
  }
  
  p = ggarrange(p1, p2, nrow = 1, widths = c(0.53, 0.47))
  
  if(log_transform){
    file =  paste0(outdir, '-log_relative_deaths_vaccine_coverage2_all', lab.fig, '.png')
  } else{
    file =  paste0(outdir, '-relative_deaths_vaccine_coverage2_all', lab.fig, '.png')
  }
  
  ggsave(p, file =file, w = 8, h = 5)
  
  return(list(p1, p2))
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

plot_relative_resurgence_vaccine2_long <- function(data_res1, outdir, labfig, selected_codes){
  
  data_res1 <- data_res1[code %in% selected_codes]
  
  prop_vac_init = unique(data_res1[, list(prop_1_init = prop_1[date == min(date)], prop_2_init = prop_2[date == min(date)]), by = 'code'])
  
  data_res = merge(data_res1, prop_vac_init, by = 'code')
  
  data_res[, `Age group` := age]

  lab = function(Age) paste0('Proportion of individuals aged ', Age, ' fully vaccinated two weeks\nbefore the beginning of Summer 2021 resurgence period')
  
  ## 18-64
  mid_code = round(length(selected_codes) / 2)
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
                          limits = range(prop_vac_init$prop_1_init),labels = scales::percent_format(accuracy = 1)) + 
    scale_fill_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)),
                         limits = range(prop_vac_init$prop_1_init),labels = scales::percent_format(accuracy = 1)) 
  
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
                          limits = range(prop_vac_init$prop_1_init),labels = scales::percent_format(accuracy = 1)) + 
    scale_fill_gradient2(high = 'darkred', low = 'cornflowerblue', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_1_init)), 
                         limits = range(prop_vac_init$prop_1_init),labels = scales::percent_format(accuracy = 1))
  
  p <- ggarrange(p1, p2, ncol = 1, common.legend = T, legend = 'bottom', heights = c(0.9, 1))
  p <- grid.arrange(p, left = 'Relative COVID-19 attributable weekly deaths')
  file =  paste0(outdir, '-relative_deaths_vaccine_coverage2', '_extended_part1', labfig, '.png')
  ggsave(p, file =file, w = 8, h = 6)
  
  ## 65+
  mid_code = round(length(selected_codes) / 2)
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
                          limits = range(prop_vac_init$prop_2_init),labels = scales::percent_format(accuracy = 1)) + 
    scale_fill_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)),
                         limits = range(prop_vac_init$prop_2_init),labels = scales::percent_format(accuracy = 1))
  
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
                          limits = range(prop_vac_init$prop_2_init),labels = scales::percent_format(accuracy = 1)) + 
    scale_fill_gradient2(low = 'lightpink', high = 'darkolivegreen', mid = 'moccasin', midpoint = mean(range(prop_vac_init$prop_2_init)), 
                         limits = range(prop_vac_init$prop_2_init),labels = scales::percent_format(accuracy = 1))
  
  p <- ggarrange(p1, p2, ncol = 1, common.legend = T, legend = 'bottom', heights = c(0.9, 1))
  p <- grid.arrange(p, left = 'Relative COVID-19 attributable weekly deaths')
  file =  paste0(outdir, '-relative_deaths_vaccine_coverage2', '_extended_part2', labfig, '.png')
  ggsave(p, file =file, w = 8, h = 6)
  
}

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

plot_PPC_relative_resurgence <- function(data_res1, data_res2, lab, outdir){
  
  data_res1[, type := 'Fit to observed data']
  data_res2[, type := 'Fit with meta-regression model']
  data_res = rbind(data_res1, data_res2, fill = T)
  
  ncol = data_res[, length(unique(code))]/5 - 1
  # data_res[, loc_label := factor(loc_label, levels = c('Florida', 'Texas', 'California', 'New York', 'Washington'))]
  
  p1 <- ggplot(data_res[age == '18-64'], aes(x = date)) + 
    geom_line(aes(y = M, col = type)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = type), alpha = 0.5) +
    facet_wrap(~loc_label, ncol = ncol, scale = 'free_y') +
    labs(y = 'Relative COVID-19 attributable weekly deaths in age group 18-64', 
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
  
  ggsave(p1, file = paste0(outdir, '-relative_deaths_vaccine_coverage_PPC_part_1_', lab, '.png'), w = 8, h = 9, limitsize = F)
  
  p1 <- p1 + scale_y_continuous(trans = 'log', breaks = base_breaks()) + labs(y = 'Relative COVID-19 attributable weekly deaths in age group 18-64 (log scale)')
  ggsave(p1, file = paste0(outdir, '-relative_deaths_vaccine_coverage_PPC_part_1_', lab, '_log.png'), w = 8, h = 9, limitsize = F)
  
  
  p2 <- ggplot(data_res[age == '65+'], aes(x = date)) + 
    geom_line(aes(y = M, col = type)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = type), alpha = 0.5) +
    facet_wrap(~loc_label, ncol = ncol) +
    labs(y = 'Relative COVID-19 attributable weekly deaths in age group 65+', 
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
  
  ggsave(p2, file = paste0(outdir, '-relative_deaths_vaccine_coverage_PPC_part_2_', lab, '.png'), w = 8, h = 9, limitsize = F)
  
  p2 <- p2 +  scale_y_continuous(trans = 'log', breaks = base_breaks()) + labs(y = 'Relative COVID-19 attributable weekly deaths in age group 65+ (log scale)')
  ggsave(p2, file = paste0(outdir, '-relative_deaths_vaccine_coverage_PPC_part_2_', lab, '_log.png'), w = 8, h = 9, limitsize = F)
  
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
  
  p = list()
  tmp1 <- subset(tmp, label_counterfactual == label_fit)
  
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
      scale_color_manual(values = 'black') + 
      scale_fill_manual(values = 'black') + 
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
  ggsave(p, file = paste0(outdir, '-predicted_weekly_deaths_vaccine_coverage_observed_', lab, '.png'), w = 6 + length(unique(data_res1$code))/6, h = 5 + 2*(length(unique(data_res1$code))/4))
  
  
  #######
  
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
  
  
  #####################################
  prop_vac_counterfactual_df <- copy(prop_vac_counterfactual)
  
  label <- function(age) paste0('Counterfactual analysis with a change in the\nvaccine coverage among individuals aged ', age)
  
  setnames(prop_vac_counterfactual_df, 'age_index', 'age_index_counterfactual')
  df <- copy(df_age_vaccination2[, .(age, age_index)])
  setnames(df, 1:2, c('age_counterfactual', 'age_index_counterfactual'))
  prop_vac_counterfactual_df <- merge(prop_vac_counterfactual_df, df, by = 'age_index_counterfactual')
  df <- unique(prop_vac_counterfactual_df[diff_value == 0, .(value_true, value, diff_value, age_counterfactual, age_index_counterfactual, state_index)])
  df[, counterfactual_index:=0]
  prop_vac_counterfactual_df <- rbind(prop_vac_counterfactual_df, df)
  
  tmp1 <- tmp[, list(max_date = max(date)), by = 'code']
  tmp1 <- merge(tmp, tmp1, by = 'code')
  tmp1 <- tmp1[date == max_date]
  
  df <- rbind(df_counterfactual, data.table(counterfactual_index = 0, label_counterfactual = label_fit))
  tmp1 <- merge(tmp1, df, by = 'label_counterfactual')
  tmp1 <- merge(tmp1, prop_vac_counterfactual_df, by = c('state_index', 'counterfactual_index'))
  
  tmp1[, age_counterfactual2 := gsub('.* aged (.+)', '\\1', label_counterfactual)]
  tmp1[, label_age_counterfactual := label(age_counterfactual2)]
  tmp1[label_counterfactual == label_fit, label_age_counterfactual := label_fit]
  
  tmp2 <- tmp1[age_counterfactual2 != '18-64 and 65+']
  tmp2 <- tmp2[age_counterfactual == age_counterfactual2 | label_counterfactual == label_fit]
  
  tmp2[, label_counterfactual := factor(label_counterfactual, 
                                        levels = c(label_fit, 
                                                   df_counterfactual[!grepl('18-64 and 65+', label_counterfactual), label_counterfactual]))]
  tmp2[, label_age_counterfactual := factor(label_age_counterfactual, 
                                            levels = c(label_fit, 
                                                       label('18-64'), label('65+')))]
  
  
  tmp2 <- tmp2[diff_value != 0 | label_counterfactual == label_fit]
  # df <- subset(tmp2, label_counterfactual == label_fit)
  # setnames(df, c('M', 'CL', 'CU'), c('M_fit', 'CL_fit', "CU_fit"))
  # tmp2 <- merge(tmp2, unique(df[, .(M_fit, CL_fit, CU_fit, code, age_index)]), by = c('code', 'age_index'))
  # 
  # 
  cols <- viridisLite::viridis(length(unique(tmp1$label_counterfactual)) , direction = -1, begin = 0.1)
  
  p <- ggplot(tmp2, aes(x = diff_value)) + 
    geom_hline(aes(yintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_vline(aes(xintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.9, width = 0, col = 'grey40') + 
    geom_point(aes(y = M, pch = label_age_counterfactual, col = loc_label)) + 
    geom_line(aes(y = M, col = loc_label), alpha = 0.5) + 
    facet_grid(.~age) +
    # scale_color_manual(values = cols) + 
    # scale_fill_manual(values = cols) + 
    # ggsci::scale_colour_npg() + 
    scale_color_jcolors(palette = "pal8") + 
    # scale_x_continuous(labels = scales::percent) +
    # scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom'
          # legend.box="vertical",
          # legend.spacing.y = unit(-0, "cm")
    ) + 
    labs(col = '', y = paste0('Change in age-specific COVID-19 attributable\nweekly deaths at the end of the resurgence period'),
         fill = '', linetype = '', 
         x = 'Change in vaccination rate', pch = '') +
    guides(col=guide_legend(nrow=4,byrow=TRUE, order =2), 
           pch = guide_legend(order=1,nrow=3,byrow=TRUE)) 
  ggsave(p, file = paste0(outdir, '-predicted_change_weekly_deaths_vaccine_coverage_', lab, '_abs.png'), w = 7.5, h = 3.5 + 2*(length(unique(tmp2$code))/4))
  
}

plot_vaccine_effects_counterfactual_change_old <- function(data_res, prop_vac_counterfactual, lab, namevar, outdir){
  
  prop_vac_counterfactual_df <- copy(prop_vac_counterfactual)
  
  label <- function(age) paste0('Counterfactual analysis with a change in the\nvaccine coverage among individuals aged ', age)
  label_higher <- function(age) paste0('Counterfactual analysis with higher vaccine coverage\namong individuals aged ', age)
  label_lower <- function(age) paste0('Counterfactual analysis with lower vaccine coverage\namong individuals aged ', age)
  
  setnames(prop_vac_counterfactual_df, 'age_index', 'age_index_counterfactual')
  df <- copy(df_age_vaccination2[, .(age, age_index)])
  setnames(df, 1:2, c('age_counterfactual', 'age_index_counterfactual'))
  prop_vac_counterfactual_df <- merge(prop_vac_counterfactual_df, df, by = 'age_index_counterfactual')
  
  tmp1 <- data_res[, list(max_date = max(date)), by = 'code']
  tmp1 <- merge(data_res, tmp1, by = 'code')
  tmp1 <- tmp1[date == max_date]
  
  tmp1 <- merge(tmp1, prop_vac_counterfactual_df, by = c('state_index', 'counterfactual_index'))
  tmp1[, age_counterfactual2 := gsub('.* aged (.+)', '\\1', label_counterfactual)]
  tmp1[, label_age_counterfactual := label(age_counterfactual2)]
  
  tmp2 <- tmp1[age_counterfactual2 != '18-64 and 65+']
  tmp2 <- tmp2[age_counterfactual == age_counterfactual2]
  
  cols <- viridisLite::viridis(length(unique(tmp1$label_counterfactual)) , direction = -1, begin = 0.1)
  
  p <- ggplot(tmp2, aes(x = diff_value)) + 
    geom_hline(aes(yintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_vline(aes(xintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.9, width = 0, col = 'grey40') + 
    geom_point(aes(y = M, col = label_counterfactual, pch = loc_label)) + 
    facet_grid(.~age) +
    # scale_color_manual(values = cols) + 
    # scale_fill_manual(values = cols) + 
    ggsci::scale_colour_npg() + 
    # scale_x_continuous(labels = scales::percent) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom'
          # legend.box="vertical",
          # legend.spacing.y = unit(-0, "cm")
    ) + 
    labs(col = '', y = paste0('Change in age-specific COVID-19 attributable\nweekly deaths at the end of the resurgence period'),
         fill = '', linetype = '', 
         x = 'Change in vaccination rate', pch = '') +
    guides(fill=guide_legend(nrow=4,byrow=TRUE, order =2), 
           col=guide_legend(nrow=4,byrow=TRUE, order =2), 
           pch = guide_legend(order=1,nrow=4,byrow=TRUE)) 
  
  if(grepl('perc', namevar)){
    p <- p + scale_y_continuous(labels = scales::percent) 
  }
  
  ggsave(p, file = paste0(outdir, '-predicted_', namevar, '_weekly_deaths_vaccine_coverage_', lab, '.png'), w = 7.5, h = 3.5 + 2*(length(unique(data_res$code))/4))
  
  
  p <- ggplot(tmp2, aes(x = diff_value)) + 
    geom_hline(aes(yintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_vline(aes(xintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.9, width = 0, col = 'grey40') + 
    geom_point(aes(y = M, col = label_counterfactual)) + 
    facet_grid(loc_label~age) +
    # scale_color_manual(values = cols) + 
    # scale_fill_manual(values = cols) + 
    ggsci::scale_colour_npg() + 
    # scale_x_continuous(labels = scales::percent) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom'
          # legend.box="vertical",
          # legend.spacing.y = unit(-0, "cm")
    ) + 
    labs(col = '', y = paste0('Change in age-specific COVID-19 attributable\nweekly deaths at the end of the resurgence period'),
         fill = '', linetype = '', 
         x = 'Change in vaccination rate', pch = '') +
    guides(fill=guide_legend(nrow=2,byrow=TRUE, order =1), 
           col=guide_legend(nrow=2,byrow=TRUE, order =1)) 
  
  if(grepl('perc', namevar)){
    p <- p + scale_y_continuous(labels = scales::percent) 
  }
  
  ggsave(p, file = paste0(outdir, '-predicted_', namevar, '_weekly_deaths_vaccine_coverage_', lab, '2.png'), w = 7, h = 5 + 2*(length(unique(data_res$code))/4))
  
}

plot_vaccine_effects_counterfactual_allages <- function(data_res1, data_res2, resurgence_dates, var, lab, outdir){
  
  label_fit <- 'Fit to observed data'
  data_res2[, label_counterfactual := label_fit]
  
  data_res1 = merge(data_res1, select(resurgence_dates, code, start_resurgence), by = 'code')
  data_res1 = data_res1[date >= start_resurgence ]
  data_res1 = select(data_res1, -start_resurgence, -counterfactual_index)
  
  tmp = rbind(data_res1, data_res2)
  
  tmp = subset(tmp, date >= as.Date('2021-01-01'))
  tmp = merge(tmp, select(resurgence_dates, code, stop_resurgence), by = 'code')
  tmp = tmp[date <= stop_resurgence]
  tmp[, label_counterfactual := factor(label_counterfactual, levels = c(label_fit, levels(df_counterfactual$label_counterfactual)))]
  
  dummy.dt = merge(resurgence_dates, df_state, by = 'code')
  dummy.dt[, text := 'Beginning of Summer 2021 resurgence period']
  
  # values_col = c('grey50', 'darkorchid4', )
  
  cols <- viridisLite::viridis(length(unique(tmp$label_counterfactual)), direction = -1, begin = 0.1)
  
  for(j in df_counterfactual$counterfactual_index){
    
    tmp1 <- subset(tmp, label_counterfactual %in% c(label_fit, df_counterfactual[counterfactual_index == j, as.character(label_counterfactual)]))
    
    p <- ggplot(tmp1, aes(x = date)) + 
      geom_line(aes(y = M, col = label_counterfactual)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = label_counterfactual), alpha = 0.5) + 
      facet_grid(loc_label~.) + 
      scale_x_date(expand = expansion(mult = c(0.05,0)), date_labels = c("%b-%y"), breaks = '1 month') + 
      theme_bw() + 
      geom_vline(data = dummy.dt, aes(xintercept = start_resurgence, linetype = text), col = 'grey50') +
      theme(strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA), 
            legend.position = 'bottom', 
            axis.text.x = element_text(angle = 70, hjust = 1),
            axis.title.x = element_blank(), 
            legend.direction = 'vertical',
            legend.box="vertical") + 
      # scale_color_manual(values = values_col) +
      # scale_fill_manual(values = values_col) + 
      scale_color_manual(values = c(cols[1], cols[j + 1])) + 
      scale_fill_manual(values = c(cols[1], cols[j + 1])) + 
      labs(col = '', y = paste0('Predicted ',var,' COVID-19 attributable weekly deaths among 18+'),
           fill = '', linetype = '') +
      scale_linetype_manual(values = 2) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE, order =1), col=guide_legend(nrow=2,byrow=TRUE, order =1), 
             linetype = guide_legend(order=2))
    
    ggsave(p, file = paste0(outdir, '-predicted_weekly_deaths_vaccine_coverage_counterfactual', j, '_', lab, 'AllAges.png'), w = 3 + length(unique(data_res1$code))/6, h = 5 + 2*(length(unique(data_res1$code))/4))
    
  }
  
  #############################################
  tmp1 <- tmp[, list(max_date = max(date)), by = 'code']
  tmp1 <- merge(tmp, tmp1, by = 'code')
  tmp1 <- tmp1[date == max_date]
  
  tmp2 <- subset(tmp1, label_counterfactual == label_fit)
  setnames(tmp2, c('M', 'CL', 'CU'), c('M_fit', 'CL_fit', "CU_fit"))
  tmp1 <- merge(tmp1, tmp2[, .(M_fit, CL_fit, CU_fit, code)], by = c('code'))
  
  p <- ggplot(tmp1, aes(x = label_counterfactual)) + 
    geom_hline(aes(yintercept=M_fit), col = cols[1], linetype = 'dashed') +
    geom_rect( aes(ymin = CL_fit, ymax = CU_fit), xmin = -Inf, xmax = Inf, fill = cols[1], alpha = 0.05) +
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.9, width = 0, col = 'grey40') + 
    geom_point(aes(y = M, col = label_counterfactual)) + 
    facet_grid(loc_label~.) +
    scale_color_manual(values = cols) + 
    scale_fill_manual(values = cols) + 
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom', 
          axis.text.x = element_blank(),
          axis.title.x = element_blank()) + 
    labs(col = '', y = paste0('Predicted ',var,' COVID-19 attributable weekly deaths\namong 18+ at the end of the resurgence period'),
         fill = '', linetype = '') +
    guides(fill=guide_legend(nrow=1 + stan_data$N_COUNTERFACTUAL,byrow=TRUE, order =1), 
           col=guide_legend(nrow=1 + stan_data$N_COUNTERFACTUAL,byrow=TRUE, order =1), 
           linetype = guide_legend(order=2)) 
  ggsave(p, file = paste0(outdir, '-predicted_weekly_deaths_vaccine_coverage_', lab, 'AllAges.png'), w = 5, h = 5 + 2*(length(unique(data_res1$code))/4))
  
  
  #####################################
  prop_vac_counterfactual_df <- copy(prop_vac_counterfactual)
  
  label <- function(age) paste0('Counterfactual analysis with a change in the\nvaccine coverage among individuals aged ', age)
  
  setnames(prop_vac_counterfactual_df, 'age_index', 'age_index_counterfactual')
  df <- copy(df_age_vaccination2[, .(age, age_index)])
  setnames(df, 1:2, c('age_counterfactual', 'age_index_counterfactual'))
  prop_vac_counterfactual_df <- merge(prop_vac_counterfactual_df, df, by = 'age_index_counterfactual')
  df <- unique(prop_vac_counterfactual_df[diff_value == 0, .(value_true, value, diff_value, age_counterfactual, age_index_counterfactual, state_index)])
  df[, counterfactual_index:=0]
  prop_vac_counterfactual_df <- rbind(prop_vac_counterfactual_df, df)
  
  tmp1 <- tmp[, list(max_date = max(date)), by = 'code']
  tmp1 <- merge(tmp, tmp1, by = 'code')
  tmp1 <- tmp1[date == max_date]
  
  df <- rbind(df_counterfactual, data.table(counterfactual_index = 0, label_counterfactual = label_fit))
  tmp1 <- merge(tmp1, df, by = 'label_counterfactual')
  tmp1 <- merge(tmp1, prop_vac_counterfactual_df, by = c('state_index', 'counterfactual_index'))
  
  tmp1[, age_counterfactual2 := gsub('.* aged (.+)', '\\1', label_counterfactual)]
  tmp1[, label_age_counterfactual := label(age_counterfactual2)]
  tmp1[label_counterfactual == label_fit, label_age_counterfactual := label_fit]
  
  tmp2 <- tmp1[age_counterfactual2 != '18-64 and 65+']
  tmp2 <- tmp2[age_counterfactual == age_counterfactual2 | label_counterfactual == label_fit]
  
  tmp2[, label_counterfactual := factor(label_counterfactual, 
                                        levels = c(label_fit, 
                                                   df_counterfactual[!grepl('18-64 and 65+', label_counterfactual), label_counterfactual]))]
  tmp2[, label_age_counterfactual := factor(label_age_counterfactual, 
                                            levels = c(label_fit, 
                                                       label('18-64'), label('65+')))]
  
  
  tmp2 <- tmp2[diff_value != 0 | label_counterfactual == label_fit]
  # df <- subset(tmp2, label_counterfactual == label_fit)
  # setnames(df, c('M', 'CL', 'CU'), c('M_fit', 'CL_fit', "CU_fit"))
  # tmp2 <- merge(tmp2, unique(df[, .(M_fit, CL_fit, CU_fit, code, age_index)]), by = c('code', 'age_index'))
  # 
  # 

}

plot_vaccine_effects_counterfactual_FL <- function(data_res1, data_res2, resurgence_dates, age_group){
  
  label_fit <- 'Fit to observed data'
  data_res2[, label_counterfactual := label_fit]
  
  data_res1 = merge(data_res1, select(resurgence_dates, code, start_resurgence), by = 'code')
  data_res1 = data_res1[date >= start_resurgence ]
  data_res1 = select(data_res1, -start_resurgence, -counterfactual_index)
  
  tmp = rbind(data_res1, data_res2)

  tmp = subset(tmp, date >= as.Date('2021-01-01'))
  tmp = merge(tmp, select(resurgence_dates, code, stop_resurgence), by = 'code')
  tmp = tmp[date <= stop_resurgence]
  tmp[, label_counterfactual := factor(label_counterfactual, levels = c(label_fit, levels(df_counterfactual$label_counterfactual)))]
  
  dummy.dt = merge(resurgence_dates, df_state, by = 'code')
  dummy.dt[, text := 'Beginning of Summer 2021 resurgence period']
  
  # values_col = c('grey50', 'darkorchid4', )
  
  tmp <- subset(tmp, code == "FL" & !grepl('and', label_counterfactual))
  lims <- tmp[, range(c(CL, CU))]
  tmp <- tmp[age == age_group & grepl(paste0(age_group, "|Fit"), label_counterfactual)]

  dummy.dt1 <- subset(dummy.dt, code == 'FL')
  
  colors <-  ggsci::pal_nejm(palette = c("default"), alpha = 1)(n = 8)[c(1, 2, 3, 4)]
  if(age_group == '65+'){
    colors <- colors[3:4]
  }

  p_FL <- ggplot(tmp, aes(x = date)) + 
    geom_line(aes(y = M, col = label_counterfactual)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = label_counterfactual), alpha = 0.15) + 
    scale_x_date(expand = expansion(mult = c(0.05,0)), date_labels = c("%b-%y"), breaks = '1 month') + 
    theme_bw() + 
    geom_vline(data = dummy.dt1, aes(xintercept = start_resurgence, linetype = text), col = 'grey50') +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'none', 
          strip.text = element_text(size= rel(1)),
          axis.text.y = element_text(size= rel(1)),
          axis.text.x = element_text(angle = 70, hjust = 1),
          axis.title.x = element_blank()
          # plot.margin = unit(c(5.5,5.5,110,5.5), "pt")
          ) + 
    scale_colour_manual(values = c('black', colors))  + 
    scale_fill_manual(values = c('black', colors)) +
    scale_linetype_manual(values = 'dashed') +
    scale_y_continuous(limits = lims) + 
    labs(col = '', y = paste0('Predicted COVID-19\ncumulative deaths among ', age_group),
         fill = '', linetype = '')+ 
    facet_grid(.~loc_label)

  return(p_FL)
}

plot_vaccine_effects_counterfactual_change_allages <- function(data_res, prop_vac_counterfactual, lab, namevar, outdir, yintercept = 0){
  
  #last data
  tmp1 <- data_res[, list(date = max(date)), by = 'code']
  tmp1 <- merge(data_res, tmp1, by = c('code', 'date'))
  
  # merge to get propo vac diff
  tmp1 <- merge(tmp1, prop_vac_counterfactual, by = c('state_index', 'counterfactual_index'))
  
  # dont keep double change
  tmp1 <- tmp1[!grepl('and', label_counterfactual)]
  
  # add one label
  label_fit <- 'Fit to observed data'
  levels_wo_and <- levels(df_counterfactual$label_counterfactual)[!grepl('and', levels(df_counterfactual$label_counterfactual))]
  tmp1[, label_counterfactual := factor(label_counterfactual, levels = c(label_fit, levels_wo_and))]
  
  colors <-  ggsci::pal_npg(palette = c("nrc"), alpha = 1)(n = 4)
  
  # keep only state wiht diff 
  tmp1 <- tmp1[diff_value != 0]
  
  p <- ggplot(tmp1, aes(x = diff_value)) + 
    geom_hline(aes(yintercept=yintercept), linetype = 'dashed', col = 'grey70') +
    geom_vline(aes(xintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.9, width = 0, col = 'grey40') + 
    geom_point(aes(y = M, col = label_counterfactual, shape = loc_label), size = 2.5) + 
    scale_colour_manual(values = c('black', colors),  drop = FALSE) + 
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom', 
          axis.title = element_text(size =rel(1.1)),
          axis.text= element_text(size =rel(1.1)),
          legend.text = element_text(size =rel(1.1))
          # legend.box="vertical",
          # legend.spacing.y = unit(-0, "cm")
    ) + 
    scale_shape_manual(values = c(15, 17, 3, 10, 11, 20, 12, 13, 14, 4)[1:length(tmp1$code)]) + 
    labs(col = '', y = paste0('Change in COVID-19 attributable deaths\namong 18+ at the end of the resurgence period'),
         fill = '', linetype = '', 
         x = 'Difference between hypothesised vaccination rate\nand observed vaccination rate', pch = '') +
    guides(fill=guide_legend(order=1,ncol=1,byrow=TRUE), 
           col=guide_legend(order=1,ncol=1,byrow=TRUE), 
           pch = guide_legend(order=2,nrow=length(unique(tmp1$code))/2,byrow=TRUE)) 
  
  if(grepl('perc', namevar)){
    p <- p + scale_y_continuous(labels = scales::percent)
    # p <- p + scale_y_log10(labels = scales::percent)
  }
  
  ggsave(p, file = paste0(outdir, '-predicted_', namevar, '_weekly_deaths_vaccine_coverage_', lab, 'AllAges.png'), w = 6.5, h = 4 + 2*6/4)
  
  return(p)
}

plot_vaccine_effects_counterfactual_change <- function(data_res, prop_vac_counterfactual, age_group){

  #last data
  tmp1 <- data_res[, list(date = max(date)), by = 'code']
  tmp1 <- merge(data_res, tmp1, by = c('code', 'date'))
  
  # merge to get propo vac diff
  tmp1 <- merge(tmp1, prop_vac_counterfactual, by = c('state_index', 'counterfactual_index'))

  # dont keep double change
  tmp1 <- tmp1[!grepl('and', label_counterfactual)]
  
  # add one label
  label_fit <- 'Fit to observed data'
  levels_wo_and <- levels(df_counterfactual$label_counterfactual)[!grepl('and', levels(df_counterfactual$label_counterfactual))]
  levels_wo_and <- levels_wo_and[grepl(paste0(age_group), levels_wo_and)]
  
  lims <- tmp1[, range(CL, CU)]
  tmp1 <- tmp1[age == age_group & grepl(paste0(age_group, "|Fit"), label_counterfactual)]
  tmp1[, label_counterfactual := factor(label_counterfactual, levels = c(label_fit, levels_wo_and))]
  
  colors <-  ggsci::pal_nejm(palette = c("default"), alpha = 1)(n = 8)[c(1, 2, 3, 4)]
  if(age_group == '65+'){
    colors = colors[3:4]
  }
  
  # keep only state wiht diff 
  tmp1 <- tmp1[diff_value != 0]

  p <- ggplot(tmp1, aes(x = diff_value)) + 
    geom_hline(aes(yintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_vline(aes(xintercept=0), linetype = 'dashed', col = 'grey70') +
    geom_errorbar(aes(ymin = CL, ymax = CU),width = 0, col = 'grey50', alpha= 0.5) + 
    geom_point(aes(y = M, x = diff_value, col = label_counterfactual), size = 0.5) +
    geom_label_repel(aes(y = M, x = diff_value, col = label_counterfactual, label = code), 
                     size = 2.5, label.size = NA, fill = NA, show.legend = FALSE) + 
    scale_colour_manual(values = c('black', colors),  drop = FALSE) + 
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom'
          # legend.box="vertical",
          # legend.spacing.y = unit(-0, "cm")
    ) + 
    labs(col = '', y = paste0('Change in predicted COVID-19\ncumulative deaths among ', age_group),
         fill = '', linetype = '', 
         x = 'Difference between counterfactual vaccination rate\nand observed vaccination rate') +
    guides(col=guide_legend(order=1, override.aes = list(size = 3), byrow = T, nrow = 3)) + 
    scale_y_continuous(labels = scales::percent, limits= lims) 

  return(list(p, lims))
}

plot_vaccine_effects_counterfactual_panel <- function(E_pdeaths_counterfactual_resurgence_cumulative, E_pdeaths_predict_resurgence_cumulative, 
                                                      perc_E_pdeaths_counterfactual, resurgence_dates, prop_vac_counterfactual, outdir){
  age_groups = c('18-64', '65+')
  p_FL <- list(); p_all <- list(); p_all_log = list()
  for(i in seq_along(age_groups)){
    age_group = age_groups[i]
    p_FL[[i]] <- plot_vaccine_effects_counterfactual_FL(E_pdeaths_counterfactual_resurgence_cumulative, E_pdeaths_predict_resurgence_cumulative, resurgence_dates, age_group) + 
      ggtitle(paste0('Counterfactual scenarios with varying\nvaccination rate in ', age_group, '\n'))+
      theme(plot.title = element_text(hjust = 0.5))
    tmp <-  plot_vaccine_effects_counterfactual_change(perc_E_pdeaths_counterfactual, prop_vac_counterfactual, age_group) 
    lims <- tmp[[2]]
    p_all[[i]] <- tmp[[1]] + coord_cartesian(ylim = c(lims[1], 5))
    p_all_log[[i]] <-  tmp[[1]]  + scale_y_continuous(labels = scales::percent, limits = lims, trans = 'pseudo_log') 
    
    if(age_group == '18-64'){
      p_FL[[i]] <- ggarrange(p_FL[[i]], labels = c('A'), label.y = 0.85  )
      p_all[[i]] <- ggarrange(p_all[[i]], labels = c('B'))
      p_all_log[[i]] <- ggarrange(p_all_log[[i]], labels = c('B'))
    }
  }
  pp <- grid.arrange(grobs = c(p_FL, p_all), layout_matrix= rbind(c(1, NA, 2, NA), c(3, 3, 4, 4)), heights = c(0.4, 0.6), widths = c(0.4, 0.1, 0.4, 0.1))
  ggsave(pp, file = paste0(outdir, '-predicted_weekly_deaths_vaccine_coverage_counterfactual_panel_plot.png'), w = 10, h = 9)
  
  pp <- grid.arrange(grobs = c(p_FL, p_all_log), layout_matrix= rbind(c(1, NA, 2, NA), c(3, 3, 4, 4)), heights = c(0.4, 0.6), widths = c(0.4, 0.1, 0.4, 0.1))
  ggsave(pp, file = paste0(outdir, '-predicted_weekly_deaths_vaccine_coverage_counterfactual_panel_plot_log.png'), w = 10, h = 9)
  
}


plot_forest_plot <- function(tmp, outdir){
  
  tmp1 <- tmp[!grepl('vacc', variable) & group == 'baseline']
  tmp1[grepl('\\["18-64', variable) | (grepl('18-64"\\]', variable) & !grepl('\\["65', variable)), age_group := 'among 18-64']
  tmp1[grepl('\\["65', variable) | (grepl('65\\+"\\]', variable) & !grepl('\\["18-64', variable)), age_group := 'among 65+']
  p1 <- ggplot(tmp1, aes(y = variable)) + 
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
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("Non-vaccine parameters on baseline")
  
  tmp1 <- tmp[!grepl('vacc', variable) & group == 'slope']
  p2 <- ggplot(tmp1, aes(y = variable)) + 
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
          plot.title = element_text(hjust = 0.5)) +
    ggtitle("Non-vaccine parameters on slope")
  
  tmp1 <- tmp[grepl('vacc', variable) & group == 'baseline']
  tmp1[, type := 'Direct\nvaccine\neffects']
  tmp1[grepl('vacc-cross', variable), type := 'Indirect\nvaccine\neffects']
  p3 <- ggplot(tmp1, aes(y = variable)) + 
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
    ggtitle("Vaccine parameters on baseline")
  
  if(length(unique(tmp1$variable))>2){
    p3 <- p3 + facet_grid(type~., scales= "free", space="free") 
  }
  
  tmp1 <- tmp[grepl('vacc', variable) & group == 'slope']
  tmp1[, type := 'Direct\nvaccine\neffects']
  tmp1[grepl('vacc-cross', variable), type := 'Indirect\nvaccine\neffects']
  p4 <- ggplot(tmp1, aes(y = variable)) + 
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
    ggtitle("Vaccine parameters on slope")
  
  if(length(unique(tmp1$variable))>2){
    p4 <- p4 + facet_grid(type~., scales= "free", space="free") 
  }
  
  p <- grid.arrange(p1, p2, p3, p4, layout_matrix = rbind(c(1, 1, 3), c(1, 1, 4), c(1, 1, NA), c(NA, 2, NA)), heights = c(0.25, 0.25, 0.35, 0.15), widths = c(0.07, 0.43, 0.5))
  ggsave(p, file = paste0(outdir, '-forest_plot.png'), w = 8, h = 7)
  
  return(p)
}

plot_contribution_vaccine <- function(contribution, vaccine_data_pop, lab, outdir){
  
  # delay = 7*2
  tmp <- vaccine_data_pop[, list(mindate = min(date[prop > 0.05])), by = 'code']
  tmp[, type:= '']
  tmp1 <- merge(tmp, unique(contribution[, .(code, loc_label)]), by = 'code')
  tmp1 <- tmp1[code %in% unique(contribution$code)]

  p1 = ggplot(contribution, aes(x = date)) + 
    geom_vline(data = tmp1, aes(xintercept = mindate, alpha = type), col = 'grey50', size=1) +
    geom_line(aes(y = M_median, col = as.factor(age), linetype = '')) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.5) + 
    geom_line(aes(y = M, col = as.factor(age))) + 
    scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0,1)) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.direction='vertical') +
    labs(y = paste0("Estimated contribution to COVID-19 weekly deaths"), 
         col = "Age group", fill = "Age group", x = '', linetype = 'National median', alpha = 'Start of vaccination\ncampaign') + 
    scale_shape_manual(values = 4)+ 
    scale_alpha_manual(values = c(1, 0.7)) + 
    scale_linetype_manual(values = 'dashed') + 
    scale_x_date(expand = c(0,0), date_labels = "%b-%y")  + 
    guides(col = guide_legend(order = 1), fill = guide_legend(order = 1), 
           linetype = guide_legend(order = 2),
          alpha = guide_legend(order = 3))  +
    ggsci::scale_color_jco() + 
    ggsci::scale_fill_jco()
  
  if(length(unique(contribution$code)) > 10){
    p1 <- p1 +     facet_wrap(~loc_label, ncol= 5) + 
      theme(axis.text.x = element_text(angle = 70, hjust = 1))
    ggsave(p1, file = paste0(outdir, '-contribution_vaccine_coverage_', lab, '.png'), w = 9, h = 12)
    ggsave(p1, file = paste0(outdir, '-contribution_vaccine_coverage_', lab, '.pdf'), w = 9, h = 12)
    
  }else if(length(unique(contribution$code)) > 3){
    p1 <- p1 +  facet_wrap(~loc_label, nrow = length(unique(contribution$code)) )
    ggsave(p1, file = paste0(outdir, '-contribution_vaccine_coverage_', lab, '.png'), w = 6, h = 9)
    
  }else{
    
    p1 <- p1 +  facet_wrap(~loc_label, nrow = length(unique(contribution$code))) + 
      labs(y = paste0("Estimated contribution to\nCOVID-19 weekly deaths")) +
      theme(axis.text.x = element_text(angle = 70, hjust = 1)) 
    ggsave(p1, file = paste0(outdir, '-contribution_vaccine_coverage_', lab, '.png'), w = 5, h = 5)
    ggsave(p1, file = paste0(outdir, '-contribution_vaccine_coverage_', lab, '.pdf'), w = 5, h = 5)
    
  }
  

}


plot_relative_mortality_all_states <- function(mortality_rateJ21, nyt_data, outdir){
  
  tmp <- merge(mortality_rateJ21, nyt_data[, .(STATE, SHARE_DEATHS)], by.x = 'loc_label', by.y = 'STATE')
  tmp <- tmp[!is.na(SHARE_DEATHS)]
  
  tmp <- tmp[!age %in% c('0-24', '25-54', '85+')]
  tmp[, `Age group`:=age]
  
  p <- ggplot(tmp, aes(y = M_rel, x = SHARE_DEATHS, col = loc_label)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = CL_rel, ymax = CU_rel)) +
    labs(x = "Share of COVID-19 deaths linked to long term care facilities", 
         y = 'Predicted COVID-19 attributable\nmortality rate ratio in 85+ relative to 55-84', 
         col = '') + 
    # facet_grid(.~`Age group`, label = 'label_both') + 
    theme_bw() + 
    theme(legend.position = 'bottom', 
          strip.text = element_text(size = rel(0.8)),
          axis.title = element_text(size = rel(0.85)),
          axis.text = element_text(size = rel(0.8)),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) + 
    scale_color_viridis_d(option = 'C', end = 0.9) + 
    geom_text(aes(label = code, y = CU_rel + 0.5, col = loc_label), size = 2.5) + 
    guides(col = 'none') + 
    scale_x_continuous(labels = scales::percent)
  ggsave(p, file = paste0(outdir, paste0('-MortalityRateRelative_CareHome.png')), w = 7, h = 5)
  
  return(p)
}

