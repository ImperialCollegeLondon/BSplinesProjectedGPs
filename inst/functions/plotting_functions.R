compare_CDC_JHU_DoH_error_plot = function(CDC_data, JHUData, scrapedData, var.weekly.deaths.CDC, outdir, Code = NULL)
{
  # find errors 
  CDCdata = CDC_data[, list(cumulative_deaths.CDC = sum(na.omit( get(var.weekly.deaths.CDC) ))), by = c('code', 'date')]
  
  # take cumulative deaths
  CDCdata[, cumulative_deaths.CDC := cumsum(cumulative_deaths.CDC), by = c('code', 'code')]
  
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  JHUData[, date := as.Date(date)]
  JHUData = subset(JHUData, date <= max(CDCdata$date))
  
  scrapedData = select(as.data.table(scrapedData), code, date, cum.deaths, age)
  scrapedData = scrapedData[, list(cum.deaths = sum(cum.deaths)), by = c('code', 'date')]
  scrapedData[, date := as.Date(date)]
  scrapedData = subset(scrapedData, date <= max(CDCdata$date))

  # plot
  JHUData[, source := 'JHU']
  CDCdata[, source := 'CDC']
  scrapedData[, source := 'DoH']
  setnames(CDCdata, 'cumulative_deaths.CDC', 'cumulative_deaths')
  setnames(scrapedData, 'cum.deaths', 'cumulative_deaths')
  
  tmp2 = rbind(rbind(JHUData, CDCdata),scrapedData)

  if(!dir.exists( basename(outdir) )){
    dir.create( basename(outdir), recursive = T )
  }
  
  n_code = length(unique(tmp2$code))
  
  p = ggplot(tmp2, aes(x = date, y = cumulative_deaths, col = source)) + 
    geom_line() +
    facet_wrap(~code, nrow = length(unique(tmp2$code)), scale = 'free') + 
    theme_bw()  +
    labs(x = '', y = 'Cumulative COVID-19 deaths', col = '') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90)) 
  ggsave(p, file = paste0(outdir, '-comparison_JHU_DoH_CDC.pdf'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
  
  if(!is.null(Code)){
    tmp = subset(tmp2, code %in% Code)
    
    p = ggplot(tmp, aes(x = date, y = cumulative_deaths, col = source)) + 
      geom_line() +
      theme_bw() + 
      facet_wrap(~code) +
      scale_color_viridis_d(option = "B", direction = -1, end = 0.8) +
      labs(x = '', y = 'Cumulative COVID-19 \ndeaths', col = '') + 
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
      theme(legend.position = 'bottom',
            axis.text.x = element_text(angle = 90)) 
    ggsave(p, file = paste0(outdir, '-comparison_JHU_DoH_CDC_Code.pdf'), w = 7, h = 6, limitsize = F)
    
  }
}

plot_data = function(deathByAge, outdir, Code = NULL)
{
  
  p = ggplot(deathByAge, aes(x = date, y = age)) + 
    geom_raster(aes(fill = weekly.deaths )) + 
    facet_wrap(~loc_label,ncol = 6) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'pseudo_log', breaks = c(10, 100,1000)) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0), breaks = unique(deathByAge$age)[rep(c(T,F), length(unique(deathByAge$age))/2)]) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.key = element_blank(), 
          strip.background = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = 'Retrievable weekly\nCOVID-19 attributable deaths')
  ggsave(p, file = paste0(outdir, '-deathByAge.png'), w = 10, h = 12)
  
  df = data.frame(date = min(deathByAge$date), age = deathByAge$age[1],
                  dummy = c('Missing\nweekly COVID-19\nattributable deaths', 'Non-retrievable\nweekly COVID-19\nattributable deaths'))
  

  tmp = subset(deathByAge, code %in%  Code)
  p <- ggplot(tmp, aes(x = date, y = age)) +
    geom_raster(aes(fill = weekly.deaths )) +
    theme_bw() +
    facet_wrap(~loc_label,ncol = 4) + 
    scale_fill_viridis_c(trans = 'sqrt',  na.value="grey70", breaks = c(0, 100, 1000,2500),) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    scale_y_discrete(expand = c(0,0)) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank(),
          # axis.text.y =  element_text(size = rel(0.8)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(),
          legend.title = element_text(size = rel(0.85)),
          # legend.text = element_text(size = rel(1)),
          strip.background = element_blank(),
          panel.background = element_rect(fill = '#FCB360', colour = 'red')) +
    geom_point(data = df, aes(color = dummy), shape = 15, size = 0) +
    labs(x = '', y = 'Age group', fill = 'Retrievable\nweekly COVID-19\nattributable deaths', col = '') +
    scale_color_manual(values = c('#FCB360', "grey70")) +
    guides(color = guide_legend(override.aes = list(size=4)), 
           fill = guide_colorbar(title.position = "right")) 
  ggsave(p, file = paste0(outdir, '-deathByAge_selected_all_states.png'), w = 8.5, h = 3*(length(Code)/4))
  
  tmp = subset(deathByAge, code %in%  c('CA', 'FL', 'NY', 'TX'))
  p <- ggplot(tmp, aes(x = date, y = age)) +
    geom_raster(aes(fill = weekly.deaths )) +
    theme_bw() +
    facet_wrap(.~loc_label, nrow = 2) + 
    scale_fill_viridis_c(trans = 'sqrt',  na.value="grey70", breaks = c(0, 100, 1000,2500),) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    scale_y_discrete(expand = c(0,0)) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)),
          axis.text.y = element_text(size = rel(1)),
          axis.title.x = element_blank(),
          legend.justification = c(0,1),
          # axis.text.y =  element_text(size = rel(0.8)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_blank(),
          legend.title = element_text(size = rel(0.85)),
          axis.title.y = element_text(size = rel(1.1)),
          strip.text = element_text(size = rel(0.9)),
          # legend.text = element_text(size = rel(1)),
          strip.background = element_blank(),
          panel.background = element_rect(fill = '#FCB360', colour = 'red')) +
    geom_point(data = df, aes(color = dummy), shape = 15, size = 0) +
    labs(x = '', y = 'Age group', fill = 'Retrievable\nweekly COVID-19\nattributable deaths', col = '') +
    scale_color_manual(values = c('#FCB360', "grey70")) +
    guides(color = guide_legend(override.aes = list(size=4)),
           fill = guide_colorbar(title.position = "right", barheight = 0.7))
  ggsave(p, file = paste0(outdir, '-deathByAge_selected_states.png'), w = 7, h = 5)
  
  p1 = ggplot(deathByAge, aes(x = date, y = age)) + 
    geom_raster(aes(fill = min.sum.weekly.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30),
                         limits = c(0,max(c(na.omit(deathByAge$sum.weekly.deaths),na.omit(deathByAge$max.sum.weekly.deaths))))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = "Value sum of covid-19 weekly deaths", title = 'lower bound')
  
  p2 = ggplot(deathByAge, aes(x = date, y = age)) + 
    geom_raster(aes(fill = max.sum.weekly.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30), 
                         limits = c(0,max(c(na.omit(deathByAge$sum.weekly.deaths),na.omit(deathByAge$max.sum.weekly.deaths))))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = "Value sum of covid-19 weekly deaths", title = 'upper bound')
  
  p3 = ggplot(deathByAge, aes(x = date, y = age)) + 
    geom_raster(aes(fill = sum.weekly.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30), 
                         limits = c(0,max(c(na.omit(deathByAge$sum.weekly.deaths), na.omit(deathByAge$max.sum.weekly.deaths))))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = "Value sum of covid-19 weekly deaths", title = 'exact sum')
  
  p_all = ggpubr::ggarrange(p1, p2, p3,  common.legend = T, legend = 'bottom')
  ggsave(p_all, file = paste0(outdir, '-deathByAge_boundaries.png'), w = 20, h = 25)
  
  return(p)
}
