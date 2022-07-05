compare_cumulative_CDC_JHU_DoH_error_plot = function(deathByAge, JHUData, scrapedData, outdir)
{
  
  # find initial cumulative deatha 
  file = file.path(dirname(indir), paste0('misc/data/CDC_data_2020-05-06.csv'))
  tmp <- as.data.table(read.csv(file))
  tmp <- tmp[grepl("Total", State)]
  tmp[, loc_label := gsub('(.+) Total', '\\1', State)]
  tmp1 <- tmp[loc_label %in% deathByAge[!(loc_label %in% 'New York'), unique(loc_label)]]
  tmp1 <- tmp1[, .(loc_label, COVID.19.Deaths)]
  tmp <- tmp[loc_label %in% c('New York', 'New York City'), list(COVID.19.Deaths = sum(COVID.19.Deaths), loc_label = 'New York')]
  tmp <- rbind(tmp,tmp1)
  
  # aggregate deaths by age across ages and add initial cumulative deaths
  CDCdata = deathByAge[, list(weekly_deaths.CDC = sum(na.omit( weekly.deaths ))), by = c('loc_label', 'date')]
  CDCdata = CDCdata[order(loc_label, date)]
  CDCdata[, cumulative_deaths := cumsum(weekly_deaths.CDC), by = c('loc_label')]
  CDCdata <- merge(tmp, CDCdata, by = 'loc_label')
  CDCdata[, cumulative_deaths := cumulative_deaths + COVID.19.Deaths]
  CDCdata <- CDCdata[, .(loc_label, date, cumulative_deaths)]
  CDCdata[, source := 'Reported by the CDC']

  # aggregate JHU by weeks
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  JHUData[, date := as.Date(date)]
  JHUData = subset(JHUData, date <= max(CDCdata$date))
  JHUData <- merge(JHUData, unique(deathByAge[, .(code, loc_label)]), by = 'code')
  JHUData <- JHUData[, .(loc_label, date, cumulative_deaths)]
  JHUData[, source := 'Reported by JHU']
  
  # prepare scraped data
  scrapedData = select(as.data.table(scrapedData), code, date, cum.deaths, age)
  scrapedData = scrapedData[, list(cumulative_deaths = sum(cum.deaths)), by = c('code', 'date')]
  scrapedData[, date := as.Date(date)]
  scrapedData = subset(scrapedData, date <= max(CDCdata$date))
  scrapedData <- merge(scrapedData, unique(deathByAge[, .(code, loc_label)]), by = 'code')
  scrapedData <- scrapedData[, .(loc_label, date, cumulative_deaths)]
  scrapedData[, source := 'DoH']
  
  # merge
  tmp2 = rbind(rbind(JHUData, CDCdata),scrapedData)

  if(!dir.exists( basename(outdir) )){
    dir.create( basename(outdir), recursive = T )
  }
  
  n_code = length(unique(tmp2$loc_label))
  
  p1 = ggplot(tmp2, aes(x = date, y = cumulative_deaths, col = source)) + 
    geom_line() +
    facet_wrap(~loc_label, nrow = length(unique(tmp2$loc_label)), scale = 'free') + 
    theme_bw()  +
    labs(x = '', y = 'Cumulative COVID-19 deaths', col = '') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90)) 
  ggsave(p1, file = paste0(outdir, '-comparison_JHU_DoH_CDC.pdf'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
  
  cols <- pal_npg("nrc")(10)
  tmp = subset(tmp2, loc_label %in% deathByAge[code %in% c('CA', 'FL', 'TX', 'NY'), unique(loc_label)])
  tmp <- tmp[source != 'DoH']
  tmp <- tmp[date >= min(CDCdata$date)]
  tmp[, source := factor(source, levels = c('Reported by the CDC', 'Reported by JHU'))]
  p2 = ggplot(tmp, aes(x = date, y = cumulative_deaths, col = source)) + 
    geom_line() +
    theme_bw() + 
    facet_wrap(~loc_label, scale = 'free_y', ncol = 1) +
    scale_color_viridis_d(option = "B", direction = 1, end = 0.8) +
    # scale_color_manual(values = cols[5:6]) +
    labs( y = 'Cumulative all-ages COVID-19 deaths', col = '') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 45, vjust =0.5), 
          axis.title.x = element_blank(), 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          panel.grid.minor = element_blank(), 
          legend.title = element_text(size = rel(0.85)),
          axis.title.y = element_text(size = rel(1.1)),
          axis.text.y = element_text(size = rel(1)),
          strip.text = element_text(size = rel(0.9)),) 
  ggsave(p2, file = paste0(outdir, '-comparison_JHU_DoH_CDC_Code.pdf'), w = 7, h = 6, limitsize = F)
  
  p <- p  + facet_wrap(.~loc_label, nrow = 4) + theme(legend.box = 'vertical')
  p <- ggarrange(p, labels = 'A')
  p2 <- ggarrange(p2, labels = 'B')
  pa <- grid.arrange(p, p2, ncol = 2, layout_matrix = rbind(c(1, 2), c(1, NA)), heights = c(0.91, 0.09))
  ggsave(pa, file = paste0(outdir, '-deathByAge_CDCJHU_panel_selected_states.png'), w = 7.7, h = 8.5)
  
}

compare_weekly_CDC_JHU_DoH_error_plot = function(deathByAge, JHUData, scrapedData, outdir)
{
  
  # find initial cumulative deatha 
  file = file.path(dirname(indir), paste0('misc/data/CDC_data_2020-05-06.csv'))
  tmp <- as.data.table(read.csv(file))
  tmp <- tmp[grepl("Total", State)]
  tmp[, loc_label := gsub('(.+) Total', '\\1', State)]
  tmp1 <- tmp[loc_label %in% deathByAge[!(loc_label %in% 'New York'), unique(loc_label)]]
  tmp1 <- tmp1[, .(loc_label, COVID.19.Deaths)]
  tmp <- tmp[loc_label %in% c('New York', 'New York City'), list(COVID.19.Deaths = sum(COVID.19.Deaths), loc_label = 'New York')]
  tmp <- rbind(tmp,tmp1)
  
  # aggregate deaths by age across ages and add initial cumulative deaths
  CDCdata = deathByAge[, list(weekly_deaths.CDC = sum(na.omit( weekly.deaths ))), by = c('loc_label', 'date')]
  CDCdata = CDCdata[order(loc_label, date)]
  CDCdata[, cumulative_deaths := cumsum(weekly_deaths.CDC), by = c('loc_label')]
  CDCdata <- merge(tmp, CDCdata, by = 'loc_label')
  CDCdata[, cumulative_deaths := cumulative_deaths + COVID.19.Deaths]
  CDCdata <- CDCdata[, .(loc_label, date, cumulative_deaths)]
  CDCdata[, source := 'Reported by the CDC']
  
  # aggregate JHU by weeks
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  JHUData[, date := as.Date(date)]
  JHUData = subset(JHUData, date <= max(CDCdata$date))
  JHUData <- merge(JHUData, unique(deathByAge[, .(code, loc_label)]), by = 'code')
  JHUData <- JHUData[, .(loc_label, date, cumulative_deaths)]
  JHUData[, source := 'Reported by JHU']
  
  # prepare scraped data
  scrapedData = select(as.data.table(scrapedData), code, date, cum.deaths, age)
  scrapedData = scrapedData[, list(cumulative_deaths = sum(cum.deaths)), by = c('code', 'date')]
  scrapedData[, date := as.Date(date)]
  scrapedData = subset(scrapedData, date <= max(CDCdata$date))
  scrapedData <- merge(scrapedData, unique(deathByAge[, .(code, loc_label)]), by = 'code')
  scrapedData <- scrapedData[, .(loc_label, date, cumulative_deaths)]
  scrapedData[, source := 'DoH']
  
  # merge
  tmp2 = rbind(rbind(JHUData, CDCdata),scrapedData)
  
  if(!dir.exists( basename(outdir) )){
    dir.create( basename(outdir), recursive = T )
  }
  
  n_code = length(unique(tmp2$loc_label))
  
  p1 = ggplot(tmp2, aes(x = date, y = cumulative_deaths, col = source)) + 
    geom_line() +
    facet_wrap(~loc_label, nrow = length(unique(tmp2$loc_label)), scale = 'free') + 
    theme_bw()  +
    labs(x = '', y = 'Cumulative COVID-19 deaths', col = '') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90)) 
  ggsave(p1, file = paste0(outdir, '-comparison_JHU_DoH_CDC.pdf'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
  
  cols <- pal_npg("nrc")(10)
  tmp = subset(tmp2, loc_label %in% deathByAge[code %in% c('CA', 'FL', 'TX', 'NY'), unique(loc_label)])
  tmp <- tmp[source != 'DoH']
  tmp <- tmp[date >= min(CDCdata$date)]
  tmp[, source := factor(source, levels = c('Reported by the CDC', 'Reported by JHU'))]
  p2 = ggplot(tmp, aes(x = date, y = cumulative_deaths, col = source)) + 
    geom_line() +
    theme_bw() + 
    facet_wrap(~loc_label, scale = 'free_y', ncol = 1) +
    scale_color_viridis_d(option = "B", direction = 1, end = 0.8) +
    # scale_color_manual(values = cols[5:6]) +
    labs( y = 'Cumulative all-ages COVID-19 deaths', col = '') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 45, vjust =0.5), 
          axis.title.x = element_blank(), 
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          panel.grid.minor = element_blank(), 
          legend.title = element_text(size = rel(0.85)),
          axis.title.y = element_text(size = rel(1.1)),
          axis.text.y = element_text(size = rel(1)),
          strip.text = element_text(size = rel(0.9)),) 
  ggsave(p2, file = paste0(outdir, '-comparison_JHU_DoH_CDC_Code.pdf'), w = 7, h = 6, limitsize = F)
  
  p <- p  + facet_wrap(.~loc_label, nrow = 4) + theme(legend.box = 'vertical')
  p <- ggarrange(p, labels = 'A')
  p2 <- ggarrange(p2, labels = 'B')
  pa <- grid.arrange(p, p2, ncol = 2, layout_matrix = rbind(c(1, 2), c(1, NA)), heights = c(0.91, 0.09))
  ggsave(pa, file = paste0(outdir, '-deathByAge_CDCJHU_panel_selected_states.png'), w = 7.7, h = 8.5)
  
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


plot_vaccine_data = function(deathByAge, vaccine_data, pop_data, Code, outdir){
  
  tmp = merge(vaccine_data, select(df_age_continuous, -age_index), by = 'age')
  tmp[, age_index := which(df_age_vaccination$age_from <= age_from & df_age_vaccination$age_to >= age_to), by = c('age_from', 'age_to')]
  tmp = merge(select(tmp, -age_from, -age_to, -age), df_age_vaccination, by = 'age_index')
  tmp = tmp[, list(prop = unique(prop), pop = sum(pop)), by = c('age', 'code', 'date', 'loc_label', 'age_index')]
  
  ggplot(tmp, aes(date, prop)) +
    geom_line(aes(col = age)) + 
    facet_wrap(~loc_label) + 
    theme_bw()+ 
    labs(x = '', y = 'Proportion of vaccinated', col = 'Age group') + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white")) 
  ggsave(paste0(outdir, '-proportion_vaccine_age_code.png'), w = 9, h = 8)
  
  tmp[, `Age group` := age]
  p2 <- ggplot(subset(tmp, code %in% c('CA', 'FL', 'NY', 'TX') & age_index >2), aes(date, prop)) +
    geom_line(aes(col = loc_label)) + 
    facet_wrap(~`Age group`, label = 'label_both') + 
    theme_bw()+ 
    labs( y = 'Proportion of fully vaccinated individuals', col = '') + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom',
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          axis.title.x = element_blank()) + 
    scale_color_jcolors('pal8') + 
    scale_y_continuous(labels = scales::percent)+ 
    scale_x_date(date_labels = c("%b-%y"), breaks = '2 months') 
  ggsave(paste0(outdir, '-proportion_vaccine_age_code_selected_states.png'), w = 6, h = 3.75)
  
  p <- p  + facet_wrap(.~loc_label, nrow = 4) 
  p <- ggarrange(p, labels = 'A')
  p2 <- p2 + facet_wrap(.~`Age group`, nrow = 2, label = 'label_both') 
  p2 <- ggarrange(p2, labels = 'B')
  pa <- grid.arrange(p, p2, layout_matrix = rbind(c(1, 2), c(1, NA)), heights = c(0.75, 0.25))
  ggsave(pa, file = paste0(outdir, '-deathByAge_proportion_vaccine_panel_selected_states.png'), w = 7.7, h = 8.5)
  
  
  if(length(Code) > 6){
    
    mid_code = length(Code) / 2
    
    p1 <- ggplot(subset(tmp, code %in% Code[1:mid_code] & age_index >2), aes(date, prop)) +
      geom_hline(yintercept = 0.5, linetype = 'dashed', col = 'grey50') + 
      geom_line(aes(col = loc_label)) + 
      theme_bw()+ 
      labs( y = '', col = '') + 
      theme(legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"),
            panel.border = element_rect(colour = "black", fill = NA), 
            legend.position = 'bottom',
            axis.text.x = element_blank(), 
            axis.title.x = element_blank()) + 
      scale_color_jcolors('pal8') + 
      facet_grid(`Age group`~loc_label) +
      scale_y_continuous(labels = scales::percent)+ 
      scale_x_date(date_labels = c("%b-%y"), breaks = '2 months') 
    
    
    p2 <- ggplot(subset(tmp, code %in% Code[(mid_code + 1):(mid_code*2)] & age_index >2), aes(date, prop)) +
      geom_hline(yintercept = 0.5, linetype = 'dashed', col = 'grey50') + 
      geom_line(aes(col = loc_label)) + 
      theme_bw()+ 
      labs( y = '', col = '') + 
      theme(legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"),
            panel.border = element_rect(colour = "black", fill = NA), 
            legend.position = 'bottom',
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
            axis.title.x = element_blank()) + 
      scale_color_jcolors('pal8') + 
      facet_grid(`Age group`~loc_label) +
      scale_y_continuous(labels = scales::percent)+ 
      scale_x_date(date_labels = c("%b-%y"), breaks = '2 months') 
    
    p <- ggarrange(p1, p2, nrow = 2,  common.legend = T, legend = 'bottom')
    grid.arrange(p, left = 'Proportion of fully vaccinated individuals')
    ggsave(paste0(outdir, '-proportion_vaccine_age_code_all_states.png'), w = 8, h = 6)
    
    
  } else{
    ggplot(subset(tmp, code %in% Code & age_index >2), aes(date, prop)) +
      geom_line(aes(col = loc_label)) + 
      theme_bw()+ 
      labs( y = 'Proportion of fully vaccinated individuals', col = '') + 
      theme(legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"),
            panel.border = element_rect(colour = "black", fill = NA), 
            legend.position = 'bottom',
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
            axis.title.x = element_blank()) + 
      scale_color_jcolors('pal8') + 
      facet_grid(`Age group`~loc_label) +
      scale_y_continuous(labels = scales::percent)+ 
      scale_x_date(date_labels = c("%b-%y"), breaks = '2 months') 
    ggsave(paste0(outdir, '-proportion_vaccine_age_code_all_states.png'), w = 8, h = 3)
  }
  
  
  # population count
  pop_data[, age_index := which(df_age_vaccination$age_from <= age & df_age_vaccination$age_to >= age), by = 'age']
  df_pop_data = pop_data[, list(pop_adj = sum(pop)), by = c('code', 'age_index')]
  
  # deaths
  df_age_close_vaccination = copy(df_age_reporting)
  df_age_close_vaccination[, age_index :=  c(1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4)]
  delay = 7 * 2
  tmp1 = merge(deathByAge, df_age_close_vaccination, by = c('age_from', 'age_to', 'age'))
  tmp1 = tmp1[, list(weekly.deaths = sum(na.omit(weekly.deaths))), by = c('code', 'date', 'loc_label', 'age_index')]
  tmp2 = tmp1[, list(total_deaths = sum(na.omit(weekly.deaths))), by = c('date', 'code')]
  tmp1 = merge(tmp1, tmp2, by = c('date', 'code'))
  tmp1 = merge(tmp1, df_pop_data, by = c('code', 'age_index'))
  
  tmp1[, prop_deaths := weekly.deaths / total_deaths]
  tmp1[, date := date - delay ]
  
  tmp = merge(tmp, tmp1, c('code', 'date', 'loc_label', 'age_index'))
  tmp = subset(tmp, code %in% Code )
  # summary(lm(prop_deaths ~ prop*age  + age, data = subset(tmp, code == 'TX')))
  
  tmp2 = tmp[, list(max_date = max(date[prop == 0])), by = c('age', 'loc_label')]
  tmp = merge(tmp, tmp2, by = c('age', 'loc_label'))
  tmp = tmp[date >= max_date,]
  
  tmp[, avg_pop := mean(pop_adj), by = 'age_index']
  tmp[, weekly.deaths_adj := weekly.deaths / pop_adj * avg_pop]
  
  
  ggplot(tmp, aes(x = prop, y = prop_deaths, col = date)) + 
    geom_smooth( col = 'black', size = 0.5) + 
    geom_point() + 
    facet_grid(loc_label~age) + 
    scale_color_viridis_c(trans = "date") + 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white")) +
    labs(x = 'Proportion of vaccinated', y = 'Contribution to weekly death 2 weeks later' )
  ggsave(paste0(outdir, '-proportion_vaccine_contribution_deaths.png'), w = 8, h = 6)
  
  ggplot(tmp, aes(x = prop, y = weekly.deaths, col = date)) +
    geom_smooth(col = 'black', size = 0.5) + 
    geom_point() + 
    facet_grid(loc_label~age) + 
    labs(x = 'Proportion of vaccinated', y = 'Weekly deaths 2 weeks later' ) +
    # scale_y_log10() + 
    scale_color_viridis_c(trans = "date") + 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white")) 
  ggsave(paste0(outdir, '-proportion_vaccine_abs_deaths.png'), w = 8, h = 6)
  
  ggplot(tmp, aes(x = prop, y = weekly.deaths_adj, col = date)) +
    geom_smooth(col = 'black', size = 0.5) + 
    geom_point() + 
    facet_grid(loc_label~age) + 
    labs(x = 'Proportion of vaccinated', y = 'Weekly deaths 2 weeks later\nadjusted for population count' ) +
    # scale_y_log10() + 
    scale_color_viridis_c(trans = "date") + 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white")) 
  ggsave(paste0(outdir, '-proportion_vaccine_adj_abs_deaths.png'), w = 8, h = 6)
  
  ggplot(tmp, aes(x = date, y = prop, col = loc_label)) + 
    geom_line() + 
    facet_grid(~age) +
    labs(y = 'Proportion of vaccinated', x= '') + 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white"))
  ggsave(paste0(outdir, '-proportion_vaccine.png'), w = 6, h = 4)
  
  ggplot(tmp, aes(x = date, y = weekly.deaths, col = loc_label)) + 
    geom_line() + 
    facet_grid(~age) +
    labs(y = 'Proportion of vaccinated', x= '') + 
    theme_bw() + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="black", fill="white"))
  ggsave(paste0(outdir, '-weekly_deaths_4states.png'), w = 6, h = 4)
  
}
