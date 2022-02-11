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
    facet_wrap(.~loc_label, nrow = 4) + 
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
  ggsave(p, file = paste0(outdir, '-deathByAge_selected_states.png'), w = 4, h = 8)
  
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

plot_data_one_state = function(tmp)
{
  
  range_wd = sqrt(range(na.omit(tmp$weekly.deaths)))
  digits_cut= ifelse(range_wd[2]/100 > 1, 2, 1)
  digits_cut= ifelse(range_wd[2]/10 > 1, digits_cut, 0)
  p1 = ggplot(tmp, aes(x = date, y = age)) + 
    geom_raster(aes(fill = weekly.deaths )) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'sqrt', breaks = round(seq(range_wd[2]/2, range_wd[2], length.out = 3)^2, digits = -digits_cut) ) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'left',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = '')
  
  p2 = ggplot(tmp, aes(x = date, y = age)) + 
    geom_raster(aes(fill = prop_deaths )) + 
    theme_bw() +
    scale_fill_viridis_c(option = "E") +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'left',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = '')
  
  return(list(p1, p2))
}

make_adjacency_plot = function()
{
  library(data.table)
  library(ggplot2)
  
  n = 2; m = 3
  
  N = n * m
  A = matrix(nrow = N, ncol = N, 0)
  B = matrix(nrow = n, ncol = m, 1:(n*m), byrow = T)
  
  for(i in 1:n){
    
    for(j in 1:m){
      
      #cat('\n Processing ', i, j)
      idx_row = n*(j-1) + i
      
      # top
      if(i - 1 > 0){
        idx_col = n*(j-1) + i - 1
        A[idx_row,idx_col] = 1
      }
      
      # bottom
      if(i + 1 <= n){
        idx_col = n*(j-1) + i + 1
        A[idx_row,idx_col] = 1
      }
      
      # left
      if(j - 1 > 0){
        idx_col = n*(j-2) + i
        A[idx_row,idx_col] = 1
      }
      
      # right
      if(j + 1 <= m){
        idx_col = n*j + i
        A[idx_row,idx_col] = 1
      }

    }
  }
  
  tmp = as.data.table( reshape2::melt(B) )
  tmp[, idx := paste0('(',Var1, ',', Var2,')')]
  ggplot(tmp, aes(x = Var2, y = Var1, label = idx)) +
    theme_minimal() +
    labs(y = 'rows', x = 'columns', title = 'Original matrix B') + 
    scale_y_reverse(breaks = 1:m)  +
    scale_x_continuous(breaks = 1:n) +
    geom_raster(aes(fill = value), fill = 'white') + 
    geom_rect(aes(ymin = Var1 - 0.5, ymax = Var1 + 0.5, xmin = Var2 - 0.5, xmax = Var2 + 0.5), colour = "grey50", fill = NA) +
    geom_text()  +
    theme(axis.text = element_blank(),
          plot.title = element_text(vjust = -2, hjust = 0.5), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(file = '~/Box\ Sync/2021/CDC/example_B.png', w = 3, h = 2)
  
  
  tmp1 = as.data.table( reshape2::melt(A) )
  ggplot(tmp1, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill = as.factor(value))) + 
    scale_fill_manual(values = c('beige',"blue")) + 
    labs(x = '', y = '', fill = '', title = 'Adjacency matrix A') + 
    scale_y_reverse(expand = c(0,0), breaks = 1:(m*n), labels = tmp$idx)  +
    scale_x_continuous(expand = c(0,0), breaks = 1:(m*n), labels = tmp$idx)  +
    theme(plot.title = element_text( hjust = 0.5)) 
  ggsave(file = '~/Box\ Sync/2021/CDC/example_Adj.png', w = 3.5, h = 3)
  
}

plot_basis_functions = function()
{
  tmp = as.data.table( reshape2::melt(stan_data$BASIS) )
  setnames(tmp, c('Var1', 'Var2'), c('basis_idx', 'age'))
  
  ggplot(tmp, aes(x = age, y= value, col = as.factor(basis_idx))) + 
    geom_line() + 
    labs(x = 'Age', y = '', col = 'basis function index') + 
    theme_bw() +
    theme(legend.position='bottom')
  ggsave(file = '~/Box\ Sync/2021/CDC/basis_functions.png', w = 6, h = 5)
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
  p <- ggplot(subset(tmp, code %in% c('CA', 'FL', 'NY', 'TX') & age_index >2), aes(date, prop)) +
    geom_line(aes(col = loc_label)) + 
    facet_wrap(~`Age group`, label = 'label_both', nrow = 2) + 
    theme_bw()+ 
    labs( y = 'Proportion of fully vaccinated individuals', col = '') + 
    theme(legend.key = element_blank(), 
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill = NA), 
          legend.position = 'bottom',
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), 
          axis.title.x = element_blank(), 
          axis.text.y = element_text(size = rel(1)), 
          axis.title.y = element_text(size = rel(1.1)), 
          strip.text = element_text(size = rel(0.9))) + 
    scale_color_jcolors('pal8') + 
    scale_y_continuous(labels = scales::percent)+ 
    scale_x_date(date_labels = c("%b-%y"), breaks = '2 months', expand = c(0,0)) 
  ggsave(p, file = paste0(outdir, '-proportion_vaccine_age_code_selected_states.png'), w = 6, h = 3.75)
  
  if(length(Code) > 6){
    
    mid_code = length(Code) / 2
    
    p1 <- ggplot(subset(tmp, code %in% Code[1:mid_code] & age_index >2), aes(date, prop)) +
      geom_hline(yintercept = 0.5, linetype = 'dashed', col = 'grey50') + 
      geom_line(aes(col = `Age group`)) + 
      theme_bw()+ 
      labs( y = '', col = '') + 
      theme(legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"),
            panel.border = element_rect(colour = "black", fill = NA), 
            legend.position = 'bottom',
            axis.text.x = element_blank(), 
            axis.title.x = element_blank()) + 
      scale_color_viridis_d(option = 'B', begin = 0.1, end = 0.9) + 
      facet_grid(`Age group`~loc_label) +
      scale_y_continuous(labels = scales::percent)+ 
      scale_x_date(date_labels = c("%b-%y"), breaks = '2 months', expand = c(0,0)) 
      
    
    p2 <- ggplot(subset(tmp, code %in% Code[(mid_code + 1):(mid_code*2)] & age_index >2), aes(date, prop)) +
      geom_hline(yintercept = 0.5, linetype = 'dashed', col = 'grey50') + 
      geom_line(aes(col = `Age group`)) + 
      theme_bw()+ 
      labs( y = '', col = '') + 
      theme(legend.key = element_blank(), 
            strip.background = element_rect(colour="white", fill="white"),
            panel.border = element_rect(colour = "black", fill = NA), 
            legend.position = 'bottom',
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
            axis.title.x = element_blank()) + 
      scale_color_viridis_d(option = 'B', begin = 0.1, end = 0.9) + 
      facet_grid(`Age group`~loc_label) +
      scale_y_continuous(labels = scales::percent)+ 
      scale_x_date(date_labels = c("%b-%y"), breaks = '2 months', expand = c(0,0)) 
    
    p_all <- ggarrange(p1, p2, nrow = 2,  common.legend = T, legend = 'none')
    grid.arrange(p_all, left = 'Proportion of fully vaccinated individuals')
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
      scale_x_date(date_labels = c("%b-%y"), breaks = '2 months', expand = c(0,0)) 
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
  
  return(p)
}

save_deathByAge_vaccineByAge_panel <- function(plot_deathByAge, plot_vaccineByAge, outdir){
  p1 <- ggarrange(plot_deathByAge, labels = c('A'), label.x = 0.03)
  p2 <- ggarrange(plot_vaccineByAge, labels = c('B'), label.x = 0.05)
  
  p <- grid.arrange(grobs = list(p1, p2), layout_matrix = rbind(c(1, 2), 
                                                                c(1, NA)), heights = c(1, 0.4), widths = c(1, 1))
  ggsave(p, file = paste0(outdir, '-deathByAge_proportion_vaccine_panel_selected_states.png'), w = 8,h = 8)
}
