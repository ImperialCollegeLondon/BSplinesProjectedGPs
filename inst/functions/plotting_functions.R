compare_CDC_JHU_Imperial_error_plot = function(CDC_data, JHU_data, scrapedData, var.daily.deaths.CDC, outdir, Code = NULL)
{
  # find errors 
  CDCdata = CDC_data[, list(cumulative_deaths.CDC = sum(na.omit( get(var.daily.deaths.CDC) ))), by = c('code', 'date')]
  
  # take cumulative deaths
  CDCdata[, cumulative_deaths.CDC := cumsum(cumulative_deaths.CDC), by = c('code', 'code')]
  
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  JHUData[, date := as.Date(date)]
  JHUData = subset(JHUData, date <= max(CDCdata$date))
  
  scrapedData = select(as.data.table(scrapedData), code, date, cum.deaths, age)
  scrapedData = scrapedData[, list(cum.deaths = sum(cum.deaths)), by = c('code', 'date')]
  scrapedData[, date := as.Date(date)]
  scrapedData = subset(scrapedData, date <= max(CDCdata$date))
  
  
  # tmp1 = merge(JHUData, CDCdata, by = c('code', 'date'))
  # tmp1[, prop_diff := abs(cumulative_deaths - cumulative_deaths.CDC) / cumulative_deaths ]
  # tmp1 = tmp1[, list(prop_diff = sum(prop_diff) / length(date)), by = c('code') ]
  
  # plot
  JHUData[, source := 'JHU']
  CDCdata[, source := 'CDC']
  scrapedData[, source := 'Imperial COVID-19 Team']
  setnames(CDCdata, 'cumulative_deaths.CDC', 'cumulative_deaths')
  setnames(scrapedData, 'cum.deaths', 'cumulative_deaths')
  
  tmp2 = rbind(rbind(JHUData, CDCdata),scrapedData)
  # tmp2 = merge(tmp2, tmp1, by = 'code')
  # tmp2[, code_2 := paste0(code, ', ', round(prop_diff*100, digits = 2), ' % error')]
  
  
  if(!dir.exists( basename(outdir) )){
    dir.create( basename(outdir), recursive = T )
  }
  
  n_code = length(unique(tmp2$code))
  
  p = ggplot(tmp2, aes(x = date, y = cumulative_deaths, col = source)) + 
    geom_line() +
    facet_wrap(~code, nrow = length(unique(tmp2$code)), scale = 'free') + 
    theme_bw()  +
    labs(x = '', y = 'Cumulative COVID-19 attributable deaths', col = '') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90)) 
  ggsave(p, file = paste0(outdir, '-comparison_JHU_Imperial_CDC.pdf'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
  
  if(!is.null(Code)){
    tmp = subset(tmp2, code == Code)
    
    p = ggplot(tmp, aes(x = date, y = cumulative_deaths, col = source)) + 
      geom_line() +
      theme_bw() + 
      scale_color_viridis_d(option = "B", direction = -1, end = 0.8) +
      labs(x = '', y = 'Cumulative COVID-19 \nattributable deaths', col = '') + 
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
      theme(legend.position = 'bottom',
            axis.text.x = element_text(angle = 90)) 
    ggsave(p, file = paste0(outdir, '-comparison_JHU_Imperial_CDC_', Code, '.pdf'), w = 5, h = 3, limitsize = F)
    
  }
}

plot_data = function(deathByAge, outdir, Code = NULL)
{
  p = ggplot(deathByAge, aes(x = date, y = age)) + 
    geom_raster(aes(fill = daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'sqrt', breaks = c(100,1000,3000)) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = 'Reported covid-19 deaths')
  ggsave(p, file = paste0(outdir, '-deathByAge.png'), w = 12, h = 10)
  
  p1 = ggplot(deathByAge, aes(x = date, y = age)) + 
    geom_raster(aes(fill = min.sum.daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30),
                         limits = c(0,max(na.omit(deathByAge$max.sum.daily.deaths)))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = "Value sum of covid-19 daily deaths", title = 'lower bound')
  
  p2 = ggplot(deathByAge, aes(x = date, y = age)) + 
    geom_raster(aes(fill = max.sum.daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30), 
                         limits = c(0,max(na.omit(deathByAge$max.sum.daily.deaths)))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = "Value sum of covid-19 daily deaths", title = 'upper bound')
  
  p3 = ggplot(deathByAge, aes(x = date, y = age)) + 
    geom_raster(aes(fill = sum.daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30), 
                         limits = c(0,max(na.omit(deathByAge$max.sum.daily.deaths)))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age group',
         fill = "Value sum of covid-19 daily deaths", title = 'exact sum')
  
  p = ggpubr::ggarrange(p1, p2, p3,  common.legend = T, legend = 'bottom')
  ggsave(p, file = paste0(outdir, '-deathByAge_boundaries.png'), w = 20, h = 25)
  
  if(!is.null(Code)){
    tmp = subset(deathByAge, code == Code)
    tmp1 = tmp[, list(total_deaths = sum(na.omit(daily.deaths))), by = 'date']
    tmp = merge(tmp, tmp1, by = 'date')
    tmp[, prop_deaths := daily.deaths / total_deaths]
      
    range_wd = sqrt(range(na.omit(tmp$daily.deaths)))
    p = ggplot(tmp, aes(x = date, y = age)) + 
      geom_raster(aes(fill = daily.deaths )) + 
      theme_bw() +
      scale_fill_viridis_c(trans = 'sqrt', breaks = round(seq(range_wd[2]/2, range_wd[2], length.out = 3)^2, digits = -2) ) +
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
      scale_y_discrete(expand = c(0,0)) + 
      theme(legend.position = 'bottom',
            axis.text.x = element_text(angle = 90), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      labs(x = '', y = 'Age group',
           fill = 'Reported covid-19 deaths')
    ggsave(p, file = paste0(outdir, '-deathByAge_', Code, '.png'), w = 5, h = 5.2)
    
    p = ggplot(tmp, aes(x = date, y = age)) + 
      geom_raster(aes(fill = prop_deaths )) + 
      theme_bw() +
      scale_fill_viridis_c(option = "E") +
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
      scale_y_discrete(expand = c(0,0)) + 
      theme(legend.position = 'bottom',
            axis.text.x = element_text(angle = 90), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      labs(x = '', y = 'Age group',
           fill = 'Share of weekly COVID-19 deaths')
    ggsave(p, file = paste0(outdir, '-deathByAge_prop_', Code, '.png'), w = 5, h = 5.2)
    
    p1 = ggplot(tmp, aes(x = date, y = age)) + 
      geom_raster(aes(fill = min.sum.daily.deaths )) + 
      theme_bw() +
      scale_fill_viridis_c(breaks = c(0,10,20,30),
                           limits = c(0,max(na.omit(tmp$max.sum.daily.deaths)))) +
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
      scale_y_discrete(expand = c(0,0)) + 
      theme(legend.position = 'bottom',
            axis.text.x = element_text(angle = 90), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = '', y = 'Age group',
           fill = "Sum of unobserved covid-19 weekly deaths", title = 'Lower bound')
    
    p2 = ggplot(tmp, aes(x = date, y = age)) + 
      geom_raster(aes(fill = max.sum.daily.deaths )) + 
      theme_bw() +
      scale_fill_viridis_c(breaks = c(0,10,20,30), 
                           limits = c(0,max(na.omit(tmp$max.sum.daily.deaths)))) +
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
      scale_y_discrete(expand = c(0,0)) + 
      theme(legend.position = 'bottom',
            axis.text.x = element_text(angle = 90), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = '', y = 'Age group',
           fill = "Sum of unobserved covid-19 weekly deaths", title = 'Upper bound')
    
    p3 = ggplot(tmp, aes(x = date, y = age)) + 
      geom_raster(aes(fill = sum.daily.deaths )) + 
      theme_bw() +
      scale_fill_viridis_c(breaks = c(0,10,20,30), 
                           limits = c(0,max(na.omit(tmp$max.sum.daily.deaths)))) +
      scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
      scale_y_discrete(expand = c(0,0)) + 
      theme(legend.position = 'bottom',
            axis.text.x = element_text(angle = 90), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = '', y = 'Age group',
           fill = "Sum of unobserved COVID-19 weekly deaths", title = 'Exact value')
    
    p = ggpubr::ggarrange(p1, p2, p3, nrow = 1, common.legend = T, legend = 'bottom')
    ggsave(p, file = paste0(outdir, '-deathByAge_boundaries_', Code, '.png'), w = 10, h = 4)
    
  }

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
