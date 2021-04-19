compare_CDC_JHU_error_plot = function(CDC_data_1, CDC_data_2 = NULL, JHU_data, var.daily.deaths.CDC, outdir)
{
  # find errors 
  CDCdata = CDC_data_1[, list(cumulative_deaths.CDC = sum(na.omit( get(var.daily.deaths.CDC) ))), by = c('code', 'date')]
  # CDCdata = subset(CDCdata, cumulative_deaths.CDC > 0)
  if(!is.null(CDC_data_2)){
    tmp = CDC_data_2[, list(cumulative_deaths.CDC = sum(na.omit( get(var.daily.deaths.CDC) ))), by = c('code', 'date')]
    # tmp = subset(tmp, cumulative_deaths.CDC > 0)
    CDCdata = rbind(CDCdata, tmp)
  }
  
  # take cumulative deaths
  CDCdata[, cumulative_deaths.CDC := cumsum(cumulative_deaths.CDC), by = c('code', 'code')]
  
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  JHUData = subset(JHUData, date <= max(CDCdata$date))
  
  tmp1 = merge(JHUData, CDCdata, by = c('code', 'date'))
  tmp1[, prop_diff := abs(cumulative_deaths - cumulative_deaths.CDC) / cumulative_deaths ]
  tmp1 = tmp1[, list(prop_diff = sum(prop_diff) / length(date)), by = c('code') ]
  
  # plot
  JHUData[, source := 'JHU']
  CDCdata[, source := 'CDC']
  setnames(CDCdata, 'cumulative_deaths.CDC', 'cumulative_deaths')
  
  tmp2 = rbind(JHUData, CDCdata)
  tmp2 = merge(tmp2, tmp1, by = 'code')
  tmp2[, code_2 := paste0(code, ', ', round(prop_diff*100, digits = 2), ' % error')]
  
  
  if(!dir.exists( basename(outdir) )){
    dir.create( basename(outdir), recursive = T )
  }
  
  n_code = length(unique(tmp2$code_2))
  
  p = ggplot(tmp2, aes(x = date, y = cumulative_deaths, col = source)) + 
    geom_line() +
    facet_wrap(~code_2, nrow = length(unique(tmp1$code)), scale = 'free') + 
    theme_bw() + 
    scale_color_viridis_d(option = "B", direction = -1, end = 0.8) 
  ggsave(p, file = paste0(outdir, '-comparison_JHU_CDC.pdf'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
}

plot_data = function(deathByAge_1, deathByAge_2, outdir)
{
  # age group first specification
  p = ggplot(deathByAge_1, aes(x = date, y = age)) + 
    geom_raster(aes(fill = daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'sqrt', breaks = c(100,300,600)) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age groups first specification',
         fill = 'Reported covid-19 deaths')
  ggsave(p, file = paste0(outdir, '-data_part1.png'), w = 12, h = 10)
  
  p1 = ggplot(deathByAge_1, aes(x = date, y = age)) + 
    geom_raster(aes(fill = min.daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30),
                         limits = c(0,max(na.omit(deathByAge_1$max.daily.deaths)))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age groups first specification',
         fill = "Bound's value covid-19 deaths", title = 'lower bound')
  
  p2 = ggplot(deathByAge_1, aes(x = date, y = age)) + 
    geom_raster(aes(fill = max.daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30), 
                         limits = c(0,max(na.omit(deathByAge_1$max.daily.deaths)))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age groups first specification',
         fill = "Bound's value covid-19 deaths", title = 'upper bound')
  
  p = ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, file = paste0(outdir, '-data_bounds_part1.png'), w = 17, h = 10)
  
  
  ## age group second specification
  p = ggplot(deathByAge_2, aes(x = date, y = age)) + 
    geom_raster(aes(fill = daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'sqrt', breaks = c(100,1000,3000)) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90)) +
    labs(x = '', y = 'Age groups second specification',
         fill = 'Reported covid-19 deaths')
  ggsave(p, file = paste0(outdir, '-data_part2.png'), w = 12, h = 10)
  
  p1 = ggplot(deathByAge_2, aes(x = date, y = age)) + 
    geom_raster(aes(fill = min.daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30),
                         limits = c(0,max(na.omit(deathByAge_2$max.daily.deaths)))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age groups first specification',
         fill = "Bound's value covid-19 deaths", title = 'lower bound')
  
  p2 = ggplot(deathByAge_2, aes(x = date, y = age)) + 
    geom_raster(aes(fill = max.daily.deaths )) + 
    facet_wrap(~loc_label) + 
    theme_bw() +
    scale_fill_viridis_c(breaks = c(0,10,20,30), 
                         limits = c(0,max(na.omit(deathByAge_2$max.daily.deaths)))) +
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age groups first specification',
         fill = "Bound's value covid-19 deaths", title = 'upper bound')
  
  p = ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, file = paste0(outdir, '-data_bounds_part2.png'), w = 17, h = 10)
  
}



make_adjacency_plot = function()
{
  n = 2; m = 3
  
  N = n * m
  A = matrix(nrow = N, ncol = N, 0)
  B = matrix(nrow = n, ncol = m, 1:(n*m), byrow = T)
  
  for(i in 1:n){
    
    for(j in 1:m){
      
      #cat('\n Processing ', i, j)
      idx_row = m*(i-1) + j
      
      if(i - 1 > 0){
        idx_col = m*(i-2) + j
        A[idx_row,idx_col] = 1 
      }
      
      if(i + 1 <= n){
        idx_col = m*i + j
        A[idx_row,idx_col] = 1 
      }
      
      if(j - 1 > 0){
        idx_col = m*(i-1) + j - 1
        A[idx_row,idx_col] = 1 
      }
      
      if(j + 1 <= m){
        idx_col = m*(i-1) + j +1
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
    scale_y_reverse(expand = c(0,0), breaks = 1:(m*n), labels = sort(tmp$idx))  +
    scale_x_continuous(expand = c(0,0), breaks = 1:(m*n), labels = sort(tmp$idx))  +
    theme(plot.title = element_text( hjust = 0.5)) 
  ggsave(file = '~/Box\ Sync/2021/CDC/example_Adj.png', w = 4.5, h = 4)
  
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
