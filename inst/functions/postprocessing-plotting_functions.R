plot_probability_deaths_age_contribution = function(fit, var_name, df_age, df_week, lab, outdir, discrete = F){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return(ggplot())
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples[[var_name]]) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps),
                       q_label=p_labs), 
              by=c('age_index', 'week_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_age$age[age_index]]
  
  if(discrete)
    tmp1[, age := factor(age, levels = df_age$age)]
  
  n_row = length(unique(tmp1$date))
  
  p = ggplot(tmp1, aes(x = age)) + 
    theme_bw() +
    labs(y = paste0("Relative contribution to ", lab), x = "Age", title = paste(Code)) + 
    facet_wrap(~date)
  
  if(discrete){
    p = p + 
      geom_point(aes(y = M)) +
      geom_errorbar(aes(ymin= CL, ymax = CU)) 
  } else {
    p = p + 
      geom_line(aes(y = M)) +
      geom_ribbon(aes(ymin= CL, ymax = CU), alpha = 0.5) 
  }
    
  ggsave(p, file = paste0(outdir, "-continuous_contribution_", var_name, '_', Code, ".png") , w= 10, h = 8, limitsize = FALSE)
  
  p = ggplot(tmp1, aes(x = date, y = age)) +
    geom_raster(aes(fill = M))  + 
    scale_fill_viridis_c(option = "E") + 
    theme(legend.position='bottom') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  
  if(discrete){
    p = p + scale_y_discrete(expand = c(0,0))  + labs(x = '', y = 'Age group', fill = 'Estimated posterior value') 
  } else {
    p = p + scale_y_continuous(expand = c(0,0))  +labs(x = '', y = 'Age', fill = 'Estimated posterior value') 
  }
  ggsave(p, file = paste0(outdir, '-continuous_contribution_allweeks_', var_name, '_', Code, '.png'), w = 5, h = 5.2)
  
}

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

plot_posterior_predictive_checks = function(data, variable, lab, outdir)
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
    labs(y = lab, x = "") + 
    facet_wrap(~date_ppp, scale = 'free_y') + 
    theme(axis.text.x = element_text(angle = 90))
  ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
  
  # posterior predictive check
  p = ggplot(data, aes(x = date)) + 
    geom_point(aes(y = M), col = "black", size = 1) +
    geom_errorbar(aes(ymin = CL, ymax = CU), width=0.3, col = "black") +
    geom_point(aes(y = get(variable)), col = "red", size = 1) + 
    theme_bw() +
    labs(y = lab, x = "") + 
    facet_wrap(~age, ncol = 3, scale = 'free_y') + 
    theme(axis.text.x = element_text(angle = 90))
  ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_2_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
  
}

plot_posterior_predictive_checks_cumulative = function(tmp, outdir)
{

  p = ggplot(tmp, aes(x = date)) + 
    geom_point(aes(y = M)) +
    geom_errorbar(aes(ymin = CL, ymax = CU)) +
    geom_ribbon(aes(ymin = min_count_censored, ymax = max_count_censored), alpha = 0.5, fill = 'red') +
    geom_point(aes(y = sum_count_censored), col = 'blue') +
    facet_wrap(~age, scale = 'free_y') + 
    labs(x = '', y = 'cumulative deaths')
  ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_cum_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
  
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

compare_CDCestimation_JHU_Imperial_plot = function(CDC_data, JHU_data, scraped_data, var.cum.deaths.CDC, outdir)
{
  # prepare JHU data
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  JHUData = subset(JHUData, code == Code)
  JHUData[, date := as.Date(date)]
  JHUData[, CL := NA]
  JHUData[, CU := NA]
  
  # prepare Imperial data
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
  scraped_data[, source := 'Imperial COVID-19 Team']
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
    labs(x = '', y = 'Cumulative COVID-19 \nattributable deaths', col = '', fill = '') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90)) 
  ggsave(p, file = paste0(outdir, '-comparison_JHU_CDC_Imperial_uncertainty_', Code, '.png'), w = 4, h = 4, limitsize = F)
}

compare_CDCestimation_JHU_error_plot = function(CDC_data, JHU_data, var.cum.deaths.CDC, outdir)
{
  # prepare JHU data
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  
  # prepare CDC estimations
  CDC_data = select(as.data.table(CDC_data), code, date, var.cum.deaths.CDC)
  setnames(CDC_data, var.cum.deaths.CDC, 'cumulative_deaths')
  
  # plot
  JHUData[, source := 'JHU']
  CDC_data[, source := 'CDC']
  
  tmp2 = rbind(JHUData, CDC_data)
  tmp2 = subset(tmp2, code %in% unique(CDC_data$code) & date <= max(CDC_data$date))
  
  n_code = length(unique(tmp2$code))
  
  p = ggplot(tmp2, aes(x = date, y = cumulative_deaths)) + 
    geom_line(aes(col = source), size = 1) +
    facet_wrap(~code, nrow = length(unique(tmp2$code)), scale = 'free') + 
    theme_bw() + 
    scale_color_viridis_d(option = "B", direction = -1, end = 0.8) + 
    scale_fill_viridis_d(option = "B", direction = -1, end = 0.8)
  ggsave(p, file = paste0(outdir, '-comparison_JHU_CDC_uncertainty_', Code, '.png'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
}

compare_CDCestimation_Imperial_age_plot = function(CDC_data, scraped_data, var.cum.deaths.CDC, outdir, overall = F)
{
  scraped_data = select(as.data.table(scraped_data), code, date, age, cum.deaths)
  scraped_data = subset(scraped_data, code == Code)
  scraped_data[, date := as.Date(date)]
  scraped_data[, CL := NA]
  scraped_data[, CU := NA]
  
  # prepare CDC estimations
  CDC_data = select(as.data.table(CDC_data), code, date, age, var.cum.deaths.CDC, CL, CU)
  setnames(CDC_data, var.cum.deaths.CDC, 'cum.deaths')
  
  # plot
  scraped_data[, source := 'Imperial COVID-19 Team']
  CDC_data[, source := 'CDC']
  
  tmp2 = rbind(scraped_data, CDC_data)
  tmp2 = subset(tmp2, code %in% unique(CDC_data$code) & date <= max(CDC_data$date))
  tmp2[, age := factor(age, levels = unique(scraped_data$age))]
  
  col = viridisLite::viridis(3, option = "B", direction = -1, end = 0.8)
  
  p = ggplot(tmp2, aes(x = date, y = cum.deaths)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = source), alpha = 0.5) +
    geom_line(aes(col = source), size = 1) +
    theme_bw() + 
    scale_color_manual(values = col[c(1,2)]) + 
    scale_fill_manual(values = col[c(1,2)]) +
    facet_wrap(~age,  scale = 'free', ncol = 2) + 
    theme(legend.position = 'bottom',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x = element_text(angle = 90)) + 
    labs(y = 'Cumulative COVID-19 attributable deaths', col = '', fill = '', x = '')
  ggsave(p, file = paste0(outdir, '-comparison_Imperial_CDC_uncertainty_', Code, '.png'), w = 4, h = 8, limitsize = F)
  
}

plot_covariance_matrix = function(fit_cum, outdir)
{
  samples = extract(fit_cum)
  
  if(is.null(stan_data$Adj)) return()
  if(!is.null(stan_data$node1)) return()
  
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

plot_posterior_plane = function(fit_cum, df_week, outdir)
{
  
  euro.levs <- as.vector(outer(c(1, 2, 5), 10^(-3:3)))   
  fit_samples = extract(fit_cum)
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  tmp1 = as.data.table( reshape2::melt( fit_samples$beta ))
  row_name = 'week_index'
  column_name = 'basis_function_index'
  if(max(tmp1$Var3) == stan_data$A) column_name = 'Age'
  setnames(tmp1, c('Var2', 'Var3'), c(row_name, column_name))
  tmp1 = tmp1[, list(q = quantile(value, prob=ps), q_label=p_labs), by=c(column_name, row_name)]
  tmp1 = dcast(tmp1, get(row_name) + get(column_name) ~ q_label, value.var = "q")
  setnames(tmp1, c('row_name', 'column_name'), c(row_name, column_name))
  tmp1 = merge(tmp1, df_week, by = row_name)
  
  ## Smooth estimate
  z <- as.matrix( dcast(tmp1, get(column_name)~get(row_name), value.var = "M")[,-1] ) 
  z.range <- range(z)
  z =  t( apply(z, 2, rev) )
  
  png(paste0(outdir, '-PlanePosterior_', Code, '.png'),width = 4, height = 4, units = 'in', res = 300, pointsize = 10)
  plot.new()
  plot.window(
    xlim = c(1 - 0.5, stan_data$W + 0.5),
    ylim = c(1 - 0.5, stan_data$num_basis + 0.5))
  image  (1:stan_data$W, 1:stan_data$num_basis,  z, zlim = z.range, useRaster = TRUE, add = TRUE, col = hcl.colors(12, "YlOrRd", rev = F))
  contour(1:stan_data$W, 1:stan_data$num_basis,  z, levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
  axis(side = 1, at = seq(1, stan_data$W, 5), labels = format(unique(tmp1$date), '%b %Y')[seq(1, stan_data$W, 5)], lwd = 0.5)
  axis(side = 2, at = seq(1, stan_data$num_basis, 2), labels = rev(seq(2, stan_data$num_basis, 2)), lwd = 0.5)
  box(lwd = 0.5)
  mtext("Date", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
  mtext("Basis function index",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
  dev.off()
  
}

plot_posterior_plane_with_data = function(deathByAge_1, deathByAge_2, outdir)
{
  library(magick)
  
  deathByAge_1 = subset(deathByAge_1, code == Code)
  deathByAge_2 = subset(deathByAge_2, code == Code)
  
  range_d = range(na.omit(c(deathByAge_1$daily.deaths, deathByAge_2$daily.deaths)))
  
  p1 = ggplot(deathByAge_1, aes(x = date, y = age)) + 
    geom_raster(aes(fill = daily.deaths )) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'sqrt', limits = range_d) +
    scale_x_date(expand = c(0,0), date_labels = c("%b %Y"), date_breaks = '1 month') + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 70, vjust  = 0.5), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age groups first specification',
         fill = 'Reported covid-19 deaths')
  
  p2 = ggplot(deathByAge_2, aes(x = date, y = age)) + 
    geom_raster(aes(fill = daily.deaths )) + 
    theme_bw() +
    scale_fill_viridis_c(trans = 'sqrt', limits = range_d) +
    scale_x_date(expand = c(0,0), date_labels = c("%b %Y"), date_breaks = '1 month') + 
    scale_y_discrete(expand = c(0,0)) + 
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 70, vjust  = 0.5), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = '', y = 'Age groups second specification',
         fill = 'Reported covid-19 deaths')
  
  p = ggpubr::ggarrange(p1, p2, nrow = 2,common.legend = T, legend = 'bottom')
  ggsave(p, file = paste0(outdir, '-panel_right_', Code, '.png'), h = 8, w = 6)
  
  panel.left <- magick::image_read(paste0(outdir, '-PlanePosterior_', Code, '.png'))
  panel.right <- magick::image_read(paste0(outdir, '-panel_right_', Code, '.png'))
  
  p = image_append(c( image_scale(panel.left, "2200"), panel.right) )

  savepdf(paste0(outdir, '-PlanePosterior_Data_', Code, '.pdf'), w = 14.5*2, h = 10*2)
  plot(p)
  dev.off()

}

savepdf <- function(fname, width=16, height=10)
{
  pdf(fname, width=width/2.54, height=height/2.54)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(0,0,0,0))
}

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

plot_probability_ratio = function(probability_ratio_table, outdir)
{
  # plot
  tmp = subset(probability_ratio_table, age %in% c('45-54', '55-64', '65-74', '75-84', '85+'))
  
  ggplot(tmp, aes(x = date)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line(aes(y = M)) + 
    geom_point(aes(y = emp.prob.ratio), col = 'darkred') + 
    scale_x_date(expand = c(0,0), date_labels = c("%b-%y")) + 
    theme_bw() + 
    geom_hline(aes(yintercept = 1)) +
    facet_wrap(~age, scales = 'free', ncol =1)+ 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle = 90)) +
    labs(x = '', y = paste0('Ratio of the share of weekly COVID-19 deaths relative to its mean before ', format(ref_date, "%d-%b-%y")))
  ggsave(file = paste0(outdir, '-ProbabilityRatio_', Code, '.png'), w = 4, h = 10)
  
  tmp = subset(probability_ratio_table, age %in% c('45-54', '55-64', '65-74', '75-84', '85+'))
  ggplot(tmp, aes(x = date)) + 
    geom_line(aes(y = M, col = age)) + 
    # geom_point(aes(y = emp.prob.ratio, col = age)) + 
    # geom_ribbon(aes(ymin = CL, ymax = CU, fill = age), alpha = 0.1) +
    theme_bw() + 
    geom_hline(aes(yintercept = 1)) +
    labs(x = '', y = paste0('Ratio of the share of weekly COVID-19 deaths relative to its mean before ', format(ref_date, "%d-%b-%y")))
  ggsave(file = paste0(outdir, '-ProbabilityRatio_elderly_', Code, '.png'), w = 8, h = 8)
}

