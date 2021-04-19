plot_continuous_age_contribution = function(fit, df_age_continuous, df_week, lab, outdir){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  if(is.null(fit)) return(ggplot())
  
  # extract samples
  fit_samples = rstan::extract(fit)
  
  tmp1 = as.data.table( reshape2::melt(fit_samples$phi) )
  setnames(tmp1, c('Var2', 'Var3'), c('age_index','week_index'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps),
                       q_label=p_labs), 
              by=c('age_index', 'week_index')]	
  tmp1 = dcast(tmp1, week_index + age_index ~ q_label, value.var = "q")
  
  tmp1 = merge(tmp1, df_week, by = 'week_index')
  tmp1[, age := df_age_continuous$age_index[age_index]]
  
  n_row = length(unique(tmp1$date))
  
  ggplot(tmp1, aes(x = date, y = age)) +
    geom_raster(aes(fill = M))  + 
    labs(x = 'Date', y = 'Age', fill = 'Estimated posterior value') + 
    scale_y_continuous(expand = c(0,0))  +
    scale_x_date(expand = c(0,0)) + 
    scale_fill_viridis_c(option = "E") + 
    theme(legend.position='bottom')
  ggsave(file = paste0(outdir, '-continuous_contribution_allweeks_', Code, '.png'), w = 5, h = 5.2)
  
  p = ggplot(tmp1, aes(x = age)) + 
    geom_line(aes(y = M)) +
    geom_ribbon(aes(ymin= CL, ymax = CU), alpha = 0.5) + 
    theme_bw() +
    labs(y = paste0("Relative contribution to ", lab), x = "Age", title = paste(Code)) + 
    facet_wrap(~date)
  ggsave(p, file = paste0(outdir, "-continuous_contribution_",Code, ".png") , w= 10, h = 8, limitsize = FALSE)
  
}

plot_convergence_diagnostics = function(fit, title, suffix, outfile)
{
  
  stopifnot(!is.null(fit))
  
  posterior <- as.array(fit)
  
  pars = list("rho", "sigma", "beta", c('aw', 'sd_aw'), c('sd_a_raw', 'a0_raw'), 'a_age', 'a0', c('gamma', 'tau', 'p'))

  for(j in 1:length(pars))
  {                     
    
    par = pars[[j]]
    
    if(any(!par %in% names(fit))) next
    
    p_trace = bayesplot::mcmc_trace(posterior, regex_pars = par) + labs(title = title) 
    p_pairs = gridExtra::arrangeGrob(bayesplot::mcmc_pairs(posterior, regex_pars = par), top = title)
    p_intervals = bayesplot::mcmc_intervals(posterior, regex_pars = par, prob = 0.95,
                                            prob_outer = 0.95) + labs(title = title)
    
    ggsave(p_trace, file = paste0(outfile, "-convergence_diagnostics-", "trace_plots_", suffix, '_', Code, "_", par,".png") , w= 10, h = 10)
    ggsave(p_pairs, file =  paste0(outfile, "-convergence_diagnostics-", "pairs_plots_",  suffix, '_',Code, "_", par,".png") , w= 15, h = 15)
    ggsave(p_intervals, file = paste0(outfile, "-convergence_diagnostics-", "intervals_plots_",  suffix, '_',Code, "_", par,".png") , w=10, h = 10)
    
  }
  
}

plot_posterior_predictive_checks = function(data_1, data_2, variable, lab, outdir)
{
  
  data_1[, PPP := paste0('inside CI: ', round(mean(na.omit(inside.CI))*100, digits = 2), '%'), by = 'date']
  data_1[, date_ppp := paste0(as.character(date), ' - ', PPP)]
  data_2[, PPP := paste0('inside CI: ', round(mean(na.omit(inside.CI))*100, digits = 2), '%'), by = 'date']
  data_2[, date_ppp := paste0(as.character(date), ' - ', PPP)]
  
  Code = unique(data_1$code)
  
  # posterior predictive check
  p1 = ggplot(data_1, aes(x = age)) + 
    geom_point(aes(y = M), col = "black", size = 1) +
    geom_errorbar(aes(ymin = CL, ymax = CU), width=0.3, col = "black") +
    geom_point(aes(y = get(variable)), col = "red", size = 1) + 
    theme_bw() +
    labs(y = lab, x = "") + 
    facet_wrap(~date_ppp, ncol = 3) + 
    theme(axis.text.x = element_text(angle = 90))
  
  p2 = ggplot(data_2, aes(x = age)) + 
    geom_point(aes(y = M), col = "black", size = 1) +
    geom_errorbar(aes(ymin = CL, ymax = CU), width=0.3, col = "black")+
    geom_point(aes(y = get(variable)), col = "red", size = 1) + 
    theme_bw() +
    labs(y = lab, x = "") + 
    facet_wrap(~date_ppp, ncol = 3) + 
    theme(axis.text.x = element_text(angle = 90))
  
  p = gridExtra::grid.arrange(p1,p2, ncol = 2, widths = c(0.8, 1))
  ggsave(p, file = paste0(outdir, "-posterior_predictive_checks_", Code,".png") , w= 15, h = 15, limitsize = FALSE)
  
}

compare_CDCestimation_JHU_error_plot_uncertainty = function(CDC_data, JHU_data, var.cum.deaths.CDC, outdir)
{
  # prepare JHU data
  JHUData = select(as.data.table(JHUData), code, date, cumulative_deaths)
  JHUData[, CL := NA]
  JHUData[, CU := NA]
  
  # prepare CDC estimations
  CDC_data = select(as.data.table(CDC_data), code, date, var.cum.deaths.CDC, CL, CU)
  setnames(CDC_data, var.cum.deaths.CDC, 'cumulative_deaths')
  
  # plot
  JHUData[, source := 'JHU']
  CDC_data[, source := 'CDC']
  
  tmp2 = rbind(JHUData, CDC_data)
  tmp2 = subset(tmp2, code %in% unique(CDC_data$code) & date <= max(CDC_data$date))
  
  n_code = length(unique(tmp2$code))
  
  p = ggplot(tmp2, aes(x = date, y = cumulative_deaths)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = source), alpha = 0.5) +
    geom_line(aes(col = source), size = 1) +
    facet_wrap(~code, nrow = length(unique(tmp2$code)), scale = 'free') + 
    theme_bw() + 
    scale_color_viridis_d(option = "B", direction = -1, end = 0.8) + 
    scale_fill_viridis_d(option = "B", direction = -1, end = 0.8)
  ggsave(p, file = paste0(outdir, '-comparison_JHU_CDC_uncertainty.png'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
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
  ggsave(p, file = paste0(outdir, '-comparison_JHU_CDC_uncertainty.png'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
}

compare_CDCestimation_scrapeddata_error_plot_uncertainty = function(CDC_data, scraped_data, var.cum.deaths.CDC, outdir)
{
  # prepare JHU data
  scraped_data = select(as.data.table(scraped_data), code, date, age, cum.deaths)
  scraped_data[, date := as.Date(date)]
  scraped_data[, CL := NA]
  scraped_data[, CU := NA]
  
  # prepare CDC estimations
  CDC_data = select(as.data.table(CDC_data), code, date, age, var.cum.deaths.CDC, CL, CU)
  setnames(CDC_data, var.cum.deaths.CDC, 'cum.deaths')
  
  # plot
  scraped_data[, source := 'DoH']
  CDC_data[, source := 'CDC']
  
  tmp2 = rbind(scraped_data, CDC_data)
  tmp2 = subset(tmp2, code %in% unique(CDC_data$code) & date <= max(CDC_data$date))
  
  n_code = length(unique(tmp2$code))
  
  p = ggplot(tmp2, aes(x = date, y = cum.deaths)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = source), alpha = 0.5) +
    geom_line(aes(col = source), size = 1) +
    facet_wrap(~age,  scale = 'free') + 
    theme_bw() + 
    scale_color_viridis_d(option = "B", direction = -1, end = 0.8) + 
    scale_fill_viridis_d(option = "B", direction = -1, end = 0.8) 
  ggsave(p, file = paste0(outdir, '-comparison_scraped_data_CDC_uncertainty.png'), w = 9, h = 2.2 * n_code + 5, limitsize = F)
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
  
  map = data.table(idx = 1:(stan_data$W * stan_data$num_basis), 
                   idx_week = rep(1:stan_data$W, each = stan_data$num_basis),
                   idx_basis = rep(1:stan_data$num_basis, stan_data$W))
  tmp1 = merge(tmp1, map, by.x = 'Var1', by.y = 'idx')
  setnames(tmp1, c('idx_week', 'idx_basis'), c('idx_week_column', 'idx_basis_column'))
  tmp1 = merge(tmp1, map, by.x = 'Var2', by.y = 'idx')
  setnames(tmp1, c('idx_week', 'idx_basis'), c('idx_week_row', 'idx_basis_row'))
  
  tmp2 = subset(tmp1, idx_week_row %in% 1 & idx_week_column %in% 1)
  ggplot(tmp2, aes(x = Var1, y = Var2)) + 
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
tmp1 = tmp1[, list( 	q= quantile(value, prob=ps),
                     q_label=p_labs),
            by=c(column_name, row_name)]
tmp1 = dcast(tmp1, get(row_name) + get(column_name) ~ q_label, value.var = "q")
setnames(tmp1, c('row_name', 'column_name'), c(row_name, column_name))
tmp1 = merge(tmp1, df_week, by = row_name)

## Smooth estimate
z <- as.matrix( dcast(tmp1, get(column_name)~get(row_name), value.var = "M")[,-1] ) 
z.range <- range(z)
z =  t( apply(z, 2, rev) )

pdf(paste0(outdir, '-PlanePosterior_', Code, '.pdf'),width=7,height=7,paper='special')
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
  ggplot(probability_ratio_table, aes(x = date)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) +
    geom_line(aes(y = M)) + 
    geom_point(aes(y = emp.prob.ratio), col = 'red') + 
    theme_bw() + 
    geom_hline(aes(yintercept = 1)) +
    facet_wrap(~age, scales = 'free')+ 
    labs(x = '', y = 'ratio probability of one additional death /n to the first week')
  ggsave(file = paste0(outdir, '-ProbabilityRatio_', Code, '.png'), w = 10, h = 10)
}

# 
# plot_prediction = function()
# {
#   source(file.path(indir, "functions", "postprocessing-summary_functions.R"))
#   
#   suffix = '_ICAR_2'
#   
#   predictive_checks_table = make_predictive_checks_table(fit_cum, "deaths_cum", df_week, df_age_reporting, tmp)
#   
#   tmp = subset(predictive_checks_table, week_index %in% c(1, stan_data$W))
#   
#   ggplot(tmp, aes(x = age)) + 
#     geom_point(aes(y = M_deaths_cum)) + 
#     geom_errorbar(aes(ymin = CL_deaths_cum, ymax = CU_deaths_cum)) + 
#     geom_point(aes(y = COVID.19.Deaths), col = 'red') + 
#     facet_wrap(~date, ncol = 1) + 
#     theme_bw() + 
#     labs(y = 'Cumulative deaths', x = '')
#   ggsave(file = paste0('~/Box\ Sync/2021/CDC/posterior_predictive_checks_AL', suffix, '.png'), w = 6, h = 5)
# }
# 
# plot_contribution = function()
# {
#   samples = extract(fit_cum)
#   
#   suffix = '_ICAR_2'
#   
#   ps <- c(0.5, 0.025, 0.975)
#   p_labs <- c('M','CL','CU')
#   tmp1 = as.data.table( reshape2::melt( samples$phi ))
#   setnames(tmp1, c('Var2', 'Var3'), c('age', 'week_index'))
#   tmp1 = tmp1[, list( 	q= quantile(value, prob=ps),
#                        q_label=p_labs), 
#               by=c('age', 'week_index')]	
#   tmp1 = dcast(tmp1, week_index + age ~ q_label, value.var = "q")
# 
#   tmp1 = merge(tmp1, df_week, by = row_name)
#   
#   ggplot(tmp1, aes(x = date, y = age)) +
#     geom_raster(aes(fill = M))  + 
#     labs(x = 'Date', y = 'Age', fill = 'Estimated posterior value') + 
#     scale_y_continuous(expand = c(0,0))  +
#     scale_x_date(expand = c(0,0)) + 
#     scale_fill_viridis_c(option = "E") + 
#     theme(legend.position='bottom')
#   ggsave(file = paste0('~/Box\ Sync/2021/CDC/continuous_contribution_AL_all_week', suffix, '.png'), w = 6, h = 6.2)
#   
#   tmp1 = subset(tmp1, week_index %in% c(1, stan_data$W))
#   
#   ggplot(tmp1, aes(x = age)) + 
#     geom_line(aes(y = M)) + 
#     geom_ribbon(aes(ymin = CL, ymax = CU), alpha  = 0.5) + 
#     facet_wrap(~date, ncol = 1) + 
#     theme_bw() + 
#     labs(y = 'Probability that one additional deaths \n falls in a', x = 'a')
#   ggsave(file = paste0('~/Box\ Sync/2021/CDC/continuous_contribution_AL', suffix, '.png'), w = 6, h = 5)
#   
# }