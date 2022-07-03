format_posterior <- function(MetrHastrw_outputs){
  
  posterior <- NULL
  vars <- c('intercept_resurgence0', 'sigma_intercept_resurgence', 'slope_resurgence0', 
            'vaccine_effect_intercept_diagonal', 'vaccine_effect_slope_diagonal', 
            'vaccine_effect_intercept_cross', 'vaccine_effect_slope_cross', 
            'sigma_r_pdeaths')
  for(i in 1:length(vars)){
    tmp <- as.data.table( do.call('rbind', MetrHastrw_outputs[[vars[i]]]) )
    M <- dim(tmp)[2]
    setnames(tmp, 1:M, paste0(vars[i], '[', 1:M, ']'))
    tmp[, iterations := as.integer(1:nrow(tmp))]
    posterior <- rbind(posterior, melt.data.table(tmp, id.vars = 'iterations'))
  }
  
  tmp <- data.table( reshape2::melt(MetrHastrw_outputs$intercept_resurgence_re) )
  setnames(tmp, 1:4, c('state_index', 'age_index', 'value', 'iterations'))
  tmp[, variable := paste0('intercept_resurgence_re[', age_index, ',', state_index, ']')]
  posterior <- rbind(posterior, tmp, fill = TRUE)
  
  max_iterations <- max(posterior$iterations)
  posterior <- posterior[iterations <= floor(max_iterations / 8) * 8 ]
  posterior[, chain := findInterval(iterations, seq(1, max_iterations,length.out = 8)), by = 'iterations']
  
  return(posterior)
}

make_convergence_diagnostics = function(res, value_name){
  
  # Convergence Rhat and effective sample size 
  n_chain = max(res$chain)
  mcmc_samples = vector(mode = 'list', length = n_chain)
  for(i in 1:n_chain) mcmc_samples[[i]] = coda::mcmc(subset(res, chain == i)[, get(value_name)])
  
  convergence_diagnostics = c(Rhat = coda::gelman.diag(mcmc.list(mcmc_samples))$psrf[1], neff = as.numeric(coda::effectiveSize(mcmc.list(mcmc_samples))))
  
  return(convergence_diagnostics)
}

save_convergence_diagnostics <- function(posterior, outdir){
  
  convergence_diagnostics <- lapply(unique(posterior$variable), function(i) make_convergence_diagnostics(posterior[variable == i], 'value') )
  names(convergence_diagnostics) = unique(posterior$variable)

  cat('\n\nConvergence Diagnostics\n')
  print(convergence_diagnostics)
  
  file = paste0(outdir, "-ConvergenceDiagnostics.rds")
  cat('\n Write file ', file, '\n')
  saveRDS(convergence_diagnostics, file=file)
  
}

make_trace_plots_parameters = function(posterior_params, outdir){ 
  
  # trace plots of the epidemiological parameter 
  tmp <- posterior_params[!grepl('_re', variable)]
  p = ggplot(tmp, aes(x = iterations, y = value, colour = chain))+
    geom_step(alpha = 0.7)+
    theme_bw() +
    xlab("Iteration") +
    facet_wrap(~variable, scales = 'free_y') +
    theme(legend.position = 'bottom')+
    guides(colour = guide_legend(nrow = 1, title.position="top", override.aes = list(size=1)))
  
  file = paste0(outdir, "-trace_plots_theta_wore.png")
  cat('\n Write file ', file, '\n')
  ggsave(p, file = file, w = 12, h = 10)
  
  # trace plots of the epidemiological parameter 
  tmp <- posterior_params[grepl('_re', variable)]
  p = ggplot(posterior_params, aes(x = iterations, y = value, colour = chain))+
    geom_step(alpha = 0.7)+
    theme_bw() +
    xlab("Iteration") +
    facet_wrap(~variable, scales = 'free_y') +
    theme(legend.position = 'bottom')+
    guides(colour = guide_legend(nrow = 1, title.position="top", override.aes = list(size=1)))
  
  file = paste0(outdir, "-trace_plots_theta_re.png")
  cat('\n Write file ', file, '\n')
  ggsave(p, file = file, w = 13, h = 10)
  
}

make_forest_plot_table <- function(summary, df_age_vaccination2, df_state, names, math_name, groups, groups_levels){
  tmp <- summary[ grepl(paste(paste0('^',names),collapse = '|'),rownames(summary)) ,]
  # tmp1 <- tmp[grepl('diagonal', rownames(tmp)),]
  # tmp <- tmp[!grepl('diagonal', rownames(tmp)),]
  
  variables <- tmp$variable; group = c()
  for(x in 1:length(variables)) group[x] = groups[which(grepl(gsub('(.+)\\[.*', '\\1', variables[x]), names))]
  for(x in 1:length(names)) variables = gsub(names[x], math_name[x], variables)
  for(x in df_age_vaccination2$age_index) variables = gsub(paste0('\\[',x), paste0('\\["', df_age_vaccination2$age[x], '"'), variables) 
  # for(x in df_age_vaccination2$age_index) variables[grepl('vac', variables) & !grepl(', ', variables)] = sub(paste0('\\[\"',df_age_vaccination2$age[x]), paste0('\\[\"', df_age_vaccination2$age[-x],', ',df_age_vaccination2$age[x]), variables[grepl('vac', variables) & !grepl(', ', variables)])
  # for(x in df_age_vaccination2$age_index) variables[grepl('vac', variables)] = gsub('(.+),.*', '\\1"\\]', variables[grepl('vac', variables)])
  
  for(x in df_state$state_index) variables[!grepl('vac', variables)] = gsub(paste0('\",', x, '\\]'), paste0(', ',df_state$loc_label[x], '"\\]'), variables[!grepl('vac', variables)]) 
  
  variables = gsub('\\+\\+', '\\+', variables)
  variables = gsub('\\+', 'plus', variables)
  
  age1 = gsub('(.+),.*', '\\1', gsub('.*\\["(.+)','\\1', variables))
  loc1 = gsub('.*, (.+)"\\]','\\1', variables)
  idx.swap = which(!grepl('\\]', age1) & !grepl('[0-9]', loc1))
  
  for(x in 1:length(idx.swap)){
    
    variables[idx.swap[x]] = gsub(age1[idx.swap[x]],paste0(age1[idx.swap[x]], '1'), variables[idx.swap[x]])
    variables[idx.swap[x]] = gsub(loc1[idx.swap[x]],age1[idx.swap[x]], variables[idx.swap[x]])
    variables[idx.swap[x]] = gsub(paste0(age1[idx.swap[x]], '1'),loc1[idx.swap[x]], variables[idx.swap[x]])
    
  }
  
  variables = gsub('plus', '+', variables)
  
  tmp <- as.data.table(tmp)
  tmp[, variable := variables]
  tmp[, group := factor(group, levels = groups_levels)]
  
  setnames(tmp, c('50%', '2.5%', '97.5%'), c('M', 'CL', "CU"))
  
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
  ggsave(p, file = paste0(outdir, '-forest_plot_mcmc.png'), w = 8, h = 7)
  
  return(p)
}

save_confidence_intervals <- function(posterior){
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('50%', '2.5%', '97.5%')
  
  tmp1 = posterior[, list( 	q= quantile(value, prob=ps, na.rm = T),
                            q_label=p_labs), by=c('variable')]	
  tmp1 = dcast(tmp1, variable~ q_label, value.var = "q")
  rownames(tmp1) <- tmp1$variable
  
  return(tmp1)
}