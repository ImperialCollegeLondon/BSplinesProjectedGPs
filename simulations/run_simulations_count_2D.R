library(mvtnorm)
library(rstan)
library(ggplot2)
library(data.table)
library(gridExtra)
library("plgp")
library(dplyr)
library(loo)
library(grid)

indir ="~/git/covid19Vaccination" # path to the repo
outdir = file.path(indir, 'simulations', 'results')
dir.create(outdir)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load functions
source(file.path(indir, 'simulations', "functions", "utils_2D.R"))
source(file.path(indir, 'inst', "functions", "stan_utility_functions.R"))

# load stan models functions
# standard GP
model_GP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'GP_count_2D.stan') )
# standard B-splines
model_BS = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'B-SPLINES_count_2D.stan') )
# bayesian P-splines
model_PS = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'P-SPLINES_count_2D.stan') )
# regularised B-splines projected GP
model_GPBS = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'GP-B-SPLINES_count_2D.stan') )


## tuning parameters

# number of samples
n <- 40
x <- seq(0,1,length=n)
X <- expand.grid(x, x)
X_index <- expand.grid(1:n, 1:n)
X_index$index = 1:nrow(X_index)
lengthscales = c(0.05, 0.1, 1)

# percentage of data included in training
prop_training = 0.4

# spline parameters 
n_knots_vec = c(10, 30, 40)
spline_degree = 3


## Run simulations
nus = c(0.1, 0.1, 0.01)
for(i in 1:length(lengthscales)){
  
  # i = 2
  cat('\n Running simulation with lengthscale ', lengthscales[i], '\n')
  
  # lab
  lab = paste0('count_lengthscale_', lengthscales[i])
  
  # simulate data 
  set.seed(28)
  y_mean = exp( generate_2DGP(X, lengthscales[i]) )
  y = find_count_2D(mean = y_mean, nu = nus[i])
  coordinates_training = X_index[sample(1:nrow(X), size = round(prop_training*nrow(X))),]
  
  # fit
  cat('\n Using a Gaussian Process \n')
  assign(paste0('GP_2D_', i), run_spatial_model_2D(x_1 = x, x_2 = x, 
                                                   coordinates_training = coordinates_training, 
                                                   y = y, y_mean = y_mean, 
                                                   lab = lab, method = 'GP', stan_model = model_GP, 
                                                   outdir = outdir, overwrite = F))
  
  for(j in 1:length(n_knots_vec)){
    #j = 1
    cat('\n Using B-splines with ',  n_knots_vec[j], 'knots \n')
    assign(paste0('BS_2D_', j, '_', i),
           run_spatial_model_2D(x_1 = x, x_2 = x, 
                                coordinates_training = coordinates_training, 
                                y = y, y_mean = y_mean, 
                                lab = lab, method = 'B-SPLINES', stan_model = model_BS, 
                                n_knots = n_knots_vec[j], spline_degree = spline_degree,
                                outdir = outdir, overwrite = F))
    
    cat('\n Using Bayesian P-splines with ',  n_knots_vec[j], 'knots \n')
    assign(paste0('PS_2D_', j, '_', i),
           run_spatial_model_2D(x_1 = x, x_2 = x, 
                                coordinates_training = coordinates_training, 
                                y = y, y_mean = y_mean, 
                                lab = lab, method = 'P-SPLINES', stan_model = model_PS, 
                                n_knots = n_knots_vec[j], spline_degree = spline_degree,
                                outdir = outdir, overwrite = F))
    
    cat('\n Using regularised B-splines projected GP with ',  n_knots_vec[j], 'knots \n')
    assign(paste0('GPBS_2D_', j, '_', i), 
           run_spatial_model_2D(x_1 = x, x_2 = x, 
                                coordinates_training = coordinates_training, 
                                y = y, y_mean = y_mean, 
                                lab = lab, method = 'GP-B-SPLINES', stan_model = model_GPBS, 
                                n_knots = n_knots_vec[j], spline_degree = spline_degree,
                                outdir = outdir, overwrite = F))
    
  }
}



## plot results
tmp = do.call('rbind', list(GP_2D_1[[1]], 
                            GPBS_2D_1_1[[1]], GPBS_2D_2_1[[1]], GPBS_2D_3_1[[1]], 
                            BS_2D_1_1[[1]], BS_2D_2_1[[1]], BS_2D_3_1[[1]],  
                            PS_2D_1_1[[1]], PS_2D_2_1[[1]], PS_2D_3_1[[1]],
                            GP_2D_2[[1]],
                            GPBS_2D_1_2[[1]], GPBS_2D_2_2[[1]], GPBS_2D_3_2[[1]],
                            BS_2D_1_2[[1]], BS_2D_2_2[[1]], BS_2D_3_2[[1]],
                            PS_2D_1_2[[1]], PS_2D_2_2[[1]], PS_2D_3_2[[1]],
                            GP_2D_3[[1]],
                            GPBS_2D_1_3[[1]], GPBS_2D_2_3[[1]], GPBS_2D_3_3[[1]],
                            BS_2D_1_3[[1]], BS_2D_2_3[[1]], BS_2D_3_3[[1]],
                            PS_2D_1_3[[1]], PS_2D_2_3[[1]], PS_2D_3_3[[1]]))

tmp[grepl('GP', method) & !grepl('SPLINES', method), method := 'Standard 2D GP']
tmp[grepl('GP-B-SPLINES', method), method := gsub('GP-B-SPLINES', 'Low-rank 2D GP', method)]
tmp[grepl('P-SPLINES', method), method := gsub('P-SPLINES', 'P-splines', method)]
tmp[grepl('B-SPLINES', method), method := gsub('B-SPLINES', 'Standard B-splines 2D surface', method)]

tmp[grepl('knots', method), n_knots := as.numeric(gsub('.*\n#knots = (.+)', '\\1', method))]
tmp[grepl('knots', method), method2 := gsub('(.+)\n#knots = .*', '\\1', method)]
tmp[!grepl('knots', method), method2 := method]
tmp[!grepl('knots', method), n_knots := 0]

variables = c('y_hat', 'f')
for(l in lengthscales){
  # l = lengthscales[2]
  for(k in 1:length(variables)){
    # k = 1
    tmp1 = subset(tmp, variable == variables[k] & lengthscale == l)

    tmp1[, method2 := factor(method2, levels = c('Standard 2D GP', 
                                                 'Standard B-splines 2D surface', 
                                                 'P-splines',
                                                 'Low-rank 2D GP'))]
    tmp1[, n_knots := factor(n_knots, levels = sort(unique(tmp1$n_knots), decreasing = F))]
    tmp1[, n_knots2:= factor(paste0(n_knots, ' knots'), levels = paste0(sort(unique(tmp1$n_knots), decreasing = F), ' knots'))]
    
    var = variables[k]
    y_var = ifelse(var == 'y_hat', 'y', 'y_mean')
    option_viridis = ifelse(var == 'y_hat', 'A', 'viridis')
    var_names = ifelse(var == 'y_hat', 'Simulated observations\n', 'Mean of the\nsimulated observations\n')
    var_names_training = ifelse(var == 'y_hat', 'Simulated observations\nincluded in the training set', 'Mean of the\nsimulated observations\nincluded in the training set')
    prediction_names = ifelse(var == 'y_hat', 'Predicted observations', 'Estimated mean of the observations')
    prediction_names_break  = ifelse(var == 'y_hat', 'Predicted observations', 'Estimated mean of the\nobservations')
    
    tmp2 = select(tmp1, x_1, x_2, y, y_mean)
    p0 = ggplot(unique(tmp2),aes(x=x_1,y=x_2)) +
      geom_raster(aes(fill=all_of(get(y_var))), interpolate = TRUE) +
      geom_contour(aes(z=all_of(get(y_var))), bins = 12, color = "gray30", 
                   size = 0.5, alpha = 0.5) +
      scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_viridis_c(option = option_viridis, limits = range(c(tmp1$M, as.vector(select(tmp1,y_var) ) )), begin = 0.1) +
      labs(x = '') + 
      ggtitle(var_names) + 
      theme(legend.position = 'none',
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_blank(),
            strip.text =  element_blank(),
            panel.spacing.x = unit(1, "lines"), 
            axis.title.y = element_blank(),
            plot.title = element_text(size = rel(1), hjust = 0.5)) 
    
    tmp2 = copy(select(tmp1, x_1, x_2, y, y_mean, y_training))
    tmp2[is.na(y_training), (y_var) := NA]
    p02 = ggplot(unique(tmp2),aes(x=x_1,y=x_2)) +
      geom_raster(aes(fill=all_of(get(y_var))), interpolate = TRUE) +
      scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_viridis_c(option = option_viridis, limits = range(c(tmp1$M, as.vector(select(tmp1,y_var) ) )), begin = 0.1) +
      labs(x = '') + 
      ggtitle(var_names_training) + 
      theme(legend.position = 'none',
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_blank(),
            strip.text =  element_blank(),
            panel.spacing.x = unit(1, "lines"), 
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            plot.title = element_text(size = rel(1), hjust = 0.5)) 
    
    tmp2 = subset(tmp1, method2 %in% c('Standard 2D GP'))
    p1 =ggplot(tmp2,aes(x=x_1,y=x_2)) +
      geom_raster(aes(fill=M), interpolate = TRUE) +
      geom_contour(aes(z=M), bins = 12, color = "gray30", 
                   size = 0.5, alpha = 0.5) +
      scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_viridis_c(option = option_viridis, limits = range(c(tmp1$M, as.vector(select(tmp1,y_var) ))), begin = 0.1) +
      labs(x = '', y= '') + 
      facet_wrap(.~method2, scale = 'free_y') +
      ggtitle(paste0('Standard 2D GP\n',prediction_names_break)) + 
      theme(legend.position = 'none',
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            strip.text =  element_blank(),
            panel.spacing.x = unit(1, "lines"), 
            plot.title = element_text(size = rel(1), hjust = 0.5)) 

    tmp2 = subset(tmp1, method2 == c('Standard B-splines 2D surface'))
    p2 = ggplot(tmp2,aes(x=x_1,y=x_2)) +
      geom_raster(aes(fill=M), interpolate = TRUE) +
      geom_contour(aes(z=M), bins = 12, color = "gray30", 
                   size = 0.5, alpha = 0.5) +
      coord_equal() +
      scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_viridis_c(option = option_viridis, limits = range(c(tmp1$M, as.vector(select(tmp1,y_var) ) )), begin = 0.1) +
      # ggtitle(paste0('Standard B-splines surface\n',prediction_names)) +
      ggtitle(paste0('Standard B-splines surface')) +
      facet_grid(method2~n_knots2) +
      theme(legend.position = 'none',
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA), 
            strip.text.y = element_blank(), 
            strip.text.x =  element_text(size = rel(1)),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            # panel.spacing.x = unit(0.6, "lines"), 
            plot.title =element_text(hjust = 0.5,size=rel(1), vjust = -1))  

    tmp2 = subset(tmp1, method2 == c('P-splines'))
    p3 = ggplot(tmp2,aes(x=x_1,y=x_2)) +
      geom_raster(aes(fill=M), interpolate = TRUE) +
      geom_contour(aes(z=M), bins = 12, color = "gray30", 
                   size = 0.5, alpha = 0.5) +
      coord_equal() +
      scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_viridis_c(option = option_viridis, limits = range(c(tmp1$M, as.vector(select(tmp1,y_var) ))), begin = 0.1) +
      # ggtitle(paste0('Bayesian P-splines\n',prediction_names)) +
      ggtitle(paste0('Bayesian P-splines')) +
      facet_grid(method2~n_knots2) +
      theme(legend.position = 'none',
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            strip.text.y = element_blank(), 
            strip.text.x =  element_text(size = rel(1.1)),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.spacing.y = unit(3, "lines"), 
            plot.title =element_text(hjust = 0.5,size=rel(1), vjust = -1))

    tmp2 = subset(tmp1, method2 == c('Low-rank 2D GP'))
    p4 = ggplot(tmp2,aes(x=x_1,y=x_2)) +
      geom_raster(aes(fill=M), interpolate = TRUE) +
      geom_contour(aes(z=M), bins = 12, color = "gray30", 
                   size = 0.5, alpha = 0.5) +
      coord_equal() +
      scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_viridis_c(option = option_viridis, limits = range(c(tmp1$M, as.vector(select(tmp1,y_var) ))), begin = 0.1) +
      ggtitle(paste0('Regularised B-splines projected 2D GP')) +
      # ggtitle(paste0('Regularised B-splines projected 2D GP\n',prediction_names)) +
      facet_grid(method2~n_knots2) +
      theme(legend.position = 'none',
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA), 
            strip.text.y = element_blank(), 
            axis.title.y = element_blank(),
            strip.text.x =  element_text(size = rel(1.1)),
            axis.title.x = element_blank(),
            panel.spacing.y = unit(3, "lines"), 
            plot.title =element_text(hjust = 0.5,size=rel(1), vjust = -1))

    p01  <- p0 + theme( plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt") )
    p021  <- p02 + theme( plot.margin=unit(c(5.5, 5.5, 0, 0), "pt") )
    p11  <- p1 + theme( plot.margin=unit(c(5.5, 5.5, 0, 0), "pt"))
    p00 = ggpubr::ggarrange(p01,p021,p11, widths = c(1.14, 1, 1), nrow =1, labels = c('A', 'B', 'C'),  font.label = list(size = 15))
    ggsave(p00, file = file.path(outdir, paste0('v2_2D_comp_count_lengthscale_', l, '_', var, '_data.png')), w = 8, h = 3.25)
    

    p31 <- p3+ theme( plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt") )
    p41 <- p4+ theme( plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt") )
    p21 <- p2+ theme( plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt") )
    
    p = ggpubr::ggarrange(p21, p31, p41, nrow = 3, labels = c('A', 'B', 'C'),  font.label = list(size = 15))
    ggsave(p, file = file.path(outdir, paste0('v2_2D_comp_count_lengthscale_', l, '_', var, '_Rsplines.png')), w = 7, h = 8.5)
    
# 
#     p = grid.arrange(p0,p02,p1, p2, p3, p4, 
#                      layout_matrix = rbind(c(NA,1,1, 2, 3), 
#                                            rep(4,5),
#                                            rep(5,5),
#                                            rep(6,5)), 
#                      widths = c(0.1,0.15, 1, 1, 1), heights = c(0.95,1, 1, 1))
# 
# 
#     ggsave(p, file = file.path(outdir, paste0('v2_2D_comp_count_lengthscale_', l, '_', var, '.png')), w = 6, h = 10)
    
    
  }
    
 }



## extract statistics for the paper

# time of execution 
time_all = list()
for(i in seq_along(lengthscales)){
  
  l = lengthscales[i]
  tmp1 = subset(tmp, variable == 'y_hat' & lengthscale == l)
  tmp1 = tmp1[, list(time = unique(time)), by = c('method', 'n_knots')]
  timeGP = tmp1[method == 'Standard 2D GP', time]
  
  tmp1 = tmp1[, timed := paste0(round((tmp1$time - timeGP) / 60)) ]
  tmp1 = tmp1[, timep := paste0(gsub(" ", "", format(round( (((tmp1$time - timeGP) / timeGP) )* 100, digits = 2), nsmall =2)), '\\%') ]
  tmp1 = tmp1[,timec := paste0(round((tmp1$time) / 60)) ]
  
  time_all[[i]] = tmp1[,c(4,5,6)]
}

tmp1 = subset(tmp, grepl('Standard 2D GP', method))
min_GP = unique(tmp1$time)
time_GP = paste0(round(min_GP / 60), ' minutes')

tmp1 = subset(tmp, grepl('Low-rank 2D GP', method))
min_GPSE_1 = unique(subset(tmp1, lengthscale == lengthscales[1] & n_knots == 30)$time)
min_GPSE_2 = unique(subset(tmp1, lengthscale == lengthscales[2] & n_knots == 30)$time)
min_GPSE_3 = unique(subset(tmp1, lengthscale == lengthscales[3] & n_knots == 30)$time)
time_GPSE = list()
time_GPSE[[1]] = paste0(round(min_GPSE_1/ 60), ' minutes')
time_GPSE[[2]] = paste0(round(min_GPSE_2/ 60), ' minutes')
time_GPSE[[3]] = paste0(round(min_GPSE_3/ 60), ' minutes')

avg_red = paste0(round(mean(1- (c(min_GPSE_1/min_GP[1], min_GPSE_2/min_GP[2], min_GPSE_3/min_GP[3] ))* 100), digits = 2) , '\\%')

saveRDS(list(time_GP, time_GPSE, avg_red, time_all), file = file.path(outdir, paste0('v2_time_execution_count.rds')))


# compare MSE
tmp[, SE := (M - y_mean)^2]
tmp[, AE := abs(M - y_mean)]
MAE1 <- subset(tmp, variable == 'f'  & lengthscale == lengthscales[1])[, list(mean = format(round(mean(AE), 2), nsmall = 2),
                                                                              sd = format(round(sd(AE), 2), nsmall = 2)), by = c('method')]
MAE2 <- subset(tmp, variable == 'f'  & lengthscale == lengthscales[2])[, list(mean = format(round(mean(AE), 2), nsmall = 2), 
                                                                              sd = format(round(sd(AE), 2), nsmall = 2)), by = c('method')]
MAE3 <- subset(tmp, variable == 'f'  & lengthscale == lengthscales[3])[, list(mean = format(round(mean(AE), 2), nsmall = 2), 
                                                                              sd = format(round(sd(AE), 2), nsmall = 2)), by = c('method')]

MAE1[order(mean)]
MAE2[order(mean)]
MAE3[order(mean)]

saveRDS(list(MAE1, MAE2, MAE3), file = file.path(outdir, paste0('v2_MAE_count.rds')))


## compare loo

# first scenario
com1 = loo_compare(loo(extract(GP_2D_1[[2]])$log_lik), 
                   loo(extract(GPBS_2D_1_1[[2]])$log_lik), loo(extract(GPBS_2D_2_1[[2]])$log_lik), loo(extract(GPBS_2D_3_1[[2]])$log_lik), 
                   loo(extract(BS_2D_1_1[[2]])$log_lik), loo(extract(BS_2D_2_1[[2]])$log_lik), loo(extract(BS_2D_3_1[[2]])$log_lik), 
                   loo(extract(PS_2D_1_1[[2]])$log_lik), loo(extract(PS_2D_2_1[[2]])$log_lik), loo(extract(PS_2D_3_1[[2]])$log_lik))

# second scenario
com2 = loo_compare(loo(extract(GP_2D_2[[2]])$log_lik), 
                   loo(extract(GPBS_2D_1_2[[2]])$log_lik), loo(extract(GPBS_2D_2_2[[2]])$log_lik), loo(extract(GPBS_2D_3_2[[2]])$log_lik), 
                   loo(extract(BS_2D_1_2[[2]])$log_lik), loo(extract(BS_2D_2_2[[2]])$log_lik), loo(extract(BS_2D_3_2[[2]])$log_lik), 
                   loo(extract(PS_2D_1_2[[2]])$log_lik), loo(extract(PS_2D_2_2[[2]])$log_lik), loo(extract(PS_2D_3_2[[2]])$log_lik))

# third scenario
com3 = loo_compare(loo(extract(GP_2D_3[[2]])$log_lik), 
                   loo(extract(GPBS_2D_1_3[[2]])$log_lik), loo(extract(GPBS_2D_2_3[[2]])$log_lik), loo(extract(GPBS_2D_3_3[[2]])$log_lik), 
                   loo(extract(BS_2D_1_3[[2]])$log_lik), loo(extract(BS_2D_2_3[[2]])$log_lik), loo(extract(BS_2D_3_3[[2]])$log_lik), 
                   loo(extract(PS_2D_1_3[[2]])$log_lik), loo(extract(PS_2D_2_3[[2]])$log_lik), loo(extract(PS_2D_3_3[[2]])$log_lik))


# extract statistics for the paper
com11 = as.data.table(com1) 
com11$model = rownames(com1)
com11[, model_index := as.numeric(gsub('model(.+)', '\\1', model))]
com11 = com11[order(model_index)]
com11 = select(com11, elpd_diff, se_diff)
com11[, elpd_diff:= gsub(" ", "", format(round(elpd_diff, digits = 2), nsmall =2)) ]
com11[, se_diff:= gsub(" ", "", format(round(se_diff, digits = 2), nsmall =2)) ]

com22  = as.data.table(com2) 
com22$model = rownames(com2)
com22[, model_index := as.numeric(gsub('model(.+)', '\\1', model))]
com22 = com22[order(model_index)]
com22 = select(com22, elpd_diff, se_diff)
com22[, elpd_diff:= gsub(" ", "", format(round(elpd_diff, digits = 2), nsmall =2)) ]
com22[, se_diff:= gsub(" ", "", format(round(se_diff, digits = 2), nsmall =2)) ]

com33  = as.data.table(com3) 
com33$model = rownames(com3)
com33[, model_index := as.numeric(gsub('model(.+)', '\\1', model))]
com33 = com33[order(model_index)]
com33 = select(com33, elpd_diff, se_diff)
com33[, elpd_diff:= gsub(" ", "", format(round(elpd_diff, digits = 2), nsmall =2)) ]
com33[, se_diff:= gsub(" ", "", format(round(se_diff, digits = 2), nsmall =2)) ]


saveRDS(list(com11, com22, com33), file = file.path(outdir, paste0('v2_LOO_comp_count_all.rds')))



