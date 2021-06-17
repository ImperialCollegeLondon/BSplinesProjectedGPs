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

model_GP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'GP-SE_2D.stan') )
model_BSGP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-GP-SE_2D.stan') )
model_BSIN = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-GP-I_2D.stan') )

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source(file.path(indir, 'simulations', "functions", "utils_2D.R"))
source(file.path(indir, 'inst', "functions", "stan_utility_functions.R"))


# number of samples
n <- 50
x <- seq(0,1,length=n)
X <- expand.grid(x, x)
sigma = 0.01
lengthscales = c(0.01, 0.1, 1)

# spline parameters 
spline_degree = 3

##
n_knots_vec = c(10, 30, 40)
lab = paste0('lengthscale_', lengthscales[1])
set.seed(23)
y = generate_2DGP(X, lengthscales[1], sigma = 0.5)
GP_2D_1 = run_GP_2D(x, x, y, lab, outdir)
BSGP_2D_1_1 = run_BSGP_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSGP_2D_2_1 = run_BSGP_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSGP_2D_3_1 = run_BSGP_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSIN_2D_1_1 = run_BSIN_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSIN_2D_2_1 = run_BSIN_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSIN_2D_3_1 = run_BSIN_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)

##
n_knots_vec = c(10, 30, 40)
lab = paste0('lengthscale_', lengthscales[2])
set.seed(24)
y = generate_2DGP(X, lengthscales[2], sigma = 0.5)
GP_2D_2 = run_GP_2D(x, x, y, lab, outdir)
BSGP_2D_1_2 = run_BSGP_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSGP_2D_2_2 = run_BSGP_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSGP_2D_3_2 = run_BSGP_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSIN_2D_1_2 = run_BSIN_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSIN_2D_2_2 = run_BSIN_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSIN_2D_3_2 = run_BSIN_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)

##
n_knots_vec = c(10, 30, 40)
lab = paste0('lengthscale_', lengthscales[3])
set.seed(23)
y = generate_2DGP(X, lengthscales[3], sigma= 0.1)
GP_2D_3 = run_GP_2D(x, x, y, lab, outdir)
BSGP_2D_1_3 = run_BSGP_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSGP_2D_2_3 = run_BSGP_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSGP_2D_3_3 = run_BSGP_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSIN_2D_1_3 = run_BSIN_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSIN_2D_2_3 = run_BSIN_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSIN_2D_3_3 = run_BSIN_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)


# plot results
tmp = do.call('rbind', list(GP_2D_1[[1]], 
                            BSGP_2D_1_1[[1]], BSGP_2D_2_1[[1]], BSGP_2D_3_1[[1]], 
                            BSIN_2D_1_1[[1]], BSIN_2D_2_1[[1]], BSIN_2D_3_1[[1]],  
                            GP_2D_2[[1]], 
                            BSGP_2D_1_2[[1]], BSGP_2D_2_2[[1]], BSGP_2D_3_2[[1]], 
                            BSIN_2D_1_2[[1]], BSIN_2D_2_2[[1]], BSIN_2D_3_2[[1]],
                            GP_2D_3[[1]], 
                            BSGP_2D_1_3[[1]], BSGP_2D_2_3[[1]], BSGP_2D_3_3[[1]], 
                            BSIN_2D_1_3[[1]], BSIN_2D_2_3[[1]], BSIN_2D_3_3[[1]]))

tmp[grepl('BS-GP-SE', method), method := gsub('BS-GP-SE', 'Low-rank 2D GP', method)]
tmp[grepl('GP-SE', method), method := gsub('GP-SE', 'Standard 2D GP', method)]
tmp[grepl('BS-GP-IN', method), method := gsub('BS-GP-IN', 'Standard B-splines 2D surface', method)]

tmp[grepl('knots', method), n_knots := as.numeric(gsub('.*\n#knots = (.+)', '\\1', method))]
tmp[grepl('knots', method), method2 := gsub('(.+)\n#knots = .*', '\\1', method)]
tmp[!grepl('knots', method), method2 := method]
tmp[!grepl('knots', method), n_knots := 0]


for(l in lengthscales){
  # l = lengthscales[1]
  tmp1 = subset(tmp, variable == 'y_hat' & lengthscale == l)
  
  tmp2 = copy(select(tmp1, x_1, x_2, mean, method2, n_knots, time, lengthscale))
  tmp1 = copy(select(tmp1, x_1, x_2, y, time, n_knots, lengthscale))
  tmp1[, time := NA_real_]
  tmp1[, method2 := 'Observations']
  setnames(tmp1, 'y', 'mean')
  tmp1 = rbind(unique(tmp1), tmp2)
  
  tmp1[, method2 := factor(method2, levels = c('Observations','Standard 2D GP', 
                                               'Standard B-splines 2D surface', 
                                               'Low-rank 2D GP'))]
  tmp1[, n_knots := factor(n_knots, levels = sort(unique(tmp1$n_knots), decreasing = F))]
  tmp1[, n_knots2:=paste0(n_knots, ' knots')]
  
  tmp2 = subset(tmp1, method2 %in% c('Observations'))
  p0 = ggplot(tmp2,aes(x=x_1,y=x_2)) +
    geom_raster(aes(fill=mean), interpolate = TRUE) +
    geom_contour(aes(z=mean), bins = 12, color = "gray30", 
                 size = 0.5, alpha = 0.5) +
    scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_viridis_c(option = "viridis", limits = range(tmp1$mean)) +
    # labs(x = expression(x[2]), y = expression(x[1])) + 
    labs(x = '', y= '') + 
    facet_wrap(.~method2, scale = 'free_y') +
    ggtitle('Observations') + 
    # facet_wrap(~ method2, nrow = 1, strip.position = "right") +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          # axis.text.x = element_text(angle = 45,hjust=1,vjust=1), 
          axis.title.x = element_blank(),
          # axis.ticks = element_blank(),
          # axis.text = element_blank(),
          strip.text =  element_blank(),
          panel.spacing.x = unit(1, "lines"), 
          plot.title = element_text(size = rel(1.1), hjust = 0.5)) 
  p0 = ggpubr::ggarrange(p0, labels = 'A', font.label = list(size = 18))
  
  tmp2 = subset(tmp1, method2 %in% c('Standard 2D GP'))
  p1 =ggplot(tmp2,aes(x=x_1,y=x_2)) +
    geom_raster(aes(fill=mean), interpolate = TRUE) +
    geom_contour(aes(z=mean), bins = 12, color = "gray30", 
                 size = 0.5, alpha = 0.5) +
    scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_viridis_c(option = "viridis", limits = range(tmp1$mean)) +
    # labs(x = expression(x[2]), y = expression(x[1])) + 
    labs(x = '', y= '') + 
    facet_wrap(.~method2, scale = 'free_y') +
    ggtitle('Standard 2D GP') + 
    # facet_wrap(~ method2, nrow = 1, strip.position = "right") +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          # axis.text.x = element_text(angle = 45,hjust=1,vjust=1), 
          axis.title.x = element_blank(),
          # axis.ticks = element_blank(),
          # axis.text = element_blank(),
          strip.text =  element_blank(),
          panel.spacing.x = unit(1, "lines"), 
          plot.title = element_text(size = rel(1.1), hjust = 0.5)) 
  p1 = ggpubr::ggarrange(p1, labels = 'B', font.label = list(size = 18))
  
  tmp2 = subset(tmp1, method2 == c('Standard B-splines 2D surface'))
  p2 = ggplot(tmp2,aes(x=x_1,y=x_2)) +
    geom_raster(aes(fill=mean), interpolate = TRUE) +
    geom_contour(aes(z=mean), bins = 12, color = "gray30", 
                 size = 0.5, alpha = 0.5) +
    coord_equal() +
    scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_viridis_c(option = "viridis", limits = range(tmp1$mean)) +
    # labs(x = expression(x[2]), y = expression(x[1])) +
    labs(y= '') + 
    # xlab(expression(decreasing %<->% 'increasing number of knots')) + 
    ggtitle('Standard B-splines 2D surface') + 
    facet_grid(method2~n_knots2) +
    # xlab(expression(decreasing %<->% 'increasing number of knots')) +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          # axis.text.x = element_text(angle = 45,hjust=1,vjust=1),
          # strip.text.x = element_text(colour = "white"), 
          strip.text.y = element_blank(), 
          strip.text.x =  element_text(size = rel(1.1)),
          # axis.ticks = element_blank(),
          # axis.text = element_blank(),
          axis.title.x = element_blank(),
          panel.spacing.y = unit(3, "lines"), 
          plot.title =element_text(hjust = 0.5,size=rel(1.1)))  
  p2 = ggpubr::ggarrange(p2, labels = 'C', font.label = list(size = 18))
  
  tmp2 = subset(tmp1, method2 == c('Low-rank 2D GP'))
  p3 = ggplot(tmp2,aes(x=x_1,y=x_2)) +
    geom_raster(aes(fill=mean), interpolate = TRUE) +
    geom_contour(aes(z=mean), bins = 12, color = "gray30", 
                 size = 0.5, alpha = 0.5) +
    coord_equal() +
    scale_x_continuous(expand=c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_viridis_c(option = "viridis", limits = range(tmp1$mean)) +
    # labs(x = expression(x[2]), y = expression(x[1])) +
    labs(y= '') + 
    # xlab(expression(decreasing %<->% 'increasing number of knots')) + 
    ggtitle('Low-rank 2D GP projected by regularised B-splines') + 
    facet_grid(method2~n_knots2) +
    # xlab(expression(decreasing %<->% 'increasing number of knots')) +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          # axis.text.x = element_text(angle = 45,hjust=1,vjust=1),
          # strip.text.x = element_text(colour = "white"), 
          strip.text.y = element_blank(), 
          strip.text.x =  element_text(size = rel(1.1)),
          # axis.ticks = element_blank(),
          # axis.text = element_blank(),
          axis.title.x = element_blank(),
          panel.spacing.y = unit(3, "lines"), 
          plot.title =element_text(hjust = 0.5,size=rel(1.1)))
  p3 = ggpubr::ggarrange(p3, labels = 'D', font.label = list(size = 18))
  
  p = grid.arrange(p0,p1, p2, p3, nrow = 5, 
                   layout_matrix = (rbind(c(NA,1,2,NA), 
                                                            c(NA, NA, NA, NA),
                                                         c(3,3,3,3),
                                                         c(NA, NA, NA,NA),
                                                         c(4,4,4,4))), widths = c(0.15, 0.4, 0.4, 0.15), 
               heights = c(0.9,0.1,1,0.1, 1))
  
  
  ggsave(p, file = file.path(outdir, paste0('2D_comp_lengthscale_', l, '.png')), w = 8, h = 10)
}

# time of execution 
time_all = list()
for(i in seq_along(lengthscales)){
  
  l = lengthscales[i]
  tmp1 = subset(tmp, variable == 'y_hat' & lengthscale == l)
  
  tmp1 = tmp1[, list(time = unique(time)), by = c('method', 'n_knots')]
  tmp1 = tmp1[, timed := paste0(round((tmp1$time - min(tmp1$time)) / 60), ' minutes') ]
  tmp1 = tmp1[, timep := paste0(gsub(" ", "", format(round( ((tmp1$time - min(tmp1$time)) / min(tmp1$time)), digits = 2), nsmall =2)), '\\%') ]

  time_all[[i]] = tmp1[,c(4,5)]
}

# GP vs GP-BS-SE
tmp1 = subset(tmp, grepl('Standard 2D GP', method))
min_GP = unique(tmp1$time)
time_GP = paste0(round(min_GP / 60), ' minutes')

tmp1 = subset(tmp, grepl('Low-rank 2D GP', method))
min_GPSE_1 = unique(subset(tmp1, lengthscale == lengthscales[1] & n_knots == 30)$time)
min_GPSE_2 = unique(subset(tmp1, lengthscale == lengthscales[2] & n_knots == 10)$time)
min_GPSE_3 = unique(subset(tmp1, lengthscale == lengthscales[3] & n_knots == 10)$time)
time_GPSE = list()
time_GPSE[[1]] = paste0(round(min_GPSE_1/ 60), ' minutes')
time_GPSE[[2]] = paste0(round(min_GPSE_2/ 60), ' minutes')
time_GPSE[[3]] = paste0(round(min_GPSE_3/ 60), ' minutes')

avg_red = paste0(round(mean((1-c(min_GPSE_1/min_GP[1], min_GPSE_2/min_GP[2], min_GPSE_3/min_GP[3] ))* 100), digits = 2) , '\\%')

saveRDS(list(time_GP, time_GPSE, avg_red, time_all), file = file.path(outdir, paste0('time_execution.rds')))


# compare loo
# first scenario
com1 = loo_compare(loo(extract(GP_2D_1[[2]])$log_lik), 
            loo(extract(BSGP_2D_1_1[[2]])$log_lik), loo(extract(BSGP_2D_2_1[[2]])$log_lik), loo(extract(BSGP_2D_3_1[[2]])$log_lik), 
            loo(extract(BSIN_2D_1_1[[2]])$log_lik), loo(extract(BSIN_2D_2_1[[2]])$log_lik), loo(extract(BSIN_2D_3_1[[2]])$log_lik))

# second scenario
com2 = loo_compare(loo(extract(GP_2D_2[[2]])$log_lik), 
            loo(extract(BSGP_2D_1_2[[2]])$log_lik), loo(extract(BSGP_2D_2_2[[2]])$log_lik), loo(extract(BSGP_2D_3_2[[2]])$log_lik), 
            loo(extract(BSIN_2D_1_2[[2]])$log_lik), loo(extract(BSIN_2D_2_2[[2]])$log_lik), loo(extract(BSIN_2D_3_2[[2]])$log_lik))

# third scenario
com3 = loo_compare(loo(extract(GP_2D_3[[2]])$log_lik), 
            loo(extract(BSGP_2D_1_3[[2]])$log_lik), loo(extract(BSGP_2D_2_3[[2]])$log_lik), loo(extract(BSGP_2D_3_3[[2]])$log_lik), 
            loo(extract(BSIN_2D_1_3[[2]])$log_lik), loo(extract(BSIN_2D_2_3[[2]])$log_lik), loo(extract(BSIN_2D_3_3[[2]])$log_lik))

com11  = as.data.table(com1) 
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


saveRDS(list(com11, com22, com33), file = file.path(outdir, paste0('LOO_comp_all.rds')))





# for text 
# LOO1 = list()
# LOO1[[1]] = round(loo_compare(loo(extract(BSGP_2D_2_1[[2]])$log_lik), loo(extract(BSIN_2D_2_1[[2]])$log_lik)), digits = 1)
# LOO1[[2]] = round(loo_compare(loo(extract(BSGP_2D_3_1[[2]])$log_lik), loo(extract(BSIN_2D_3_1[[2]])$log_lik)), digits = 1)
# 
# LOO1[[3]] = round(loo_compare(loo(extract(BSGP_2D_2_1[[2]])$log_lik), loo(extract(GP_2D_1[[2]])$log_lik)), digits = 1)
# LOO1[[4]] = round(loo_compare(loo(extract(BSGP_2D_3_1[[2]])$log_lik), loo(extract(GP_2D_1[[2]])$log_lik)), digits = 1)
# 
# LOO1[[5]] = round(loo_compare(loo(extract(BSGP_2D_1_3[[2]])$log_lik), loo(extract(GP_2D_3[[2]])$log_lik)), digits = 1)
# saveRDS(LOO1, file = file.path(outdir, paste0('LOO_1.rds')))
# 
# 
# LOO_GP = list(gsub(" ", "", format(round(loo(extract(GP_2D_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#               gsub(" ", "", format(round(loo(extract(GP_2D_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#               gsub(" ", "", format(round(loo(extract(GP_2D_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )
# 
# LOO_GP_I = list(gsub(" ", "", format(round(loo(extract(GP_I_2D_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#               gsub(" ", "", format(round(loo(extract(GP_I_2D_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#               gsub(" ", "", format(round(loo(extract(GP_I_2D_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )
# 
# LOO_BSGP_1 = list(gsub(" ", "", format(round(loo(extract(BSGP_2D_1_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#                   gsub(" ", "", format(round(loo(extract(BSGP_2D_1_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#                   gsub(" ", "", format(round(loo(extract(BSGP_2D_1_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )
# 
# LOO_BSGP_2 = list(gsub(" ", "",format(round(loo(extract(BSGP_2D_2_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#                   gsub(" ", "",format(round(loo(extract(BSGP_2D_2_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#                   gsub(" ", "",format(round(loo(extract(BSGP_2D_2_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )
# 
# LOO_BSGP_3 = list(gsub(" ", "",format(round(loo(extract(BSGP_2D_3_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#                   gsub(" ", "",format(round(loo(extract(BSGP_2D_3_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
#                   gsub(" ", "",format(round(loo(extract(BSGP_2D_3_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )
# 
# LOO_BSIN_1 = list( gsub(" ", "",format( round(loo(extract(BSIN_2D_1_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)),
#                    gsub(" ", "",format( round(loo(extract(BSIN_2D_1_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)),
#                    gsub(" ", "",format( round(loo(extract(BSIN_2D_1_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)))
# 
# LOO_BSIN_2 = list(gsub(" ", "", format( round(loo(extract(BSIN_2D_2_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2) ),
#                   gsub(" ", "", format( round(loo(extract(BSIN_2D_2_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)),
#                   gsub(" ", "", format( round(loo(extract(BSIN_2D_2_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)))
# 
# LOO_BSIN_3 = list(gsub(" ", "", format( round(loo(extract(BSIN_2D_3_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2) ),
#                   gsub(" ", "", format( round(loo(extract(BSIN_2D_3_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)),
#                   gsub(" ", "", format( round(loo(extract(BSIN_2D_3_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)))
# 
# LOO_ = list(LOO_GP_I, LOO_GP, LOO_BSGP_1, LOO_BSGP_2, LOO_BSGP_3, 
#      LOO_BSIN_1, LOO_BSIN_2, LOO_BSIN_3)
# saveRDS(LOO_, file = file.path(outdir, paste0('LOO_all.rds')))
# 
# 

# p = grid.arrange(grobs = list(p1, p2), 
#                  heights = c(0.05, 0.6,0.05,1.22, 0.05), 
#                  layout_matrix = rbind(c(NA,NA, NA), 
#                                        c(NA, 1, NA), 
#                                        c(NA,NA, NA),
#                                        c(2,2,2),
#                                        c(NA,NA, NA)), widths = c(0.,1, 0))
# 
# 
# pdf(file = file.path(outdir, paste0('2D_comp_lengthscale_', l, '.pdf')), width = 7.5,height    = 8.5)
# 
# grid.newpage()
# grid.draw(arrangeGrob(p))
# 
# grid.text(expression(decreasing %<->% 'increasing number of knots'), x = c(0.5), y = 0.365)
# grid.text(expression(decreasing %<->% 'increasing number of knots'), x = c(0.5), y = 0.04)
# 
# grid.segments(x0 = unit(0.04, "npc"), y0 = unit(0.025, "npc"), x1 = unit(0.98, "npc"), y1 = unit(0.025, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.04, "npc"), y0 = unit(0.025, "npc"), x1 = unit(0.04, "npc"), y1 = unit(0.33, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.04, "npc"), y0 = unit(0.33, "npc"), x1 = unit(0.98, "npc"), y1 = unit(0.33, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.98, "npc"), y0 = unit(0.025, "npc"), x1 = unit(0.98, "npc"), y1 = unit(0.33, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# 
# grid.segments(x0 = unit(0.04, "npc"), y0 = unit(0.35, "npc"), x1 = unit(0.04, "npc"), y1 = unit(0.653, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.04, "npc"), y0 = unit(0.35, "npc"), x1 = unit(0.98, "npc"), y1 = unit(0.35, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.04, "npc"), y0 = unit(0.653, "npc"), x1 = unit(0.98, "npc"), y1 = unit(0.653, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.98, "npc"), y0 = unit(0.35, "npc"), x1 = unit(0.98, "npc"), y1 = unit(0.653, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# 
# grid.segments(x0 = unit(0.05, "npc"), y0 = unit(0.67, "npc"), x1 = unit(0.05, "npc"), y1 = unit(0.97, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.05, "npc"), y0 = unit(0.67, "npc"), x1 = unit(0.355, "npc"), y1 = unit(0.67, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.355, "npc"), y0 = unit(0.67, "npc"), x1 = unit(0.355, "npc"), y1 = unit(0.97, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.05, "npc"), y0 = unit(0.97, "npc"), x1 = unit(0.355, "npc"), y1 = unit(0.97, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# 
# grid.segments(x0 = unit(0.365, "npc"), y0 = unit(0.67, "npc"), x1 = unit(0.365, "npc"), y1 = unit(0.97, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.365, "npc"), y0 = unit(0.67, "npc"), x1 = unit(0.665, "npc"), y1 = unit(0.67, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.665, "npc"), y0 = unit(0.67, "npc"), x1 = unit(0.665, "npc"), y1 = unit(0.97, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.365, "npc"), y0 = unit(0.97, "npc"), x1 = unit(0.665, "npc"), y1 = unit(0.97, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# 
# grid.segments(x0 = unit(0.98, "npc"), y0 = unit(0.67, "npc"), x1 = unit(0.98, "npc"), y1 = unit(0.97, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.98, "npc"), y0 = unit(0.67, "npc"), x1 = unit(0.674, "npc"), y1 = unit(0.67, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.674, "npc"), y0 = unit(0.67, "npc"), x1 = unit(0.674, "npc"), y1 = unit(0.97, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# grid.segments(x0 = unit(0.98, "npc"), y0 = unit(0.97, "npc"), x1 = unit(0.674, "npc"), y1 = unit(0.97, "npc"), default.units = "npc", gp = gpar(col = 'grey'))
# 
# grid.rect(x = unit(0.06, "npc"), y = unit(0.97, "npc"), width = unit(0.04, "npc"), height = unit(0.04, "npc"), gp = gpar(fill = 'white', col = 'white'))
# grid.rect(x = unit(0.38, "npc"), y = unit(0.97, "npc"), width = unit(0.04, "npc"), height = unit(0.04, "npc"), gp = gpar(fill = 'white', col = 'white'))
# grid.rect(x = unit(0.69, "npc"), y = unit(0.97, "npc"), width = unit(0.04, "npc"), height = unit(0.04, "npc"), gp = gpar(fill = 'white', col = 'white'))
# 
# grid.rect(x = unit(0.05, "npc"), y = unit(0.33, "npc"), width = unit(0.04, "npc"), height = unit(0.02, "npc"), gp = gpar(fill = 'white', col = 'white'))
# grid.rect(x = unit(0.04, "npc"), y = unit(0.33, "npc"), width = unit(0.02, "npc"), height = unit(0.035, "npc"), gp = gpar(fill = 'white', col = 'white'))
# grid.rect(x = unit(0.04, "npc"), y = unit(0.653, "npc"), width = unit(0.02, "npc"), height = unit(0.035, "npc"), gp = gpar(fill  = 'white', col = 'white'))
# grid.rect(x = unit(0.05, "npc"), y = unit(0.653, "npc"), width = unit(0.04, "npc"), height = unit(0.02, "npc"), gp = gpar(fill  = 'white', col = 'white'))
# 
# grid.text(c("E", 'D', 'C', 'B', 'A'), x = c(0.05, 0.05, 0.684,0.375, 0.06), y = c(0.33, 0.653, 0.97, 0.97, 0.97), 
#           gp=gpar(fontface = 'bold', fontsize = 20))
# 
# dev.off() 

