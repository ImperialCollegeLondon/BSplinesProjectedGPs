library(mvtnorm)
library(rstan)
library(ggplot2)
library(data.table)
library(gridExtra)
library("plgp")
library(dplyr)
library(loo)
library(grid)

indir ="~/git/CDC-covid19-agespecific-mortality-data" # path to the repo
outdir = file.path(indir, 'simulations', 'results')

model_GP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'GP-SE_2D.stan') )
model_GPI = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'GP-I_2D.stan') )
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
GP_I_2D_1 = run_GP_I_2D(x, x, y, lab, outdir)
BSGP_2D_1_1 = run_BSGP_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSGP_2D_2_1 = run_BSGP_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSGP_2D_3_1 = run_BSGP_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSIN_2D_1_1 = run_BSIN_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSIN_2D_2_1 = run_BSIN_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSIN_2D_3_1 = run_BSIN_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)

##
n_knots_vec = c(5, 10, 20)
lab = paste0('lengthscale_', lengthscales[2])
set.seed(24)
y = generate_2DGP(X, lengthscales[2], sigma = 0.5)
GP_2D_2 = run_GP_2D(x, x, y, lab, outdir)
GP_I_2D_2 = run_GP_I_2D(x, x, y, lab, outdir)
BSGP_2D_1_2 = run_BSGP_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSGP_2D_2_2 = run_BSGP_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSGP_2D_3_2 = run_BSGP_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSIN_2D_1_2 = run_BSIN_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSIN_2D_2_2 = run_BSIN_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSIN_2D_3_2 = run_BSIN_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)

##
n_knots_vec = c(5, 10, 20)
lab = paste0('lengthscale_', lengthscales[3])
set.seed(23)
y = generate_2DGP(X, lengthscales[3], sigma= 0.1)
GP_2D_3 = run_GP_2D(x, x, y, lab, outdir)
GP_I_2D_3 = run_GP_I_2D(x, x, y, lab, outdir)
BSGP_2D_1_3 = run_BSGP_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSGP_2D_2_3 = run_BSGP_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSGP_2D_3_3 = run_BSGP_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSIN_2D_1_3 = run_BSIN_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSIN_2D_2_3 = run_BSIN_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSIN_2D_3_3 = run_BSIN_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)


# plot results
tmp = do.call('rbind', list(GP_2D_1[[1]], GP_I_2D_1[[1]],
                            BSGP_2D_1_1[[1]], BSGP_2D_2_1[[1]], BSGP_2D_3_1[[1]], 
                            BSIN_2D_1_1[[1]], BSIN_2D_2_1[[1]], BSIN_2D_3_1[[1]],  
                            GP_2D_2[[1]], GP_I_2D_2[[1]],
                            BSGP_2D_1_2[[1]], BSGP_2D_2_2[[1]], BSGP_2D_3_2[[1]], 
                            BSIN_2D_1_2[[1]], BSIN_2D_2_2[[1]], BSIN_2D_3_2[[1]],
                            GP_2D_3[[1]], GP_I_2D_3[[1]],
                            BSGP_2D_1_3[[1]], BSGP_2D_2_3[[1]], BSGP_2D_3_3[[1]], 
                            BSIN_2D_1_3[[1]], BSIN_2D_2_3[[1]], BSIN_2D_3_3[[1]]))

unique(tmp$method)

tmp[grepl('GP-I', method), method := gsub('GP-I', 'GP-GN', method)]
tmp[grepl('GP-GNN', method), method := gsub('GP-GNN', 'GP-GN', method)]
tmp[grepl('BS-GP', method), method := gsub('BS-GP', 'GP-BS', method)]

tmp[grepl('knots', method), n_knots := as.numeric(gsub('.*\n#knots = (.+)', '\\1', method))]
tmp[grepl('knots', method), method2 := gsub('(.+)\n#knots = .*', '\\1', method)]
tmp[!grepl('knots', method), method2 := method]
tmp[!grepl('knots', method), n_knots := 0]


for(l in lengthscales){
  # l = lengthscales[1]
  tmp1 = subset(tmp, variable == 'y_hat' & lengthscale == l)
  
  tmp2 = copy(select(tmp1, x_1, x_2, M, method2, n_knots, time, lengthscale))
  tmp1 = copy(select(tmp1, x_1, x_2, y, time, n_knots, lengthscale))
  tmp1[, time := NA_real_]
  tmp1[, method2 := 'observation']
  setnames(tmp1, 'y', 'M')
  tmp1 = rbind(unique(tmp1), tmp2)
  
  tmp1[, method2 := factor(method2, levels = c('observation', 'GP-GN', 'GP-SE', 'GP-BS-GN', 'GP-BS-SE'))]
  
  tmp2 = subset(tmp1, method2 %in% c('observation',  'GP-GN', 'GP-SE' ))
  p1= ggplot(tmp2,aes(x=x_1,y=x_2)) +
    geom_raster(aes(fill=M), interpolate = TRUE) +
    geom_contour(aes(z=M), bins = 12, color = "gray30", 
                 size = 0.5, alpha = 0.5) +
    coord_equal() +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_viridis_c(option = "viridis") +
    labs(x = expression(x[2]), y = expression(x[1])) + 
    facet_grid(.~method2) +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.text.x = element_text(angle=90)) 
  
  tmp2 = subset(tmp1, !method2 %in% c('observation', 'GP-GN','GP-SE'))
  p2= ggplot(tmp2,aes(x=x_1,y=x_2)) +
    geom_raster(aes(fill=M), interpolate = TRUE) +
    geom_contour(aes(z=M), bins = 12, color = "gray30", 
                 size = 0.5, alpha = 0.5) +
    coord_equal() +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_viridis_c(option = "viridis") +
    labs(x = expression(x[2]), y = expression(x[1])) + 
    facet_grid(n_knots~method2) +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.text.x = element_text(angle=90)) 
  
  # p2 = grid.arrange(p2, right = textGrob('Number of knots', vjust = 1.5, rot = 270))
  p = grid.arrange(grobs = list(p1, p2), heights = c(0.5, 1.1), 
                   layout_matrix= rbind(c(1,NA), c(2,2)), widths = c(1, 0.1), 
                   right = textGrob('Number of knots', vjust = 1.5, rot = 270, hjust = -0.3))
  
  ggsave(p, file = file.path(outdir, paste0('2D_comp_lengthscale_', l, '.png')), w = 5, h = 7)
}

# time of execution 
tmp1 = subset(tmp, grepl('GP-SE', method))
time_GP = paste0(round(unique(tmp1$time) / 60), ' minutes')

tmp1 = subset(tmp, grepl('GP-BS-SE', method))
time_GPSE = list()
time_GPSE[[1]] = paste0(round(max(subset(tmp1, lengthscale == lengthscales[1])$time)/ 60), ' minutes')
time_GPSE[[2]] = paste0(round(max(subset(tmp1, lengthscale == lengthscales[2])$time)/ 60), ' minutes')
time_GPSE[[3]] = paste0(round(max(subset(tmp1, lengthscale == lengthscales[3])$time)/ 60), ' minutes')

saveRDS(list(time_GP, time_GPSE), file = file.path(outdir, paste0('time_execution.rds')))

loo_compare(loo(extract(GP_2D_1[[2]])$log_lik), 
            loo(extract(BSGP_2D_1_1[[2]])$log_lik), loo(extract(BSGP_2D_2_1[[2]])$log_lik), loo(extract(BSGP_2D_3_1[[2]])$log_lik))

# compare loo
# first scenario
loo_compare(loo(extract(GP_2D_1[[2]])$log_lik), loo(extract(GP_I_2D_1[[2]])$log_lik),
            loo(extract(BSGP_2D_1_1[[2]])$log_lik), loo(extract(BSGP_2D_2_1[[2]])$log_lik), loo(extract(BSGP_2D_3_1[[2]])$log_lik), 
            loo(extract(BSIN_2D_1_1[[2]])$log_lik), loo(extract(BSIN_2D_2_1[[2]])$log_lik), loo(extract(BSIN_2D_3_1[[2]])$log_lik))

# second scenario
loo_compare(loo(extract(GP_2D_2[[2]])$log_lik), loo(extract(GP_I_2D_2[[2]])$log_lik),
            loo(extract(BSGP_2D_1_2[[2]])$log_lik), loo(extract(BSGP_2D_2_2[[2]])$log_lik), loo(extract(BSGP_2D_3_2[[2]])$log_lik), 
            loo(extract(BSIN_2D_1_2[[2]])$log_lik), loo(extract(BSIN_2D_2_2[[2]])$log_lik), loo(extract(BSIN_2D_3_2[[2]])$log_lik))

# third scenario
loo_compare(loo(extract(GP_2D_3[[2]])$log_lik), loo(extract(GP_I_2D_3[[2]])$log_lik),
            loo(extract(BSGP_2D_1_3[[2]])$log_lik), loo(extract(BSGP_2D_2_3[[2]])$log_lik), loo(extract(BSGP_2D_3_3[[2]])$log_lik), 
            loo(extract(BSIN_2D_1_3[[2]])$log_lik), loo(extract(BSIN_2D_2_3[[2]])$log_lik), loo(extract(BSIN_2D_3_3[[2]])$log_lik))

# for text 
LOO1 = list()
LOO1[[1]] = round(loo_compare(loo(extract(BSGP_2D_2_1[[2]])$log_lik), loo(extract(BSIN_2D_2_1[[2]])$log_lik)), digits = 1)
LOO1[[2]] = round(loo_compare(loo(extract(BSGP_2D_3_1[[2]])$log_lik), loo(extract(BSIN_2D_3_1[[2]])$log_lik)), digits = 1)

LOO1[[3]] = round(loo_compare(loo(extract(BSGP_2D_2_1[[2]])$log_lik), loo(extract(GP_2D_1[[2]])$log_lik)), digits = 1)
LOO1[[4]] = round(loo_compare(loo(extract(BSGP_2D_3_1[[2]])$log_lik), loo(extract(GP_2D_1[[2]])$log_lik)), digits = 1)

LOO1[[5]] = round(loo_compare(loo(extract(BSGP_2D_1_3[[2]])$log_lik), loo(extract(GP_2D_3[[2]])$log_lik)), digits = 1)
saveRDS(LOO1, file = file.path(outdir, paste0('LOO_1.rds')))


LOO_GP = list(gsub(" ", "", format(round(loo(extract(GP_2D_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
              gsub(" ", "", format(round(loo(extract(GP_2D_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
              gsub(" ", "", format(round(loo(extract(GP_2D_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )

LOO_GP_I = list(gsub(" ", "", format(round(loo(extract(GP_I_2D_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
              gsub(" ", "", format(round(loo(extract(GP_I_2D_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
              gsub(" ", "", format(round(loo(extract(GP_I_2D_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )

LOO_BSGP_1 = list(gsub(" ", "", format(round(loo(extract(BSGP_2D_1_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
                  gsub(" ", "", format(round(loo(extract(BSGP_2D_1_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
                  gsub(" ", "", format(round(loo(extract(BSGP_2D_1_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )

LOO_BSGP_2 = list(gsub(" ", "",format(round(loo(extract(BSGP_2D_2_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
                  gsub(" ", "",format(round(loo(extract(BSGP_2D_2_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
                  gsub(" ", "",format(round(loo(extract(BSGP_2D_2_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )

LOO_BSGP_3 = list(gsub(" ", "",format(round(loo(extract(BSGP_2D_3_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
                  gsub(" ", "",format(round(loo(extract(BSGP_2D_3_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) ,
                  gsub(" ", "",format(round(loo(extract(BSGP_2D_3_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)) )

LOO_BSIN_1 = list( gsub(" ", "",format( round(loo(extract(BSIN_2D_1_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)),
                   gsub(" ", "",format( round(loo(extract(BSIN_2D_1_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)),
                   gsub(" ", "",format( round(loo(extract(BSIN_2D_1_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)))

LOO_BSIN_2 = list(gsub(" ", "", format( round(loo(extract(BSIN_2D_2_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2) ),
                  gsub(" ", "", format( round(loo(extract(BSIN_2D_2_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)),
                  gsub(" ", "", format( round(loo(extract(BSIN_2D_2_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)))

LOO_BSIN_3 = list(gsub(" ", "", format( round(loo(extract(BSIN_2D_3_1[[2]])$log_lik)$estimates, digits = 2), nsmall = 2) ),
                  gsub(" ", "", format( round(loo(extract(BSIN_2D_3_2[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)),
                  gsub(" ", "", format( round(loo(extract(BSIN_2D_3_3[[2]])$log_lik)$estimates, digits = 2), nsmall = 2)))

LOO_ = list(LOO_GP_I, LOO_GP, LOO_BSGP_1, LOO_BSGP_2, LOO_BSGP_3, 
     LOO_BSIN_1, LOO_BSIN_2, LOO_BSIN_3)
saveRDS(LOO_, file = file.path(outdir, paste0('LOO_all.rds')))


