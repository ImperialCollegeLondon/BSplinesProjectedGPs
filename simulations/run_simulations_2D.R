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

model_GP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'GP_2D_Kronecker.stan') )
model_BSGP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-GP-SE_2D.stan') )
model_BSCAR = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-GP-CAR_2D.stan') )
model_BSIN = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-GP-I_2D.stan') )

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source(file.path(indir, 'simulations', "functions", "utils_2D.R"))
source(file.path(indir, 'inst', "functions", "stan_utility_functions.R"))


# number of samples
n <- 20
x <- seq(0,10,length=n)
X <- expand.grid(x, x)
sigma = 0.1
lengthscales = c(2, 5, 10)

# spline parameters 
n_knots_vec = c(5, 10, 20)
spline_degree = 3

##
lab = paste0('lengthscale_', lengthscales[1])
set.seed(23)
y = generate_2DGP(X, lengthscales[1], sigma)
GP_2D_1 = run_GP_2D(x, x, y, lab, outdir)
BSGP_2D_1_1 = run_BSGP_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSGP_2D_2_1 = run_BSGP_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSGP_2D_3_1 = run_BSGP_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSCAR_2D_1_1 = run_BSCAR_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSCAR_2D_2_1 = run_BSCAR_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSCAR_2D_3_1 = run_BSCAR_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSIN_2D_1_1 = run_BSIN_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSIN_2D_2_1 = run_BSIN_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSIN_2D_3_1 = run_BSIN_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)

##
lab = paste0('lengthscale_', lengthscales[2])
set.seed(23)
y = generate_2DGP(X, lengthscales[2], sigma)
GP_2D_2 = run_GP_2D(x, x, y, lab, outdir)
BSGP_2D_1_2 = run_BSGP_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSGP_2D_2_2 = run_BSGP_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSGP_2D_3_2 = run_BSGP_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSCAR_2D_1_2 = run_BSCAR_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSCAR_2D_2_2 = run_BSCAR_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSCAR_2D_3_2 = run_BSCAR_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSIN_2D_1_2 = run_BSIN_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSIN_2D_2_2 = run_BSIN_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSIN_2D_3_2 = run_BSIN_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)

##
lab = paste0('lengthscale_', lengthscales[3])
set.seed(23)
y = generate_2DGP(X, lengthscales[3], sigma)
GP_2D_3 = run_GP_2D(x, x, y, lab, outdir)
BSGP_2D_1_3 = run_BSGP_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSGP_2D_2_3 = run_BSGP_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSGP_2D_3_3 = run_BSGP_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSCAR_2D_1_3 = run_BSCAR_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSCAR_2D_2_3 = run_BSCAR_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSCAR_2D_3_3 = run_BSCAR_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)
BSIN_2D_1_3 = run_BSIN_2D(x, x, y, lab, n_knots_vec[1], spline_degree, outdir)
BSIN_2D_2_3 = run_BSIN_2D(x, x, y, lab, n_knots_vec[2], spline_degree, outdir)
BSIN_2D_3_3 = run_BSIN_2D(x, x, y, lab, n_knots_vec[3], spline_degree, outdir)


# plot results
tmp = do.call('rbind', list(GP_2D_1[[1]], BSGP_2D_1_1[[1]], BSGP_2D_2_1[[1]], BSGP_2D_3_1[[1]], 
                            BSCAR_2D_1_1[[1]], BSCAR_2D_2_1[[1]], BSCAR_2D_3_1[[1]], 
                            BSIN_2D_1_1[[1]], BSIN_2D_2_1[[1]], BSIN_2D_3_1[[1]],  
                            GP_2D_2[[1]], BSGP_2D_1_2[[1]], BSGP_2D_2_2[[1]], BSGP_2D_3_2[[1]], 
                            BSCAR_2D_1_2[[1]], BSCAR_2D_2_2[[1]], BSCAR_2D_3_2[[1]], 
                            BSIN_2D_1_2[[1]], BSIN_2D_2_2[[1]], BSIN_2D_3_2[[1]],
                            GP_2D_3[[1]], BSGP_2D_1_3[[1]], BSGP_2D_2_3[[1]], BSGP_2D_3_3[[1]], 
                            BSCAR_2D_1_3[[1]], BSCAR_2D_2_3[[1]], BSCAR_2D_3_3[[1]], 
                            BSIN_2D_1_3[[1]], BSIN_2D_2_3[[1]], BSIN_2D_3_3[[1]]))

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
  
  tmp1[, method2 := factor(method2, levels = c('observation', 'GP-SE', 'BS-GP-IN', 'BS-GP-CAR', 'BS-GP-SE'))]
  
  tmp2 = subset(tmp1, method2 %in% c('observation', 'GP-SE'))
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
          panel.border = element_rect(colour = "black", fill = NA)) 
  
  tmp2 = subset(tmp1, !method2 %in% c('observation', 'GP-SE'))
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
          panel.border = element_rect(colour = "black", fill = NA)) 
  
  # p2 = grid.arrange(p2, right = textGrob('Number of knots', vjust = 1.5, rot = 270))
  p = grid.arrange(grobs = list(p1, p2), heights = c(0.5, 1.1), 
                   layout_matrix= rbind(c(1,NA), c(2,2)), widths = c(1, 0.1), 
                   right = textGrob('Number of knots', vjust = 1.5, rot = 270, hjust = -0.3))
  
  ggsave(p, file = file.path(outdir, paste0('2D_comp_lengthscale_', l, '.png')), w = 5, h = 7)
}

# compare loo
loo_compare(loo(extract(BSGP_2D_1_1[[2]])$log_lik), loo(extract(BSGP_2D_2_1[[2]])$log_lik))


 # BSGP_2D_3_1[[2]], 
 #            BSCAR_2D_1_1[[2]], BSCAR_2D_2_1[[2]], BSCAR_2D_3_1[[2]], 
 #            BSIN_2D_1_1[[2]], BSIN_2D_2_1[[2]], BSIN_2D_3_1[[2]])

