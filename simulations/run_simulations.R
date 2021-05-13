library(mvtnorm)
library(rstan)
library(ggplot2)
library(data.table)
library(gridExtra)

indir ="~/git/CDC-covid19-agespecific-mortality-data" # path to the repo
outdir = file.path(indir, 'simulations', 'results')
model_BSGP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-GP_1D.stan') )
model_GP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'GP_1D.stan') )

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source(file.path(indir, 'simulations', "functions", "utils.R"))
source(file.path(indir, 'inst', "functions", "stan_utility_functions.R"))

# gp parameters and noise
n = 100
set.seed(11)
x = sort(runif(n))
alpha = 1
sigma = 0.1
lengthscapes = c(0.01, 0.1, 0.25)

# spline parameters for BS-GP
n_knots_vec = c(5, 10, 20)
spline_degree = 3

# length scale = 0.01
l_1 = lengthscapes[1]
GP_1 = run_GP(x, n, alpha, l_1, sigma, outdir)
BSGP_1_1 = run_BSGP(x, n, alpha, l_1, sigma, n_knots_vec[1], spline_degree, outdir)
BSGP_1_2 = run_BSGP(x, n, alpha, l_1, sigma, n_knots_vec[2], spline_degree, outdir)
BSGP_1_3 = run_BSGP(x, n, alpha, l_1, sigma, n_knots_vec[3], spline_degree, outdir)

# length scale = 0.1
l_2 = lengthscapes[2]
GP_2 = run_GP(x, n, alpha, l_2, sigma, outdir)
BSGP_2_1 = run_BSGP(x, n, alpha, l_2, sigma, n_knots_vec[1], spline_degree, outdir)
BSGP_2_2 = run_BSGP(x, n, alpha, l_2, sigma, n_knots_vec[2], spline_degree, outdir)
BSGP_2_3 = run_BSGP(x, n, alpha, l_2, sigma, n_knots_vec[3], spline_degree, outdir)

# length scale = 0.25
l_3 = lengthscapes[3]
GP_3 = run_GP(x, n, alpha, l_3, sigma, outdir)
BSGP_3_1 = run_BSGP(x, n, alpha, l_3, sigma, n_knots_vec[1], spline_degree, outdir)
BSGP_3_2 = run_BSGP(x, n, alpha, l_3, sigma, n_knots_vec[2], spline_degree, outdir)
BSGP_3_3 = run_BSGP(x, n, alpha, l_3, sigma, n_knots_vec[3], spline_degree, outdir)

# comparative figure

tmp = do.call(rbind, list(GP_1, GP_2, GP_3, 
                          BSGP_1_1, BSGP_2_1, BSGP_3_1,
                          BSGP_1_2, BSGP_2_2, BSGP_3_2, 
                          BSGP_1_3, BSGP_2_3, BSGP_3_3))

tmp[, lengthscape_name := paste0('l = ', lengthscape)]
tmp[, method := factor(method, levels = c('GP',
                                          'BS-GP\n#knots = 5', 
                                          'BS-GP\n#knots = 10',
                                          'BS-GP\n#knots = 20'))]

tmp1 = subset(tmp, variable == 'f')
ggplot(tmp1, aes(x = x)) + 
  geom_point(aes(y = y), col = 'red') + 
  geom_line(aes(y = M)) + 
  geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
  theme_bw() + 
  labs(x = "x", y = 'y') + 
  facet_grid(method~lengthscape_name) + 
  geom_text(
    data    = tmp,
    mapping = aes(x = -Inf, y = -Inf, label = paste0('Execution time\n ', round(time/60, digits = 2), ' minutes')),
    hjust   = -1,
    vjust   = -0.1) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))
ggsave(file.path(outdir, 'simulations_comp_f.png'), w = 8, h = 8)

tmp1 = subset(tmp, variable == 'y_hat')
ggplot(tmp1, aes(x = x)) + 
  geom_line(aes(y = M)) + 
  geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
  geom_point(aes(y = y), col = 'red') + 
  theme_bw() + 
  labs(x = "x", y = 'y') + 
  facet_grid(method~lengthscape_name) + 
  geom_text(
    data    = tmp,
    mapping = aes(x = -Inf, y = -Inf, label = paste0('Execution time\n ', round(time/60, digits = 2), ' minutes')),
    hjust   = -1,
    vjust   = -0.1) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))
ggsave(file.path(outdir, 'simulations_comp_yhat.png'), w = 8, h = 8)

