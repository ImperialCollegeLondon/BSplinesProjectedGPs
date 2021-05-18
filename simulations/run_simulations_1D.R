library(mvtnorm)
library(rstan)
library(ggplot2)
library(data.table)
library(gridExtra)

indir ="~/git/CDC-covid19-agespecific-mortality-data" # path to the repo
outdir = file.path(indir, 'simulations', 'results')

model_GP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'GP_1D.stan') )
model_BSGP = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-GP-SE_1D.stan') )
model_BSCAR = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-GP-CAR_1D.stan') )
model_BSIN = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-GP-I_1D.stan') )

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source(file.path(indir, 'simulations', "functions", "utils_1D.R"))
source(file.path(indir, 'inst', "functions", "stan_utility_functions.R"))

# gp parameters and noise
n = 100
set.seed(11)
x = sort(runif(n))
alpha = 1
sigma = 0.1
lengthscapes = c(0.01, 0.1, 0.25)

# spline parameters 
n_knots_vec = c(5, 10, 20)
spline_degree = 3


# length scale = 0.01
l_1 = lengthscapes[1]
set.seed(11)
y = generate_GP(x, n, alpha, l_1, sigma)
#plot(x, y)

lab_1 = paste0('lengthscale_', l_1)
GP_1 = run_GP(x, y, lab_1, outdir)
BSGP_1_1 = run_BSGP(x, y, lab_1, n_knots_vec[1], spline_degree, outdir)
BSGP_1_2 = run_BSGP(x, y, lab_1, n_knots_vec[2], spline_degree, outdir)
BSGP_1_3 = run_BSGP(x, y, lab_1, n_knots_vec[3], spline_degree, outdir)
BSCAR_1_1 = run_BSCAR(x, y, lab_1, n_knots_vec[1], spline_degree, outdir)
BSCAR_1_2 = run_BSCAR(x, y, lab_1, n_knots_vec[2], spline_degree, outdir)
BSCAR_1_3 = run_BSCAR(x, y, lab_1, n_knots_vec[3], spline_degree, outdir)
BSIN_1_1 = run_BSIN(x, y, lab_1, n_knots_vec[1], spline_degree, outdir)
BSIN_1_2 = run_BSIN(x, y, lab_1, n_knots_vec[2], spline_degree, outdir)
BSIN_1_3 = run_BSIN(x, y, lab_1, n_knots_vec[3], spline_degree, outdir)


# length scale = 0.1
l_2 = lengthscapes[2]
set.seed(11)
y = generate_GP(x, n, alpha, l_2, sigma)
#plot(x, y)

lab_2 = paste0('lengthscale_', l_2)
GP_2 = run_GP(x, y, lab_2, outdir)
BSGP_2_1 = run_BSGP(x, y, lab_2, n_knots_vec[1], spline_degree, outdir)
BSGP_2_2 = run_BSGP(x, y, lab_2, n_knots_vec[2], spline_degree, outdir)
BSGP_2_3 = run_BSGP(x, y, lab_2, n_knots_vec[3], spline_degree, outdir)
BSCAR_2_1 = run_BSCAR(x, y, lab_2, n_knots_vec[1], spline_degree, outdir)
BSCAR_2_2 = run_BSCAR(x, y, lab_2, n_knots_vec[2], spline_degree, outdir)
BSCAR_2_3 = run_BSCAR(x, y, lab_2, n_knots_vec[3], spline_degree, outdir)
BSIN_2_1 = run_BSIN(x, y, lab_2, n_knots_vec[1], spline_degree, outdir)
BSIN_2_2 = run_BSIN(x, y, lab_2, n_knots_vec[2], spline_degree, outdir)
BSIN_2_3 = run_BSIN(x, y, lab_2, n_knots_vec[3], spline_degree, outdir)


# length scale = 0.25
l_3 = lengthscapes[3]
set.seed(11)
y = generate_GP(x, n, alpha, l_3, sigma)
#plot(x, y)

lab_3 = paste0('lengthscale_', l_3)
GP_3 = run_GP(x, y, lab_3, outdir)
BSGP_3_1 = run_BSGP(x, y, lab_3, n_knots_vec[1], spline_degree, outdir)
BSGP_3_2 = run_BSGP(x, y, lab_3, n_knots_vec[2], spline_degree, outdir)
BSGP_3_3 = run_BSGP(x, y, lab_3, n_knots_vec[3], spline_degree, outdir)
BSCAR_3_1 = run_BSCAR(x, y, lab_3, n_knots_vec[1], spline_degree, outdir)
BSCAR_3_2 = run_BSCAR(x, y, lab_3, n_knots_vec[2], spline_degree, outdir)
BSCAR_3_3 = run_BSCAR(x, y, lab_3, n_knots_vec[3], spline_degree, outdir)
BSIN_3_1 = run_BSIN(x, y, lab_3, n_knots_vec[1], spline_degree, outdir)
BSIN_3_2 = run_BSIN(x, y, lab_3, n_knots_vec[2], spline_degree, outdir)
BSIN_3_3 = run_BSIN(x, y, lab_3, n_knots_vec[3], spline_degree, outdir)



# comparative figure

tmp = do.call(rbind, list(GP_1, GP_2, GP_3, 
                          BSGP_1_1, BSGP_2_1, BSGP_3_1,
                          BSGP_1_2, BSGP_2_2, BSGP_3_2, 
                          BSGP_1_3, BSGP_2_3, BSGP_3_3,
                          BSCAR_1_1, BSCAR_2_1, BSCAR_3_1,
                          BSCAR_1_2, BSCAR_2_2, BSCAR_3_2, 
                          BSCAR_1_3, BSCAR_2_3, BSCAR_3_3,
                          BSIN_1_1, BSIN_2_1, BSIN_3_1,
                          BSIN_1_2, BSIN_2_2, BSIN_3_2, 
                          BSIN_1_3, BSIN_2_3, BSIN_3_3))

tmp[, lengthscape_name := paste0('l = ', lengthscape)]
tmp[, method := factor(method, levels = c('GP',
                                          'BS-GP\n#knots = 5', 
                                          'BS-GP\n#knots = 10',
                                          'BS-GP\n#knots = 20',
                                          'BS-CAR\n#knots = 5', 
                                          'BS-CAR\n#knots = 10',
                                          'BS-CAR\n#knots = 20',
                                          'BS-IN\n#knots = 5', 
                                          'BS-IN\n#knots = 10',
                                          'BS-IN\n#knots = 20'))]

tmp1 = subset(tmp, variable == 'f')
ggplot(tmp1, aes(x = x)) + 
  geom_point(aes(y = y), col = 'red', size = 0.1) + 
  geom_line(aes(y = M)) + 
  geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
  theme_bw() + 
  labs(x = "x", y = 'y') + 
  facet_grid(method~lengthscape_name) + 
  geom_text(
    data    = tmp,
    mapping = aes(x = -Inf, y = -Inf, label = paste0('Execution time\n ', round(time/60, digits = 2), ' minutes')),
    hjust   = -1,
    vjust   = -0.1, 
    size = 2) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))
ggsave(file.path(outdir, 'simulations_comp_f.png'), w = 5, h = 10)

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


selected_methods = c('GP', 'BS-GP\n#knots = 20','BS-CAR\n#knots = 20','BS-IN\n#knots = 20')
selected_lengthscape = c(0.1, 0.25)
tmp2 = subset(tmp1, method %in%  selected_methods &
       lengthscape %in% selected_lengthscape)

tmp2[, lengthscape := factor(lengthscape, levels = selected_lengthscape)]
tmp2[, lengthscape_name := paste0('l = ', lengthscape)]
tmp2[, method := factor(method, levels = selected_methods)]

ggplot(tmp2, aes(x = x)) + 
  geom_point(aes(y = y), col = 'red', size = 0.1) + 
  geom_line(aes(y = M)) + 
  geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
  theme_bw() + 
  labs(x = "x", y = 'y') + 
  facet_grid(method~lengthscape_name) + 
  geom_text(
    data    = tmp2,
    mapping = aes(x = -Inf, y = -Inf, label = paste0('Execution time\n ', round(time/60, digits = 2), ' minutes')),
    hjust   = -1,
    vjust   = -0.1, 
    size = 2) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))
ggsave(file.path(outdir, 'simulations_comp_f_zoom.png'), w = 5, h = 6)
