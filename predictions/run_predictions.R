library(rstan)
library(ggplot2)
library(data.table)
library(gridExtra)
library(dplyr)
library(loo)
library(grid)

indir = "~/git/CDC-covid19-agespecific-mortality-data" # path to the repo
outdir = file.path(indir, 'predictions', 'results')
model = 'BS-GP-I'

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load data
training = read.csv(file.path(indir, 'predictions', 'data', 'training.csv'))
test = read.csv(file.path(indir, 'predictions', 'data', 'prediction.csv'))

# load functions
source(file.path(indir, 'inst', "functions", "stan_utility_functions.R"))

# compile stan models
if(model == 'BS-GP-SE')
  model_stan = rstan::stan_model( file.path(indir, 'predictions', 'stan-models', 'GP-BS-SE_2D.stan') )
if(model == 'BS-GP-I')
  model_stan = rstan::stan_model( file.path(indir, 'predictions', 'stan-models', 'GP-BS-I_2D.stan') )


ps <- c(0.5, 0.025, 0.975)
p_labs = c('M', 'CL', 'CU', 'mean')

# tune 
n_knots_x = 100
n_knots_y = 100
spline_degree = 3

# run GP-BS-SE

# find knots
x = unique(sort(c(training$x, test$x)))
x_grid = data.table(x = seq(min(x) - min(diff(x)), max(x) + min(diff(x)), min(diff(x))))
x_grid[, x_index := 1:nrow(x_grid)]
n = nrow(x_grid)
knots_rows = x_grid$x_index[seq(1, nrow(x_grid), length.out = n_knots_x)]
num_basis_rows = length(knots_rows) + spline_degree - 1
IDX_BASIS_ROWS = 1:num_basis_rows
BASIS_ROWS = bsplines(x_grid$x_index, knots_rows, spline_degree)

y = unique(sort(c(training$y, test$y)))
y_grid = data.table(y = seq(min(y) - min(diff(y)), max(y) + min(diff(y)), min(diff(y))))
y_grid[, y_index := 1:nrow(y_grid)]
m = nrow(y_grid)
knots_columns = y_grid$y_index[seq(1, nrow(y_grid), length.out = n_knots_y)]
num_basis_columns = length(knots_columns) + spline_degree - 1
IDX_BASIS_COLUMNS = 1:num_basis_columns
BASIS_COLUMNS = bsplines(y_grid$y_index, knots_columns, spline_degree)

# find coordinates with observation
grid = as.data.table( expand.grid(x_index = x_grid$x_index, 
                                  y_index = y_grid$y_index) )
grid = merge(grid, x_grid, by = 'x_index')
grid = merge(grid, y_grid, by = 'y_index')

# training coordinates
training = as.data.table(training)
training[, x_index := which(x < x_grid$x[-1] &  x >= x_grid$x[-nrow(x_grid)]), by = 'x']
training[, y_index := which(y < y_grid$y[-1] &  y >= y_grid$y[-nrow(y_grid)]), by = 'y']
value.coordinates = select(training, x_index, y_index)

# test coordinates
test = as.data.table(test)
test[, x_index := which(x < x_grid$x[-1] &  x >= x_grid$x[-nrow(x_grid)]), by = 'x']
test[, y_index := which(y < y_grid$y[-1] &  y >= y_grid$y[-nrow(y_grid)]), by = 'y']
value.coordinates_test = select(test, x_index, y_index)


stan_data = list(n = n, m = m, N = nrow(training),
                 y = training$obs, 
                 coordinates = value.coordinates,
                 num_basis_rows = num_basis_rows, num_basis_columns = num_basis_columns,
                 BASIS_ROWS = BASIS_ROWS, BASIS_COLUMNS = BASIS_COLUMNS,
                 IDX_BASIS_ROWS = IDX_BASIS_ROWS, IDX_BASIS_COLUMNS = IDX_BASIS_COLUMNS)
fit <- rstan::sampling(model_stan,data=stan_data,iter=1000,warmup=200,chains=3, 
                       control = list(max_treedepth = 15, adapt_delta = 0.99))
saveRDS(fit, file.path(outdir, paste0('2D_', model, '_nknots_', n_knots_x, '.rds')))


samples = extract(fit)

tmp1 = as.data.table( reshape2::melt(samples$y_hat))
setnames(tmp1, 2:3, c('x_index', 'y_index'))
tmp1 = tmp1[, list(q = c(quantile(value, probs = ps, na.rm = T), mean(value)),
                   q_label = p_labs),
            by = c('x_index', 'y_index')]
tmp1 = dcast.data.table(tmp1, x_index + y_index ~ q_label, value.var = 'q')
tmp2 = merge(tmp1, test, by = c('x_index', 'y_index'))
tmp2[, MSE := (mean- obs)^2]

cat('MSE is ', mean(tmp2$MSE), '\n')
cat('Time of execution is ', max(apply(get_elapsed_time(fit), 1, sum))/60/60, 'hours \n')


tmp1 = merge(tmp1, x_grid, by = 'x_index')
tmp1 = merge(tmp1, y_grid, by = 'y_index')

tmpp = merge(tmp1, select(test, x_index, y_index), by = c('x_index', 'y_index'))


p0 = ggplot(tmpp, aes(x=x,y=y)) +
  geom_raster(aes(fill=M)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(low = 'green3', high = 'lightpink', mid = 'khaki1', 
                       limits = range(training$obs)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(fill = '') +
  ggtitle('Estimated surface given training data')

p2 = ggplot(test,aes(x=x,y=y)) +
  geom_raster(aes(fill=obs)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(low = 'green3', high = 'lightpink', mid = 'khaki1') + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(fill = '')+
  ggtitle('Complete data')

p3 = ggplot(training,aes(x=x,y=y)) +
  geom_point(aes(col=obs), size = 2) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_color_gradient2(low = 'green3', high = 'lightpink', mid = 'khaki1') + 
  theme_bw() + 
  labs(col = '')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('Training data') 

p = ggpubr::ggarrange(p2,p3,p0,nrow =1,common.legend = T,  legend = 'bottom')
ggsave(p, file= file.path(outdir, paste0('2D_', model, '_nknots_', n_knots_x, '.png')), w = 11, h = 7)
