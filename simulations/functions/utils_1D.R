generate_GP = function(x, n, alpha, l, sigma){
  d = abs(outer(x,x,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
  
  Sigma_SE = alpha^2 * exp(-d^2/(2*l^2)) # squared exponential kernel
  y = mvtnorm::rmvnorm(1,sigma=Sigma_SE) + rnorm(n, 0,sigma)
  
  return(as.vector(y))
}

generate_TPSB = function(x, n, j, sigma){
  
  y = TPSB(x,j) + rnorm(n, 0,sigma)
  
  return(y)
}

TPSB = function(x,j){
  sqrt(x*(1-x)) * sin( (2*pi* (1+ 2^((9-4*j)/5) )) / (x + 2^((9-4*j)/5) ) )
}

run_BSGP = function(x, y, lab, n_knots, spline_degree, outdir)
{

  if(length(x) != length(y)) stop('x and y should have the same length')
  
  knots = x[seq(1, length(x), length.out = n_knots)]
  num_basis = length(knots) + spline_degree - 1
  IDX_BASIS = 1:num_basis
  BASIS = bsplines(x, knots, spline_degree)
  
  stan_data = list(N = length(x), x = x, y = y, 
                   num_basis = num_basis, BASIS = BASIS, IDX_BASIS = IDX_BASIS)
  fit_BSGP <- rstan::sampling(model_BSGP,data=stan_data,iter=1000,warmup=200,chains=3, 
                                  control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit_BSGP, file.path(outdir, paste0('BS-GP_', lab, '_nknots_', n_knots, '.rds')))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  samples_BSGP = extract(fit_BSGP)
  
  tmp = as.data.table( reshape2::melt(samples_BSGP$f) )
  setnames(tmp, 2, c('observations_idx'))
  tmp = tmp[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('observations_idx')]	
  tmp = dcast(tmp, observations_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  tmp1 = as.data.table( reshape2::melt(samples_BSGP$y_hat) )
  setnames(tmp1, 2, c('observations_idx'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('observations_idx')]	
  tmp1 = dcast(tmp1, observations_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  tmp = merge(tmp, data.table(x = x, y = y, observations_idx = 1:length(x)))
  
  tmp[, method := paste0('BS-GP\n#knots = ', n_knots)]
  
  if(grepl('lengthscale', lab))
    tmp[, lengthscape := as.numeric(gsub('lengthscale_(.+)', '\\1', lab))]
  
  time = sum(get_elapsed_time(fit_BSGP))
  tmp[, time := time]

  return(tmp)
}

run_GP = function(x, y, lab, outdir){
  
  if(length(x) != length(y)) stop('x and y should have the same length')
  
  stan_data = list(N = length(x), x = x, y = as.vector(y))
  fit_GP <- rstan::sampling(model_GP,data=stan_data,iter=1000,warmup=200,chains=3, 
                                control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit_GP, file.path(outdir, paste0('GP_', lab, '.rds')))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  samples_GP = extract(fit_GP)
  
  tmp = as.data.table( reshape2::melt(samples_GP$f) )
  setnames(tmp, 2, c('observations_idx'))
  tmp = tmp[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('observations_idx')]	
  tmp = dcast(tmp, observations_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  tmp1 = as.data.table( reshape2::melt(samples_GP$y_hat) )
  setnames(tmp1, 2, c('observations_idx'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('observations_idx')]	
  tmp1 = dcast(tmp1, observations_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  tmp = merge(tmp, data.table(x = x, y = as.vector(y), observations_idx = 1:length(x)))

  if(grepl('lengthscale', lab))
    tmp[, lengthscape := as.numeric(gsub('lengthscale_(.+)', '\\1', lab))]
  
  tmp[, method := paste0('GP')]
  
  time = sum(get_elapsed_time(fit_GP))
  tmp[, time := time]

  return(tmp)
}

run_BSCAR = function(x, y, lab, n_knots, spline_degree, outdir)
{
  
  if(length(x) != length(y)) stop('x and y should have the same length')
  
  knots = x[seq(1, length(x), length.out = n_knots)]
  num_basis = length(knots) + spline_degree - 1
  IDX_BASIS = 1:num_basis
  BASIS = bsplines(x, knots, spline_degree)
  
  Adj = matrix(ncol = num_basis, nrow = num_basis, 0)
  Adj[row(Adj) == (col(Adj) - 1)] = 1
  Adj[row(Adj) == (col(Adj) + 1)] = 1
  Adj_n = sum(Adj)/2
  
  stan_data = list(N = length(x), x = x, y = y, 
                   num_basis = num_basis, BASIS = BASIS, IDX_BASIS = IDX_BASIS, 
                   Adj = Adj, Adj_n = Adj_n)
  
  model_BSCAR = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-CAR_1D.stan') )
  fit_BSCAR <- rstan::sampling(model_BSCAR,data=stan_data,iter=2000,warmup=200,chains=3, 
                               control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit_BSCAR, file.path(outdir, paste0('BS-CAR_', lab, '_nknots_', n_knots, '.rds')))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  samples_BSCAR = extract(fit_BSCAR)
  
  tmp = as.data.table( reshape2::melt(samples_BSCAR$f) )
  setnames(tmp, 2, c('observations_idx'))
  tmp = tmp[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('observations_idx')]	
  tmp = dcast(tmp, observations_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  tmp1 = as.data.table( reshape2::melt(samples_BSCAR$y_hat) )
  setnames(tmp1, 2, c('observations_idx'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('observations_idx')]	
  tmp1 = dcast(tmp1, observations_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  tmp = merge(tmp, data.table(x = x, y = as.vector(y), observations_idx = 1:length(x)))
  
  if(grepl('lengthscale', lab))
    tmp[, lengthscape := as.numeric(gsub('lengthscale_(.+)', '\\1', lab))]
  
  tmp[, method := paste0('BS-CAR\n#knots = ', n_knots)]
  
  time = sum(get_elapsed_time(fit_BSCAR))
  tmp[, time := time]
  
  return(tmp)
}

run_BSIN = function(x, y, lab, n_knots, spline_degree, outdir)
{
  
  if(length(x) != length(y)) stop('x and y should have the same length')
  
  knots = x[seq(1, length(x), length.out = n_knots)]
  num_basis = length(knots) + spline_degree - 1
  IDX_BASIS = 1:num_basis
  BASIS = bsplines(x, knots, spline_degree)

  stan_data = list(N = length(x), x = x, y = y, 
                   num_basis = num_basis, BASIS = BASIS)
  
  model_BSIN = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'BS-IN_1D.stan') )
  fit_BSIN <- rstan::sampling(model_BSIN,data=stan_data,iter=2000,warmup=200,chains=3, 
                               control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit_BSIN, file.path(outdir, paste0('BS-IN_', lab, '_nknots_', n_knots, '.rds')))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  samples_BSIN = extract(fit_BSIN)
  
  tmp = as.data.table( reshape2::melt(samples_BSIN$f) )
  setnames(tmp, 2, c('observations_idx'))
  tmp = tmp[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('observations_idx')]	
  tmp = dcast(tmp, observations_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  tmp1 = as.data.table( reshape2::melt(samples_BSIN$y_hat) )
  setnames(tmp1, 2, c('observations_idx'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('observations_idx')]	
  tmp1 = dcast(tmp1, observations_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  tmp = merge(tmp, data.table(x = x, y = y, observations_idx = 1:length(x)))
  
  if(grepl('lengthscale', lab))
    tmp[, lengthscape := as.numeric(gsub('lengthscale_(.+)', '\\1', lab))]
  
  tmp[, method := paste0('BS-IN\n#knots = ', n_knots)]
  
  time = sum(get_elapsed_time(fit_BSIN))
  tmp[, time := time]
  
  return(tmp)
}
