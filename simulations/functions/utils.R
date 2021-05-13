generate_GP = function(x, n, alpha, l, sigma){
  d = abs(outer(x,x,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
  
  Sigma_SE = alpha^2 * exp(-d^2/(2*l^2)) # squared exponential kernel
  y = mvtnorm::rmvnorm(1,sigma=Sigma_SE) + rnorm(n, 0,sigma)
  
  return(as.vector(y))
}

run_BSGP = function(x, n, alpha, l, sigma, n_knots, spline_degree, outdir)
{
  set.seed(11)
  y = generate_GP(x, n, alpha, l, sigma)
  # plot(x, y)
  
  knots = x[seq(1, n, length.out = n_knots)]
  num_basis = length(knots) + spline_degree - 1
  IDX_BASIS = 1:num_basis
  BASIS = bsplines(x, knots, spline_degree)
  
  stan_data = list(N = n, x = x, y = y, 
                   num_basis = num_basis, BASIS = BASIS, IDX_BASIS = IDX_BASIS)
  fit_BSGP <- rstan::sampling(model_BSGP,data=stan_data,iter=1000,warmup=200,chains=3, 
                                  control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit_BSGP, file.path(outdir, paste0('BS-GP_lengthscale_', l, '_nknots_', n_knots, '.rds')))
  
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
  
  tmp = merge(tmp, data.table(x = x, y = as.vector(y), observations_idx = 1:n))
  tmp[, lengthscape := l]
  tmp[, method := paste0('BS-GP\n#knots = ', n_knots)]
  
  time = sum(get_elapsed_time(fit_BSGP))
  tmp[, time := time]

  return(tmp)
}

run_GP = function(x, n, alpha, l, sigma, outdir){
  
  set.seed(11)
  y = generate_GP(x, n, alpha, l = l, sigma)
  # plot(x,y)
  
  stan_data = list(N = n, x = x, y = as.vector(y))
  fit_GP <- rstan::sampling(model_GP,data=stan_data,iter=1000,warmup=200,chains=3, 
                                control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit_GP, file.path(outdir, paste0('GP_lengthscale_', l, '.rds')))
  
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
  
  tmp = merge(tmp, data.table(x = x, y = as.vector(y), observations_idx = 1:n))

  tmp[, lengthscape := l]
  tmp[, method := paste0('GP')]
  
  time = sum(get_elapsed_time(fit_GP))
  tmp[, time := time]

  return(tmp)
}

