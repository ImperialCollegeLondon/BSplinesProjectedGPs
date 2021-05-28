
rbf_D <- function(X,l=1, eps = sqrt(.Machine$double.eps) ){
  # squared exponential kernel
  D <- plgp::distance(X)
  Sigma <- exp(-D/l)^2 + diag(eps, nrow(X))
}

generate_2DGP = function(X, l, sigma){
  Sigma <- rbf_D(X,l=l)
  y <- MASS::mvrnorm(1,rep(0,dim(Sigma)[1]), Sigma) + rnorm(dim(Sigma)[1], 0, sigma)
  return(y)
}

run_BSGP_2D = function(x_1, x_2, y, lab, n_knots, spline_degree, outdir){
  
  knots_rows = x_1[seq(1, length(x_1), length.out = n_knots)]
  num_basis_rows = length(knots_rows) + spline_degree - 1
  IDX_BASIS_ROWS = 1:num_basis_rows
  BASIS_ROWS = bsplines(x_1, knots_rows, spline_degree)
  
  knots_columns = x_2[seq(1, length(x_2), length.out = n_knots)]
  num_basis_columns = length(knots_columns) + spline_degree - 1
  IDX_BASIS_COLUMNS = 1:num_basis_columns
  BASIS_COLUMNS = bsplines(x_2, knots_columns, spline_degree)
  
  stan_data = list(n = length(x_1), m = length(x_2), x_1 = x_1, x_2 = x_2, 
                   y = matrix(y, nrow = length(x_1), ncol = length(x_2), byrow = F), 
                   num_basis_rows = num_basis_rows, num_basis_columns = num_basis_columns,
                   BASIS_ROWS = BASIS_ROWS, BASIS_COLUMNS = BASIS_COLUMNS,
                   IDX_BASIS_ROWS = IDX_BASIS_ROWS, IDX_BASIS_COLUMNS = IDX_BASIS_COLUMNS)
  
  fit <- rstan::sampling(model_BSGP,data=stan_data,iter=1000,warmup=200,chains=3, 
                              control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit, file.path(outdir, paste0('2D_BS-GP_', lab, '_nknots_', n_knots, '.rds')))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  samples = extract(fit)
  
  tmp = as.data.table( reshape2::melt(samples$f) )
  setnames(tmp, 2:3, c('rows_idx', 'column_idx'))
  tmp = tmp[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('rows_idx', 'column_idx')]	
  tmp = dcast(tmp, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  tmp1 = as.data.table( reshape2::melt(samples$y_hat) )
  setnames(tmp1, 2:3, c('rows_idx', 'column_idx'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('rows_idx', 'column_idx')]	
  tmp1 = dcast(tmp1, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  tmp1 = data.table(y=y, expand.grid(x_1=x_1,x_2 = x_2), expand.grid(rows_idx = 1:length(x_1), column_idx = 1:length(x_2)))
  tmp = merge(tmp, tmp1, by = c('rows_idx', 'column_idx'))
  
  tmp[, method := paste0('BS-GP-SE\n#knots = ', n_knots)]
  
  if(grepl('lengthscale', lab))
    tmp[, lengthscale := as.numeric(gsub('lengthscale_(.+)', '\\1', lab))]
  
  time = sum(get_elapsed_time(fit))
  tmp[, time := time]
  
  return(list(tmp, fit))
}

run_BSCAR_2D = function(x_1, x_2, y, lab, n_knots, spline_degree, outdir){
  
  knots_rows = x_1[seq(1, length(x_1), length.out = n_knots)]
  num_basis_rows = length(knots_rows) + spline_degree - 1
  IDX_BASIS_ROWS = 1:num_basis_rows
  BASIS_ROWS = bsplines(x_1, knots_rows, spline_degree)
  
  knots_columns = x_2[seq(1, length(x_2), length.out = n_knots)]
  num_basis_columns = length(knots_columns) + spline_degree - 1
  IDX_BASIS_COLUMNS = 1:num_basis_columns
  BASIS_COLUMNS = bsplines(x_2, knots_columns, spline_degree)
  
  N = num_basis_rows * num_basis_columns
  A = find_adjacency_matrix(num_basis_rows, num_basis_columns)
  Adj_n = sum(A) / 2
  
  stan_data = list(n = length(x_1), m = length(x_2), x_1 = x_1, x_2 = x_2, 
                   y = matrix(y, nrow = length(x_1), ncol = length(x_2), byrow = F), 
                   num_basis_rows = num_basis_rows, num_basis_columns = num_basis_columns,
                   BASIS_ROWS = BASIS_ROWS, BASIS_COLUMNS = BASIS_COLUMNS, 
                   N = N, Adj = A, Adj_n = Adj_n)
  
  fit <- rstan::sampling(model_BSCAR,data=stan_data,iter=1000,warmup=200,chains=3, 
                              control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit, file.path(outdir, paste0('2D_BS-CAR_', lab, '_nknots_', n_knots, '.rds')))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  samples = extract(fit)
  
  tmp = as.data.table( reshape2::melt(samples$f) )
  setnames(tmp, 2:3, c('rows_idx', 'column_idx'))
  tmp = tmp[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('rows_idx', 'column_idx')]	
  tmp = dcast(tmp, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  tmp1 = as.data.table( reshape2::melt(samples$y_hat) )
  setnames(tmp1, 2:3, c('rows_idx', 'column_idx'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('rows_idx', 'column_idx')]	
  tmp1 = dcast(tmp1, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  tmp1 = data.table(y=y, expand.grid(x_1=x_1,x_2 = x_2), expand.grid(rows_idx = 1:length(x_1), column_idx = 1:length(x_2)))
  tmp = merge(tmp, tmp1, by = c('rows_idx', 'column_idx'))
  
  tmp[, method := paste0('BS-GP-CAR\n#knots = ', n_knots)]
  
  if(grepl('lengthscale', lab))
    tmp[, lengthscale := as.numeric(gsub('lengthscale_(.+)', '\\1', lab))]
  
  time = sum(get_elapsed_time(fit))
  tmp[, time := time]
  
  return(list(tmp, fit))
}

run_BSIN_2D = function(x_1, x_2, y, lab, n_knots, spline_degree, outdir){
  
  knots_rows = x_1[seq(1, length(x_1), length.out = n_knots)]
  num_basis_rows = length(knots_rows) + spline_degree - 1
  IDX_BASIS_ROWS = 1:num_basis_rows
  BASIS_ROWS = bsplines(x_1, knots_rows, spline_degree)
  
  knots_columns = x_2[seq(1, length(x_2), length.out = n_knots)]
  num_basis_columns = length(knots_columns) + spline_degree - 1
  IDX_BASIS_COLUMNS = 1:num_basis_columns
  BASIS_COLUMNS = bsplines(x_2, knots_columns, spline_degree)
  
  stan_data = list(n = length(x_1), m = length(x_2), x_1 = x_1, x_2 = x_2, 
                   y = matrix(y, nrow = length(x_1), ncol = length(x_2), byrow = F), 
                   num_basis_rows = num_basis_rows, num_basis_columns = num_basis_columns,
                   BASIS_ROWS = BASIS_ROWS, BASIS_COLUMNS = BASIS_COLUMNS)
  
  fit <- rstan::sampling(model_BSIN,data=stan_data,iter=1000,warmup=200,chains=3, 
                         control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit, file.path(outdir, paste0('2D_BS-IN_', lab, '_nknots_', n_knots, '.rds')))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  samples = extract(fit)
  
  tmp = as.data.table( reshape2::melt(samples$f) )
  setnames(tmp, 2:3, c('rows_idx', 'column_idx'))
  tmp = tmp[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('rows_idx', 'column_idx')]	
  tmp = dcast(tmp, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  tmp1 = as.data.table( reshape2::melt(samples$y_hat) )
  setnames(tmp1, 2:3, c('rows_idx', 'column_idx'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('rows_idx', 'column_idx')]	
  tmp1 = dcast(tmp1, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  tmp1 = data.table(y=y, expand.grid(x_1=x_1,x_2 = x_2), expand.grid(rows_idx = 1:length(x_1), column_idx = 1:length(x_2)))
  tmp = merge(tmp, tmp1, by = c('rows_idx', 'column_idx'))
  
  tmp[, method := paste0('BS-GP-IN\n#knots = ', n_knots)]
  
  if(grepl('lengthscale', lab))
    tmp[, lengthscale := as.numeric(gsub('lengthscale_(.+)', '\\1', lab))]
  
  time = sum(get_elapsed_time(fit))
  tmp[, time := time]
  
  return(list(tmp, fit))
}

run_GP_I_2D = function(x_1, x_2, y, lab, outdir){
  
  stan_data = list(n = length(x_1), m = length(x_2), x_1 = x_1, x_2 = x_2, 
                   y = matrix(y, nrow = length(x_1), ncol = length(x_2), byrow = F))
  
  fit <- rstan::sampling(model_GPI,data=stan_data,iter=1000,warmup=200,chains=3, 
                         control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit, file.path(outdir, paste0('2D_GP-I_', lab, '.rds')))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  samples = extract(fit)
  
  tmp = as.data.table( reshape2::melt(samples$f) )
  setnames(tmp, 2:3, c('rows_idx', 'column_idx'))
  tmp = tmp[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('rows_idx', 'column_idx')]	
  tmp = dcast(tmp, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  tmp1 = as.data.table( reshape2::melt(samples$y_hat) )
  setnames(tmp1, 2:3, c('rows_idx', 'column_idx'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('rows_idx', 'column_idx')]	
  tmp1 = dcast(tmp1, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  tmp1 = data.table(y=y, expand.grid(x_1=x_1,x_2 = x_2), expand.grid(rows_idx = 1:length(x_1), column_idx = 1:length(x_2)))
  tmp = merge(tmp, tmp1, by = c('rows_idx', 'column_idx'))
  
  tmp[, method := paste0('GP-I')]
  
  if(grepl('lengthscale', lab))
    tmp[, lengthscale := as.numeric(gsub('lengthscale_(.+)', '\\1', lab))]
  
  time = sum(get_elapsed_time(fit))
  tmp[, time := time]
  
  return(list(tmp, fit))
}

run_GP_2D = function(x_1, x_2, y, lab, outdir){
  
  stan_data = list(n = length(x_1), m = length(x_2), x_1 = x_1, x_2 = x_2, 
                   y = matrix(y, nrow = length(x_1), ncol = length(x_2), byrow = F))
  
  fit <- rstan::sampling(model_GP,data=stan_data,iter=1000,warmup=200,chains=3, 
                         control = list(max_treedepth = 15, adapt_delta = 0.99))
  saveRDS(fit, file.path(outdir, paste0('2D_GP_', lab, '.rds')))
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  samples = extract(fit)
  
  tmp = as.data.table( reshape2::melt(samples$f) )
  setnames(tmp, 2:3, c('rows_idx', 'column_idx'))
  tmp = tmp[, list( 	q= quantile(value, prob=ps, na.rm = T),
                     q_label=p_labs), 
            by=c('rows_idx', 'column_idx')]	
  tmp = dcast(tmp, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  tmp1 = as.data.table( reshape2::melt(samples$y_hat) )
  setnames(tmp1, 2:3, c('rows_idx', 'column_idx'))
  tmp1 = tmp1[, list( 	q= quantile(value, prob=ps, na.rm = T),
                       q_label=p_labs), 
              by=c('rows_idx', 'column_idx')]	
  tmp1 = dcast(tmp1, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  tmp1 = data.table(y=y, expand.grid(x_1=x_1,x_2 = x_2), expand.grid(rows_idx = 1:length(x_1), column_idx = 1:length(x_2)))
  tmp = merge(tmp, tmp1, by = c('rows_idx', 'column_idx'))
  
  tmp[, method := paste0('GP-SE')]
  
  if(grepl('lengthscale', lab))
    tmp[, lengthscale := as.numeric(gsub('lengthscale_(.+)', '\\1', lab))]
  
  time = sum(get_elapsed_time(fit))
  tmp[, time := time]
  
  return(list(tmp, fit))
}
