distance <- function (X) 
{
  X = as.matrix(X)
  
  D <- sapply(1:nrow(X), function(i) apply( (matrix(X[i,], ncol = ncol(X), nrow = nrow(X), byrow = T) - X)^2 , 1, sum)  )

  return(D)
}

rbf_D <- function(X,l=1, eps = sqrt(.Machine$double.eps) ){
  # squared exponential kernel
  D <- distance(X)
  Sigma <- exp(-D/l)^2 + diag(eps, nrow(X))
}

generate_2DGP = function(X, l){
  Sigma <- rbf_D(X,l=l)
  y <- MASS::mvrnorm(n = 1, mu = rep(0,dim(Sigma)[1]), Sigma = Sigma) 
  return(y)
}

find_count_2D = function(mean, nu){
  
  stan_generating_model = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'rnbinom.stan') )
  stan_data = list(mu = mean, n = length(mean), nu = nu)
  fit <- rstan::sampling(stan_generating_model,data=stan_data,iter=1, chains = 1, algorithm = "Fixed_param", verbose =F)
  y = as.vector(extract(fit)$y)
  
  return(y)
  
}

simulate_data <- function(X, l, nu, indir){
  
  file = file.path(indir, 'simulations', 'data', paste0('simulated_data_', lab, '.rds'))
  
  if(file.exists(file)){
    tmp = readRDS(file)
    y <<- tmp[[1]]
    y_mean <<- tmp[[2]]
  }else{
    y_mean <<- exp( generate_2DGP(X, lengthscales[i]) )
    y <<- find_count_2D(mean = y_mean, nu = nus[i])
    saveRDS(list(y, y_mean), file)
  }

}



run_spatial_model_2D = function(x_1, x_2, coordinates_training, y, y_mean, lab, method, stan_model, outdir, overwrite = F, 
                                n_knots = NULL, spline_degree = NULL){
  
  if(grepl('SPLINE', method)){
    knots_rows = x_1[seq(1, length(x_1), length.out = n_knots)]
    num_basis_rows = length(knots_rows) + spline_degree - 1
    IDX_BASIS_ROWS = 1:num_basis_rows
    BASIS_ROWS = bsplines(x_1, knots_rows, spline_degree)
  
    knots_columns = x_2[seq(1, length(x_2), length.out = n_knots)]
    num_basis_columns = length(knots_columns) + spline_degree - 1
    IDX_BASIS_COLUMNS = 1:num_basis_columns
    BASIS_COLUMNS = bsplines(x_2, knots_columns, spline_degree)

    stan_data = list(n = length(x_1), m = length(x_2), 
                   N = nrow(coordinates_training), y = y[coordinates_training$index], 
                   coordinates = coordinates_training[, 1:2],
                   num_basis_rows = num_basis_rows, num_basis_columns = num_basis_columns,
                   BASIS_ROWS = BASIS_ROWS, BASIS_COLUMNS = BASIS_COLUMNS,
                   IDX_BASIS_ROWS = IDX_BASIS_ROWS, IDX_BASIS_COLUMNS = IDX_BASIS_COLUMNS)
  } else{
    stan_data = list(n = length(x_1), m = length(x_2), x_1 = x_1, x_2 = x_2, 
                     N = nrow(coordinates_training), y = y[coordinates_training$index], 
                     coordinates = coordinates_training[, 1:2])
  }

  if(grepl('P-SPLINE', method)){
    K = num_basis_rows * num_basis_columns
    A = find_adjacency_matrix(num_basis_rows, num_basis_columns)
    tmp = subset(reshape2::melt( A ), value == 1)
    
    stan_data = c(stan_data, list(K = K, node1 = tmp$Var2, node2 = tmp$Var1, N_edges = nrow(tmp)))
  }
  
  file = file.path(outdir, paste0(method, '_', lab, '_nknots_', n_knots, '.rds'))
  file2 = file.path(outdir, paste0(method, '_', lab, '_nknots_', n_knots, '_summary.rds'))
  
  if(!file.exists(file) | overwrite){
    fit <- rstan::sampling(stan_model,data=stan_data,iter=2000,warmup=500,chains=3,
                          control = list(max_treedepth = 15, adapt_delta = 0.99))
    saveRDS(fit, file)
  } else{
    fit = readRDS(file)
  }
  
  if(file.exists(file2) & !overwrite){
    tmp = readRDS(file2)
    return(list(tmp, fit))
  }
  
  # predictions
  samples = rstan::extract(fit)
  stan_data_predictions <- list(M = dim(samples$alpha)[1], 
                                n = dim(samples$alpha)[2], 
                                m = dim(samples$alpha)[3], 
                                alpha = samples$alpha, 
                                theta = samples$theta)
  
  stan_prediction_model = rstan::stan_model( file.path(indir, 'simulations', 'stan-models', 'rnbinom_matrix.stan') )
  pred <- rstan::sampling(stan_prediction_model,
                         data=stan_data_predictions,iter=1, chains = 1, algorithm = "Fixed_param", verbose =F)
  

  # summarize
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp = as.data.table( reshape2::melt(samples$f) )
  setnames(tmp, 2:3, c('rows_idx', 'column_idx'))
  tmp = tmp[, list( 	q= c(quantile(value, prob=ps, na.rm = T), mean(value)),
                     q_label=c(p_labs, 'mean')), 
            by=c('rows_idx', 'column_idx')]	
  tmp = dcast(tmp, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp[, variable := 'f']
  
  samples.pred <- rstan::extract(pred, pars = 'y_hat')
  tmp1 = as.data.table( reshape2::melt(samples.pred$y_hat[1,,,]) )
  setnames(tmp1, 2:3, c('rows_idx', 'column_idx'))
  tmp1 = tmp1[, list( 	q= c(quantile(value, prob=ps, na.rm = T), mean(value)),
                       q_label=c(p_labs, 'mean')), 
              by=c('rows_idx', 'column_idx')]	
  tmp1 = dcast(tmp1, rows_idx + column_idx ~ q_label, value.var = "q")
  tmp1[, variable := 'y_hat']
  
  tmp = rbind(tmp, tmp1)
  
  y_training = y
  y_training[!1:length(y) %in% coordinates_training$index] = NA
  tmp1 = data.table(y=y, 
                    y_training = y_training,
                    y_mean = y_mean,
                    expand.grid(x_1=x_1,x_2 = x_2), expand.grid(rows_idx = 1:length(x_1), column_idx = 1:length(x_2)))
  tmp = merge(tmp, tmp1, by = c('rows_idx', 'column_idx'))
  
  tmp[, method := paste0(method, '\n#knots = ', n_knots)]
  
  if(grepl('lengthscale', lab))
    tmp[, lengthscale := as.numeric(gsub('.*lengthscale_(.+)', '\\1', lab))]
  
  time = sum(get_elapsed_time(fit))
  tmp[, time := time]
  
  neff = summary(fit)$summary[, 9]
  tmp[, min_neff := min(neff)]
  tmp[, max_neff := max(neff)]

  saveRDS(tmp, file2)
  return(list(tmp, fit))
}
