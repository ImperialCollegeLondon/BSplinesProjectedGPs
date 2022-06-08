prepare_stan_data = function(deathByAge, loc_name, ref_date, last_date_previous_spec = NULL)
  {
  
  # create map of original age groups 
  df_state_age_strata = unique(select(deathByAge, age_from, age_to, age))
  df_state_age_strata[, age_index := 1:nrow(df_state_age_strata)]
  df_state_age_strata[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_state_age_strata[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  
  # number of age groups 
  B = nrow(df_state_age_strata)
  
  # map week index
  last_date_previous_spec <<- max(deathByAge$date)
  WEEKS = seq.Date(min(deathByAge$date), max(deathByAge$date), by = 'week')
  W = length(WEEKS)
  IDX_WEEKS_OBSERVED = which(WEEKS %in% unique(deathByAge$date))
  W_OBSERVED = length(IDX_WEEKS_OBSERVED)
  
  IDX_WEEKS_NON_OBSERVED = which(!WEEKS %in% unique(deathByAge$date))
  IDX_WEEKS_OBSERVED_REPEATED = IDX_WEEKS_OBSERVED
  if(length(IDX_WEEKS_NON_OBSERVED) > 0){
    IDX_WEEKS_OBSERVED_REPEATED = insert.at(IDX_WEEKS_OBSERVED, min(IDX_WEEKS_NON_OBSERVED - 1), rep(min(IDX_WEEKS_NON_OBSERVED) - 1,length(IDX_WEEKS_NON_OBSERVED)) )
    IDX_WEEKS_OBSERVED_REPEATED[(max(IDX_WEEKS_NON_OBSERVED)+1):W] = IDX_WEEKS_OBSERVED_REPEATED[(max(IDX_WEEKS_NON_OBSERVED)+1):W] -length(IDX_WEEKS_NON_OBSERVED)
  }
  
  stopifnot(length(IDX_WEEKS_OBSERVED_REPEATED) == W)
  
  df_week <<- data.table(week_index = 1:W, date = WEEKS)
  if(W == 0) return(NULL)
  
  # ref date
  W_ref_index = 6*4
  w_ref_index = sort(rev(subset(df_week, date <= ref_date )$week_index)[1:W_ref_index])
  
  # create list
  N_idx_non_missing = vector(mode = 'list', length = length(loc_name))
  idx_non_missing =vector(mode = 'list', length = length(loc_name))
  deaths = vector(mode = 'list', length = length(loc_name))
  N_missing = vector(mode = 'list', length = length(loc_name))
  min_count_censored = vector(mode = 'list', length = length(loc_name))
  max_count_censored = vector(mode = 'list', length = length(loc_name))
  sum_count_censored = vector(mode = 'list', length = length(loc_name))
  start_or_end_period = vector(mode = 'list', length = length(loc_name))
  age_missing = vector(mode = 'list', length = length(loc_name))
  idx_weeks_missing_min = vector(mode = 'list', length = length(loc_name))
  idx_weeks_missing_max = vector(mode = 'list', length = length(loc_name))
  inv_sum_deaths = vector(mode = 'list', length = length(loc_name))
  
  data = NULL
  
  for(m in 1:length(loc_name)){
    tmp = subset(deathByAge, loc_label == loc_name[m])
    tmp = tmp[order(date, age_from)]
    data = rbind(data, tmp)

    # some checks
    stopifnot(all(tmp$age_from <= tmp$age_to))
    # if(!is.null(last_date_previous_spec))
      # stopifnot(last_date_previous_spec == min(tmp$date)-7)
    
    if(!all(WEEKS %in% unique(tmp$date)))
      cat('\n Missing weeks: ', as.character(WEEKS[which(!WEEKS %in% unique(tmp$date))]), '\n')
    
    # smooth weekly deaths on the spike
    date.spike = as.Date("2021-01-23")
    ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}
    tmp[, smooth.weekly.deaths := ma(weekly.deaths, 5), by = c('age')]
    tmp[is.na(smooth.weekly.deaths), smooth.weekly.deaths := weekly.deaths, by = c('age')]
    tmp[date == date.spike, weekly.deaths := as.integer(smooth.weekly.deaths)]
    tmp = select(tmp, - smooth.weekly.deaths)
    
    # create map of original age groups without NA 
    N_idx_non_missing_m = vector(mode = 'integer', length = W_OBSERVED)
    idx_non_missing_m = matrix(nrow = B, ncol = W_OBSERVED, 0)
    deaths_m = matrix(nrow = B, ncol = W_OBSERVED, 0)
    
    for(w in 1:W_OBSERVED){
      
      Week = sort(unique(tmp$date))[w]
      
      tmp1 = subset(tmp, date == Week & !is.na( weekly.deaths ))
      df_non_missing = unique(select(tmp1, age_from, age_to, age))
      
      # number of non missing age category 
      N_idx_non_missing_m[w] = nrow(df_non_missing)
      
      # index non missing
      .idx_non_missing = which(df_state_age_strata$age %in% df_non_missing$age)
      idx_non_missing_m[,w] = c(.idx_non_missing, rep(-1, B - length(.idx_non_missing)))
      
      # deaths
      tmp1 = copy(tmp)
      tmp1[is.na(weekly.deaths), weekly.deaths := -1]
      deaths_m = reshape2::dcast(tmp1, age ~ date, value.var = 'weekly.deaths')[,-1]
    }
    
    # create map for series with missing data
    N_missing_m = 0
    min_count_censored_m = c(); max_count_censored_m = c(); sum_count_censored_m = c(); 
    start_or_end_period_m = c(); age_missing_m = c()
    idx_weeks_missing = list(); N_weeks_missing = list()
    
    for( a in 1:length(unique(tmp$age)) )
    {
      Age = unique(tmp$age)[a]
      tmp1 = subset(tmp, age == Age)
      
      .idx_missing = which(is.na(tmp1$weekly.deaths))
      .idx_missing_full = which(df_week$date %in% tmp1[is.na(weekly.deaths)]$date)
      
      if(length(.idx_missing) == 0) next
      
      stopifnot(all(min(.idx_missing):max(.idx_missing) %in% .idx_missing))
      
      N_missing_m = N_missing_m + 1
      N_weeks_missing[[N_missing_m]] = length(.idx_missing_full)
      idx_weeks_missing[[N_missing_m]] = .idx_missing_full
      
      age_missing_m[N_missing_m] = a
      
      if(is.na(unique(tmp1$sum.weekly.deaths[.idx_missing])))
      {
        stopifnot(length(unique(tmp1$min.sum.weekly.deaths[.idx_missing])) == 1)
        min_count_censored_m[N_missing_m] = unique(tmp1$min.sum.weekly.deaths[.idx_missing])
        max_count_censored_m[N_missing_m] = unique(tmp1$max.sum.weekly.deaths[.idx_missing])
        sum_count_censored_m[N_missing_m] = -1
        start_or_end_period_m[N_missing_m] = 1
      } else {
        sum_count_censored_m[N_missing_m] = unique(tmp1$sum.weekly.deaths[.idx_missing])
        min_count_censored_m[N_missing_m] = max_count_censored_m[N_missing_m] = -1
        start_or_end_period_m[N_missing_m] = 0
      }
      
    }
    
    
    for(n in 1:N_missing_m){
      
      if( any(c(min(IDX_WEEKS_NON_OBSERVED) - 1, max(IDX_WEEKS_NON_OBSERVED) + 1) %in% idx_weeks_missing[[n]]) ){
        
        idx_m1 = which(idx_weeks_missing[[n]] == min(IDX_WEEKS_NON_OBSERVED) - 1)
        idx_p1 = which(idx_weeks_missing[[n]] == max(IDX_WEEKS_NON_OBSERVED) + 1)
        
        if(length(idx_m1) == 1 & length(idx_p1) == 0){
          idx_weeks_missing[[n]] = c(idx_weeks_missing[[n]], IDX_WEEKS_NON_OBSERVED)
          N_weeks_missing[[n]] = N_weeks_missing[[n]] + length(IDX_WEEKS_NON_OBSERVED)
        } else if(length(idx_m1) == 0 & length(idx_p1) == 1){
          idx_weeks_missing[[n]] = c(IDX_WEEKS_NON_OBSERVED, idx_weeks_missing[[n]])
          N_weeks_missing[[n]] = N_weeks_missing[[n]] + length(IDX_WEEKS_NON_OBSERVED)
        } else {
          idx_weeks_missing[[n]] = insert.at(idx_weeks_missing[[n]], idx_p1 - 1, IDX_WEEKS_NON_OBSERVED)
          N_weeks_missing[[n]] = N_weeks_missing[[n]] + length(IDX_WEEKS_NON_OBSERVED)
        }
        
      }
    }
    
    idx_weeks_missing_min_m = list(); idx_weeks_missing_max_m = list();
    for(n in 1:N_missing_m){
      idx_weeks_missing_min_m[[n]] = min(idx_weeks_missing[[n]])
      idx_weeks_missing_max_m[[n]] = max(idx_weeks_missing[[n]])
    }
    
    
    # inverse sum of deaths
    inv_sum_deaths_m = vector(mode = 'double', length = W_OBSERVED)
    for(w in 1:W_OBSERVED){
      
      deaths_w = sum(deaths_m[idx_non_missing_m[1:N_idx_non_missing_m[w],w],w])
      
      if(deaths_w == 0){
        inv_sum_deaths_m[w] = 1
        next
      }
      
      inv_sum_deaths_m[w] = 1 / deaths_w
    }
    
    
    N_idx_non_missing[[m]] = N_idx_non_missing_m
    idx_non_missing[[m]] = idx_non_missing_m
    deaths[[m]] = as.matrix(deaths_m)
    N_missing[[m]] = N_missing_m
    min_count_censored[[m]] = min_count_censored_m
    max_count_censored[[m]] = max_count_censored_m
    sum_count_censored[[m]] = sum_count_censored_m
    start_or_end_period[[m]] = start_or_end_period_m
    age_missing[[m]] = age_missing_m
    idx_weeks_missing_min[[m]] = unlist(idx_weeks_missing_min_m)
    idx_weeks_missing_max[[m]] = unlist(idx_weeks_missing_max_m)
    inv_sum_deaths[[m]] = inv_sum_deaths_m
    
  }

  tmp <<- merge(data, df_state, by = c('loc_label', 'code'))
  
  # create stan data list
  stan_data <- list()
  
  # age bands
  stan_data = c(stan_data, 
                list(W = W,
                     W_OBSERVED = W_OBSERVED,
                     W_NOT_OBSERVED = length(IDX_WEEKS_NON_OBSERVED),
                     IDX_WEEKS = 1:W,
                     IDX_WEEKS_OBSERVED = IDX_WEEKS_OBSERVED,
                     IDX_WEEKS_OBSERVED_REPEATED = IDX_WEEKS_OBSERVED_REPEATED,
                     W_ref_index = W_ref_index,
                     w_ref_index = w_ref_index,
                     A = nrow(df_age_continuous),
                     age = df_age_continuous$age,
                     B = B, 
                     age_from_state_age_strata = df_state_age_strata$age_from_index,
                     age_to_state_age_strata = df_state_age_strata$age_to_index,
                     N_idx_non_missing = N_idx_non_missing,
                     idx_non_missing = idx_non_missing,
                     N_missing = as.array(unlist(N_missing)),
                     start_or_end_period = start_or_end_period,
                     age_missing = age_missing,
                     idx_weeks_missing_min = idx_weeks_missing_min,
                     idx_weeks_missing_max = idx_weeks_missing_max,
                     sum_count_censored = sum_count_censored,
                     min_count_censored = min_count_censored,
                     max_count_censored = max_count_censored,
                     deaths = deaths,
                     M = length(loc_name),
                     inv_sum_deaths = inv_sum_deaths
                ))

  stan_data$N_missing = unlist(stan_data$N_missing)
  
  variables = c('start_or_end_period', 'age_missing', 'idx_weeks_missing_min', 'idx_weeks_missing_max',
                'sum_count_censored', 'min_count_censored', 'max_count_censored')
  for(var in variables){
    
    n_max = max(unlist(lapply(stan_data[[var]], length)))
    stan_data[[var]] = lapply(stan_data[[var]], function(x) c(x, rep(-1, n_max - length(x))))

  }
  
  return(stan_data)
}

add_JHU_data <- function(stan_data, df_week, Code){
  
  JHUData = as.data.table(JHUData)
  JHUData = subset(JHUData, code %in% Code)
  JHUData[, date := as.Date(date)]
  tmp2 = JHUData[, list(daily_deaths = sum(daily_deaths)), by = c('date', 'code')]
  tmp2 = tmp2[order(code, date)]
  tmp2[, cum.death := cumsum(daily_deaths)]
  tmp2 = unique(subset(tmp2, date %in% c(df_week$date, max(df_week$date)+7)))
  tmp2[, weekly.deaths := c(diff(cum.death),NA), by = 'code']
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = 'date')
  tmp2 = select(tmp2, code, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))
  tmp2[weekly.deaths<0, weekly.deaths := 0]
  tmp2 = tmp2[order(code)]
  
  stan_data[['deaths_JHU']] = as.matrix(reshape2::dcast(tmp2, code ~ week_index, value.var = 'weekly.deaths')[,-1])
  
  return(stan_data)
}

add_2D_splines_stan_data = function(stan_data, spline_degree = 3, n_knots_rows = 8, n_knots_columns = 8, 
                                    knots_rows = NULL, knots_columns = NULL)
{
  
  if(is.null(knots_rows)){
    knots_rows = stan_data$age[seq(1, length(stan_data$age), length.out = n_knots_rows)] 
  }
  if(is.null(knots_columns)){
    knots_columns = (1:stan_data$W)[seq(1, stan_data$W, length.out = n_knots_columns)]
  }
  
  cat('Knots on the age dimension ', knots_rows, '\n')
  cat('Knots on the week dimension ',  knots_columns, '\n')
  cat('Number of knots on the age dimension: ', length(knots_rows), '\n')
  cat('Number of knots on the week dimension: ', length(knots_columns),'\n')
  
  stan_data$num_basis_rows = length(knots_rows) + spline_degree - 1
  stan_data$num_basis_columns = length(knots_columns) + spline_degree - 1
  
  stan_data$IDX_BASIS_ROWS = 1:stan_data$num_basis_rows
  stan_data$IDX_BASIS_COLUMNS = 1:stan_data$num_basis_columns
  
  stan_data$BASIS_ROWS = bsplines(stan_data$age, knots_rows, spline_degree)
  stan_data$BASIS_COLUMNS = bsplines(1:stan_data$W, knots_columns, spline_degree)
  
  stopifnot(all( apply(stan_data$BASIS_ROWS, 1, sum) > 0  ))
  stopifnot(all( apply(stan_data$BASIS_COLUMNS, 1, sum) > 0  ))

  return(stan_data)
}

bspline = function(x, k, order, intervals)
  {
  
  if(order == 1){
    return(x >= intervals[k] & x < intervals[k+1])
  }
  
  w1 = 0; w2 = 0
  
  if(intervals[k] != intervals[k+order-1])
    w1 = (x - intervals[k]) / (intervals[k+order-1] - intervals[k])
  if(intervals[k+1] != intervals[k+order])
    w2 = 1 - (x - intervals[k+1]) / (intervals[k+order] - intervals[k+1])
  
  spline = w1 * bspline(x, k, order - 1, intervals) +
    w2 * bspline(x, k+1, order - 1, intervals)
  
  return(spline)
}

find_intervals = function(knots, degree, repeating = T)
  {
  
  K = length(knots)
  
  intervals = vector(mode = 'double', length = 2*degree + K)
  
  # support of knots
  intervals[(degree+1):(degree+K)] = knots
  
  # extreme
  if(repeating)
  {
    intervals[1:degree] = min(knots)
    intervals[(degree+K+1):(2*degree+K)] = max(knots)
  } else {
    gamma = 0.1
    intervals[1:degree] = min(knots) - gamma*degree:1
    intervals[(degree+K+1):(2*degree+K)] = max(knots) + gamma*1:degree
  }

  return(intervals)
}

bsplines = function(data, knots, degree)
 {
  K = length(knots)
  num_basis = K + degree - 1
  
  intervals = find_intervals(knots, degree)
  
  m = matrix(nrow = num_basis, ncol = length(data), 0)
  
  for(k in 1:num_basis)
  {
    m[k,] = bspline(data, k, degree + 1, intervals) 
  }
  
  m[num_basis,length(data)] = 1
  
  return(m)
}

add_adjacency_matrix_stan_data = function(stan_data, n, m)
  {

  N = n * m
  
  A = find_adjacency_matrix(n,m)

  stan_data$K = N
  stan_data$Adj = A
  stan_data$Adj_n = sum(A) / 2

  return(stan_data)
}

find_adjacency_matrix = function(n,m){
  N = n * m
  
  A = matrix(nrow = N, ncol = N, 0)
  
  for(i in 1:n){
    
    for(j in 1:m){
      
      #cat('\n Processing ', i, j)
      idx_row = n*(j-1) + i
      
      # top
      if(i - 1 > 0){
        idx_col = n*(j-1) + i - 1
        A[idx_row,idx_col] = 1
      }
      
      # bottom
      if(i + 1 <= n){
        idx_col = n*(j-1) + i + 1
        A[idx_row,idx_col] = 1
      }
      
      # left
      if(j - 1 > 0){
        idx_col = n*(j-2) + i
        A[idx_row,idx_col] = 1
      }
      
      # right
      if(j + 1 <= m){
        idx_col = n*j + i
        A[idx_row,idx_col] = 1
      }
      
    }
  }
  return(A)
}

add_nodes_stan_data = function(stan_data)
 {
  tmp = reshape2::melt( stan_data$Adj )
  tmp = subset(tmp, value == 1)
  
  stan_data$node1 = tmp$Var2
  stan_data$node2 = tmp$Var1
  stan_data$N_edges = nrow(tmp)
  
  return(stan_data)
}

add_prior_parameters_lambda = function(stan_data, distribution = 'exponential')
{

  lambda_prior_parameters = vector(mode = 'list', length = length(stan_data$inv_sum_deaths))
  
  for(m in 1:length(stan_data$inv_sum_deaths)){
    
    tmp1 = data.table(weekly_deaths = 1/stan_data$inv_sum_deaths[[m]])
    set(tmp1, NULL, 'diff_sum_deaths', c(NA, tmp1[, diff(weekly_deaths)]))
    tmp1[, rel_diff := ifelse(weekly_deaths, diff_sum_deaths / weekly_deaths, weekly_deaths)]
    tmp1[, diff_sum_deaths_abs := abs(diff_sum_deaths)]
    tmp1[, source := 'CDC']
    
    if(0)
    {
      ggplot(tmp1, aes(x = diff_sum_deaths_abs, y = weekly_deaths, col = source)) + 
        geom_point() 
    }
    
    lin_fit = lm(diff_sum_deaths_abs~weekly_deaths-1, data = tmp1)
    
    mean = tmp1$weekly_deaths
    sd = predict(lin_fit, tmp1) * 2
    
    if( any(sd == 0) ) sd[sd ==0] = 0.01
    if( any(mean == 0) ) mean[mean ==0] = 0.01
    
    if(distribution == 'gamma'){
      alpha = mean^2 / (sd^2)
      beta = mean / (sd^2) 
      lambda_prior_parameters[[m]] = rbind(alpha, beta)
    }
    if(distribution == 'log_normal'){
      mu = 2 * log(mean) - 1/2 * log(sd^2 + mean^2)
      sigma = sqrt(-2*log(mean) + log(sd^2 + mean^2))
      lambda_prior_parameters[[m]] = rbind(mu, sigma)
    }
    if(distribution == 'exponential'){
      lambda_prior_parameters[[m]] = mean
    }
    
  }
  
  stan_data[['lambda_prior_parameters']] = lambda_prior_parameters
  
  return(stan_data)
}

insert.at <- function(a, pos, ...){
  dots <- list(...)
  stopifnot(length(dots)==length(pos))
  result <- vector("list",2*length(pos)+1)
  result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
  result[c(FALSE,TRUE)] <- dots
  unlist(result)
}

add_vaccine_prop = function(stan_data, df_week, Code, vaccine_data, resurgence_dates){
  
  delay = 2*7
  
  age_index_min = df_age_vaccination[age == '18-64']$age_index
  
  vaccine_data[, age_index := which(df_age_vaccination$age_from <= age & df_age_vaccination$age_to >= age), by = 'age']
  tmp = vaccine_data[age_index >= age_index_min, list(prop = unique(prop)), by = c('code', 'date', 'loc_label', 'age_index')]
  tmp[, date := date + delay]
  tmp = merge(tmp, resurgence_dates, by = 'code')
  tmp = tmp[code %in% Code]
  tmp = tmp[date %in% df_week$date]
  
  # resurgence period
  tmp = tmp[date >= start_resurgence & date <= stop_resurgence]
  tmp[, idx.resurgence.date := 1:length(date), by = c('age_index', 'loc_label')]
  
  # set proportion
  tmp1 = tmp[, list(pre_prop = prop[date == start_resurgence]), by = c('loc_label', 'age_index')]
  tmp = merge(tmp, tmp1,  by = c('loc_label', 'age_index'))
  tmp[, prop := prop - pre_prop]
  tmp = tmp[order(code)]
  
  prop_vac = list(); prop_vac_start = list(); 
  N_COUNTERFACTUAL = 6; prop_vac_start_counterfactual = list()
  for(i in 1:length(unique(tmp$age_index))){
    
    tmp1 = subset(tmp, age_index == unique(tmp$age_index)[i])
    
    prop_vac[[i]] = as.matrix( reshape2::dcast(tmp1, idx.resurgence.date ~ code, value.var = 'prop')[,-1] )
    prop_vac_start[[i]] = unique(select(tmp1, code, pre_prop))$pre_prop
    prop_vac_start_counterfactual[[i]] = prop_vac_start[[i]]
    
    if(i == 1){
      values_max = rep(max(prop_vac_start_counterfactual[[i]]), length(prop_vac_start_counterfactual[[i]]))
      values_min = rep(min(prop_vac_start_counterfactual[[i]]), length(prop_vac_start_counterfactual[[i]]))
      prop_vac_start_counterfactual[[i]] = matrix(nrow = N_COUNTERFACTUAL, ncol = length(values_max), 
                                                  c(values_max, prop_vac_start[[i]], values_max,
                                                    values_min, prop_vac_start[[i]], values_min), byrow = T)
      
    }
    
    if(i == 2){
      values_max = rep(max(prop_vac_start_counterfactual[[i]]), length(prop_vac_start_counterfactual[[i]]))
      values_min = rep(min(prop_vac_start_counterfactual[[i]]), length(prop_vac_start_counterfactual[[i]]))
      prop_vac_start_counterfactual[[i]] = matrix(nrow = N_COUNTERFACTUAL, ncol = length(values_max), 
                                                  c(prop_vac_start[[i]], values_max, values_max,
                                                    prop_vac_start[[i]], values_min, values_min), byrow = T)
      
    }
  }
  
  stan_data[['age_from_vac_age_strata']] = df_age_vaccination[age_index >= age_index_min]$age_from
  stan_data[['age_to_vac_age_strata']] = df_age_vaccination[age_index >= age_index_min]$age_to
  stan_data[['C']] = nrow(df_age_vaccination[age_index >= age_index_min])
  
  stan_data[['prop_vac']] = prop_vac
  stan_data[['prop_vac_start']] = prop_vac_start
  
  stan_data[['N_COUNTERFACTUAL']] = N_COUNTERFACTUAL
  stan_data[['prop_vac_start_counterfactual']] = prop_vac_start_counterfactual
  
  cat('Minimum vaccine coverage among 18-54', min(stan_data[['prop_vac_start']][[1]]), ' in ', Code[which.min(stan_data[['prop_vac_start']][[1]])], '\n')
  cat('Maximum vaccine coverage among 18-54', max(stan_data[['prop_vac_start']][[1]]), ' in ', Code[which.max(stan_data[['prop_vac_start']][[1]])], '\n\n')
  cat('Minimum vaccine coverage among 65+', min(stan_data[['prop_vac_start']][[2]]), ' in ', Code[which.min(stan_data[['prop_vac_start']][[2]])], '\n')
  cat('Maximum vaccine coverage among 65+', max(stan_data[['prop_vac_start']][[2]]), ' in ', Code[which.max(stan_data[['prop_vac_start']][[2]])], '\n\n')
  
  # substract mean to remove correlations between intercept and vaccine effects parameters
  means <- sapply(stan_data[['prop_vac_start']], mean)
  stan_data[['prop_vac_start']] <- lapply(1:stan_data[['C']], function(x) stan_data[['prop_vac_start']][[x]] - means[x])
  stan_data[['prop_vac_start_counterfactual']] <- lapply(1:stan_data[['C']], function(x) stan_data[['prop_vac_start_counterfactual']][[x]] - means[x])
  
  return(stan_data)
}

find_resurgence_dates <- function(JHUData, deathByAge, Code){
  
  ma <- function(x, n = 5){as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))}
  
  JHUData = as.data.table(JHUData)
  JHUData = subset(JHUData, code %in% Code)
  JHUData[, date := as.Date(date)]
  tmp2 = JHUData[, list(daily_deaths = sum(daily_deaths)), by = c('date', 'code')]
  tmp2 = tmp2[order(code, date)]
  tmp2[, cum.death := cumsum(daily_deaths)]
  tmp2 = unique(subset(tmp2, date %in% c(df_week$date, max(df_week$date)+7)))
  tmp2[, weekly.deaths := c(diff(cum.death),NA), by = 'code']
  tmp2 = unique(subset(tmp2, date %in% df_week$date))
  tmp2 = merge(tmp2, df_week, by = 'date')
  tmp2 = select(tmp2, code, week_index, weekly.deaths)
  tmp2 = subset(tmp2, !is.na(weekly.deaths))
  tmp2[weekly.deaths<0, weekly.deaths := 0]
  tmp2[, dummy := 1:nrow(tmp2)]
  tmp2[, smooth.weekly.deaths := ma(weekly.deaths, 4), by = c('code')]
  tmp2[, diff.smooth.weekly.deaths := c(NA, diff(smooth.weekly.deaths)), by = c('code')]
  tmp2[, lag.smooth.weekly.deaths := (shift(smooth.weekly.deaths)), by = c('code')]
  tmp2[, change.smooth.weekly.deaths := diff.smooth.weekly.deaths / lag.smooth.weekly.deaths]
  
  # find start resurgence
  tmp2 = merge(tmp2, df_week, by = 'week_index')
  tmp3 <- tmp2[change.smooth.weekly.deaths > 0.05 & date >= as.Date('2021-07-01'), list(start_resurgence = min(date) ), by = c('code')]
  # tmp2[is.na(diff.smooth.weekly.deaths), number.positive.days.ahead := NA_real_, by = c('code')]
  # tmp2[!is.na(diff.smooth.weekly.deaths), number.positive.days.ahead := c(lag_neg(diff.smooth.weekly.deaths), NA_real_), by = c('code')]
  # tmp3 <- tmp2[date >= as.Date('2021-07-01'), list(start_resurgence = date[which.max(number.positive.days.ahead)] ), by = c('code')]
  
  # subset(tmp2, code == 'CA')
  # find stop resurgence
  # tmp4 = merge(tmp3, tmp4, by = 'code')
  # max_resurgence_period = tmp4[, min(stop_resurgence - start_resurgence)] / 7
  max_resurgence_period = (as.Date("2021-09-25") - tmp3[, max(start_resurgence)] )/ 7
  # max_resurgence_period = min(tmp2[date >= as.Date('2021-07-01'), list(resurgence_period = max(na.omit(number.positive.days.ahead))), by = c('code')]$resurgence_period)
  tmp3[, stop_resurgence := start_resurgence + 7*max_resurgence_period]
  
  stopifnot(max(tmp3$stop_resurgence) <= max(deathByAge$date))
  
  
  if(0){ #plot
    ggplot(subset(tmp2, date >= as.Date('2021-04-01')), aes(x = date)) + 
      geom_line(aes(y = weekly.deaths), col = 'red')+ 
      geom_line(aes(y = smooth.weekly.deaths), col = 'blue') + 
      facet_wrap(~code, nrow = length(Code), scale = 'free') + 
      geom_vline(data = tmp3, aes(xintercept = start_resurgence), linetype = 'dashed') + 
      geom_vline(data = tmp3, aes(xintercept = stop_resurgence), linetype = 'dashed') + 
      theme_bw()
    ggsave('~/Downloads/file.png', h= 30, w =5) 
  }
  
  tmp3 <- tmp3[order(code)]
  
  return(tmp3)
}

add_resurgence_period = function(stan_data, df_week, resurgence_dates){
  
  resurgence_dates = resurgence_dates[order(code)]
  
  stan_data[['w_start_resurgence']] = sapply(resurgence_dates$start_resurgence, function(x) df_week[date == x]$week_index)
  stan_data[['w_stop_resurgence']] = sapply(resurgence_dates$stop_resurgence, function(x) df_week[date == x]$week_index)
  
  T = sapply(1:length( stan_data[['w_start_resurgence']]), function(x)  stan_data[['w_stop_resurgence']][[x]] - stan_data[['w_start_resurgence']][[x]] + 1)
  stan_data[['T']] = unique(unlist(T))
  stopifnot(length(stan_data[['T']]) == 1)
  
  stan_data[['week_indices_resurgence']] = 0:(stan_data[['T']]-1)
  
  return(stan_data)
}
