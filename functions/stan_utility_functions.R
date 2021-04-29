prepare_stan_data = function(deathByAge, loc_name, ref_date, last_date_previous_spec = NULL)
  {
  
  tmp = subset(deathByAge, loc_label == loc_name)
  tmp <<- tmp[order(date, age_from)]
  
  # some checks
  stopifnot(all(tmp$age_from <= tmp$age_to))
  if(!is.null(last_date_previous_spec))
    stopifnot(last_date_previous_spec == min(tmp$date)-7)
  
  # create map of original age groups 
  df_state_age_strata = unique(select(tmp, age_from, age_to, age))
  df_state_age_strata[, age_index := 1:nrow(df_state_age_strata)]
  df_state_age_strata[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age"]
  df_state_age_strata[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age"]
  
  # number of age groups 
  B = nrow(df_state_age_strata)
  
  # map week index
  last_date_previous_spec <<- max(tmp$date)
  WEEKS = seq.Date(min(tmp$date), max(tmp$date), by = 'week')
  W = length(WEEKS)
  IDX_WEEKS_OBSERVED = which(WEEKS %in% unique(tmp$date))
  W_OBSERVED = length(IDX_WEEKS_OBSERVED)
  
  IDX_WEEKS_NON_OBSERVED = which(!WEEKS %in% unique(tmp$date))
  IDX_WEEKS_OBSERVED_REPEATED = IDX_WEEKS_OBSERVED
  if(length(IDX_WEEKS_NON_OBSERVED) > 0){
    IDX_WEEKS_OBSERVED_REPEATED = insert.at(IDX_WEEKS_OBSERVED, IDX_WEEKS_NON_OBSERVED - 1, IDX_WEEKS_NON_OBSERVED - 1)
    IDX_WEEKS_OBSERVED_REPEATED[(IDX_WEEKS_NON_OBSERVED+1):W] = IDX_WEEKS_OBSERVED_REPEATED[(IDX_WEEKS_NON_OBSERVED+1):W] -1
  }
    
  stopifnot(length(IDX_WEEKS_OBSERVED_REPEATED) == W)
  
  if(!all(WEEKS %in% unique(tmp$date)))
    cat('\n Missing weeks: ', as.character(unique(tmp$date)[which(!WEEKS %in% unique(tmp$date))]), '\n')
  
  df_week <<- data.table(week_index = 1:W, date = WEEKS)
  if(W == 0) return(NULL)
  
  # ref date
  w_ref_index = NA
  if(min(df_week$date) <= ref_date & max(df_week$date) >= ref_date)
    w_ref_index = subset(df_week, date == ref_date )$week_index
  
  # create map of original age groups without NA 
  N_idx_non_missing = vector(mode = 'integer', length = W_OBSERVED)
  idx_non_missing = matrix(nrow = B, ncol = W_OBSERVED, 0)
  deaths = matrix(nrow = B, ncol = W_OBSERVED, 0)
  
  for(w in 1:W_OBSERVED){
    
    Week = sort(unique(tmp$date))[w]
    
    tmp1 = subset(tmp, date == Week & !is.na( daily.deaths ))
    df_non_missing = unique(select(tmp1, age_from, age_to, age))
    
    # number of non missing age category 
    N_idx_non_missing[w] = nrow(df_non_missing)
    
    # index non missing
    .idx_non_missing = which(df_state_age_strata$age %in% df_non_missing$age)
    idx_non_missing[,w] = c(.idx_non_missing, rep(-1, B - length(.idx_non_missing)))

    # deaths
    tmp1 = copy(tmp)
    tmp1[is.na(daily.deaths), daily.deaths := -1]
    deaths = reshape2::dcast(tmp1, age ~ date, value.var = 'daily.deaths')[,-1]
  }
  
  # create map for series with missing data
  N_missing = 0
  min_count_censored = c(); max_count_censored = c(); sum_count_censored = c(); 
  start_or_end_period = c(); age_missing = c()
  idx_weeks_missing = list(); N_weeks_missing = c()
  
  for( a in 1:length(unique(tmp$age)) )
  {
    Age = unique(tmp$age)[a]
    tmp1 = subset(tmp, age == Age)
    
    .idx_missing = which(is.na(tmp1$daily.deaths))
    .idx_missing_full = which(df_week$date %in% tmp1[is.na(daily.deaths)]$date)
    
    if(length(.idx_missing) == 0) next
    
    stopifnot(all(min(.idx_missing):max(.idx_missing) %in% .idx_missing))

    N_missing = N_missing + 1
    N_weeks_missing[N_missing] = length(.idx_missing_full)
    idx_weeks_missing[[N_missing]] = .idx_missing_full
    
    age_missing[N_missing] = a
    
    if(is.na(unique(tmp1$sum.daily.deaths[.idx_missing])))
    {
      stopifnot(length(unique(tmp1$min.sum.daily.deaths[.idx_missing])) == 1)
      min_count_censored[N_missing] = unique(tmp1$min.sum.daily.deaths[.idx_missing])
      max_count_censored[N_missing] = unique(tmp1$max.sum.daily.deaths[.idx_missing])
      sum_count_censored[N_missing] = -1
      start_or_end_period[N_missing] = 1
    } else {
      sum_count_censored[N_missing] = unique(tmp1$sum.daily.deaths[.idx_missing])
      min_count_censored[N_missing] = max_count_censored[N_missing] = -1
      start_or_end_period[N_missing] = 0
    }

  }
  
  for(n in 1:N_missing){
    if( any(c(IDX_WEEKS_NON_OBSERVED - 1, IDX_WEEKS_NON_OBSERVED + 1) %in% idx_weeks_missing[[n]]) ){
      
      idx_m1 = which(idx_weeks_missing[[n]] == IDX_WEEKS_NON_OBSERVED - 1)
      idx_p1 = which(idx_weeks_missing[[n]] == IDX_WEEKS_NON_OBSERVED + 1)
      
      if(length(idx_m1) == 1 & length(idx_p1) == 0){
        idx_weeks_missing[[n]] = c(idx_weeks_missing[[n]], IDX_WEEKS_NON_OBSERVED)
        N_weeks_missing[[n]] = N_weeks_missing[[n]] + 1
      } else if(length(idx_m1) == 0 & length(idx_p1) == 1){
        idx_weeks_missing[[n]] = c(IDX_WEEKS_NON_OBSERVED, idx_weeks_missing[[n]])
        N_weeks_missing[[n]] = N_weeks_missing[[n]] + 1
      } else {
        idx_weeks_missing[[n]] = insert.at(idx_weeks_missing[[n]], idx_p1 - 1, IDX_WEEKS_NON_OBSERVED)
        N_weeks_missing[[n]] = N_weeks_missing[[n]] + 1
      }

    }
  }
      
  for(n in 1:N_missing){
    idx_weeks_missing[[n]] = c(idx_weeks_missing[[n]], rep(-1, max(N_weeks_missing) - N_weeks_missing[n]) )
  }
  idx_weeks_missing = matrix(unlist(idx_weeks_missing), nrow = max(N_weeks_missing) , ncol = N_missing)
  
  # inverse sum of deaths
  sum_deaths = vector(mode = 'double', length = W_OBSERVED)
  inv_sum_deaths = vector(mode = 'double', length = W_OBSERVED)
  for(w in 1:W_OBSERVED){
     
    deaths_w = sum(deaths[idx_non_missing[1:N_idx_non_missing[w],w],w])
    
    sum_deaths[w] = deaths_w
    
    if(deaths_w == 0){
      inv_sum_deaths[w] = 1
      next
    }
      
    inv_sum_deaths[w] = 1 / deaths_w
  }

  # create stan data list
  stan_data <- list()
  
  # age bands
  stan_data = c(stan_data, 
                list(W = W,
                     W_OBSERVED = W_OBSERVED,
                     IDX_WEEKS = 1:W,
                     IDX_WEEKS_OBSERVED = IDX_WEEKS_OBSERVED,
                     IDX_WEEKS_OBSERVED_REPEATED = IDX_WEEKS_OBSERVED_REPEATED,
                     w_ref_index = w_ref_index,
                     A = nrow(df_age_continuous),
                     age = df_age_continuous$age,
                     B = B, 
                     age_from_state_age_strata = df_state_age_strata$age_from_index,
                     age_to_state_age_strata = df_state_age_strata$age_to_index,
                     N_idx_non_missing = N_idx_non_missing,
                     idx_non_missing = idx_non_missing,
                     N_missing = N_missing,
                     start_or_end_period = start_or_end_period,
                     age_missing = age_missing,
                     N_weeks_missing = N_weeks_missing,
                     idx_weeks_missing = idx_weeks_missing,
                     sum_count_censored = sum_count_censored,
                     min_count_censored = min_count_censored,
                     max_count_censored = max_count_censored,
                     deaths = deaths,
                     inv_sum_deaths = inv_sum_deaths,
                     sum_deaths = sum_deaths
                ))

  return(stan_data)
}

add_splines_stan_data = function(stan_data, spline_degree = 3, n_knots = 8)
{
  
  knots = stan_data$age[seq(1, length(stan_data$age), length.out = n_knots)]

  stan_data$num_basis = length(knots) + spline_degree - 1

  stan_data$IDX_BASIS = 1:stan_data$num_basis
  
  stan_data$BASIS = bsplines(stan_data$age, knots, spline_degree)
  
  stopifnot(all( apply(stan_data$BASIS, 1, sum) > 0  ))
  # B <- t(splines::bs(age, knots=age, degree=spline_degree, intercept = T)) 
  return(stan_data)
}

add_diff_matrix = function(stan_data, n, m, order = 2)
{
  library(Matrix)
  
  In <- Diagonal(n)   
  Im <- Diagonal(m) 
  
  # rows difference 
  D0 <- diff(In, diff = order)  
  stan_data$D1 <- as.matrix( kronecker(Im, D0) ) 
  
  # cols difference 
  D0 <- diff(Im, diff = order) 
  stan_data$D2 <- as.matrix( kronecker(D0, In) )
  
  stan_data$D1_N = nrow(stan_data$D1)
  stan_data$D2_N = nrow(stan_data$D2)
  
  return(stan_data)
}

bspline = function(x, k, degree, intervals)
  {
  
  if(degree == 1){
    return(x >= intervals[k] & x < intervals[k+1])
  }
  
  w1 = 0; w2 = 0
  
  if(intervals[k] != intervals[k+degree-1])
    w1 = (x - intervals[k]) / (intervals[k+degree-1] - intervals[k])
  if(intervals[k+1] != intervals[k+degree])
    w2 = 1 - (x - intervals[k+1]) / (intervals[k+degree] - intervals[k+1])
  
  spline = w1 * bspline(x, k, degree - 1, intervals) +
    w2 * bspline(x, k+1, degree - 1, intervals)
  
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

      # # top-left diagonal
      # if(i - 1 > 0 & j - 1 > 0){
      #   idx_col = n*(j-2) + i - 1
      #   A[idx_row,idx_col] = 1
      # }
      #
      # # bottom-right diagonal
      # if(i + 1 <= n & j + 1 <= m){
      #   idx_col = n*j + i + 1
      #   A[idx_row,idx_col] = 1
      # }
      #
      # # top-right diagonal
      # if(i - 1 > 0 & j + 1 <= m){
      #   idx_col = n*(j-2) + i + 1
      #   A[idx_row,idx_col] = 1
      # }
      #
      # # bottom-left diagonal
      # if(i + 1  <= n & j - 1 > 0){
      #   idx_col = n*j + i - 1
      #   A[idx_row,idx_col] = 1
      # }
    }
  }

  stan_data$N = N
  stan_data$Adj = A
  stan_data$Adj_n = sum(A) / 2

  return(stan_data)
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


add_prior_parameters_lambda = function(stan_data, distribution = 'gamma')
{

  tmp1 = data.table(weekly_deaths = stan_data$sum_deaths)
  set(tmp1, NULL, 'weekly_deaths', stan_data$sum_deaths)
  set(tmp1, NULL, 'diff_sum_deaths', c(NA, tmp1[, diff(weekly_deaths)]))
  tmp1[, rel_diff := diff_sum_deaths / weekly_deaths]
  tmp1[, diff_sum_deaths_abs := abs(diff_sum_deaths)]
  tmp1[, source := 'CDC']
  
  if(0)
  {
    ggplot(tmp1, aes(x = diff_sum_deaths_abs, y = weekly_deaths, col = source)) + 
      geom_point() 
  }
  
  lin_fit = lm(weekly_deaths~diff_sum_deaths_abs-1, data = tmp1)
  
  mean = stan_data$sum_deaths
  sd = stan_data$sum_deaths / (lin_fit$coefficients*2)
  
  if(distribution == 'gamma'){
    alpha = mean^2 / (sd^2)
    beta = mean / (sd^2) 
    stan_data$lambda_prior_parameters = rbind(alpha, beta)
  }
  if(distribution == 'log_normal'){
    mu = 2 * log(mean) - 1/2 * log(sd^2 + mean^2)
    sigma = sqrt(-2*log(mean) + log(sd^2 + mean^2))
    stan_data$lambda_prior_parameters = rbind(mu, sigma)
  }
  
  return(stan_data)
}
# 
# add_adjacency_matrix_stan_data = function(stan_data, n, m){
# 
#   N = n * m
#   A = matrix(nrow = N, ncol = N, 0)
# 
#   for(i in 1:n){
# 
#     for(j in 1:m){
# 
#       #cat('\n Processing ', i, j)
#       idx_row = m*(i-1) + j
# 
#       if(i - 1 > 0){
#         idx_col = m*(i-2) + j
#         A[idx_row,idx_col] = 1
#       }
# 
#       if(i + 1 <= n){
#         idx_col = m*i + j
#         A[idx_row,idx_col] = 1
#       }
# 
#       if(j - 1 > 0){
#         idx_col = m*(i-1) + j - 1
#         A[idx_row,idx_col] = 1
#       }
# 
#       if(j + 1 <= m){
#         idx_col = m*(i-1) + j +1
#         A[idx_row,idx_col] = 1
#       }
# 
#     }
#   }
# 
#   stan_data$N = N
#   stan_data$Adj = A
#   stan_data$Adj_n = sum(A) / 2
# 
#   return(stan_data)
# }
# 
insert.at <- function(a, pos, ...){
  dots <- list(...)
  stopifnot(length(dots)==length(pos))
  result <- vector("list",2*length(pos)+1)
  result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
  result[c(FALSE,TRUE)] <- dots
  unlist(result)
}
