#
# Metropolis-Hastings Monte Carlo
#

MetrHastrw <- cmpfun(function(iter, warmup, r_pdeaths, data, adapt = T){
  
  # iter = 1000; warmup = 100
  N = iter #+ warmup + 1
  
  # SAVE DATA OBJECTS IN ENV TO AVOID EXCESSIVE FCT INPUTS
  w_stop_resurgence <<- data$w_stop_resurgence
  w_start_resurgence <<- data$w_start_resurgence; 
  age_from_vac_age_strata <<- data$age_from_vac_age_strata
  age_to_vac_age_strata <<- data$age_to_vac_age_strata
  prop_vac_start <<- t(do.call('rbind', data$prop_vac_start))
  N_COUNTERFACTUAL <<- data$N_COUNTERFACTUAL
  prop_vac_start_counterfactual <<- data$prop_vac_start_counterfactual
  week_indices_resurgence <- data$week_indices_resurgence
  
  #######################
  ### STARTING VALUES ###
  #######################
  
  params = list()
  
  params$intercept_resurgence0 = rep(list(numeric(C)), N)
  params$intercept_resurgence0[[1]] = rnorm(C, 0, 5)
  
  params$sigma_intercept_resurgence = rep(list(numeric(C)), N)
  params$sigma_intercept_resurgence[[1]] = extraDistr::rhcauchy(C, sigma = 1)
  
  params$intercept_resurgence_re = rep(list(matrix(nrow = M, ncol = C, 0)), N)
  params$intercept_resurgence_re[[1]][,1] = rnorm(M, 0, params$sigma_intercept_resurgence[[1]])
  params$intercept_resurgence_re[[1]][,2] = rnorm(M, 0, params$sigma_intercept_resurgence[[1]])
  
  params$slope_resurgence0 = rep(list(numeric(C)), N)
  params$slope_resurgence0[[1]] = rnorm(C, 0, 5)
  
  params$vaccine_effect_intercept_diagonal = rep(list(numeric(C)), N)
  params$vaccine_effect_intercept_diagonal[[1]] = rnorm(C, 0, 2.5)
  params$vaccine_effect_slope_diagonal = rep(list(numeric(C)), N)
  params$vaccine_effect_slope_diagonal[[1]] = rnorm(C, 0, 2.5)

  params$vaccine_effect_intercept_cross = rep(list(numeric(C)), N)
  params$vaccine_effect_intercept_cross[[1]] = rnorm(C, 0, 2.5)
  params$vaccine_effect_slope_cross = rep(list(numeric(C)), N)
  params$vaccine_effect_slope_cross[[1]] = rnorm(C, 0, 2.5)
  
  params$sigma_r_pdeaths = rep(list(numeric(M)), N)
  params$sigma_r_pdeaths[[1]] = extraDistr::rhcauchy(M, sigma = 1)

  # likelihood
  params$log_lik <- numeric(N)

  # fix sigma for proposal
  s_intercept_resurgence0 <- rep(5, C)
  s_slope_resurgence0 <- rep(5, C)
  s_sigma_intercept_resurgence <- rep(1, C)
  s_intercept_resurgence_re <- matrix(nrow = M, ncol = C, 1)
  s_vaccine_effect_intercept_diagonal <- rep(2.5, C)
  s_vaccine_effect_slope_diagonal <- rep(2.5, C)
  s_vaccine_effect_intercept_cross <- rep(2.5, C)
  s_vaccine_effect_slope_cross <- rep(2.5, C)
  s_sigma_r_pdeaths <- rep(1, M)
  
  # parameters
  parameters = list(intercept_resurgence0 = params$intercept_resurgence0[[1]], 
                    intercept_resurgence_re = params$intercept_resurgence_re[[1]], 
                    slope_resurgence0 = params$slope_resurgence0[[1]], 
                    vaccine_effect_intercept_diagonal = params$vaccine_effect_intercept_diagonal[[1]], 
                    vaccine_effect_slope_diagonal = params$vaccine_effect_slope_diagonal[[1]], 
                    vaccine_effect_intercept_cross = params$vaccine_effect_intercept_cross[[1]], 
                    vaccine_effect_slope_cross = params$vaccine_effect_slope_cross[[1]], 
                    sigma_r_pdeaths = params$sigma_r_pdeaths[[1]],
                    sigma_intercept_resurgence = params$sigma_intercept_resurgence[[1]])
  params$log_lik[[1]] <- log_likelihood(M, C, TT, r_pdeaths[[1]], prop_vac_start, week_indices_resurgence, parameters) 
  
  #############
  ### SAMPLE ##
  #############  
  
  for(n in 2:N){
    
    ## intercept_resurgence0 ##
    res <- update_parameters_size_C(n, C, 'intercept_resurgence0', s_intercept_resurgence0, parameters, logprior_intercept_resurgence0)
    parameters$intercept_resurgence0 <- res[[1]]
    s_intercept_resurgence0 <- res[[2]]
    
    ### intercept_resurgence_re ###
    parameters_star <- parameters

    for(c in 1:C){
      
      #propose
      u_sigma_intercept_resurgence <- s_sigma_intercept_resurgence[c]*extraDistr::dtnorm(1, 0, 1, a = 0)
      parameters_star$sigma_intercept_resurgence[c] <- parameters$sigma_intercept_resurgence[c] + u_sigma_intercept_resurgence
      
      u_s_intercept_resurgence_re <- rep(0, M)
      
      for(m in 1:M){
        
        #propose
        u_s_intercept_resurgence_re[m] <- s_intercept_resurgence_re[m,c]*rnorm(1)
        parameters_star$intercept_resurgence_re[m,c] <- parameters$intercept_resurgence_re[m,c] + u_s_intercept_resurgence_re[m]
      }
      
      #accept-reject
      acceptance_prob = min(1, exp(
        log_likelihood(M, C, TT, r_pdeaths[[n]], 
                       prop_vac_start, week_indices_resurgence, parameters) - 
          log_likelihood(M, C, TT, r_pdeaths[[n-1]], 
                         prop_vac_start, week_indices_resurgence, parameters_star) + 
          logprior_sigma_intercept_resurgence(parameters_star$sigma_intercept_resurgence[c]) - logprior_sigma_intercept_resurgence(parameters$sigma_intercept_resurgence[c]) + 
          logprior_intercept_resurgence_re(parameters_star$intercept_resurgence_re[,c], parameters_star$sigma_intercept_resurgence[c]) - logprior_intercept_resurgence_re(parameters$intercept_resurgence_re[,c], parameters$sigma_intercept_resurgence[c]) +
          proposal_sigma_intercept_resurgence(parameters$sigma_intercept_resurgence[c], parameters_star$sigma_intercept_resurgence[c], s_sigma_intercept_resurgence[c]) -  proposal_sigma_intercept_resurgence(parameters_star$sigma_intercept_resurgence[c], parameters$sigma_intercept_resurgence[c], s_sigma_intercept_resurgence[c])
        ))
      if(runif(1) <= acceptance_prob){
        parameters$sigma_intercept_resurgence[c] <- parameters_star$sigma_intercept_resurgence[c]
        parameters$intercept_resurgence_re[,c] <- parameters_star$intercept_resurgence_re[,c]
      }
      if(n < warmup+1 & adapt == T){#adaptative move
        s_sigma_intercept_resurgence[c] = ramcmc::adapt_S(S = s_sigma_intercept_resurgence[c], u = u_sigma_intercept_resurgence, current = acceptance_prob, n = n-1)
        
        for(m in 1:M){
          s_intercept_resurgence_re[m,c] = ramcmc::adapt_S(S = s_intercept_resurgence_re[m,c], u = u_s_intercept_resurgence_re[m], current = acceptance_prob, n = n-1)
          }
      }
    }
    
    ## slope_resurgence0 ##
    res <- update_parameters_size_C(n, C, 'slope_resurgence0', s_slope_resurgence0, parameters, logprior_slope_resurgence0)
    parameters$slope_resurgence0 <- res[[1]]
    s_slope_resurgence0 <- res[[2]]

    ## vaccine_effects ##
    res <- update_parameters_size_C(n, C, 'vaccine_effect_intercept_diagonal', s_vaccine_effect_intercept_diagonal, parameters, logprior_vaccine_effects)
    parameters$vaccine_effect_intercept_diagonal <- res[[1]]
    s_vaccine_effect_intercept_diagonal <- res[[2]]
    
    res <- update_parameters_size_C(n, C, 'vaccine_effect_slope_diagonal', s_vaccine_effect_slope_diagonal, parameters, logprior_vaccine_effects)
    parameters$vaccine_effect_slope_diagonal <- res[[1]]
    s_vaccine_effect_slope_diagonal <- res[[2]]
    
    res <- update_parameters_size_C(n, C, 'vaccine_effect_intercept_cross', s_vaccine_effect_intercept_cross, parameters, logprior_vaccine_effects)
    parameters$vaccine_effect_intercept_cross <- res[[1]]
    s_vaccine_effect_intercept_cross <- res[[2]]
    
    res <- update_parameters_size_C(n, C, 'vaccine_effect_slope_cross', s_vaccine_effect_slope_cross, parameters, logprior_vaccine_effects)
    parameters$vaccine_effect_slope_cross <- res[[1]]
    s_vaccine_effect_slope_cross <- res[[2]]
    
    ## sigma deaths ##
    res <- update_parameters_size_C(n, M, 'sigma_r_pdeaths', s_sigma_r_pdeaths, parameters, logprior_sigma_intercept_resurgence, 
                                    function(x) extraDistr::dtnorm(x, 0, 1, a = 0), proposal_sigma_intercept_resurgence)
    parameters$sigma_r_pdeaths <- res[[1]]
    s_sigma_r_pdeaths <- res[[2]]
    
    # update
    params$intercept_resurgence0[[n]] <- parameters$intercept_resurgence0
    params$intercept_resurgence_re[[n]] <- parameters$intercept_resurgence_re
    params$slope_resurgence0[[n]] <- parameters$slope_resurgence0
    params$vaccine_effect_intercept_diagonal[[n]] <- parameters$vaccine_effect_intercept_diagonal
    params$vaccine_effect_slope_diagonal[[n]] <- parameters$vaccine_effect_slope_diagonal
    params$vaccine_effect_intercept_cross[[n]] <- parameters$vaccine_effect_intercept_cross
    params$vaccine_effect_slope_cross[[n]] <- parameters$vaccine_effect_slope_cross
    params$sigma_r_pdeaths[[n]] <- parameters$sigma_r_pdeaths
    params$sigma_intercept_resurgence[[n]] <- parameters$sigma_intercept_resurgence
    params$log_lik[[n]] <- log_likelihood(M, C, TT, r_pdeaths[[n]], prop_vac_start, week_indices_resurgence, parameters) 
    
    # compute generated quantities
    
  }
    
    
  ######################
  ### REMOVE BURN-IN ###
  ######################
  
  params$intercept_resurgence0 = params$intercept_resurgence0[-c(1:(warmup+1))]
  params$intercept_resurgence_re = params$intercept_resurgence_re[-c(1:(warmup+1))]
  params$slope_resurgence0 = params$slope_resurgence0[-c(1:(warmup+1))]
  params$vaccine_effect_intercept_diagonal = params$vaccine_effect_intercept_diagonal[-c(1:(warmup+1))]
  params$vaccine_effect_slope_diagonal = params$vaccine_effect_slope_diagonal[-c(1:(warmup+1))]
  params$vaccine_effect_intercept_cross = params$vaccine_effect_intercept_cross[-c(1:(warmup+1))]
  params$vaccine_effect_slope_cross = params$vaccine_effect_slope_cross[-c(1:(warmup+1))]
  params$sigma_r_pdeaths = params$sigma_r_pdeaths[-c(1:(warmup+1))]
  params$sigma_intercept_resurgence = params$sigma_intercept_resurgence[-c(1:(warmup+1))]
  params$log_lik = params$log_lik[-c(1:(warmup+1))]
  
  
  return(params)
  
})
