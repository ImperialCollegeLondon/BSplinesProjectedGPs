get_expected_resurgence_deaths <- function(C, M, TT,
                               prop_vac_start, week_indices_resurgence,
                               intercept_resurgence0, intercept_resurgence_re, 
                               slope_resurgence0, 
                               vaccine_effect_intercept_diagonal, vaccine_effect_slope_diagonal,
                               vaccine_effect_intercept_cross, vaccine_effect_slope_cross){
  
  intercept_resurgence = matrix(nrow = M, ncol = C, 0)
  slope_resurgence = matrix(nrow = M, ncol = C, 0)

  xi = array(0, dim = c(M, C, TT))
  log_xi = array(0, dim = c(M, C, TT))
  
  for(c in 1:C){
    intercept_resurgence[,c] = rep(intercept_resurgence0[c], M) + intercept_resurgence_re[,c];
    slope_resurgence[,c] = rep(slope_resurgence0[c], M) ;
    
    for(c_prime in 1:C){

      if(c_prime == c){
        intercept_resurgence[,c] = intercept_resurgence[,c] + (prop_vac_start[,c] * rep(vaccine_effect_intercept_diagonal[c], M));
        slope_resurgence[,c] =  slope_resurgence[,c] + prop_vac_start[,c] * rep(vaccine_effect_slope_diagonal[c], M);
      } else{
        intercept_resurgence[,c] = intercept_resurgence[,c] + (prop_vac_start[,c_prime] * rep(vaccine_effect_intercept_cross[c], M));
        slope_resurgence[,c] =  slope_resurgence[,c] + prop_vac_start[,c_prime] * rep(vaccine_effect_slope_cross[c], M);
      }
    }
    
    log_xi[,c,] = matrix(nrow = M, ncol = TT, intercept_resurgence[,c], byrow = F) + matrix(nrow = M, ncol = TT, week_indices_resurgence, byrow =T) * matrix(nrow = M, ncol = TT, slope_resurgence[,c], byrow = F);
    xi[,c,] = exp(log_xi[,c,]);
  }

  
  return(xi)
}

log_likelihood <- function(M, C, TT, r_pdeaths, 
                           prop_vac_start, week_indices_resurgence,
                           parameters){
  
  intercept_resurgence0 <- parameters$intercept_resurgence0
  intercept_resurgence_re <- parameters$intercept_resurgence_re
  slope_resurgence0 <- parameters$slope_resurgence0
  vaccine_effect_intercept_diagonal <- parameters$vaccine_effect_intercept_diagonal
  vaccine_effect_slope_diagonal <- parameters$vaccine_effect_slope_diagonal
  vaccine_effect_intercept_cross <- parameters$vaccine_effect_intercept_cross
  vaccine_effect_slope_cross <- parameters$vaccine_effect_slope_cross
  sigma_r_pdeaths <- parameters$sigma_r_pdeaths
  
  xi = get_expected_resurgence_deaths (C, M, TT, 
                                       prop_vac_start, week_indices_resurgence,
                                       intercept_resurgence0, intercept_resurgence_re, 
                                       slope_resurgence0, 
                                       vaccine_effect_intercept_diagonal, vaccine_effect_slope_diagonal,
                                       vaccine_effect_intercept_cross, vaccine_effect_slope_cross)
  
  log_lik <- 0
  
  for(m in 1:M){
    for(c in 1:C){
      
      shape = (xi[m, c, 1:TT])^2 / sigma_r_pdeaths[m]^2
      scale = sigma_r_pdeaths[m]^2 / xi[m, c, 1:TT]
      log_lik = log_lik + sum(dgamma(r_pdeaths[m, c, 1:TT], shape, scale, log = T))

    }
  }

  return(log_lik)
}

logprior_intercept_resurgence0 <- function(intercept_resurgence0){
  dnorm(intercept_resurgence0, 0, 5, log = T)
}
logprior_slope_resurgence0 <- function(slope_resurgence0){
  dnorm(slope_resurgence0, 0, 5, log = T)
}

logprior_vaccine_effects <- function(slope_resurgence0){
  dnorm(slope_resurgence0, 0, 2.5, log = T)
}

logprior_sigma_intercept_resurgence <- function(sigma_intercept_resurgence){
  dcauchy(sigma_intercept_resurgence, 0, 1, log = T)
}

logprior_intercept_resurgence_re <- function(intercept_resurgence_re, sigma_intercept_resurgence){
  sum(dnorm(intercept_resurgence_re, 0, sigma_intercept_resurgence, log = T))
}

proposal_sigma_intercept_resurgence <- function(proposal_sigma_intercept_star, proposal_sigma_intercept, sigma_proposal_sigma_intercept){
  extraDistr::dtnorm(proposal_sigma_intercept_star, proposal_sigma_intercept, sigma_proposal_sigma_intercept, a = 0, log = T)
}

update_parameters_size_C <- function(n, D, theta, s_theta, parameters, log_prior, proposal = rnorm, log_proposal = NULL){
  
  ### intercept_resurgence0 ###
  theta_star <- parameters[[theta]]
  
  for(c in 1:D){
    
    #propose
    u <- s_theta[c]*proposal(1)
    theta_star[c] <- theta_star[c] + u
    parameters_star <- parameters
    parameters_star[[theta]][c] <- theta_star[c]
    
    #accept-reject
    ratio = log_likelihood(M, C, TT, r_pdeaths[[n]], 
                           prop_vac_start, week_indices_resurgence, parameters_star) - 
      log_likelihood(M, C, TT, r_pdeaths[[n-1]], 
                     prop_vac_start, week_indices_resurgence, parameters) + 
      log_prior(theta_star[c]) - log_prior(parameters[[theta]][c])
    
    if(!is.null(log_proposal)){
      ratio = ratio + 
        log_proposal(parameters[[theta]][c], parameters_star[[theta]][c], s_theta[c]) - log_proposal(parameters_star[[theta]][c], parameters[[theta]][c], s_theta[c]) 
    }
    
    acceptance_prob = min(1, exp(ratio))
    
    if(runif(1) <= acceptance_prob){
      parameters[[theta]][c] <- theta_star[c]
    }
    if(n < warmup+1 & adapt == T){#adaptative move
      s_theta[c] = ramcmc::adapt_S(S = s_theta[c], u = u, current = acceptance_prob, n = n-1)
    }
  }
  
  return(list(parameters[[theta]], s_theta))
}


