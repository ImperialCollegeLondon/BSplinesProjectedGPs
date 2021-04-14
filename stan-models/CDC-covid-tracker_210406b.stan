data{
  int<lower=0> W; // number of weeks
  int<lower=0> A; // continuous age
  int<lower=0> B; // original age bands
  int<lower=0,upper=B> N_idx_non_missing[W];
  int<lower=0,upper=B> N_idx_missing[W];
  int<lower=-1,upper=B> idx_non_missing[B,W]; // indices non-missing deaths
  int<lower=-1,upper=B> idx_missing[B,W]; // indices missing deaths
  vector[A] age; // age continuous
  int deaths[B,W]; // cumulative deaths in age band b at time n
  int age_from_state_age_strata[B]; // age from of age band b
  int age_to_state_age_strata[B];// age to of age band b
  int min_count_censored[B,W]; // range of the censored data
  int max_count_censored[B,W]; // range of the censored data
  
  //splines
  int num_basis;
  matrix[num_basis, A] BASIS; 
}

transformed data {
  real covariates_w[W];
  real covariates_k[num_basis];
  
  for(w in 1:W) 
    covariates_w[w] = w;
    
  for(i in 1:num_basis)
    covariates_k[i] = i;
}

parameters {
  real gamma[W]; 
  vector<lower=0>[W] nu;
  real<lower=0> lambda[W];
  
  vector[W] delta_w;  
  vector[num_basis] delta_k;  
  real<lower=0> rho_w;             // length-scale
  real<lower=0> sigma_w;           // amplitude
  real<lower=0> rho_k;             // length-scale
  real<lower=0> sigma_k;           // amplitude
  
  simplex[2] weights_mixture;
}

transformed parameters {
  vector<lower=0>[W] theta = nu ./ (1 + nu);
  matrix[A,W] phi;
  matrix[A,W] alpha;
  matrix[B,W] alpha_reduced;
  vector[num_basis] zeta; 
  vector[W] eta; 
  matrix[W,num_basis] beta; 
  
  {
    matrix[W, W] cov =  cov_exp_quad(covariates_w, sigma_w, rho_w) + diag_matrix(rep_vector(1e-10, W));
    matrix[W, W] L_cov = cholesky_decompose(cov);
    eta = L_cov * delta_w;
  }

  {
    matrix[num_basis, num_basis] cov =  cov_exp_quad(covariates_k, sigma_k, rho_k) + diag_matrix(rep_vector(1e-10, num_basis));
    matrix[num_basis, num_basis] L_cov = cholesky_decompose(cov);
    zeta = L_cov * delta_k;
  }
  
  for(w in 1:W)
  {
    for( i in 1:num_basis)
      beta[w,i] = weights_mixture[1] * eta[w] + weights_mixture[2] * zeta[i];
    
    phi[:,w] = softmax(gamma[w] + to_vector(beta[w,:]*BASIS) ); 
    
    alpha[:,w] = phi[:,w] * lambda[w] / nu[w];

    for(b in 1:B){
      alpha_reduced[b,w] = sum(alpha[age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
    }
  }

}

model {
  nu ~ exponential(1);
  gamma ~ normal(0, 1); 
  delta_w ~ normal(0, 1);
  delta_k ~ normal(0, 1);
  
  rho_w ~ inv_gamma(5, 5);
  sigma_w ~ std_normal();
  rho_k ~ inv_gamma(5, 5);
  sigma_k ~ std_normal();
  
  for(w in 1:W){

    lambda[w] ~ exponential(1.0 / sum(deaths[idx_non_missing[1:N_idx_non_missing[w],w],w]));

    target += neg_binomial_lpmf(deaths[idx_non_missing[1:N_idx_non_missing[w],w],w] | alpha_reduced[idx_non_missing[1:N_idx_non_missing[w],w], w] , theta[w] );
  
    if(N_idx_missing[w] > 0){

      for(n in 1:N_idx_missing[w])
        for(i in min_count_censored[idx_missing[n,w],w]:max_count_censored[idx_missing[n,w],w])
          target += neg_binomial_lpmf( i | alpha_reduced[idx_missing[n,w], w] , theta[w] ) ;
      }
        
          
    }

  }

generated quantities {
  real log_lik[W];
  // int deaths_predict[A,W];
  // int deaths_predict_state_age_strata_non_missing[B,W] = rep_array(0, B, W);
  int deaths_predict_state_age_strata[B,W];

  for(w in 1:W){
    log_lik[w] = neg_binomial_lpmf(deaths[idx_non_missing[1:N_idx_non_missing[w],w],w] | alpha_reduced[idx_non_missing[1:N_idx_non_missing[w],w], w] , theta[w] );
    
    if(N_idx_missing[w] > 0){
      
      for(n in 1:N_idx_missing[w])
        for(i in min_count_censored[idx_missing[n,w],w]:max_count_censored[idx_missing[n,w],w])
          log_lik[w] += neg_binomial_lpmf( i | alpha_reduced[idx_missing[n,w], w] , theta[w] ) ;
    }

    
    // deaths_predict[:,w] = neg_binomial_rng(alpha[:,w], theta[w]);
    // deaths_predict_state_age_strata_non_missing[idx_non_missing[1:N_idx_non_missing[w],w],w] = neg_binomial_rng(alpha_reduced[idx_non_missing[1:N_idx_non_missing[w],w], w], theta[w]);
    deaths_predict_state_age_strata[:,w] = neg_binomial_rng(alpha_reduced[:,w], theta[w]);
  }

}


