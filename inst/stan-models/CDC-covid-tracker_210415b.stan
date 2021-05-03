data{
  int<lower=0> W1; // number of weeks with the first age band specification
  int<lower=0> W2; // number of weeks with the second age band specification
  int<lower=0> W; // number of weeks in total (W1 + W2)
  int<lower=0,upper=W1> W_OBSERVED1; // number of weeks observed with the first age band specification
  int<lower=0,upper=W2> W_OBSERVED2; // number of weeks observed with the secon age band specification
  int<lower=0,upper=W> W_OBSERVED; // number of weeks observed 
  int<lower=1, upper=W1> IDX_WEEKS_OBSERVED1[W_OBSERVED1]; // index of the weeks observed with the first age band specification
  int<lower=1, upper=W2> IDX_WEEKS_OBSERVED2[W_OBSERVED2]; // index of the weeks observed with the first age band specification
  int<lower=1, upper=W1> IDX_WEEKS_OBSERVED_REPEATED1[W1]; // index of the weeks observed where missing is equal to the previous one the first age band specification
  int<lower=1, upper=W2> IDX_WEEKS_OBSERVED_REPEATED2[W2]; // index of the weeks observed  where missing is equal to the previous one with the first age band specification
  int<lower=0,upper=W> w_ref_index; // week index to compare the death prob
  int<lower=0> A; // continuous age
  int<lower=0> B1; // first age band specification
  int<lower=0> B2; // second age band specification
  int<lower=0,upper=max(B1,B2)> N_idx_non_missing[W_OBSERVED];
  int<lower=0,upper=max(B1,B2)> N_idx_missing[W_OBSERVED];
  int<lower=-1,upper=B1> idx_non_missing_1[B1,W_OBSERVED1]; // indices non-missing deaths for W1
  int<lower=-1,upper=B1> idx_missing_1[B1,W_OBSERVED1]; // indices missing deaths for W1
  int<lower=-1,upper=B2> idx_non_missing_2[B2,W_OBSERVED2]; // indices non-missing deaths for W2
  int<lower=-1,upper=B2> idx_missing_2[B2,W_OBSERVED2]; // indices missing deaths for W2
  vector[A] age; // age continuous
  real inv_sum_deaths[W_OBSERVED]; // inverse sum of deaths
  int deaths_1[B1,W_OBSERVED1]; // cumulative deaths in age band b at time n
  int deaths_2[B2,W_OBSERVED2]; // cumulative deaths in age band b at time n
  int age_from_state_age_strata_1[B1]; // age from of age band b
  int age_to_state_age_strata_1[B1];// age to of age band b
  int min_count_censored_1[B1,W_OBSERVED1]; // range of the censored data
  int max_count_censored_1[B1,W_OBSERVED1]; // range of the censored data
  int age_from_state_age_strata_2[B2]; // age from of age band b
  int age_to_state_age_strata_2[B2];// age to of age band b
  int min_count_censored_2[B2,W_OBSERVED2]; // range of the censored data
  int max_count_censored_2[B2,W_OBSERVED2]; // range of the censored data
  
  //splines
  int num_basis;
  matrix[num_basis, A] BASIS; 
  
  // ICAR model
  int N; // W * num_basis
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];
  int<lower=1, upper=N> node2[N_edges];
}

parameters {
  vector[N] beta_raw; 
  vector<lower=0>[W-1] nu_raw;
  vector<lower=0>[W-1] lambda_raw;
}

transformed parameters {
  vector<lower=0>[W] nu = append_row(nu_raw[IDX_WEEKS_OBSERVED_REPEATED1], nu_raw[IDX_WEEKS_OBSERVED_REPEATED2]);
  vector<lower=0>[W] lambda = append_row(lambda_raw[IDX_WEEKS_OBSERVED_REPEATED1], lambda_raw[IDX_WEEKS_OBSERVED_REPEATED2]);
  vector<lower=0>[W] theta = nu ./ (1 + nu);
  matrix[A,W] phi;
  matrix[B2,W] phi_reduced;
  matrix[A,W] alpha;
  matrix[B1,W1] alpha_reduced_1;
  matrix[B2,W2] alpha_reduced_2;
  matrix[W,num_basis] beta = to_matrix(beta_raw, W, num_basis); 

  for(w in 1:W)
  {
    
    phi[:,w] = softmax( to_vector(beta[w,:]*BASIS) ); 
    
    alpha[:,w] = phi[:,w] * lambda[w] / nu[w];
    
  }
  
  for(w in 1:W1){
    for(b in 1:B1){
      alpha_reduced_1[b,w] = sum(alpha[age_from_state_age_strata_1[b]:age_to_state_age_strata_1[b], w]);
    }
  }
    
  for(w in 1:W2){
    for(b in 1:B2){
      alpha_reduced_2[b,w] = sum(alpha[age_from_state_age_strata_2[b]:age_to_state_age_strata_2[b], W1 + w]);
    }
  }
  
  for(w in 1:W){
    for(b in 1:B2){
      phi_reduced[b,w] = sum(phi[age_from_state_age_strata_2[b]:age_to_state_age_strata_2[b], w]);
    }
  }

}


model {
  nu_raw ~ exponential(1);
  
  target += -0.5 * dot_self(beta_raw[node1] - beta_raw[node2]);
  // soft sum-to-zero constraint on phi)
  sum(beta_raw) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)

  for(w in 1:W_OBSERVED1){
    
    lambda_raw[w] ~ exponential( inv_sum_deaths[w]);

    target += neg_binomial_lpmf(deaths_1[idx_non_missing_1[1:N_idx_non_missing[w],w],w] | 
                                alpha_reduced_1[idx_non_missing_1[1:N_idx_non_missing[w],w], IDX_WEEKS_OBSERVED1[w]] , theta[IDX_WEEKS_OBSERVED1[w]] );
  
    if(N_idx_missing[w] > 0){
      for(n in 1:N_idx_missing[w])
        for(i in min_count_censored_1[idx_missing_1[n,w],w]:max_count_censored_1[idx_missing_1[n,w],w])
          target += neg_binomial_lpmf( i | alpha_reduced_1[idx_missing_1[n,w], IDX_WEEKS_OBSERVED1[w]] , theta[IDX_WEEKS_OBSERVED1[w]] ) ;
    }
        
  }

  for(w in 1:W_OBSERVED2){
    int w_cum = w + W_OBSERVED1;
    
    lambda_raw[w_cum] ~ exponential( inv_sum_deaths[w]);

    target += neg_binomial_lpmf(deaths_2[idx_non_missing_2[1:N_idx_non_missing[w_cum],w],w] | 
                                alpha_reduced_2[idx_non_missing_2[1:N_idx_non_missing[w_cum],w], IDX_WEEKS_OBSERVED2[w]] , theta[IDX_WEEKS_OBSERVED2[w] + W1] );
  
    if(N_idx_missing[w_cum] > 0){
      for(n in 1:N_idx_missing[w_cum])
        for(i in min_count_censored_2[idx_missing_2[n,w],w]:max_count_censored_2[idx_missing_2[n,w],w])
          target += neg_binomial_lpmf( i | alpha_reduced_2[idx_missing_2[n,w], IDX_WEEKS_OBSERVED2[w]] , theta[w_cum] ) ;
    }
        
  }
  
  }

generated quantities {
  real log_lik = 0;
  int deaths_predict[A,W];
  int deaths_predict_state_age_strata_1[B1,W1];
  int deaths_predict_state_age_strata_2[B2,W2];
  matrix[A,W] probability_ratio;
  matrix[B2,W] probability_ratio_age_strata;

  for(w in 1:W){

    // phi ratio
    probability_ratio[:,w] = phi[:,w] ./ phi[:,w_ref_index];
    probability_ratio_age_strata[:,w] = phi_reduced[:,w] ./ phi_reduced[:,w_ref_index];
    
    // predict deaths
    deaths_predict[:,w] = neg_binomial_rng(alpha[:,w], theta[w]);
    
  }
  
  for(w in 1:W1){
    deaths_predict_state_age_strata_1[:,w] = neg_binomial_rng(alpha_reduced_1[:,w], theta[w]);
  }
  
  for(w in 1:W2){
    int w_cum = w + W1;
    deaths_predict_state_age_strata_2[:,w] = neg_binomial_rng(alpha_reduced_2[:,w], theta[w_cum]);
  }
  
    for(w in 1:W_OBSERVED1){

    log_lik += neg_binomial_lpmf(deaths_1[idx_non_missing_1[1:N_idx_non_missing[w],w],w] | 
                                alpha_reduced_1[idx_non_missing_1[1:N_idx_non_missing[w],w], IDX_WEEKS_OBSERVED1[w]] , theta[IDX_WEEKS_OBSERVED1[w]] );
  
    if(N_idx_missing[w] > 0){
      for(n in 1:N_idx_missing[w])
        for(i in min_count_censored_1[idx_missing_1[n,w],w]:max_count_censored_1[idx_missing_1[n,w],w])
          log_lik += neg_binomial_lpmf( i | alpha_reduced_1[idx_missing_1[n,w], IDX_WEEKS_OBSERVED1[w]] , theta[IDX_WEEKS_OBSERVED1[w]] ) ;
    }
        
  }

  for(w in 1:W_OBSERVED2){
    int w_cum = w + W_OBSERVED1;

    log_lik += neg_binomial_lpmf(deaths_2[idx_non_missing_2[1:N_idx_non_missing[w_cum],w],w] | 
                                alpha_reduced_2[idx_non_missing_2[1:N_idx_non_missing[w_cum],w], IDX_WEEKS_OBSERVED2[w]] , theta[IDX_WEEKS_OBSERVED2[w] + W1] );
  
    if(N_idx_missing[w_cum] > 0){
      for(n in 1:N_idx_missing[w_cum])
        for(i in min_count_censored_2[idx_missing_2[n,w],w]:max_count_censored_2[idx_missing_2[n,w],w])
          log_lik += neg_binomial_lpmf( i | alpha_reduced_2[idx_missing_2[n,w], IDX_WEEKS_OBSERVED2[w]] , theta[w_cum] ) ;
    }
        
  }

}




