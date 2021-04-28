data{
  int<lower=0> W; // number of weeks 
  int<lower=0,upper=W> W_OBSERVED; // number of weeks observed 
  int<lower=1, upper=W> IDX_WEEKS_OBSERVED[W_OBSERVED]; // index of the weeks observed 
  int<lower=1, upper=W> IDX_WEEKS_OBSERVED_REPEATED[W]; // index of the weeks observed where missing is equal to the previous one 
  int<lower=0,upper=W> w_ref_index; // week index to compare the death prob
  int<lower=0> A; // continuous age
  int<lower=0> B; // first age band specification
  int<lower=0,upper=B> N_idx_non_missing[W_OBSERVED];
  int<lower=-1,upper=B> idx_non_missing[B,W_OBSERVED]; // indices non-missing deaths for W
  vector[A] age; // age continuous
  real inv_sum_deaths[W_OBSERVED]; // inverse sum of deaths
  int deaths[B,W_OBSERVED]; // daily deaths in age band b at time t
  int age_from_state_age_strata[B]; // age from of age band b
  int age_to_state_age_strata[B];// age to of age band b
  
  // missing death count
  int<lower=0> N_missing; // number of missing series 
  int<lower=0, upper=1> start_or_end_period[N_missing]; // is the serie cut by the end of the period
  int<lower=1,upper=B> age_missing[N_missing]; // age index with the missing death count in the serie
  int<lower=1,upper=W> N_weeks_missing[N_missing]; // numbers weeks missing in each serie
  int<lower=-1,upper=W> idx_weeks_missing[max(N_weeks_missing),N_missing]; // index of weeks missing in the series
  int<lower=-1> sum_count_censored[N_missing]; // sum of the missing deaths over the serie if it ends before the period
  int min_count_censored[N_missing]; // range of the censored data if it ends after the period
  int max_count_censored[N_missing]; // range of the censored data if it ends after the period
  
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
  real<lower=0> nu;
  vector<lower=0>[W-1] lambda_raw;
}

transformed parameters {
  vector<lower=0>[W] lambda = lambda_raw[IDX_WEEKS_OBSERVED_REPEATED];
  real<lower=0> theta = nu / (1 + nu);
  matrix[A,W] phi;
  matrix[A,W] alpha;
  matrix[B,W] phi_reduced;
  matrix[B,W] alpha_reduced;
  matrix[W,num_basis] beta = to_matrix(beta_raw, W, num_basis); 

  for(w in 1:W)
  {
    
    phi[:,w] = softmax( to_vector(beta[w,:]*BASIS) ); 
    
    alpha[:,w] = phi[:,w] * lambda[w] / nu;
    
  }
  
  for(w in 1:W){
    for(b in 1:B){
      alpha_reduced[b,w] = sum(alpha[age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
      phi_reduced[b,w] = sum(phi[age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
    }
  }

}

model {
  
  nu ~ exponential(1);
  
  target += -0.5 * dot_self(beta_raw[node1] - beta_raw[node2]);
  // soft sum-to-zero constraint on phi)
  sum(beta_raw) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)

  for(w in 1:W_OBSERVED){
    
    lambda_raw[w] ~ exponential( inv_sum_deaths[w] );

    target += neg_binomial_lpmf(deaths[idx_non_missing[1:N_idx_non_missing[w],w],w] | 
                                alpha_reduced[idx_non_missing[1:N_idx_non_missing[w],w], IDX_WEEKS_OBSERVED[w]] , theta );
  
        
  }
  
  for(n in 1:N_missing){
    if(!start_or_end_period[n])
    {
       target += neg_binomial_lpmf( sum_count_censored[n] | sum(alpha_reduced[age_missing[n], idx_weeks_missing[1:N_weeks_missing[n], n] ]) , theta ) ;
       

    } else {
       for(i in min_count_censored[n]:max_count_censored[n])
          target += neg_binomial_lpmf( i | sum(alpha_reduced[age_missing[n], idx_weeks_missing[1:N_weeks_missing[n], n]] ) , theta ) ;
    }
  }
      
}

generated quantities {
  real log_lik = 0;
  int deaths_predict[A,W];
  int deaths_predict_state_age_strata[B,W];
  matrix[A,W] probability_ratio;
  matrix[B,W] probability_ratio_age_strata;

  for(w in 1:W){

    // phi ratio
    probability_ratio[:,w] = phi[:,w] ./ phi[:,w_ref_index];
    probability_ratio_age_strata[:,w] = phi_reduced[:,w] ./ phi_reduced[:,w_ref_index];
    
    // predict deaths
    deaths_predict[:,w] = neg_binomial_rng(alpha[:,w], theta );
    deaths_predict_state_age_strata[:,w] = neg_binomial_rng(alpha_reduced[:,w], theta );
  }
  
  
  for(w in 1:W_OBSERVED){

    log_lik += neg_binomial_lpmf(deaths[idx_non_missing[1:N_idx_non_missing[w],w],w] | 
                                alpha_reduced[idx_non_missing[1:N_idx_non_missing[w],w], IDX_WEEKS_OBSERVED[w]] , theta );
  
        
  }
  
  for(n in 1:N_missing){
    if(!start_or_end_period[n])
    {
       log_lik += neg_binomial_lpmf( sum_count_censored[n] | sum(alpha_reduced[age_missing[n], idx_weeks_missing[1:N_weeks_missing[n], n] ]) , theta ) ;

    } else {
       for(i in min_count_censored[n]:max_count_censored[n])
          log_lik += neg_binomial_lpmf( i | sum(alpha_reduced[age_missing[n], idx_weeks_missing[1:N_weeks_missing[n], n] ]) , theta ) ;
    }
  }

}



