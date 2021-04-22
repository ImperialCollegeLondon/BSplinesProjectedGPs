functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}

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
  int age_from_state_age_strata_2[B2]; // age from of age band b
  int age_to_state_age_strata_2[B2];// age to of age band b
  
  // missing death count
  int<lower=0> N_missing_1; // number of missing series 
  int<lower=0, upper=1> start_or_end_period_1[N_missing_1]; // is the serie cut by the end of the period
  int<lower=1,upper=B1> age_missing_1[N_missing_1]; // ageindex with the missing death count in the serie
  int<lower=1,upper=W1> N_weeks_missing_1[N_missing_1]; // numbers weeks missing in each serie
  int<lower=-1,upper=W1> idx_weeks_missing_1[max(N_weeks_missing_1),N_missing_1]; // number of weeks missing in the series
  int<lower=-1> sum_missing_deaths_1[N_missing_1]; // summ of the missing deaths over the serie if it ends before the period
  int min_count_censored_1[N_missing_1]; // range of the censored data if it ends after the period
  int max_count_censored_1[N_missing_1]; // range of the censored data if it ends after the period
  int<lower=0> N_missing_2; // number of missing series 
  int<lower=0, upper=1> start_or_end_period_2[N_missing_2]; // is the serie cut by the end of the period
  int<lower=1,upper=B2> age_missing_2[N_missing_2]; // ageindex with the missing death count in the serie
  int<lower=1,upper=W2> N_weeks_missing_2[N_missing_2]; // numbers weeks missing in each serie
  int<lower=-1,upper=W2> idx_weeks_missing_2[max(N_weeks_missing_2),N_missing_2]; // number of weeks missing in the series
  int<lower=-1> sum_missing_deaths_2[N_missing_2]; // summ of the missing deaths over the serie if it ends before the period
  int min_count_censored_2[N_missing_2]; // range of the censored data if it ends after the period
  int max_count_censored_2[N_missing_2]; // range of the censored data if it ends after the period
  
  //splines
  int num_basis;
  matrix[num_basis, A] BASIS; 
  
  // CAR model
  int N; // W * num_basis
  matrix<lower = 0, upper = 1>[N, N] Adj; // adjacency matrix
  int Adj_n;                // number of adjacent region pairs
}

transformed data {
  int Adj_sparse[Adj_n, 2];   // adjacency pairs
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[N] egv;       // eigenvalues of invsqrtD * Adj * invsqrtD
  
  { // generate sparse representation for Ajd
  int counter;
  counter = 1;
  // loop over upper triangular part of Adj to identify neighbor pairs
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (Adj[i, j] == 1) {
          Adj_sparse[counter, 1] = i;
          Adj_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:N) D_sparse[i] = sum(Adj[i]);
  {
    vector[N] invsqrtD;  
    for (i in 1:N) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    egv = eigenvalues_sym(quad_form(Adj, diag_matrix(invsqrtD)));
  }
}

parameters {
  vector[N] beta_raw; 
  real<lower=0> nu;
  vector<lower=0>[W-1] lambda_raw;
  real<lower = 0, upper = 1> p;
  real<lower = 0> inv_tau;
}

transformed parameters {
  real<lower = 0> tau = 1.0 / inv_tau;
  vector<lower=0>[W] lambda = append_row(lambda_raw[IDX_WEEKS_OBSERVED_REPEATED1], lambda_raw[IDX_WEEKS_OBSERVED_REPEATED2]);
  real<lower=0> theta = nu / (1 + nu);
  matrix[A,W] phi;
  matrix[B2,W] phi_reduced;
  matrix[A,W] alpha;
  matrix[B1,W1] alpha_reduced_1;
  matrix[B2,W2] alpha_reduced_2;
  matrix[W,num_basis] beta = to_matrix(beta_raw, W, num_basis); 

  for(w in 1:W)
  {
    
    phi[:,w] = softmax( to_vector(beta[w,:]*BASIS) ); 
    
    alpha[:,w] = phi[:,w] * lambda[w] / nu;
    
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
  nu ~ exponential(1);
  
  inv_tau ~ exponential(50);
  
  beta_raw ~ sparse_car(tau, p, Adj_sparse, D_sparse, egv, N, Adj_n);

  for(w in 1:W_OBSERVED1){
    
    lambda_raw[w] ~ exponential( inv_sum_deaths[w]);

    target += neg_binomial_lpmf(deaths_1[idx_non_missing_1[1:N_idx_non_missing[w],w],w] | 
                                alpha_reduced_1[idx_non_missing_1[1:N_idx_non_missing[w],w], IDX_WEEKS_OBSERVED1[w]] , theta );
  
        
  }

  for(w in 1:W_OBSERVED2){
    int w_cum = w + W_OBSERVED1;
    
    lambda_raw[w_cum] ~ exponential( inv_sum_deaths[w]);

    target += neg_binomial_lpmf(deaths_2[idx_non_missing_2[1:N_idx_non_missing[w_cum],w],w] | 
                                alpha_reduced_2[idx_non_missing_2[1:N_idx_non_missing[w_cum],w], IDX_WEEKS_OBSERVED2[w]] , theta );
        
  }
  
  for(n in 1:N_missing_1){
    if(!start_or_end_period_1[n])
    {
       target += neg_binomial_lpmf( sum_missing_deaths_1[n] | sum(alpha_reduced_1[age_missing_1[n],idx_weeks_missing_1[1:N_weeks_missing_1[n], n]]) , theta ) ;

    } else {
       for(i in min_count_censored_1[n]:max_count_censored_1[n])
          target += neg_binomial_lpmf( i | sum(alpha_reduced_1[age_missing_1[n],idx_weeks_missing_1[1:N_weeks_missing_1[n], n]]) , theta ) ;
    }
  }
      
  
  for(n in 1:N_missing_2){
    
    if(!start_or_end_period_2[n])
    {
       target += neg_binomial_lpmf( sum_missing_deaths_2[n] | sum(alpha_reduced_2[age_missing_2[n],idx_weeks_missing_2[1:N_weeks_missing_2[n], n]]) , theta ) ;

    } else {
       for(i in min_count_censored_2[n]:max_count_censored_2[n])
          target += neg_binomial_lpmf( i | sum(alpha_reduced_2[age_missing_2[n],idx_weeks_missing_2[1:N_weeks_missing_2[n], n]]) , theta ) ;
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
    deaths_predict[:,w] = neg_binomial_rng(alpha[:,w], theta );
    
  }
  
  for(w in 1:W1){
    deaths_predict_state_age_strata_1[:,w] = neg_binomial_rng(alpha_reduced_1[:,w], theta );
  }
  
  for(w in 1:W2){
    deaths_predict_state_age_strata_2[:,w] = neg_binomial_rng(alpha_reduced_2[:,w], theta );
  }
  
    for(w in 1:W_OBSERVED1){

    log_lik += neg_binomial_lpmf(deaths_1[idx_non_missing_1[1:N_idx_non_missing[w],w],w] | 
                                alpha_reduced_1[idx_non_missing_1[1:N_idx_non_missing[w],w], IDX_WEEKS_OBSERVED1[w]] , theta );
        
  }

  for(w in 1:W_OBSERVED2){
    int w_cum = w + W_OBSERVED1;

    log_lik += neg_binomial_lpmf(deaths_2[idx_non_missing_2[1:N_idx_non_missing[w_cum],w],w] | 
                                alpha_reduced_2[idx_non_missing_2[1:N_idx_non_missing[w_cum],w], IDX_WEEKS_OBSERVED2[w]] , theta );
        
  }
  
    for(n in 1:N_missing_1){
    if(!start_or_end_period_1[n])
    {
       log_lik += neg_binomial_lpmf( sum_missing_deaths_1[n] | sum(alpha_reduced_1[age_missing_1[n],idx_weeks_missing_1[1:N_weeks_missing_1[n], n]]) , theta ) ;

    } else {
       for(i in min_count_censored_1[n]:max_count_censored_1[n])
          log_lik += neg_binomial_lpmf( i | sum(alpha_reduced_1[age_missing_1[n],idx_weeks_missing_1[1:N_weeks_missing_1[n], n]]) , theta ) ;
    }
  }
      
  
  for(n in 1:N_missing_2){
    
    if(!start_or_end_period_2[n])
    {
       log_lik += neg_binomial_lpmf( sum_missing_deaths_2[n] | sum(alpha_reduced_2[age_missing_2[n],idx_weeks_missing_2[1:N_weeks_missing_2[n], n]]) , theta ) ;

    } else {
       for(i in min_count_censored_2[n]:max_count_censored_2[n])
          log_lik += neg_binomial_lpmf( i | sum(alpha_reduced_2[age_missing_2[n],idx_weeks_missing_2[1:N_weeks_missing_2[n], n]]) , theta ) ;
    }
  }

}





