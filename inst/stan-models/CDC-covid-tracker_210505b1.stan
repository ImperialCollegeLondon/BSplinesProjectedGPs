functions {
    matrix kron_mvprod(matrix A, matrix B, matrix V) 
    {
        return transpose(A*transpose(B*V));
    }
    
  matrix gp(int T, int D, real[] time, real[] delay,
            real delta0,
            real alpha_gp1_t, real alpha_gp1_d, 
            vector rho_gp1_t_dist, vector rho_gp1_d_dist,
            matrix z1)
  {
    
    matrix[T,D] GP1;//long range
    
    matrix[T, T] K1_t;
    matrix[T, T] L_K1_t;
    
    matrix[D, D] K1_d;
    matrix[D, D] L_K1_d;
    
    real sq_alpha1_t = square(alpha_gp1_t);
    real sq_alpha1_d = square(alpha_gp1_d);
    
    real K1diag_t = sq_alpha1_t + delta0;
    real K1diag_d = sq_alpha1_d + delta0;
    
 //time   
    for (i in 1:T) 
        {
            K1_t[i, i] = K1diag_t;
            for (j in (i + 1):T) 
                {
                    K1_t[i, j] = sq_alpha1_t
                            * exp(-0.5 * dot_self((time[i] - time[j]) ./ rho_gp1_t_dist)) ; //long range sq exp                         
                    K1_t[j, i] = K1_t[i, j];
                }
        }
    K1_t[T, T] = K1diag_t;
    
 //delay   
    for (i in 1:D) 
        {
            K1_d[i, i] = K1diag_d;
            for (j in (i + 1):D) 
                {
                    K1_d[i, j] = sq_alpha1_d
                            * exp(-0.5 * dot_self((delay[i] - delay[j]) ./ rho_gp1_d_dist)) ; //long range sq exp
                    K1_d[j, i] = K1_d[i, j];
                }
        }
    K1_d[D, D] = K1diag_d;

    L_K1_t = cholesky_decompose(K1_t);
    L_K1_d = cholesky_decompose(K1_d);
    
    GP1 = kron_mvprod(L_K1_d, L_K1_t, z1);

    return(GP1);
  }
}

data{
  int<lower=0> W; // number of weeks 
  int<lower=0,upper=W> W_OBSERVED; // number of weeks observed 
  int<lower=0,upper=W> W_NOT_OBSERVED; // number of weeks not observed 
  int<lower=1, upper=W> IDX_WEEKS_OBSERVED[W_OBSERVED]; // index of the weeks observed 
  int<lower=1, upper=W> IDX_WEEKS_OBSERVED_REPEATED[W]; // index of the weeks observed where missing is equal to the previous one 
  int<lower=0,upper=W> W_ref_index; // number of index to compare the death prob
  int<lower=0,upper=W> w_ref_index[W_ref_index]; // week index to compare the death prob
  int<lower=0> A; // continuous age
  int<lower=0> B; // first age band specification
  int<lower=0,upper=B> N_idx_non_missing[W_OBSERVED];
  int<lower=-1,upper=B> idx_non_missing[B,W_OBSERVED]; // indices non-missing deaths for W
  real age[A]; // age continuous
  real inv_sum_deaths[W_OBSERVED]; // inverse sum of deaths
  matrix[2,W_OBSERVED] lambda_prior_parameters; // parameters of the prior distribution of lambda
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
  
  // GP
  real IDX_WEEKS[W];
  real IDX_BASIS[num_basis];
}

transformed data
{   
    real delta0 = 1e-9;  
    int N_log_lik = 0;
    
  for(w in 1:W_OBSERVED){
    for(i in idx_non_missing[1:N_idx_non_missing[w],w]){
      N_log_lik += 1;
    }
  }
  for(n in 1:N_missing){
    if(!start_or_end_period[n])
    {
       N_log_lik += 1;

    } else {
       for(i in min_count_censored[n]:max_count_censored[n])
          N_log_lik += 1;
    }
  }
}

parameters {
  real<lower=0> nu;
  vector<lower=0>[W-W_NOT_OBSERVED] lambda_raw;
  real<lower=0> alpha_gp1_t;
  real<lower=0> alpha_gp1_d;
  matrix[W,num_basis] z1;
  vector<lower=0>[1] rho_gp1_t_dist; 
  vector<lower=0>[1] rho_gp1_d_dist;
}

transformed parameters {
  vector<lower=0>[W] lambda = lambda_raw[IDX_WEEKS_OBSERVED_REPEATED];
  real<lower=0> theta = nu / (1 + nu);
  matrix[A,W] phi;
  matrix[A,W] alpha;
  matrix[B,W] phi_reduced;
  matrix[B,W] alpha_reduced;
  vector[N_missing] alpha_reduced_missing;
  matrix[W,num_basis] beta = gp(W, num_basis, IDX_WEEKS, IDX_BASIS,
                              delta0,
                              alpha_gp1_t, alpha_gp1_d, 
                              rho_gp1_t_dist,  rho_gp1_d_dist,
                              z1); 

  for(w in 1:W)
  {
    
    phi[:,w] = softmax( to_vector( beta[w,:]*BASIS ) ); 
    
    alpha[:,w] = phi[:,w] * lambda[w] / nu ;
    
  }
  
  for(w in 1:W){
    for(b in 1:B){
      alpha_reduced[b,w] = sum(alpha[age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
      phi_reduced[b,w] = sum(phi[age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
    }
  }
  
  for(n in 1:N_missing){
    alpha_reduced_missing[n] = sum(alpha_reduced[age_missing[n],  idx_weeks_missing[1:N_weeks_missing[n], n] ]);
  }

}

model {
  
  nu ~ exponential(1);
  lambda_raw ~ gamma( lambda_prior_parameters[1,:],lambda_prior_parameters[2,:]);
  
  alpha_gp1_t ~ normal(0,1);
  alpha_gp1_d ~ normal(0,1);
  
  rho_gp1_t_dist ~ inv_gamma(5, 5);
  rho_gp1_d_dist ~ inv_gamma(5, 5);

    for(t in 1:W)
    {
        for(d in 1:num_basis)
            {
                {
                   z1[t,d] ~ normal(0,0.1);
                }
            }
    }
    
  for(w in 1:W_OBSERVED){
  

    target += neg_binomial_lpmf(deaths[idx_non_missing[1:N_idx_non_missing[w],w],w] |
                                alpha_reduced[idx_non_missing[1:N_idx_non_missing[w],w], IDX_WEEKS_OBSERVED[w]] , theta );
  
        
  }
  
  for(n in 1:N_missing){
    if(!start_or_end_period[n])
    {
       
      target += neg_binomial_lpmf( sum_count_censored[n] | alpha_reduced_missing[n], theta ) ;
    } 
    else {
       for(i in min_count_censored[n]:max_count_censored[n])
          target += neg_binomial_lpmf( i | alpha_reduced_missing[n] , theta ) ;
    }
  }

}

generated quantities {
  real log_lik[N_log_lik];
  int deaths_predict[A,W];
  int deaths_predict_state_age_strata[B,W];
  matrix[A,W] probability_ratio;
  matrix[B,W] probability_ratio_age_strata;

  for(w in 1:W){

    // phi ratio
    probability_ratio[:,w] = phi[:,w] ./ (phi[:,w_ref_index] * rep_vector(1.0 / W_ref_index, W_ref_index));
    probability_ratio_age_strata[:,w] = phi_reduced[:,w] ./ (phi_reduced[:,w_ref_index] * rep_vector(1.0 / W_ref_index, W_ref_index));
    
    // predict deaths
    deaths_predict[:,w] = neg_binomial_rng(alpha[:,w], theta );
    deaths_predict_state_age_strata[:,w] = neg_binomial_rng(alpha_reduced[:,w], theta );
  }
  
  
  {
    int idx_log_lik = 0;
    for(w in 1:W_OBSERVED){
      for(i in idx_non_missing[1:N_idx_non_missing[w],w]){
        idx_log_lik += 1; 
        log_lik[idx_log_lik] = neg_binomial_lpmf(deaths[i,w] | alpha_reduced[i, IDX_WEEKS_OBSERVED[w]] , theta );
      }
    }
    for(n in 1:N_missing){
      if(!start_or_end_period[n])
      {
      idx_log_lik += 1; 
       log_lik[idx_log_lik] = neg_binomial_lpmf( sum_count_censored[n] |  alpha_reduced_missing[n] , theta ) ;
      } else {
       for(i in min_count_censored[n]:max_count_censored[n]){
          idx_log_lik += 1; 
          log_lik[idx_log_lik] = neg_binomial_lpmf( i |  alpha_reduced_missing[n], theta ) ;
       }
      }
    }
  }

}



