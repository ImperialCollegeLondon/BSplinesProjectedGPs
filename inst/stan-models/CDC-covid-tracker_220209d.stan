functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2, real inv_tau_squared){
    return - inv_tau_squared * 0.5 * dot_self(phi[node1] - phi[node2])
      + (N-1) * 0.5 * log(inv_tau_squared)
      + normal_lpdf(sum(phi) | 0, 0.001 * N);
  }
}

data{
  int<lower=1> W; // number of weeks 
  int<lower=1> M; // number of states 
  int<lower=0> A; // continuous age
  int<lower=0> B; // first age band specification
  int<lower=0,upper=W> W_OBSERVED; // number of weeks observed 
  int<lower=0,upper=W> W_NOT_OBSERVED; // number of weeks not observed 
  int<lower=1, upper=W> IDX_WEEKS_OBSERVED[W_OBSERVED]; // index of the weeks observed 
  int<lower=1, upper=W> IDX_WEEKS_OBSERVED_REPEATED[W]; // index of the weeks observed where missing is equal to the previous one 
  int<lower=0,upper=B> N_idx_non_missing[M,W_OBSERVED];
  int<lower=-1,upper=B> idx_non_missing[M,B,W_OBSERVED]; // indices non-missing deaths for W
  real age[A]; // age continuous
  real inv_sum_deaths[M,W_OBSERVED]; // inverse sum of deaths
  matrix[2,W_OBSERVED] lambda_prior_parameters[M]; // parameters of the prior distribution of lambda
  int deaths[M,B,W_OBSERVED]; // daily deaths in age band b at time t
  int age_from_state_age_strata[B]; // age from of age band b
  int age_to_state_age_strata[B];// age to of age band b
  
  // missing death count
  int<lower=0> N_missing[M]; // number of missing series 
  int<lower=-1, upper=1> start_or_end_period[M,max(N_missing)]; // is the serie cut by the end of the period
  int<lower=-1,upper=B> age_missing[M,max(N_missing)]; // age index with the missing death count in the serie
  // int<lower=1,upper=W> N_weeks_missing[M,max(N_missing)]; // numbers weeks missing in each serie
  int<lower=-1,upper=W> idx_weeks_missing_min[M,max(N_missing)]; // index of weeks missing in the series
  int<lower=-1,upper=W> idx_weeks_missing_max[M,max(N_missing)]; // index of weeks missing in the series
  int<lower=-1> sum_count_censored[M,max(N_missing)]; // sum of the missing deaths over the serie if it ends before the period
  int min_count_censored[M,max(N_missing)]; // range of the censored data if it ends after the period
  int max_count_censored[M,max(N_missing)]; // range of the censored data if it ends after the period

  //splines
  int num_basis_rows;
  int num_basis_columns;
  matrix[num_basis_rows, A] BASIS_ROWS; 
  matrix[num_basis_columns, W] BASIS_COLUMNS; 
  
  // ICAR model
  int<lower=1> K;
  int<lower=1> N_edges;
  int<lower=1, upper=K> node1[N_edges];
  int<lower=1, upper=K> node2[N_edges];
}

transformed data
{   
    int N_log_lik = 0;
    
  for(m in 1:M){
    for(w in 1:W_OBSERVED){
      for(i in idx_non_missing[m][1:N_idx_non_missing[m][w],w]){
        N_log_lik += 1;
      }
    }
 
   for(n in 1:N_missing[m]){
     if(!start_or_end_period[m][n])
     {
        N_log_lik += 1;
 
     } else {
        for(i in min_count_censored[m][n]:max_count_censored[m][n])
           N_log_lik += 1;
     }
    }
  }

}

parameters {
  real<lower=0> nu_inverse[M];
  vector<lower=0>[W-W_NOT_OBSERVED] lambda_raw[M];
  vector[K] beta_raw[M]; 
  real<lower=0> tau[M];
}

transformed parameters {
  real<lower=0> inv_tau_squared[M];
  vector<lower=0>[W] lambda[M];
  real<lower=0> nu[M];
  matrix[A,W] phi[M];
  matrix[A,W] alpha[M];
  matrix[B,W] phi_reduced[M];
  matrix[B,W] alpha_reduced[M];
  matrix[num_basis_rows,num_basis_columns] beta[M]; 
  matrix[A, W] f[M];

  for(m in 1:M){
    lambda[m] = lambda_raw[m][IDX_WEEKS_OBSERVED_REPEATED];
    nu[m] = (1/nu_inverse[m]);

    inv_tau_squared[m] = 1/(tau[m]^2);

    beta[m] = to_matrix(beta_raw[m], num_basis_rows,num_basis_columns); 

    f[m] = (BASIS_ROWS') * beta[m] * BASIS_COLUMNS;
    
    for(w in 1:W){
      phi[m][:,w] = softmax( f[m][:,w] ); 
      alpha[m][:,w] = phi[m][:,w] * lambda[m][w] / nu[m] ;
      
      for(b in 1:B){
        alpha_reduced[m][b,w] = sum(alpha[m][age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
        phi_reduced[m][b,w] = sum(phi[m][age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
      }
    }
  }
}

model {
  
  nu_inverse ~ normal(0,5);

  tau ~ cauchy(0,1);

  for(m in 1:M){
    
    beta_raw[m] ~ icar_normal_lpdf(K, node1, node2, inv_tau_squared[m]);
    lambda_raw[m] ~ gamma( lambda_prior_parameters[m][1,:],lambda_prior_parameters[m][2,:]);

    // Note on the neg bin parametrisation related to the paper:
    // mean neg_binomial_lpmf is alpha_reduced / nu_inverse = alpha_reduced * nu
    // var neg_binomial_lpmf is alpha_reduced / nu_inverse^2 * (nu_inverse + 1) = alpha_reduced * nu (1 + nu)
    
    for(w in 1:W_OBSERVED){
      
      int indx_w_m[N_idx_non_missing[m][w]] = idx_non_missing[m][1:N_idx_non_missing[m][w],w];
      target += neg_binomial_lpmf(deaths[m][indx_w_m,w] | alpha_reduced[m][indx_w_m, IDX_WEEKS_OBSERVED[w]] , nu_inverse[m] );
        
    }

    for(n in 1:N_missing[m]){
      real alpha_reduced_missing = sum(alpha_reduced[m][ age_missing[m][n], idx_weeks_missing_min[m][n]:idx_weeks_missing_max[m][n] ]);

      if(!start_or_end_period[m][n])
      {
        target += neg_binomial_lpmf( sum_count_censored[m][n] | alpha_reduced_missing , nu_inverse[m] ) ;
      } 
      else {
       for(i in min_count_censored[m][n]:max_count_censored[m][n])
          target += neg_binomial_lpmf( i | alpha_reduced_missing , nu_inverse[m] ) ;
      }
    } 
  }

}

generated quantities {
  real log_lik[N_log_lik];
  int deaths_predict[M,A,W];
  int deaths_predict_state_age_strata[M,B,W];
  
    // predction and log likelihood for the weekly age model
  {
    int idx_log_lik = 0;
        
    for(m in 1:M){


        for(w in 1:W){
            // predict deaths
            deaths_predict[m,:,w] = neg_binomial_rng(alpha[m,:,w], nu_inverse[m] );
            deaths_predict_state_age_strata[m,:,w] = neg_binomial_rng(alpha_reduced[m,:,w], nu_inverse[m] );
        }

        for(w in 1:W_OBSERVED){
            for(i in idx_non_missing[m][1:N_idx_non_missing[m][w],w]){
                idx_log_lik += 1;
                log_lik[idx_log_lik] = neg_binomial_lpmf(deaths[m,i,w] | alpha_reduced[m,i, IDX_WEEKS_OBSERVED[w]] , nu_inverse[m] );
            }
        }

        for(n in 1:N_missing[m]){
            real alpha_reduced_missing = sum(alpha_reduced[m][ age_missing[m][n], idx_weeks_missing_min[m][n]:idx_weeks_missing_max[m][n] ]);

            if(!start_or_end_period[m][n])
                {
                idx_log_lik += 1;
                log_lik[idx_log_lik] = neg_binomial_lpmf( sum_count_censored[m][n] | alpha_reduced_missing , nu_inverse[m] ) ;
            } else {
                for(i in min_count_censored[m][n]:max_count_censored[m][n]){
                    idx_log_lik += 1;
                    log_lik[idx_log_lik] = neg_binomial_lpmf( i | alpha_reduced_missing , nu_inverse[m] ) ;
                }
            }
        }
    }
  }
    
}


