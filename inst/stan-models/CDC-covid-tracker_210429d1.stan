functions {
  vector row_sums(int[,] X) 
  {
      int n_rows = dims(X)[1]; // this will give number of rows
      vector[n_rows] s;
      matrix [dims(X)[1], dims(X)[2]] mat_X;
      mat_X = to_matrix(X);
      for (i in 1:n_rows) s[i] = sum(row(mat_X, i));
      return s;
  }
  

    matrix kron_mvprod(matrix A, matrix B, matrix V) 
    {
        return transpose(A*transpose(B*V));
    }
  
    matrix calculate_eigenvalues(vector A, vector B, int n1, int n2, real sigma2) 
    {
        matrix[n1,n2] e;
            for(i in 1:n1) 
            {
                for(j in 1:n2) 
                {
                    e[i,j] = (A[i]*B[j] + sigma2);
                }
            }
        return(e);
    }
    
  matrix gp(int T, int D, real[] time, real[] delay,
            real delta0, real delta1, real delta2, 
            real alpha_gp1_t, real alpha_gp2_t, real alpha_gp1_d, real alpha_gp2_d,
            vector rho_gp1_t_dist, vector rho_gp2_t_dist, vector rho_gp1_d_dist, vector rho_gp2_d_dist,
            matrix z1, matrix z2)
  {
    matrix[T,D] lambda;
    matrix[T,D] lambda_gp1;
    matrix[T,D] lambda_gp2;
    
    matrix[T,D] GP1;//long range
    matrix[T,D] GP2;//short range
    
    matrix[T, T] K1_t;
    matrix[T, T] K2_t;
    matrix[T, T] L_K1_t;
    matrix[T, T] L_K2_t;
    
    matrix[D, D] K1_d;
    matrix[D, D] K2_d;
    matrix[D, D] L_K1_d;
    matrix[D, D] L_K2_d;
    
    real sq_alpha1_t = square(alpha_gp1_t);
    real sq_alpha2_t = square(alpha_gp2_t);
    real sq_alpha1_d = square(alpha_gp1_d);
    real sq_alpha2_d = square(alpha_gp2_d);
    
    real K1diag_t = sq_alpha1_t + delta0 + delta1;
    real K2diag_t = sq_alpha2_t + delta0 ;
    
    real K1diag_d = sq_alpha1_d + delta0;
    real K2diag_d = sq_alpha2_d + delta0 + delta2;
    
 //time   
    for (i in 1:T) 
        {
            K1_t[i, i] = K1diag_t;
            K2_t[i, i] = K2diag_t;
            for (j in (i + 1):T) 
                {
                    K1_t[i, j] = sq_alpha1_t
                            * exp(-0.5 * dot_self((time[i] - time[j]) ./ rho_gp1_t_dist)) ; //long range sq exp                         
                    K2_t[i, j] = sq_alpha2_t
                            * exp( -0.5 * dot_self((time[i] - time[j]) ./ rho_gp2_t_dist));// shortrange  sq exp
                    K1_t[j, i] = K1_t[i, j];
                    K2_t[j, i] = K2_t[i, j];
                }
        }
    K1_t[T, T] = K1diag_t;
    K2_t[T, T] = K2diag_t;
    
 //delay   
    for (i in 1:D) 
        {
            K1_d[i, i] = K1diag_d;
            K2_d[i, i] = K2diag_d;
            for (j in (i + 1):D) 
                {
                    K1_d[i, j] = sq_alpha1_d
                            * exp(-0.5 * dot_self((delay[i] - delay[j]) ./ rho_gp1_d_dist)) ; //long range sq exp
                    K2_d[i, j] = sq_alpha2_d
                            * exp( -0.5 * dot_self((delay[i] - delay[j]) ./ rho_gp2_d_dist));// shortrange  sq exp
                    K1_d[j, i] = K1_d[i, j];
                    K2_d[j, i] = K2_d[i, j];
                }
        }
    K1_d[D, D] = K1diag_d;
    K2_d[D, D] = K2diag_d;

    L_K1_t = cholesky_decompose(K1_t);
    L_K1_d = cholesky_decompose(K1_d);
    L_K2_t = cholesky_decompose(K2_t);
    L_K2_d = cholesky_decompose(K2_d);
    
    GP1 = kron_mvprod(L_K1_d, L_K1_t, z1);
    GP2 = kron_mvprod(L_K2_d, L_K2_t, z2);

    for (t in 1:T)
        {
            for(d in 1:D)
                {
                    lambda_gp1[t,d] = exp(GP1[t,d]);
                    lambda_gp2[t,d] = exp(GP2[t,d]);
                    lambda[t,d] = lambda_gp1[t,d] * lambda_gp2[t,d];
                }
        }
    return(lambda);
  }
}

data{
  int<lower=0> W; // number of weeks 
  int<lower=0,upper=W> W_OBSERVED; // number of weeks observed 
  int<lower=0,upper=W> W_NOT_OBSERVED; // number of weeks not observed 
  int<lower=1, upper=W> IDX_WEEKS_OBSERVED[W_OBSERVED]; // index of the weeks observed 
  int<lower=1, upper=W> IDX_WEEKS_OBSERVED_REPEATED[W]; // index of the weeks observed where missing is equal to the previous one 
  int<lower=0,upper=W> w_ref_index; // week index to compare the death prob
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
    
    vector[1] rho_gp1_t_dist; 
    vector[1] rho_gp2_t_dist;
    vector[1] rho_gp1_d_dist;
    vector[1] rho_gp2_d_dist;

    rho_gp1_t_dist[1] = W;
    rho_gp2_t_dist[1] = 1;
    rho_gp1_d_dist[1] = num_basis;
    rho_gp2_d_dist[1] = 1;
}

parameters {
  real<lower=0> nu;
  vector<lower=0>[W-W_NOT_OBSERVED] lambda_raw;
  real<lower=0> alpha_gp1_t;
  real<lower=0> alpha_gp2_t;
  real<lower=0> alpha_gp1_d;
  real<lower=0> alpha_gp2_d;
  real<lower=0> delta1;
  real<lower=0> delta2;
  matrix[W,num_basis] z1;
  matrix[W,num_basis] z2;
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
                              delta0, delta1, delta2, 
                              alpha_gp1_t, alpha_gp2_t, alpha_gp1_d, alpha_gp2_d,
                              rho_gp1_t_dist, rho_gp2_t_dist, rho_gp1_d_dist, rho_gp2_d_dist,
                              z1, z2); 

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
  
  delta1 ~ normal(0,1e-7);
  delta2 ~ normal(0,1e-7);
  
  alpha_gp1_t ~ normal(W,5);
  alpha_gp1_d ~ normal(0,1);
  
  alpha_gp2_t ~ normal(num_basis,5);
  alpha_gp2_d ~ normal(0,1);

    for(t in 1:W)
    {
        for(d in 1:num_basis)
            {
                {
                   z1[t,d] ~ normal(0,0.1);
                   z2[t,d] ~ normal(0,0.1);
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
  real log_lik = 0;
  int deaths_predict[A,W];
  int deaths_predict_state_age_strata[B,W];
  matrix[A,W] probability_ratio;
  matrix[B,W] probability_ratio_age_strata;

  for(w in 1:W){

    // phi ratio
    probability_ratio[:,w] = phi[:,w] ./ (phi[:,1:w_ref_index] * rep_vector(1.0 / w_ref_index, w_ref_index));
    probability_ratio_age_strata[:,w] = phi_reduced[:,w] ./ (phi_reduced[:,1:w_ref_index] * rep_vector(1.0 / w_ref_index, w_ref_index));
    
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
       log_lik += neg_binomial_lpmf( sum_count_censored[n] |  alpha_reduced_missing[n] , theta ) ;

    } else {
       for(i in min_count_censored[n]:max_count_censored[n])
          log_lik += neg_binomial_lpmf( i |  alpha_reduced_missing[n], theta ) ;
    }
  }

}




