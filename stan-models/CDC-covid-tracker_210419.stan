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
  real IDX_WEEKS[W];
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
  real IDX_BASIS[num_basis];
  matrix[num_basis, A] BASIS; 
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
  vector<lower=0>[W-1] nu_raw;
  vector<lower=0>[W-1] lambda_raw;
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
  vector<lower=0>[W] nu = append_row(nu_raw[IDX_WEEKS_OBSERVED_REPEATED1], nu_raw[IDX_WEEKS_OBSERVED_REPEATED2]);
  vector<lower=0>[W] lambda = append_row(lambda_raw[IDX_WEEKS_OBSERVED_REPEATED1], lambda_raw[IDX_WEEKS_OBSERVED_REPEATED2]);
  vector<lower=0>[W] theta = nu ./ (1 + nu);
  matrix[A,W] phi;
  matrix[B2,W] phi_reduced;
  matrix[A,W] alpha;
  matrix[B1,W1] alpha_reduced_1;
  matrix[B2,W2] alpha_reduced_2;
  matrix[W,num_basis] beta = gp(W, num_basis, IDX_WEEKS, IDX_BASIS,
                                delta0, delta1, delta2, 
                                alpha_gp1_t, alpha_gp2_t, alpha_gp1_d, alpha_gp2_d,
                                rho_gp1_t_dist, rho_gp2_t_dist, rho_gp1_d_dist, rho_gp2_d_dist,
                                z1, z2); 

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


