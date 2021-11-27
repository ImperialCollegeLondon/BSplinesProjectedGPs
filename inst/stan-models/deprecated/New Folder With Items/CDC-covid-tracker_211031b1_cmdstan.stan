functions {
  matrix kron_mvprod(matrix A, matrix B, matrix V) 
  {
      return transpose(A*transpose(B*V));
  }
    
  matrix gp(int N_rows, int N_columns, real[] rows_idx, real[] columns_index,
            real delta0,
            real alpha_gp1, real alpha_gp2, 
            real rho_gp1, real rho_gp2,
            matrix z1)
  {
    
    matrix[N_rows,N_columns] GP;
    
    matrix[N_rows, N_rows] K1;
    matrix[N_rows, N_rows] L_K1;
    
    matrix[N_columns, N_columns] K2;
    matrix[N_columns, N_columns] L_K2;
    
    K1 = cov_exp_quad(rows_idx, alpha_gp1, rho_gp1) + diag_matrix(rep_vector(delta0, N_rows));
    K2 = cov_exp_quad(columns_index, alpha_gp2, rho_gp2) + diag_matrix(rep_vector(delta0, N_columns));

    L_K1 = cholesky_decompose(K1);
    L_K2 = cholesky_decompose(K2);
    
    GP = kron_mvprod(L_K2, L_K1, z1);

    return(GP);
  }
  
  matrix low_rank_surface(int num_basis_rows, int num_basis_columns, 
                          real[] IDX_BASIS_ROWS, real[] IDX_BASIS_COLUMNS, 
                          real delta0, real alpha_gp1, real alpha_gp2, 
                          real rho_gp1, real rho_gp2, matrix z1)
  {
      matrix[num_basis_rows,num_basis_columns] beta = gp(num_basis_rows, num_basis_columns, IDX_BASIS_ROWS, IDX_BASIS_COLUMNS, delta0,
              alpha_gp1, alpha_gp2, rho_gp1, rho_gp2, z1);
              
      return(beta);
    }
  
  matrix age_week_surface(int A, int W, 
                          matrix BASIS_ROWS, matrix BASIS_COLUMNS, 
                          matrix beta)
  {
      matrix[A, W] f = (BASIS_ROWS') * beta * BASIS_COLUMNS;
      
      return(f);
    }
  
  vector total_deaths(int W, vector lambda_raw, int[] IDX_WEEKS_OBSERVED_REPEATED)
  {
      vector[W] lambda = lambda_raw[IDX_WEEKS_OBSERVED_REPEATED];
      
      return(lambda);
    }
  
  matrix scale_contribution_to_deaths(int A, int W, matrix f)
  {
      matrix[A,W] phi; 
      
      for(w in 1:W){
        phi[:,w] = softmax( f[:,w] ); 
      }
      
      return(phi);
    }

  matrix alpha_coefficients(int A, int W, matrix phi, vector lambda, real nu)
  {
      matrix[A,W] alpha;
      for(w in 1:W)
        alpha[:,w] = phi[:,w] * lambda[w] / nu ;
      
      return(alpha);
    }
  
  matrix alpha_coefficients_reduced(int B, int W, 
                                    matrix alpha,
                                    int[] age_from_state_age_strata, int[] age_to_state_age_strata)
  {
      matrix[B,W] alpha_reduced;
      
      for(w in 1:W){
        for(b in 1:B){
          alpha_reduced[b,w] = sum(alpha[age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
        }
      }
      
      return(alpha_reduced);
    }
  
  matrix scale_contribution_to_deaths_reduced(int B, int W, 
                                    matrix phi,
                                    int[] age_from_state_age_strata, int[] age_to_state_age_strata)
  {
      matrix[B,W] phi_reduced;
      
      for(w in 1:W){
        for(b in 1:B){
          phi_reduced[b,w] = sum(phi[age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
        }
      }
      
      return(phi_reduced);
    }
  
  matrix E_deaths_vaccination_age_groups(int W, int C, 
                                         matrix phi_reduced_vac, row_vector deaths_JHU)
  {
      matrix[C, W] E_pdeaths;
      
      for(w in 1:W){
        for(c in 1:C){
          E_pdeaths[c,w] = phi_reduced_vac[c,w] * deaths_JHU[w];
        }
      }
      return(E_pdeaths);
    }
  
  matrix resurgence_deaths_vaccination_age_groups(int W, int T, int C,
                                                  matrix E_pdeaths, 
                                                  int w_start_resurgence, int w_stop_resurgence)
  {
      matrix[C,T] r_pdeaths;
      
      for(w in 1:W){
        for(c in 1:C){
          if(w >= w_start_resurgence && w <= w_stop_resurgence){
            r_pdeaths[c,w - w_start_resurgence + 1] = E_pdeaths[c,w] / max(E_pdeaths[c,1:(w_start_resurgence - 1)]) ;
          }
        }
      }
      return(r_pdeaths);
    }
    
    
  real countries_log_dens(int[,,] deaths_slice,
                          int start,
                          int end,
                          int A, 
                          int B,
                          int C,
                          int W, 
                          int T,
                          int[,] start_or_end_period,
                          int W_OBSERVED,
                          int[] IDX_WEEKS_OBSERVED,
                          int[] IDX_WEEKS_OBSERVED_REPEATED,
                          int[,,] idx_non_missing,
                          int[,] N_idx_non_missing,
                          int[] N_missing,
                          int[,] age_missing,
                          int[,] idx_weeks_missing_min,
                          int[,] idx_weeks_missing_max,
                          int[,] sum_count_censored,
                          int[,] min_count_censored,
                          int[,] max_count_censored,
                          int[] age_from_state_age_strata, 
                          int[] age_to_state_age_strata,
                          int[] age_from_vac_age_strata, 
                          int[] age_to_vac_age_strata,
                          int[] w_start_resurgence,
                          int[] w_stop_resurgence,
                          int num_basis_rows, 
                          int num_basis_columns, 
                          real[] IDX_BASIS_ROWS, 
                          real[] IDX_BASIS_COLUMNS, 
                          matrix BASIS_ROWS, 
                          matrix BASIS_COLUMNS, 
                          matrix deaths_JHU,
                          real delta0, 
                          real[] alpha_gp1, 
                          real[] alpha_gp2, 
                          real[] rho_gp1, 
                          real[] rho_gp2, 
                          matrix[] z1,
                          vector[] lambda_raw, 
                          real[] nu_unscaled,
                          matrix[] xi, 
                          real[] sigma_r_pdeaths)
  {
        
    real lpmf = 0.0;
    int M_slice = end - start + 1;
    int index_country_slice; 
    int min_age_slice;
    int max_age_slice;
    vector[W] lambda;
    matrix[num_basis_rows,num_basis_columns] beta;
    matrix[A, W] f;
    matrix[A,W] phi; 
    matrix[A,W] alpha;
    matrix[B,W] alpha_reduced;
    matrix[B,W] phi_reduced;
    matrix[C,W] phi_reduced_vac;
    matrix[C,W] E_pdeaths;
    matrix[C,T] r_pdeaths;
    matrix[C,T] log_r_pdeaths;

    for(m_slice in 1:M_slice) 
    {
      int m = m_slice + start - 1;
      real nu = (1/nu_unscaled[m])^2;
      real theta = (1 / nu);

      // transformed parameters
      lambda = total_deaths(W, lambda_raw[m], IDX_WEEKS_OBSERVED_REPEATED);

      beta = low_rank_surface(num_basis_rows, num_basis_columns, IDX_BASIS_ROWS, IDX_BASIS_COLUMNS, delta0, alpha_gp1[m], alpha_gp2[m], 
                              rho_gp1[m], rho_gp2[m], z1[m]);
      
      f = age_week_surface(A, W, BASIS_ROWS, BASIS_COLUMNS, beta);
      phi = scale_contribution_to_deaths(A, W, f);
      alpha = alpha_coefficients(A, W, phi, lambda, nu);
      alpha_reduced = alpha_coefficients_reduced(B, W, alpha, age_from_state_age_strata, age_to_state_age_strata);
      phi_reduced = scale_contribution_to_deaths_reduced(B, W, phi, age_from_state_age_strata, age_to_state_age_strata);

      phi_reduced_vac = scale_contribution_to_deaths_reduced(C, W, phi, age_from_vac_age_strata, age_to_vac_age_strata);
      E_pdeaths = E_deaths_vaccination_age_groups(W, C, phi_reduced_vac, deaths_JHU[m,]);
      r_pdeaths = resurgence_deaths_vaccination_age_groups(W, T, C, E_pdeaths, w_start_resurgence[m], w_stop_resurgence[m]);
      log_r_pdeaths = log(r_pdeaths);

      // likelihood 
       for(w in 1:W_OBSERVED){
          int indx_w_m[N_idx_non_missing[m][w]] = idx_non_missing[m][1:N_idx_non_missing[m][w],w];
          lpmf += neg_binomial_lpmf(deaths_slice[m_slice, indx_w_m, w] | alpha_reduced[indx_w_m, IDX_WEEKS_OBSERVED[w]] , theta );
      }

      for(n in 1:N_missing[m]){
          real alpha_reduced_missing = sum(alpha_reduced[ age_missing[m][n], idx_weeks_missing_min[m][n]:idx_weeks_missing_max[m][n] ]);

          if(!start_or_end_period[m][n])
          {
           lpmf += neg_binomial_lpmf( sum_count_censored[m][n] | alpha_reduced_missing , theta ) ;
          } 
          else {
              for(i in min_count_censored[m][n]:max_count_censored[m][n])
                  lpmf += neg_binomial_lpmf( i | alpha_reduced_missing , theta ) ;
          }
      }    
      
      for(c in 1:C){
        lpmf += gamma_lpdf(r_pdeaths[c,:] | square(xi[c][:,m] / sigma_r_pdeaths[m]), xi[c][:,m] ./ rep_vector(square(sigma_r_pdeaths[m]), T) );
      }
    }
    return(lpmf);
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
  
  // JHU data
  matrix[M,W] deaths_JHU; // deaths reported by JHU

  //splines
  int num_basis_rows;
  int num_basis_columns;
  matrix[num_basis_rows, A] BASIS_ROWS; 
  matrix[num_basis_columns, W] BASIS_COLUMNS; 
  
  // GP
  real IDX_BASIS_ROWS[num_basis_rows];
  real IDX_BASIS_COLUMNS[num_basis_columns];

  // vaccine effect
  int<lower=1,upper=W> w_stop_resurgence[M]; // index of the week when Summer 2021 resurgences stop
  int<lower=1,upper=W> w_start_resurgence[M]; // index of the week when Summer 2021 resurgences starts
  int<lower=1> T; // number of weeks during the Summer 2021 resurgences (w_stop_resurgence-w_start_resurgence+1)
  vector[T] week_indices_resurgence; // 0:(T-1)
  int<lower=1> C; // number of age groups in vaccination data
  int age_from_vac_age_strata[C]; // age from of age band c
  int age_to_vac_age_strata[C];// age to of age band c
  int<lower=1,upper=C> c_counterfactual; // age group in counterfactual scenario
  matrix[T,M] prop_vac[C]; // proportion of vaccinated individuals
  row_vector[M] prop_vac_start[C]; // pre-resurgence proportion of vaccinated individuals
  
}

transformed data
{   
    real delta0 = 1e-9;  
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
  real<lower=0> nu_unscaled[M];
  vector<lower=0>[W-W_NOT_OBSERVED] lambda_raw[M];
  matrix[num_basis_rows,num_basis_columns] z1[M];
  real<lower=0> alpha_gp1[M];
  real<lower=0> alpha_gp2[M];
  real<lower=0> rho_gp1[M]; 
  real<lower=0> rho_gp2[M];

  real intercept_resurgence0[C];
  row_vector[M] intercept_resurgence_re[C];
  real<lower=0> sigma_intercept_resurgence[C];
  real slope_resurgence0[C];
  real vaccine_effect_intercept[C,C];
  real vaccine_effect_slope[C,C];
  real<lower=0> sigma_r_pdeaths[M];
}

transformed parameters {
  row_vector[M] intercept_resurgence[C];
  row_vector[M] slope_resurgence[C];
  matrix[T, M] xi[C];
  matrix[T, M] log_xi[C];
  
  for(c in 1:C){
    intercept_resurgence[c] = rep_row_vector(intercept_resurgence0[c], M) + intercept_resurgence_re[c];
    slope_resurgence[c] = rep_row_vector(slope_resurgence0[c], M) ;

    for(c_prime in 1:C){
      intercept_resurgence[c] += prop_vac_start[c_prime] .* rep_row_vector(vaccine_effect_intercept[c,c_prime], M);
      slope_resurgence[c] += prop_vac_start[c_prime] .* rep_row_vector(vaccine_effect_slope[c,c_prime], M);
    }

    log_xi[c] = rep_matrix(intercept_resurgence[c], T) + rep_matrix(week_indices_resurgence, M) .* rep_matrix(slope_resurgence[c], T);
    xi[c] = exp(log_xi[c]);
  }
  
}

model {
  
  nu_unscaled ~ normal(0,1);

  alpha_gp1 ~ cauchy(0,1);
  alpha_gp2 ~ cauchy(0,1);
  rho_gp1 ~ inv_gamma(5, 5);
  rho_gp2 ~ inv_gamma(5, 5);

  intercept_resurgence0 ~ normal(0,0.5);
  sigma_intercept_resurgence ~ cauchy(0,1);
  slope_resurgence0 ~ normal(0,0.5);
  sigma_r_pdeaths ~ cauchy(0,1);

  for(c in 1:C){
    intercept_resurgence_re[c] ~ normal(0,sigma_intercept_resurgence[c]);

    vaccine_effect_intercept[c,:] ~ normal(0,0.5);
    vaccine_effect_slope[c,:] ~ normal(0,0.5);
  }

  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
      z1[:,i,j] ~ normal(0,1);
    }
  }
  

  for(m in 1:M){
    
    lambda_raw[m] ~ gamma( lambda_prior_parameters[m][1,:],lambda_prior_parameters[m][2,:]);
    
  }

  // rstan version
  target += countries_log_dens(deaths, 1, M,
  // cmdstan version
  // target += reduce_sum(countries_log_dens, deaths, 1,
                       A, 
                       B,
                       C,
                       W, 
                       T,
                       start_or_end_period,
                       W_OBSERVED,
                       IDX_WEEKS_OBSERVED,
                       IDX_WEEKS_OBSERVED_REPEATED,
                       idx_non_missing,
                       N_idx_non_missing,
                       N_missing,
                       age_missing,
                       idx_weeks_missing_min,
                       idx_weeks_missing_max,
                       sum_count_censored,
                       min_count_censored,
                       max_count_censored,
                       age_from_state_age_strata, 
                       age_to_state_age_strata,
                       age_from_vac_age_strata, 
                       age_to_vac_age_strata,
                       w_start_resurgence,
                       w_stop_resurgence,
                       num_basis_rows, 
                       num_basis_columns, 
                       IDX_BASIS_ROWS, 
                       IDX_BASIS_COLUMNS, 
                       BASIS_ROWS, 
                       BASIS_COLUMNS, 
                       deaths_JHU,
                       delta0, 
                       alpha_gp1, 
                       alpha_gp2, 
                       rho_gp1, 
                       rho_gp2, 
                       z1,
                       lambda_raw, 
                       nu_unscaled,
                       xi, 
                       sigma_r_pdeaths);


}







