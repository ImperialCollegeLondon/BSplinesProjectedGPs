functions {
    matrix kron_mvprod(matrix A, matrix B, matrix V) 
    {
        return transpose(A*transpose(B*V));
    }
    
  matrix gp(int N_rows, int N_columns, real[] rows_idx, real[] columns_index,
            real delta0,
            real zeta_gp, 
            real rho_gp1, real rho_gp2,
            matrix z1)
  {
    
    matrix[N_rows,N_columns] GP;
    
    matrix[N_rows, N_rows] K1;
    matrix[N_rows, N_rows] L_K1;
    
    matrix[N_columns, N_columns] K2;
    matrix[N_columns, N_columns] L_K2;
    
    K1 = cov_exp_quad(rows_idx, sqrt(zeta_gp), rho_gp1) + diag_matrix(rep_vector(delta0, N_rows));
    K2 = cov_exp_quad(columns_index, sqrt(zeta_gp), rho_gp2) + diag_matrix(rep_vector(delta0, N_columns));

    L_K1 = cholesky_decompose(K1);
    L_K2 = cholesky_decompose(K2);
    
    GP = kron_mvprod(L_K2, L_K1, z1);

    return(GP);
  }
  
  matrix get_2D_surface(int A, int W, 
                        int num_basis_rows, int num_basis_columns, matrix BASIS_ROWS, matrix BASIS_COLUMNS, 
                        real[] IDX_BASIS_ROWS, real[] IDX_BASIS_COLUMNS, 
                        real zeta_gp, real gamma_gp1, real gamma_gp2, matrix z1, real delta0){
                          
    matrix[num_basis_rows,num_basis_columns] beta = gp(num_basis_rows, num_basis_columns, IDX_BASIS_ROWS, IDX_BASIS_COLUMNS, delta0,
              zeta_gp, gamma_gp1,  gamma_gp2, z1);
    matrix[A, W] f  = (BASIS_ROWS') * beta * BASIS_COLUMNS;
    
    return(f);
  }
  
  matrix get_age_profile(int A, int W, matrix f){
    matrix[A,W] phi;
    
    for(w in 1:W){
      phi[:,w] = softmax( f[:,w] ); 
    }
    
    return(phi);         
   }
  
  matrix get_dirichlet_parameter(int A, int W, matrix phi, vector lambda, real nu){
    matrix[A,W] alpha = phi .* rep_matrix(to_row_vector(lambda ./ rep_vector(nu, W)), A);
    
    return(alpha);
  }
  
  matrix get_reduced_matrix(int B, int W, matrix matrix_to_reduce, matrix state_age_strata){
    
    matrix[B,W] matrix_reduced = state_age_strata * matrix_to_reduce;
    
     return(matrix_reduced);
   }
   
   matrix get_expected_deaths(int C, int W, row_vector deaths_JHU, matrix phi_reduced_vac){
     matrix[C,W] E_pdeaths = phi_reduced_vac .* rep_matrix(deaths_JHU, C);
     
     return(E_pdeaths);
   }
   
   matrix get_resurgence_deaths(int C, int T, matrix E_pdeaths, int w_start_resurgence, int w_stop_resurgence){
     
     vector[C] max_Epdeaths;
     matrix[C,T] r_pdeaths = rep_matrix(1.0, C, T);
     
     for(c in 1:C){
        max_Epdeaths[c] = max(E_pdeaths[c,1:(w_start_resurgence - 1)]);
     }
     
     r_pdeaths = E_pdeaths[:,w_start_resurgence:w_stop_resurgence] ./ rep_matrix(max_Epdeaths, T);
     
    return(r_pdeaths);
   }
   
   matrix get_expected_resurgence_deaths(int C, int T, vector prop_vac_start, vector week_indices_resurgence,
          vector intercept_resurgence0, vector intercept_resurgence_re, vector vaccine_effect_intercept_diagonal,
          vector slope_resurgence0, vector vaccine_effect_slope_diagonal, 
          vector vaccine_effect_intercept_cross, vector vaccine_effect_slope_cross){

    vector[C] intercept_resurgence = intercept_resurgence0 + intercept_resurgence_re ;
    vector[C] slope_resurgence = slope_resurgence0;
    matrix[C,T] log_xi;
    matrix[C,T] xi;
     
    for(c in 1:C){
        for(c_prime in 1:C){
          if(c_prime == c){
            intercept_resurgence[c] += prop_vac_start[c] .* vaccine_effect_intercept_diagonal[c];
            slope_resurgence[c] += prop_vac_start[c] .* vaccine_effect_slope_diagonal[c];
          } else{
            intercept_resurgence[c] += prop_vac_start[c_prime] .* vaccine_effect_intercept_cross[c];
            slope_resurgence[c] += prop_vac_start[c_prime] .* vaccine_effect_slope_cross[c];
          }
      }
    }

    log_xi = rep_matrix(intercept_resurgence, T) + (rep_matrix(to_row_vector(week_indices_resurgence),C) .* rep_matrix(slope_resurgence, T));
    xi = exp(log_xi);
    
    return(xi);
  }
  
  real my_neg_binomial_lpmf(int y, real alpha, real beta){
    return(lchoose(y + alpha - 1, alpha - 1) + alpha * (log(beta) - log(beta + 1)) + y * (log(1) - log(beta + 1)));
  }
  
  real countries_log_dens(int[,] sum_count_censored, 
                          int start,
                          int end,
                          // data
                          int[,,] deaths,
                          row_vector[] lambda_prior_parameters,
                          int W_OBSERVED,
                          int[,] N_idx_non_missing,
                          int[,,] idx_non_missing,
                          int[] IDX_WEEKS_OBSERVED,
                          int[] IDX_WEEKS_OBSERVED_REPEATED,
                          int[] N_missing,
                          int[,] idx_weeks_missing_min,
                          int[,] idx_weeks_missing_max,
                          int[,] age_missing,
                          int[,] start_or_end_period,
                          int[,] min_count_censored,
                          int[,] max_count_censored,
                          int A,
                          int B,
                          int W,
                          int T,
                          int C,
                          row_vector[] prop_vac_start,
                          vector week_indices_resurgence,
                          int[] w_start_resurgence, 
                          int[] w_stop_resurgence,
                          int num_basis_rows, 
                          int num_basis_columns, 
                          matrix BASIS_ROWS, 
                          matrix BASIS_COLUMNS, 
                          real[] IDX_BASIS_ROWS, 
                          real[] IDX_BASIS_COLUMNS, 
                          real delta0,
                          matrix state_age_strata,
                          matrix vac_age_strata,
                          matrix deaths_JHU,
                          // parameters 
                           real[] zeta_gp, 
                          real[] gamma_gp1, 
                          real[] gamma_gp2, 
                          matrix[] z1, 
                          vector[] lambda_raw,
                          real[] nu,
                          real[] sigma_r_pdeaths,
                          real[] intercept_resurgence0, 
                          row_vector[] intercept_resurgence_re, 
                          real[] vaccine_effect_intercept_diagonal,
                          real[] slope_resurgence0, 
                          real[] vaccine_effect_slope_diagonal,
                          real[] vaccine_effect_intercept_cross, 
                          real[] vaccine_effect_slope_cross
                          ){
    real lpmf = 0.0;
    int M_slice = end - start + 1;
    
    for(m_slice in 1:M_slice) {
      int m = m_slice + start - 1;
      vector[W] lambda;
      real nu_inverse;
      matrix[A,W] phi;
      matrix[A,W] alpha;
      matrix[B,W] phi_reduced;
      matrix[C,W] phi_reduced_vac;
      matrix[B,W] alpha_reduced;
      matrix[A, W] f;
      matrix[C,W] E_pdeaths;
      matrix[C,T] r_pdeaths = rep_matrix(1.0, C, T);
      matrix[C,T] log_r_pdeaths = rep_matrix(0.0, C, T);
      matrix[C,T] xi;
  
      lambda = lambda_raw[m][IDX_WEEKS_OBSERVED_REPEATED];
      nu_inverse = (1/nu[m]);
        
      f = get_2D_surface(A, W, num_basis_rows, num_basis_columns, BASIS_ROWS, BASIS_COLUMNS,
                          IDX_BASIS_ROWS, IDX_BASIS_COLUMNS,
                          zeta_gp[m], gamma_gp1[m], gamma_gp2[m], z1[m], delta0);
      phi  =  get_age_profile(A,W,f);
      alpha = get_dirichlet_parameter(A,W,phi, lambda, nu[m]);
      phi_reduced = get_reduced_matrix(B, W, phi, state_age_strata);
      alpha_reduced = get_reduced_matrix(B, W, alpha, state_age_strata);
      phi_reduced_vac = get_reduced_matrix(C, W, phi, vac_age_strata);
      
      E_pdeaths = get_expected_deaths(C,W,to_row_vector(deaths_JHU[m,]), phi_reduced_vac);
      r_pdeaths = get_resurgence_deaths(C, T,E_pdeaths,w_start_resurgence[m], w_stop_resurgence[m]);
      log_r_pdeaths = log(r_pdeaths);
       
      xi = get_expected_resurgence_deaths(C, T, to_vector(prop_vac_start[:,m]), week_indices_resurgence,
             to_vector(intercept_resurgence0), to_vector(intercept_resurgence_re[:,m]), to_vector(vaccine_effect_intercept_diagonal),
            to_vector(slope_resurgence0), to_vector(vaccine_effect_slope_diagonal), 
            to_vector(vaccine_effect_intercept_cross), to_vector(vaccine_effect_slope_cross));

      lpmf +=  exponential_lpdf(lambda_raw[m] | rep_row_vector(1.0, W_OBSERVED) ./ lambda_prior_parameters[m] );
      
      // Note on the neg bin parametrisation related to the paper:
      // mean neg_binomial_lpmf is alpha_reduced / nu_inverse = alpha_reduced * nu
      // var neg_binomial_lpmf is alpha_reduced / nu_inverse^2 * (nu_inverse + 1) = alpha_reduced * nu (1 + nu)
    
      for(w in 1:W_OBSERVED){
      
      int indx_w_m[N_idx_non_missing[m][w]] = idx_non_missing[m][1:N_idx_non_missing[m][w],w];
      lpmf += neg_binomial_lpmf(deaths[m][indx_w_m,w] | alpha_reduced[indx_w_m, IDX_WEEKS_OBSERVED[w]] , nu_inverse );
        
      }
      
      for(n in 1:N_missing[m]){
      real alpha_reduced_missing = sum(alpha_reduced[ age_missing[m][n], idx_weeks_missing_min[m][n]:idx_weeks_missing_max[m][n] ]);

      if(!start_or_end_period[m][n])
      {
        lpmf += neg_binomial_lpmf( sum_count_censored[m][n] | alpha_reduced_missing , nu_inverse ) ;
      } 
      else {
       for(i in min_count_censored[m][n]:max_count_censored[m][n]){
         // if(i == 1){ // bug cmdstan when argument is 1
           lpmf += my_neg_binomial_lpmf( i | alpha_reduced_missing , nu_inverse ) ;
         // }else{
         //   lpmf += neg_binomial_lpmf( i | alpha_reduced_missing , nu_inverse ) ;
         // }
          
       }
      }
    } 

    for(c in 1:C){
      lpmf += normal_lpdf(log_r_pdeaths[c,1:T] | log(xi[c,1:T]), sigma_r_pdeaths[m]  );
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
  row_vector[W_OBSERVED] lambda_prior_parameters[M]; // parameters of the prior distribution of lambda
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
  row_vector[M] prop_vac_start[C]; // pre-resurgence proportion of vaccinated individuals
  
  // counterfactual
  int<lower=1> N_COUNTERFACTUAL; //number of counterfactual analysis
  matrix[N_COUNTERFACTUAL,M] prop_vac_start_counterfactual[C]; // pre-resurgence proportion of vaccinated individuals
  
}

transformed data
{   
    real delta0 = 1e-9;  
    int N_log_lik = 0;
    matrix[B, A] state_age_strata = rep_matrix(0.0, B, A);
    matrix[C, A] vac_age_strata = rep_matrix(0.0, C, A);
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
  
  for(b in 1:B){
   state_age_strata[b, age_from_state_age_strata[b]:age_to_state_age_strata[b]] = rep_row_vector(1.0, age_to_state_age_strata[b] - age_from_state_age_strata[b] + 1);
  }
  
  for(c in 1:C){
    vac_age_strata[c, age_from_vac_age_strata[c]:age_to_vac_age_strata[c]] = rep_row_vector(1.0, age_to_vac_age_strata[c] - age_from_vac_age_strata[c] + 1);
  }

}

parameters {
  real<lower=0> nu[M];
  vector<lower=0>[W-W_NOT_OBSERVED] lambda_raw[M];
  matrix[num_basis_rows,num_basis_columns] z1[M];
  real<lower=0> zeta_gp[M];
  real<lower=0> gamma_gp1[M]; 
  real<lower=0> gamma_gp2[M];

  real intercept_resurgence0[C];
  row_vector[M] intercept_resurgence_re[C];
  real<lower=0> sigma_intercept_resurgence;
  real slope_resurgence0[C];
  real vaccine_effect_intercept_diagonal[C];
  real vaccine_effect_slope_diagonal[C];
  real<lower=0> sigma_r_pdeaths[M];
  real vaccine_effect_intercept_cross[C];
  real vaccine_effect_slope_cross[C];
}

model {
  
  target += exponential_lpdf(nu | 1);

  target += cauchy_lpdf(zeta_gp | 0,1);
  target += inv_gamma_lpdf(gamma_gp1 | 2, 2);
  target += inv_gamma_lpdf(gamma_gp2 | 2, 2);

  target += normal_lpdf(intercept_resurgence0 | 0,5);
  target += normal_lpdf(slope_resurgence0 | 0,5);
  
  target += normal_lpdf(vaccine_effect_intercept_diagonal | 0,2.5);
  target += normal_lpdf(vaccine_effect_slope_diagonal | 0,2.5);
    
  target += normal_lpdf(vaccine_effect_intercept_cross | 0,2.5);
  target += normal_lpdf(vaccine_effect_slope_cross | 0,2.5);
  
  target += cauchy_lpdf(sigma_intercept_resurgence | 0,1);
  target += cauchy_lpdf(sigma_r_pdeaths | 0,1);

  for(c in 1:C){
    target += normal_lpdf(intercept_resurgence_re[c] | 0,sigma_intercept_resurgence);
  }

  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
      target += normal_lpdf(z1[:,i,j] | 0,1);
    }
  }
  

  // rstan version
  // target += countries_log_dens(sum_count_censored, 1, M,
  // cmdstan version
  target += reduce_sum(countries_log_dens, sum_count_censored, 1,
                              deaths,
                              lambda_prior_parameters,
                              W_OBSERVED,
                              N_idx_non_missing,
                              idx_non_missing,
                              IDX_WEEKS_OBSERVED,
                              IDX_WEEKS_OBSERVED_REPEATED,
                              N_missing,
                              idx_weeks_missing_min,
                              idx_weeks_missing_max,
                              age_missing,
                              start_or_end_period,
                              min_count_censored,
                              max_count_censored,
                              A,
                              B,
                              W,
                              T,
                              C,
                              prop_vac_start,
                              week_indices_resurgence,
                              w_start_resurgence,
                              w_stop_resurgence,
                              num_basis_rows,
                              num_basis_columns,
                              BASIS_ROWS,
                              BASIS_COLUMNS,
                              IDX_BASIS_ROWS,
                              IDX_BASIS_COLUMNS,
                              delta0,
                              state_age_strata,
                              vac_age_strata,
                              deaths_JHU,
                              zeta_gp,
                              gamma_gp1,
                              gamma_gp2,
                              z1,
                              lambda_raw,
                              nu,
                              sigma_r_pdeaths,
                              intercept_resurgence0,
                              intercept_resurgence_re,
                              vaccine_effect_intercept_diagonal,
                              slope_resurgence0,
                              vaccine_effect_slope_diagonal,
                              vaccine_effect_intercept_cross, 
                              vaccine_effect_slope_cross
                          );

}

