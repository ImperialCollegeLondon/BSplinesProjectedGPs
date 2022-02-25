functions {
    matrix kron_mvprod(matrix A, matrix B, matrix V) 
    {
        return transpose(A*transpose(B*V));
    }
    
  matrix gp(int N_rows, int N_columns, real[] rows_idx, real[] columns_index,
            real delta0,
            real alpha_gp, 
            real rho_gp1, real rho_gp2,
            matrix z1)
  {
    
    matrix[N_rows,N_columns] GP;
    
    matrix[N_rows, N_rows] K1;
    matrix[N_rows, N_rows] L_K1;
    
    matrix[N_columns, N_columns] K2;
    matrix[N_columns, N_columns] L_K2;
    
    K1 = cov_exp_quad(rows_idx, alpha_gp, rho_gp1) + diag_matrix(rep_vector(delta0, N_rows));
    K2 = cov_exp_quad(columns_index, alpha_gp, rho_gp2) + diag_matrix(rep_vector(delta0, N_columns));

    L_K1 = cholesky_decompose(K1);
    L_K2 = cholesky_decompose(K2);
    
    GP = kron_mvprod(L_K2, L_K1, z1);

    return(GP);
  }
}

data {
  int<lower=1> n;
  int<lower=1> m;
  real x_1[n];
  real x_2[m];
  int<lower=1> N;
  int coordinates[N,2];
  int y[N];
}

transformed data {
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho_1;
  real<lower=0> rho_2;
  real<lower=0> alpha_gp;
  matrix[n,m] eta;
  real<lower=0> nu;
}
transformed parameters {
  real<lower=0> nu_inverse = (1/nu);
  matrix[n,m] f = exp(gp(n, m, x_1, x_2,
                              delta,
                              alpha_gp, 
                              rho_1,  rho_2,
                              eta));
  matrix[n,m] alpha = f / nu;
}

model {
  nu ~ exponential(1);
    
  rho_1 ~ inv_gamma(2,2);
  rho_2 ~ inv_gamma(2,2);
  alpha_gp ~ cauchy(0,1);
  
  for(i in 1:n){
    for(j in 1:m){
        eta[i,j] ~ std_normal();
    }
  }
  
  for(k in 1:N){
      y[k] ~ neg_binomial(alpha[coordinates[k,1],coordinates[k,2]], nu_inverse);
  }

  
}

generated quantities {
 real log_lik[N];
  
  for(k in 1:N)
    log_lik[k] = neg_binomial_lpmf(y[k] | alpha[coordinates[k,1],coordinates[k,2]], nu_inverse);

}




