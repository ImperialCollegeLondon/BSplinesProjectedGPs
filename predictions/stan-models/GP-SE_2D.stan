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
}

data {
  // training
  int<lower=1> n;
  int coordinates[n,2];
  real x_1[n];
  real x_2[n];
  vector[n] y;
  
  // test
  int<lower=1> n_pred;
  real x_1_pred[n_pred];
  real x_2_pred[n_pred];
}

transformed data {
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho_1;
  real<lower=0> rho_2;
  real<lower=0> alpha_1;
  real<lower=0> alpha_2;
  matrix[n,n] eta;
  real<lower=0> sigma;
}
transformed parameters {
  matrix[n,n] f = gp(n, n, x_1, x_2,
                              delta,
                              alpha_1, alpha_2, 
                              rho_1,  rho_2,
                              eta);
}

model {
  rho_1 ~ inv_gamma(5, 5);
  rho_2 ~ inv_gamma(5, 5);
  alpha_1 ~ std_normal();
  alpha_2 ~ std_normal();
  sigma ~ std_normal();
  
  for(i in 1:n){
    for(j in 1:n){
        eta[i,j] ~ std_normal();
    }
  }
  
    for(k in 1:n){
        y[k] ~ normal(f[coordinates[k,1],coordinates[k,2]], sigma);
        //print("y[k]", y[k]);
        //print("f[coordinates[k,1],coordinates[k,2]]", f[coordinates[k,1],coordinates[k,2]]);
  }


  
}

generated quantities {
  matrix[n,n] y_hat;
  matrix[n,n] eta_pred;
  matrix[n,n] f_pred;
  matrix[n_pred,n_pred] y_pred;
  real log_lik[n];
  
  {
    for(i in 1:n){
        log_lik[i] = normal_lpdf(y[i] | f[coordinates[i,1],coordinates[i,2]], sigma);
        
      for(j in 1:n){
      y_hat[i,j] = normal_rng(f[i,j], sigma);
      }
    }
  }

  for(i_pred in 1:n_pred){
    for(j_pred in 1:n_pred){
      eta_pred[i_pred,j_pred] = normal_rng(0, 1);
    }
  }
  
  f_pred = gp(n_pred, n_pred, x_1_pred, x_2_pred, delta, alpha_1, alpha_2, rho_1, rho_2, eta_pred);
  for(i_pred in 1:n_pred){
    for(j_pred in 1:n_pred){
      y_pred[i_pred,j_pred] = normal_rng(f_pred[i_pred,j_pred], sigma);
    }
  }
}

