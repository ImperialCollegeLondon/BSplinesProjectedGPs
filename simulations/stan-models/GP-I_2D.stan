functions {
    matrix kron_mvprod(matrix A, matrix B, matrix V) 
    {
        return transpose(A*transpose(B*V));
    }
    
  matrix gp(int N_rows, int N_columns, 
            real delta0,
            real alpha_gp, 
            matrix z1)
  {
    
    matrix[N_rows,N_columns] GP;
    
    matrix[N_rows, N_rows] K1 = alpha_gp^2 * diag_matrix(rep_vector(1.0, N_rows));
    matrix[N_rows, N_rows] L_K1;
    
    matrix[N_columns, N_columns] K2 = alpha_gp^2 *  diag_matrix(rep_vector(1.0, N_columns));
    matrix[N_columns, N_columns] L_K2;
    
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
  matrix[n,m] y;
}

transformed data {
  real delta = 1e-9;
}

parameters {
  real<lower=0> alpha_gp;
  matrix[n,m] eta;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[n,m] f = gp(n, m, 
                              delta,
                              alpha_gp, 
                              eta);
}

model {
  alpha_gp ~ std_normal();
  sigma ~ std_normal();
  
  for(i in 1:n){
    for(j in 1:m){
        eta[i,j] ~ std_normal();
        y[i,j] ~ normal(f[i,j], sigma);
    }
  }

  
}

generated quantities {
  matrix[n,m] y_hat;
  real log_lik[n*m];
  
  {
    int counter = 1;
    for(i in 1:n){
      for(j in 1:m){
      y_hat[i,j] = normal_rng(f[i,j], sigma);
      log_lik[counter] = normal_lpdf(y[i,j] | f[i,j], sigma);
      counter = counter + 1;
      }
    }
  }
}


