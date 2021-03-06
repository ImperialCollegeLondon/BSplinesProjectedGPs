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
  int<lower=1> n; // number of rows
  int<lower=1> m; // number of columns
  int<lower=1> N; // number of entries observed
  int coordinates[N,2]; // coordinate of entries observed
  vector[N] y; // data on entries observed
  
  //splines
  int num_basis_rows; // number of B-Splines basis functions rows 
  int num_basis_columns; // number of B-Splines basis functions columns 
  matrix[num_basis_rows, n] BASIS_ROWS; // B-splines basis functions on the rows
  matrix[num_basis_columns, m] BASIS_COLUMNS; // B-splines basis functions on the columns
  
  // GP
  real IDX_BASIS_ROWS[num_basis_rows]; // index of the B-splines basis functions rows
  real IDX_BASIS_COLUMNS[num_basis_columns]; // index of the B-splines basis functions columns
}

transformed data {
  real delta = 1e-9;
}

parameters {
  real<lower=0> rho_1; // length scale rows
  real<lower=0> rho_2; // length scale columns
  real<lower=0> alpha_gp; // output variance
  matrix[num_basis_rows,num_basis_columns] eta; // GP variables
  real<lower=0> sigma; // observational noise 
}

transformed parameters {
  matrix[num_basis_rows,num_basis_columns] beta = gp(num_basis_rows, num_basis_columns, 
                                                     IDX_BASIS_ROWS, IDX_BASIS_COLUMNS,
                                                     delta,
                                                     alpha_gp,
                                                     rho_1,  rho_2,
                                                     eta); 
                                                     
  matrix[n,m] f = (BASIS_ROWS') * beta * BASIS_COLUMNS;
}

model {
  rho_1 ~ inv_gamma(2,2);
  rho_2 ~ inv_gamma(2,2);
  alpha_gp ~ cauchy(0,1);
  
  sigma ~ cauchy(0,1);
  
  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
        eta[i,j] ~ std_normal();
    }
  }
  
  for(k in 1:N){
        y[k] ~ normal(f[coordinates[k,1],coordinates[k,2]], sigma);
  }

}

generated quantities {
 matrix[n,m] y_hat;
 real log_lik[N];

  {
    for(k in 1:N)
      log_lik[k] = normal_lpdf(y[k] | f[coordinates[k,1],coordinates[k,2]], sigma);

    for(i in 1:n){
      for(j in 1:m){
      y_hat[i,j] = normal_rng(f[i,j], sigma);
      }
    }
  }


}


