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
    
    matrix[N_rows, N_rows] K1 = alpha_gp * diag_matrix(rep_vector(1.0, N_rows));
    matrix[N_rows, N_rows] L_K1;
    
    matrix[N_columns, N_columns] K2 = alpha_gp *  diag_matrix(rep_vector(1.0, N_columns));
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
  int<lower=1> N;
  int coordinates[N,2];
  int y[N];
  
  //splines
  int num_basis_rows;
  int num_basis_columns;
  matrix[num_basis_rows, n] BASIS_ROWS; 
  matrix[num_basis_columns, m] BASIS_COLUMNS; 
}

transformed data {
  real delta = 1e-9;
}

parameters {
  real<lower=0> alpha_gp;
  matrix[num_basis_rows,num_basis_columns] eta;
  real<lower=0> nu_unscaled;
}

transformed parameters {
  real<lower=0> nu = (1/nu_unscaled)^2;
  real<lower=0> theta = (1 / nu);
  matrix[num_basis_rows,num_basis_columns] beta = gp(num_basis_rows, num_basis_columns, 
                                                     delta,
                                                     alpha_gp, 
                                                     eta); 
   matrix[n,m] f = exp( (BASIS_ROWS') * beta * BASIS_COLUMNS );
   matrix[n,m] alpha = f / nu;
}

model {
  alpha_gp ~ cauchy(0,1);
  nu_unscaled ~ normal(0, 1);
  
  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
        eta[i,j] ~ std_normal();
    }
  }
  
  for(k in 1:N){
      y[k] ~ neg_binomial(alpha[coordinates[k,1],coordinates[k,2]], theta);
  }

}

generated quantities {
 real log_lik[N];
  for(k in 1:N)
    log_lik[k] = neg_binomial_lpmf(y[k] | alpha[coordinates[k,1],coordinates[k,2]], theta);
}






