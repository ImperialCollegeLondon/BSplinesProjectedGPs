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

parameters {
  matrix[num_basis_rows,num_basis_columns] beta;
  real<lower=0> nu_unscaled;
}

transformed parameters {
  real<lower=0> nu = (1/nu_unscaled)^2;
  real<lower=0> theta = (1 / nu);
   matrix[n,m] f = exp( (BASIS_ROWS') * beta * BASIS_COLUMNS );
   matrix[n,m] alpha = f / nu;
}

model {
  nu_unscaled ~ normal(0, 1);
  
  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
        beta[i,j] ~ std_normal();
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






