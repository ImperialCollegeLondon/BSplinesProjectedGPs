data {
  // training
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> N;
  int coordinates[N,2];
  vector[N] y;
  
  //splines
  int num_basis_rows;
  int num_basis_columns;
  matrix[num_basis_rows, n] BASIS_ROWS; 
  matrix[num_basis_columns, m] BASIS_COLUMNS; 
}

parameters {
  matrix[num_basis_rows,num_basis_columns] beta;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[n,m] f = (BASIS_ROWS') * beta * BASIS_COLUMNS;
}

model {
  sigma ~ cauchy(0,1);
  
  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
        beta[i,j] ~ normal(0, 5);
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


