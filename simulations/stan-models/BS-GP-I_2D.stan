data {
  int<lower=1> n;
  int<lower=1> m;
  real x_1[n];
  real x_2[m];
  matrix[n,m] y;
  
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
  sigma ~ std_normal();
  
  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){        
      beta[i,j] ~ std_normal();
      }
  }
  
  
  for(i in 1:n){
    for(j in 1:m){
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






