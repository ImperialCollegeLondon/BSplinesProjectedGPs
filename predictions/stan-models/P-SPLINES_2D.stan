functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2){
    return -0.5 * dot_self(phi[node1] - phi[node2])
      + normal_lpdf(sum(phi) | 0, 0.001 * N);
  }
}

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
  
  // ICAR model
  int<lower=1> K;
  int<lower=1> N_edges;
  int<lower=1, upper=K> node1[N_edges];
  int<lower=1, upper=K> node2[N_edges];
}

parameters {
  real<lower=0> sigma;
  vector[K] beta_raw; 
}

transformed parameters {
  matrix[num_basis_rows,num_basis_columns] beta = to_matrix(beta_raw, num_basis_rows,num_basis_columns); 
                                                     
  matrix[n,m] f = (BASIS_ROWS') * beta * BASIS_COLUMNS;
}

model {
  sigma ~ cauchy(0,1);

  beta_raw ~ icar_normal_lpdf(K, node1, node2);
  
  for(i in 1:N){
        y[i] ~ normal(f[coordinates[i,1],coordinates[i,2]], sigma);
  }

}

generated quantities {
 matrix[n,m] y_hat;
 real log_lik[N];

  {
    for(i in 1:N)
      log_lik[i] = normal_lpdf(y[i] | f[coordinates[i,1],coordinates[i,2]], sigma);

    for(i in 1:n){
      for(j in 1:m){
      y_hat[i,j] = normal_rng(f[i,j], sigma);
      }
    }
  }


}


