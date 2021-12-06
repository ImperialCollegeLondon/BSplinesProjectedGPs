functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2, real inv_tau_squared){
    return - inv_tau_squared * 0.5 * dot_self(phi[node1] - phi[node2])
      + (N-1) * 0.5 * log(inv_tau_squared)
      + normal_lpdf(sum(phi) | 0, 0.001 * N);
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
  
  // ICAR model
  int<lower=1> K;
  int<lower=1> N_edges;
  int<lower=1, upper=K> node1[N_edges];
  int<lower=1, upper=K> node2[N_edges];
}

parameters {
  vector[K] beta_raw; 
  real<lower=0> nu_unscaled;
  real<lower=0> tau;
}

transformed parameters {
  real<lower=0> inv_tau_squared = 1/(tau^2);
  real<lower=0> nu = (1/nu_unscaled)^2;
  real<lower=0> theta = (1 / nu);
  matrix[num_basis_rows,num_basis_columns] beta = to_matrix(beta_raw, num_basis_rows,num_basis_columns); 
  matrix[n,m] f = exp( (BASIS_ROWS') * beta * BASIS_COLUMNS );
  matrix[n,m] alpha = f / nu;
}

model {
  nu_unscaled ~ exponential(1);
  
  beta_raw ~ icar_normal_lpdf(K, node1, node2, inv_tau_squared);
  
  for(k in 1:N){
      y[k] ~ neg_binomial(alpha[coordinates[k,1],coordinates[k,2]], theta);
  }
}

generated quantities {
 real log_lik[N];
  
  for(k in 1:N)
    log_lik[k] = neg_binomial_lpmf(y[k] | alpha[coordinates[k,1],coordinates[k,2]], theta);
}






