functions {
  /**
  * From https://github.com/mbjoseph/CARstan
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}

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
  
  // CAR model
  int N; // n *m
  matrix<lower = 0, upper = 1>[N,N] Adj; // adjacency matrix
  int Adj_n;                // number of adjacent region pairs
}

transformed data {
  int Adj_sparse[Adj_n, 2];   // adjacency pairs
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[N] egv;       // eigenvalues of invsqrtD * Adj * invsqrtD
   { // generate sparse representation for Ajd
  int counter;
  counter = 1;
  // loop over upper triangular part of Adj to identify neighbor pairs
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (Adj[i, j] == 1) {
          Adj_sparse[counter, 1] = i;
          Adj_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:N) D_sparse[i] = sum(Adj[i]);
  {
    vector[N] invsqrtD;  
    for (i in 1:N) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    egv = eigenvalues_sym(quad_form(Adj, diag_matrix(invsqrtD)));
  }
}

parameters {
  real<lower=0> sigma;
  real<lower = 0> tau;
  real<lower = 0, upper = 1> p;
  vector[N] beta_raw; 
}

transformed parameters {
  matrix[num_basis_rows,num_basis_columns] beta = to_matrix(beta_raw, num_basis_rows,num_basis_columns); 
  matrix[n,m] f = (BASIS_ROWS') * beta * BASIS_COLUMNS;
}

model {
  sigma ~ std_normal();

  target += exponential_lpdf(tau | 1);
  
  beta_raw ~ sparse_car(tau, p, Adj_sparse, D_sparse, egv, N, Adj_n);
  
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




