data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
  
  //splines
  int num_basis;
  matrix[num_basis, N] BASIS; 
  
  // GP
  real IDX_BASIS[num_basis];
}

transformed data {
  real delta = 1e-9;
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[num_basis] eta;
}

transformed parameters {
  vector[N] f;
  row_vector[num_basis] GP;
  {
    matrix[num_basis, num_basis] L_K;
    matrix[num_basis, num_basis] K = cov_exp_quad(IDX_BASIS, alpha, rho);

    // diagonal elements
    for (n in 1:num_basis)
      K[n, n] = K[n, n] + delta;

    L_K = cholesky_decompose(K);
    
    GP = to_row_vector( L_K * eta );
    f = to_vector( GP*BASIS );
  }
}

model {
  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  eta ~ std_normal();

  y ~ normal(f, sigma);
}

generated quantities {
  real y_hat[N] = normal_rng(f, sigma);
}
