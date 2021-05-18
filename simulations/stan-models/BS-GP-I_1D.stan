data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
  
  //splines
  int num_basis;
  matrix[num_basis, N] BASIS; 
}
transformed data {
  real delta = 1e-9;
}
parameters {
  row_vector[num_basis] beta; 
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] f;
  f = to_vector( beta*BASIS );
}

model {
  sigma ~ std_normal();
  beta ~ std_normal();

  y ~ normal(f, sigma);
}

generated quantities {
  real y_hat[N] = normal_rng(f, sigma);
}

