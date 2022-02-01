functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2, real inv_tau_squared){
    return - inv_tau_squared * 0.5 * dot_self(phi[node1] - phi[node2])
      + (N-1) * 0.5 * log(inv_tau_squared)
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
  
  // ICAR model
  int<lower=1> K;
  int<lower=1> N_edges;
  int<lower=1, upper=K> node1[N_edges];
  int<lower=1, upper=K> node2[N_edges];
}

parameters {
  real<lower=0> sigma;
  real<lower=0> tau;
  vector[K] f_raw; 
}

transformed parameters {
  real<lower=0> inv_tau_squared = 1/(tau^2);
  
  matrix[n,m] f = to_matrix(f_raw, n,m);
}

model {
  sigma ~ cauchy(0,1);
  tau ~ cauchy(0,1);
  
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


