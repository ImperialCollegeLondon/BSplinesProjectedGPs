data {
  int<lower=1> n;
  int<lower=1> m;
  matrix[n,m] y;
}

parameters {
  real<lower=0> alpha_gp;
  matrix[n,m] f;
  real<lower=0> sigma;
}

model {
  alpha_gp ~ cauchy(0,1);
  sigma ~ cauchy(0,1);
  
  for(i in 1:n){
    for(j in 1:m){
        f[i,j] ~ normal(0, alpha_gp);
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


