// prametrization E[y] = mu, Var[y] = mu (1 + nu)

data {
  int<lower=1> M; // number of samples
  int<lower=1> n;
  int<lower=1> m;
  matrix[n,m] alpha[M]; // mean
  real<lower=0> theta[M]; // over dispersion
}

generated quantities {
  int y_hat[M, n, m];
  
  for(k in 1:M){
    for(i in 1:n){
      y_hat[k,i,] = neg_binomial_rng(alpha[k][i,], theta[k]);
    }
  }
}
