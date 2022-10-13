// prametrization E[y] = mu, Var[y] = mu (1 + nu)

data {
  int<lower=0> n; // number of samples
  vector[n] mu; // mean
  real<lower=0> nu; // over dispersion
}

transformed data {
  vector[n] alpha = mu / nu;
  real nu_inverse = 1 / nu;
}

generated quantities {
  int y[n] = neg_binomial_rng(alpha, nu_inverse);
}

