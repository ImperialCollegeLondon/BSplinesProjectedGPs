functions {
  
	//GP covariance function
	vector gp(vector[] x, real sdgp, real lscale, vector zgp) { 
		matrix[size(x), size(x)] cov;
		cov = cov_exp_quad(x, sdgp, lscale) + diag_matrix(rep_vector(1e-10, size(x)));

		return cholesky_decompose(cov) * zgp;
	}
	
}

data{
  int<lower=0> M; // number of regions
  int<lower=0> W; // number of weeks
  int<lower=0> A; // continuous age
  int<lower=0> B; // original age bands
  int<lower=0,upper=B> N_idx_non_missing[M,W];
  int<lower=0,upper=B> N_idx_missing[M,W];
  int<lower=-1,upper=B> idx_non_missing[M,B,W]; // indices non-missing deaths
  int<lower=-1,upper=B> idx_missing[M,B,W]; // indices missing deaths
  vector[1] age[A]; // age continuous
  int deaths[M,B,W]; // cumulative deaths in age band b at time n
  int age_from_state_age_strata[B]; // age from of age band b
  int age_to_state_age_strata[B];// age to of age band b
  int range_censored[2]; // range of the censored data
}

transformed data{
	vector[A] age_1; 				//product in linear component

	for(i in 1:A){
		age_1[i] = age[i][1]; 
		}
}

parameters {
  real<lower=0> rho[M,W];
  real<lower=0> sigma[M,W];
  vector[A] eta[M,W];
  vector[W] nu[M];
  real<lower=0> lambda[M,W];
  real c0[M];	
	real c1[M];
}

transformed parameters {
  vector<lower=0>[W] theta[M];
  matrix[A,W] phi[M];
  matrix[A,W] alpha[M];
  matrix[B,W] alpha_reduced[M];
  
  for(m in 1:M){
    
    theta[m] = nu[m] ./ (1 + nu[m]);
    
    for(w in 1:W)
    {
      phi[m][:,w] = softmax( c0[m] + c1[m]*age_1 + gp(age, sigma[m,w], rho[m,w], eta[m,w]) ); 
  
      alpha[m][:,w] = phi[m][:,w] * lambda[m,w] / nu[m,w];
      
      for(b in 1:B){
        alpha_reduced[m][b,w] = sum(alpha[m][age_from_state_age_strata[b]:age_to_state_age_strata[b], w]);
      }
      
    }
  }


}

model {
  c0 ~ normal(0,5);
	c1 ~ normal(0,5);
	  
  for(m in 1:M){
    nu[m,:] ~ exponential(1);
    rho[m,:] ~ inv_gamma(5, 5);
    sigma[m,:] ~ std_normal();
	  
    for(w in 1:W){
      eta[m,w] ~ std_normal();
      lambda[m,w] ~ exponential(1.0 / sum(deaths[m,idx_non_missing[m,1:N_idx_non_missing[m,w],w],w]));
      
      target += neg_binomial_lpmf(deaths[m,idx_non_missing[m,1:N_idx_non_missing[m,w],w],w] | alpha_reduced[m][idx_non_missing[m,1:N_idx_non_missing[m,w],w], w] , theta[m][w] );
    
      for(i in range_censored[1]:range_censored[2])
        target += neg_binomial_lpmf(i| alpha_reduced[m][idx_missing[m,1:N_idx_missing[m,w],w], w] , theta[m][w] ) ;
    }
  }



}

generated quantities {
  real log_lik[M,W];
  int deaths_predict[M,A,W];
  // int deaths_predict_state_age_strata_non_missing[B,W] = rep_array(0, B, W);
  int deaths_predict_state_age_strata[M,B,W];

  for(m in 1:M){
    for(w in 1:W){
      log_lik[m,w] = neg_binomial_lpmf(deaths[m,idx_non_missing[m][1:N_idx_non_missing[m,w],w],w] | 
                                      alpha_reduced[m][idx_non_missing[m][1:N_idx_non_missing[m,w],w],w] , theta[m][w] );
      for(i in range_censored[1]:range_censored[2])
        log_lik[m,w] += neg_binomial_lpmf(i | alpha_reduced[m][idx_missing[m,:N_idx_missing[m,w],w], w] , theta[m][w] );
      
      deaths_predict[m,:,w] = neg_binomial_rng(alpha[m][:,w], theta[m][w]);
      // deaths_predict_state_age_strata_non_missing[idx_non_missing[1:N_idx_non_missing[w],w],w] = neg_binomial_rng(alpha_reduced[idx_non_missing[1:N_idx_non_missing[w],w], w], theta[w]);
      deaths_predict_state_age_strata[m,:,w] = neg_binomial_rng(alpha_reduced[m][:,w], theta[m][w]);
    }
  }

}


