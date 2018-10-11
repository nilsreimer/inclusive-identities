data {
  int<lower = 1> N_row;
  int<lower = 1> N_col;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N_row];
  vector[N_col] Y[N_row];
  // k-fold cross-validation
  int<lower = 1> N_t;
  int<lower = 1> N_h;
  int<lower = 1, upper = N_row> ii_t[N_t];
  int<lower = 1, upper = N_row> ii_h[N_h];
}
parameters {
  vector[N_col] b_0;
  vector[N_col] b_j[J];
  vector<lower = 0>[N_col] sigma;
  vector<lower = 0>[N_col] sigma_j;
  cholesky_factor_corr[N_col] Lcorr;
}
transformed parameters {
  vector[N_col] Mu[N_row];
      
  // Regression
  for (n in 1:N_row)
    for (m in 1:N_col)
      Mu[n,m] = b_0[m] + b_j[jj[n],m]*sigma_j[m]; 
}
model {
  // Priors
  for (m in 1:N_col) {
    b_0[m] ~ normal(0, 1);
    b_j[,m] ~ normal(0, 1);
    sigma[m] ~ cauchy(0, 1);
    sigma_j[m] ~ cauchy(0, 1);
  }
  Lcorr ~ lkj_corr_cholesky(1);
  
  // Likelihood
  for (i in 1:N_t)
    Y[ii_t[i],] ~ multi_normal_cholesky(Mu[ii_t[i],], diag_pre_multiply(sigma, Lcorr));
}
generated quantities {
  vector[N_t] log_lik_t;
  vector[N_h] log_lik_h;
  for (i in 1:N_t)
    log_lik_t[i] = multi_normal_cholesky_lpdf(Y[ii_t[i],] | Mu[ii_t[i],], diag_pre_multiply(sigma, Lcorr));  
  for (i in 1:N_h)
    log_lik_h[i] = multi_normal_cholesky_lpdf(Y[ii_h[i],] | Mu[ii_h[i],], diag_pre_multiply(sigma, Lcorr));
}
