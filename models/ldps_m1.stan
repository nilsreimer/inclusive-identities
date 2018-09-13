data {
  // Numbers
  int<lower = 1> N_row;
  int<lower = 1> N_col;
  int<lower = 0> N_obs;
  int<lower = 0> N_mis;
  int<lower = 0> N_all_obs;
  // Indices
  int<lower = 1, upper = (N_obs + N_mis)> ii_obs[N_obs];
  int<lower = 1, upper = (N_obs + N_mis)> ii_mis[N_mis];
  int<lower = 1, upper = N_row> nn_obs[N_all_obs];
  // Vectors
  int<lower = 0, upper = 1> x_scst[N_row];
  int<lower = 0, upper = 1> x_obc[N_row];
  vector[N_obs] y_obs;
}
parameters {
  vector[N_col] b_0;
  vector[N_col] b_scst;
  vector[N_col] b_obc;
  vector<lower = 0>[N_col] sigma;
  cholesky_factor_corr[N_col] Lcorr;
  
  // Missing: Parameters
  vector[N_mis] y_mis;
}
transformed parameters {
  vector[N_col * N_row] y;
  vector[N_col] Y[N_row];
  vector[N_col] Mu[N_row];
  
  // Missing: Merging
  y[ii_obs] = y_obs;
  y[ii_mis] = y_mis;
  for (n in 1:N_row) 
    for (m in 1:N_col)
      Y[n, m] = y[(m - 1) * N_row + n];
      
  // Missing: Regression
  for (n in 1:N_row) {
    for (m in 1:N_col) {
      Mu[n,m] = b_0[m] + x_scst[n]*b_scst[m] + x_obc[n]*b_obc[m];
    }
  }
}
model {
  // Priors
  b_0 ~ normal(0, 1);
  b_scst ~ normal(0, 1);
  b_obc ~ normal(0, 1);
  sigma ~ cauchy(0, 1);
  Lcorr ~ lkj_corr_cholesky(1);
  
  // Likelihood
  Y ~ multi_normal_cholesky(Mu, diag_pre_multiply(sigma, Lcorr));
}
generated quantities {
  matrix[N_col, N_col] Rho = multiply_lower_tri_self_transpose(Lcorr);
  vector[N_all_obs] log_lik;
  for (n in 1:N_all_obs) 
    log_lik[n] = multi_normal_cholesky_lpdf(Y[nn_obs[n],] | Mu[nn_obs[n],], diag_pre_multiply(sigma, Lcorr));

}
