data {
  // Numbers
  int<lower = 0> N_row;
  int<lower = 0> N_col;
  int<lower = 0> N_obs;
  int<lower = 0> N_mis;
  // Indices
  int<lower = 1, upper = (N_obs + N_mis)> ii_obs[N_obs];
  int<lower = 1, upper = (N_obs + N_mis)> ii_mis[N_mis];
  int<lower = 1, upper = N_col> ii_item[N_obs];
  // Summaries
  vector[N_col] x_mean;
  vector[N_col] x_sd;
  // Observations
  vector[N_obs] x_obs;
  int<lower = 0, upper = 1> x_obc[N_row];
  int<lower = 0, upper = 1> x_scst[N_row];
}
transformed data {
  vector[N_obs] x_std;
  for (i in 1:N_obs)
    x_std[i] = ( x_obs[i] - x_mean[ii_item[i]] ) / x_sd[ii_item[i]];
}
parameters {
  vector[N_col] b_0;
  vector[N_col] b_obc;
  vector[N_col] b_scst;
  vector<lower = 0>[N_col] x_sigma;
  cholesky_factor_corr[N_col] Lcorr;
  
  // Missing: Parameters
  vector[N_mis] x_mis;
}
transformed parameters {
  vector[N_col * N_row] x;
  vector[N_col] X[N_row];
  vector[N_col] Mu[N_row];
  
  // Missing: Merging
  x[ii_obs] = x_std;
  x[ii_mis] = x_mis;
  for (n in 1:N_row) 
    for (m in 1:N_col)
      X[n, m] = x[(m - 1) * N_row + n];
      
  // Missing: Regression
  for (n in 1:N_row)
    for (m in 1:N_col)
      Mu[n,m] = b_0[m] + b_obc[m]*x_obc[n] + b_scst[m]*x_scst[n]; 
}
model {
  // Priors
  b_0 ~ normal(0, 1);
  b_obc ~ normal(0, 1);
  b_scst ~ normal(0, 1);
  x_sigma ~ cauchy(0, 1);
  
  // Missing: Imputation
  X ~ multi_normal_cholesky(Mu, diag_pre_multiply(x_sigma, Lcorr));
  Lcorr ~ lkj_corr_cholesky(1);
}
generated quantities {
  matrix[N_col, N_col] Rho = multiply_lower_tri_self_transpose(Lcorr);
  vector[N_col * N_row] x_imp;
  
  // Missing: Rescale
  for (n in 1:N_row)
    for (m in 1:N_col)
      x_imp[(m - 1) * N_row + n] = X[n, m] * x_sd[m] + x_mean[m];
}
