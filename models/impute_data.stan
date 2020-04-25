data {
  int<lower = 0> N_row;
  int<lower = 0> N_col;
  int<lower = 0> N_obs;
  int<lower = 0> N_mis;
  int<lower = 1, upper = (N_obs + N_mis)> ii_obs[N_obs];
  int<lower = 1, upper = (N_obs + N_mis)> ii_mis[N_mis];
  vector[N_obs] x_obs;
  int<lower = 0> K;
  int<lower = 1, upper = 3> kk[N_row];
}
parameters {
  vector[N_col] x_mu[K];
  vector<lower = 0>[N_col] x_sigma[K];
  cholesky_factor_corr[N_col] L_Omega[K];
  vector[N_mis] x_mis; 
}
transformed parameters {
  vector[N_obs + N_mis] x_imp;
  vector[N_col] X[N_row];
  
  // merge missing and observed responses
  x_imp[ii_obs] = x_obs;
  x_imp[ii_mis] = x_mis;
  for (m in 1:N_col)
    for (n in 1:N_row) 
      X[n, m] = x_imp[(m - 1) * N_row + n];
}
model {
  matrix[N_col, N_col] L_Sigma[K];
  for (k in 1:K) {
    x_mu[k] ~ normal(0, 1);
    L_Omega[k] ~ lkj_corr_cholesky(4);
    x_sigma[k] ~ cauchy(1, 1);
    L_Sigma[k] = diag_pre_multiply(x_sigma[k], L_Omega[k]);
  }
  for (n in 1:N_row)
    X[n] ~ multi_normal_cholesky(x_mu[kk[n]], L_Sigma[kk[n]]);
}
generated quantities {
  matrix[N_col, N_col] Rho[K];
  for (k in 1:K) 
    Rho[k] = multiply_lower_tri_self_transpose(L_Omega[k]);
}
