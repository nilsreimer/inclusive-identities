data {
  // Numbers
  int<lower = 0> N_row;
  int<lower = 0> N_col;
  int<lower = 0> N_obs;
  int<lower = 0> N_mis;
  // Indices
  int<lower = 1, upper = (N_obs + N_mis)> ii_obs[N_obs];
  int<lower = 1, upper = (N_obs + N_mis)> ii_mis[N_mis];
  // Observations
  vector[N_obs] x_obs;
}
parameters {
  vector[N_col] x_mu;
  vector<lower = 0>[N_col] x_sigma;
  cholesky_factor_corr[N_col] Lcorr;
  
  // Missing: Parameters
  vector[N_mis] x_mis;
}
transformed parameters {
  vector[N_col * N_row] x;
  vector[N_col] X[N_row];
  
  // Missing: Merging
  x[ii_obs] = x_obs;
  x[ii_mis] = x_mis;
  for (n in 1:N_row) 
    for (m in 1:N_col)
      X[n, m] = x[(m - 1) * N_row + n];
}
model {
  // Priors
  x_mu ~ normal(0, 1);
  x_sigma ~ cauchy(0, 1);
  
  // Missing: Imputation
  X ~ multi_normal_cholesky(x_mu, diag_pre_multiply(x_sigma, Lcorr));
  Lcorr ~ lkj_corr_cholesky(1);
}
generated quantities {
  matrix[N_col, N_col] Rho;
  Rho = multiply_lower_tri_self_transpose(Lcorr);
}
