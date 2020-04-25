data {
  int<lower = 1> K;
  int<lower = 0> N;
  vector[K] y[N];
  vector[N] x_gm;
  vector[N] x_obc;
  vector[N] x_scst;
  vector[N] x_of_gm;
  vector[N] x_of_obc;
  vector[N] x_of_scst;
  vector[N] x_of_muslim;
  vector[N] x_nc_gm;
  vector[N] x_nc_obc;
  vector[N] x_nc_scst;
  vector[N] x_nc_muslim;
}
parameters {
  vector[K] b_0;
  vector[3] b_gm_of;
  vector[3] b_gm_nc;
  vector[K] b_obc;
  vector[4] b_obc_of;
  vector[4] b_obc_nc;
  vector[K] b_scst;
  vector[4] b_scst_of;
  vector[4] b_scst_nc;
  cholesky_factor_corr[K] L_Omega;
  vector<lower = 0>[K] L_sigma;
}
transformed parameters {
  matrix[K, K] L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  vector[K] mu[N];
  for (n in 1:N) {
    mu[n, 1] = b_0[1] + b_obc[1] * x_obc[n] + b_scst[1] * x_scst[n]; // ig
    mu[n, 2] = b_0[2] + (b_obc[2] + b_obc_of[1] * x_of_scst[n] + b_obc_nc[1] * x_nc_scst[n]) * x_obc[n] + (b_scst[2] + b_scst_of[1] * x_of_gm[n] + b_scst_nc[1] * x_nc_gm[n]) * x_scst[n] + (b_gm_of[1] * x_of_scst[n] + b_gm_nc[1] * x_nc_scst[n]) * x_gm[n]; // scst
    mu[n, 3] = b_0[3] + (b_obc[3] + b_obc_of[2] * x_of_gm[n] + b_obc_nc[2] * x_nc_gm[n]) * x_obc[n] + (b_scst[3] + b_scst_of[2] * x_of_obc[n] + b_scst_nc[2] * x_nc_obc[n]) * x_scst[n] + (b_gm_of[2] * x_of_obc[n] + b_gm_nc[2] * x_nc_obc[n]) * x_gm[n]; // obc
    mu[n, 4] = b_0[4] + (b_obc[4] + b_obc_of[3] * x_of_gm[n] + b_obc_nc[3] * x_nc_gm[n]) * x_obc[n] + (b_scst[4] + b_scst_of[3] * x_of_gm[n] + b_scst_nc[3] * x_nc_gm[n]) * x_scst[n]; // gm
    mu[n, 5] = b_0[5] + b_obc[5] * x_obc[n] + b_scst[5] * x_scst[n]; // hindu  
    mu[n, 6] = b_0[6] + (b_obc[6] + b_obc_of[4] * x_of_muslim[n] + b_obc_nc[4] * x_nc_muslim[n]) * x_obc[n] + (b_scst[6] + b_scst_of[4] * x_of_gm[n] + b_scst_nc[4] * x_nc_gm[n]) * x_scst[n] + (b_gm_of[3] * x_of_muslim[n] + b_gm_nc[3] * x_nc_muslim[n]) * x_gm[n]; // muslim
  }
}
model {
  b_0 ~ normal(0, 1);  
  b_gm_of ~ normal(0, 1);
  b_gm_nc ~ normal(0, 1);
  b_obc ~ normal(0, 1);
  b_obc_of ~ normal(0, 1);
  b_obc_nc ~ normal(0, 1);
  b_scst ~ normal(0, 1);
  b_scst_of ~ normal(0, 1);
  b_scst_nc ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 1);
  y ~ multi_normal_cholesky(mu, L_Sigma);
}
generated quantities {
  matrix[K, K] Rho = multiply_lower_tri_self_transpose(L_Omega);
  vector[N] log_lik;
  vector[K] R2;
  for (n in 1:N)
    log_lik[n] = multi_normal_cholesky_lpdf(y[n] | mu[n], L_Sigma);
  for (j in 1:K)
    R2[j] = variance(mu[, j])/(variance(mu[, j]) + L_sigma[j]^2);
}
