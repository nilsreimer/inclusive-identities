data {
  int<lower = 1> K;
  int<lower = 0> N;
  vector[K] y[N];
}
parameters {
  vector[K] b_0;
  cholesky_factor_corr[K] L_Omega;
  vector<lower = 0>[K] L_sigma;
}
transformed parameters {
  matrix[K, K] L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  vector[K] mu[N];
  for (n in 1:N)
    mu[n] = b_0;
}
model {
  b_0 ~ normal(0, 1);  
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
