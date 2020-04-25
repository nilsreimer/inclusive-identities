data {
  int<lower = 1> K;
  int<lower = 0> N;
  vector[K] y[N];
  vector[N] x_obc;
  vector[N] x_scst;
}
parameters {
  vector[K] b_0;
  vector[K] b_obc;
  vector[K] b_scst;
  cholesky_factor_corr[K] L_Omega;
  vector<lower = 0>[K] L_sigma;
}
transformed parameters {
  matrix[K, K] L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  vector[K] mu[N];
  for (n in 1:N) {
    mu[n, 1] = b_0[1] + b_obc[1] * x_obc[n] + b_scst[1] * x_scst[n];  
    mu[n, 2] = b_0[2] + b_obc[2] * x_obc[n] + b_scst[2] * x_scst[n];  
    mu[n, 3] = b_0[3] + b_obc[3] * x_obc[n] + b_scst[3] * x_scst[n];  
  }
}
model {
  b_0 ~ normal(0, 1);  
  b_obc ~ normal(0, 1);  
  b_scst ~ normal(0, 1);  
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
