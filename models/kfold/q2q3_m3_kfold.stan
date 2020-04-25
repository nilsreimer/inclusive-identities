data {
  int<lower = 1> N_row;
  int<lower = 1> N_col;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N_row];
  int<lower = 1, upper = K> kk[N_row];
  vector[N_col] Y[N_row];
  int<lower = 0, upper = 1> x_scst[N_row];
  int<lower = 0, upper = 1> x_obc[N_row];
  // k-fold cross-validation
  int<lower = 1> N_t;
  int<lower = 1> N_h;
  int<lower = 1, upper = N_row> ii_t[N_t];
  int<lower = 1, upper = N_row> ii_h[N_h];
}
parameters {
  vector[N_col] b_0;
  vector[N_col] b_scst;
  vector[N_col] b_obc;
  matrix[N_col, J] b_j_z;  
  matrix[K, N_col] b_k_z;
  matrix[4, N_col] b_scst_k_z; 
  matrix[4, N_col] b_obc_k_z; 
  vector<lower = 0>[N_col] sigma;
  vector<lower = 0>[N_col] sigma_j;
  vector<lower = 0>[N_col] sigma_k;
  vector<lower = 0>[N_col] sigma_scst_k;
  vector<lower = 0>[N_col] sigma_obc_k;
  cholesky_factor_corr[N_col] Lcorr;
  cholesky_factor_corr[N_col] Lcorr_j;
}
transformed parameters {
  vector[N_col] Mu[N_row];
  matrix[N_col, N_col] L_Sigma = diag_pre_multiply(sigma, Lcorr);
  
  // varying effects, non-centred Cholesky parameterization
  matrix[J, N_col] b_j = (diag_pre_multiply(sigma_j, Lcorr_j) * b_j_z)';

  // regression
  for (m in 1:N_col) {
    for (n in 1:N_row) {
      Mu[n, m] = b_0[m] + b_j[jj[n], m] + b_k_z[kk[n], m]*sigma_k[m];
      if (kk[n] <= 4) {
        Mu[n, m] += x_scst[n]*(b_scst[m] + b_scst_k_z[kk[n], m]*sigma_scst_k[m]) + x_obc[n]*(b_obc[m] + b_obc_k_z[kk[n], m]*sigma_obc_k[m]);
      }
    }
  }
}
model {
  b_0 ~ normal(0, 1);
  b_scst ~ normal(0, 1);
  b_obc ~ normal(0, 1);
  to_vector(b_j_z) ~ normal(0, 1);
  to_vector(b_k_z) ~ normal(0, 1);
  to_vector(b_scst_k_z) ~ normal(0, 1);
  to_vector(b_obc_k_z) ~ normal(0, 1);
  sigma ~ cauchy(0, 1);
  sigma_j ~ cauchy(0, 1);
  sigma_k ~ cauchy(0, 1);
  sigma_scst_k ~ cauchy(0, 1);
  sigma_obc_k ~ cauchy(0, 1);
  Lcorr ~ lkj_corr_cholesky(4);
  Lcorr_j ~ lkj_corr_cholesky(4);

  // likelihood
  Y[ii_t] ~ multi_normal_cholesky(Mu[ii_t], L_Sigma);
}
generated quantities {
  vector[N_h] log_lik;
  for (i in 1:N_h)
    log_lik[i] = multi_normal_cholesky_lpdf(Y[ii_h[i]] | Mu[ii_h[i]], L_Sigma);
}
