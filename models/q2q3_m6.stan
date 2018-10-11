data {
  // Numbers
  int<lower = 1> N_row;
  int<lower = 1> N_col;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 0> N_obs;
  int<lower = 0> N_mis;
  int<lower = 0> N_all_obs;
  // Indices
  int<lower = 1, upper = J> jj[N_row];
  int<lower = 1, upper = J> kk[N_row];
  int<lower = 1, upper = (N_obs + N_mis)> ii_obs[N_obs];
  int<lower = 1, upper = (N_obs + N_mis)> ii_mis[N_mis];
  int<lower = 1, upper = N_row> nn_obs[N_all_obs];
  // Vectors
  int<lower = 0, upper = 1> x_scst[N_row];
  int<lower = 0, upper = 1> x_obc[N_row];
  int<lower = 0, upper = 1> x_q1[N_row];
  vector[N_obs] y_obs;
}
parameters {
  vector[N_col] b_0;
  vector[N_col] b_j[J];
  vector[N_col] b_k[K];
  vector[N_col] b_scst_k[4];
  vector[N_col] b_obc_k[4];
  vector[N_col] b_q1[K];
  vector[N_col] b_q1_scst[4];
  vector[N_col] b_q1_obc[4];
  vector<lower = 0>[N_col] sigma;
  vector<lower = 0>[N_col] sigma_j;
  vector<lower = 0>[N_col] sigma_k;
  vector<lower = 0>[N_col] sigma_scst_k;
  vector<lower = 0>[N_col] sigma_obc_k;
  vector<lower = 0>[N_col] sigma_q1_k;
  vector<lower = 0>[N_col] sigma_q1_scst_k;
  vector<lower = 0>[N_col] sigma_q1_obc_k;
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
      
  // Regression
  for (n in 1:N_row)
    for (m in 1:N_col) {
      Mu[n,m] = b_0[m] + b_j[jj[n],m]*sigma_j[m] + b_k[kk[n],m]*sigma_k[m] + b_q1[kk[n],m]*x_q1[n]*sigma_q1_k[m];
      if (kk[n] <= 4) Mu[n,m] = Mu[n,m] + x_scst[n]*(b_scst_k[kk[n],m]*sigma_scst_k[m] + b_q1_scst[kk[n],m]*x_q1[n]*sigma_q1_scst_k[m]) + x_obc[n]*(b_obc_k[kk[n],m]*sigma_obc_k[m] + b_q1_obc[kk[n],m]*x_q1[n]*sigma_q1_obc_k[m]);
    }
}
model {
  // Priors
  for (m in 1:N_col) {
    b_0[m] ~ normal(0, 1);
    b_j[,m] ~ normal(0, 1);
    b_k[,m] ~ normal(0, 1);
    b_scst_k[,m] ~ normal(0, 1);
    b_obc_k[,m] ~ normal(0, 1);
    b_q1[,m] ~ normal(0, 1);
    b_q1_scst[,m] ~ normal(0, 1);
    b_q1_obc[,m] ~ normal(0, 1);
    sigma[m] ~ cauchy(0, 1);
    sigma_j[m] ~ cauchy(0, 1);
    sigma_k[m] ~ cauchy(0, 1);
    sigma_scst_k[m] ~ cauchy(0, 1);
    sigma_obc_k[m] ~ cauchy(0, 1);
    sigma_q1_k[m] ~ cauchy(0, 1);
    sigma_q1_scst_k[m] ~ cauchy(0, 1);
    sigma_q1_obc_k[m] ~ cauchy(0, 1);
  }
  Lcorr ~ lkj_corr_cholesky(1);

  // Likelihood
  Y ~ multi_normal_cholesky(Mu, diag_pre_multiply(sigma, Lcorr));
}
generated quantities {
  matrix[N_col, N_col] Rho = multiply_lower_tri_self_transpose(Lcorr);
  vector[N_all_obs] log_lik;
  vector[N_col] R2;
  
  for (n in 1:N_all_obs) 
    log_lik[n] = multi_normal_cholesky_lpdf(Y[nn_obs[n],] | Mu[nn_obs[n],], diag_pre_multiply(sigma, Lcorr));
    
  for (m in 1:N_col)
    R2[m] = variance(Mu[,m]) / ( variance(Mu[,m]) + variance(to_vector(Y[,m]) - to_vector(Mu[,m])) );
}
