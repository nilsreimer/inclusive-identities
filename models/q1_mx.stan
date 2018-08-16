data {
  // Numbers
  int<lower = 0> N_row;
  int<lower = 0> N_col;
  int<lower = 0> N_obs;
  int<lower = 0> N_mis;
  int<lower = 1> J;
  int<lower = 1> K;
  // Indices
  int<lower = 1, upper = (N_obs + N_mis)> ii_obs[N_obs];
  int<lower = 1, upper = (N_obs + N_mis)> ii_mis[N_mis];
  int<lower = 1, upper = J> jj[N_row];
  int<lower = 1, upper = K> kk[N_row];
  // Vectors
  int<lower = 0, upper = 1> x_scst[N_row];
  int<lower = 0, upper = 1> x_og[N_row];
  vector[N_obs] x_obs;
  vector[N_mis] x_imp_m;
  vector<lower = 0>[N_mis] x_imp_sd;
  int<lower = 0, upper = 1> y[N_row];
}
parameters {
  real b_0;
  real b_cq;
  real b_pc;
  real b_nc;
  real b_of;
  vector[J] b_j;
  vector[K] b_k;
  vector[4] b_scst_k;
  real<lower = 0> sigma_j;
  real<lower = 0> sigma_k;
  real<lower = 0> sigma_scst_k;
  
  // Missing: Parameters
  vector[N_mis] x_imp;
}
transformed parameters {
  vector[N_col * N_row] x;
  vector[N_col] X[N_row];
  vector[N_row] alpha;
  
  // Missing: Merging
  x[ii_obs] = x_obs;
  x[ii_mis] = x_imp;
  for (n in 1:N_row) 
    for (m in 1:N_col)
      X[n, m] = x[(m - 1) * N_row + n];
  
  // Regression
  for (i in 1:N_row) {
    alpha[i] = b_0 + b_j[jj[i]]*sigma_j + b_k[kk[i]]*sigma_k;
    if (kk[i] <= 4) {
      alpha[i] = alpha[i] + x_scst[i]*(b_scst_k[kk[i]]*sigma_scst_k) + x_og[i]*(b_cq*X[i,1] + b_pc*X[i,2] + b_nc*X[i,3] + b_of*X[i,4]);
    }
  }
}
model {
  b_0 ~ student_t(2.5, 0, 1);
  b_cq ~ student_t(2.5, 0, 1);
  b_pc ~ student_t(2.5, 0, 1);
  b_nc ~ student_t(2.5, 0, 1);
  b_of ~ student_t(2.5, 0, 1);
  b_j ~ normal(0, 1);
  b_k ~ normal(0, 1);
  b_scst_k ~ normal(0, 1);
  sigma_j ~ cauchy(0, 1);
  sigma_k ~ cauchy(0, 1);
  sigma_scst_k ~ cauchy(0, 1);
  
  // Missing: Imputation
  x_imp ~ normal(x_imp_m, x_imp_sd);
  
  for (i in 1:N_row)
    y[i] ~ bernoulli_logit(alpha[i]);
}
generated quantities {
  vector[3] l_k[K];
  vector[N_row] log_lik;

  for (k in 1:K) {
    if (k <= 4) {
      l_k[k,1] = b_0 + b_k[k]*sigma_k;
      l_k[k,2] = b_0 + b_k[k]*sigma_k;
      l_k[k,3] = b_0 + b_k[k]*sigma_k + b_scst_k[k]*sigma_scst_k;      
    } else {
      l_k[k,1] = b_0 + b_k[k]*sigma_k;
      l_k[k,2] = b_0 + b_k[k]*sigma_k;
      l_k[k,3] = b_0 + b_k[k]*sigma_k;
    }
  }

  for (i in 1:N_row)
    log_lik[i] = bernoulli_logit_lpmf(y[i] | alpha[i]);
}
