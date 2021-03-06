data {
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N];
  int<lower = 1, upper = K> kk[N];
  int<lower = 0, upper = 1> x_scst[N];
  int<lower = 0, upper = 1> x_obc[N];
  int<lower = 0, upper = 1> x_og[N];
  int<lower = 0, upper = 1> y[N];
  // missing data
  int<lower = 0> N_obs;
  int<lower = 0> N_mis;
  int<lower = 1, upper = (N_obs + N_mis)> ii_obs[N_obs];
  int<lower = 1, upper = (N_obs + N_mis)> ii_mis[N_mis];
  vector[N_obs] x_obs;
  vector[N_mis] x_mu;
  vector<lower = 0>[N_mis] x_sigma;
  // k-fold cross-validation
  int<lower = 1> N_t;
  int<lower = 1> N_h;
  int<lower = 1, upper = N> ii_t[N_t];
  int<lower = 1, upper = N> ii_h[N_h];
}
parameters {
  real b_0;
  real b_nc;
  real b_of;
  vector[J] b_j;
  vector[K] b_k;
  real b_scst;
  real b_obc;
  vector[4] b_scst_k;
  vector[4] b_obc_k;
  real<lower = 0> sigma_j;
  real<lower = 0> sigma_k;
  real<lower = 0> sigma_scst_k;
  real<lower = 0> sigma_obc_k;
  vector[N_mis] x_mis;
}
transformed parameters {
  vector[N] alpha;
  vector[N_obs + N_mis] x_imp;
  vector[N] x_nc;
  vector[N] x_of;

  // merge missing and observed responses
  x_imp[ii_obs] = x_obs;
  x_imp[ii_mis] = x_mis;
  x_nc = to_vector(x_imp[(2*N + 1):(3*N)]);
  x_of = to_vector(x_imp[(3*N + 1):(4*N)]);

  // regression
  for (i in 1:N) {
    alpha[i] = b_0 + b_j[jj[i]]*sigma_j + b_k[kk[i]]*sigma_k;
    if (kk[i] <= 4) {
      alpha[i] += x_scst[i]*(b_scst + b_scst_k[kk[i]]*sigma_scst_k) + x_obc[i]*(b_obc + b_obc_k[kk[i]]*sigma_obc_k) + x_og[i]*(b_nc*x_nc[i] + b_of*x_of[i]);
    }
  }
}
model {
  b_0 ~ student_t(2.5, 0, 1);
  b_nc ~ student_t(2.5, 0, 1);
  b_of ~ student_t(2.5, 0, 1);
  b_j ~ normal(0, 1);
  b_k ~ normal(0, 1);
  b_scst ~ student_t(2.5, 0, 1);
  b_scst_k ~ normal(0, 1);
  b_obc ~ student_t(2.5, 0, 1);
  b_obc_k ~ normal(0, 1);
  sigma_j ~ cauchy(0, 1);
  sigma_k ~ cauchy(0, 1);
  sigma_scst_k ~ cauchy(0, 1);
  sigma_obc_k ~ cauchy(0, 1);
  
  // impute missing responses
  x_mis ~ normal(x_mu, x_sigma);
  
  // likelihood
  for (i in ii_t)
    y[i] ~ bernoulli_logit(alpha[i]);
}
generated quantities {
  vector[N_h] log_lik;
  for (i in 1:N_h)
    log_lik[i] = bernoulli_logit_lpmf(y[ii_h[i]] | alpha[ii_h[i]]);
}
