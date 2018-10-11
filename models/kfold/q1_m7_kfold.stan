data {
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N];
  int<lower = 1, upper = K> kk[N];
  int<lower = 0, upper = 1> x_scst[N];
  int<lower = 0, upper = 1> x_gs[N];
  vector[N] x_sdo_d;
  vector[N] x_sdo_e;
  int<lower = 0, upper = 1> y[N];
  // k-fold cross-validation
  int<lower = 1> N_t;
  int<lower = 1> N_h;
  int<lower = 1, upper = N> ii_t[N_t];
  int<lower = 1, upper = N> ii_h[N_h];
}
parameters {
  real b_0;
  real b_sdo_d;
  real b_sdo_e;
  vector[J] b_j;
  vector[K] b_k;
  vector[4] b_scst_k;
  real<lower = 0> sigma_j;
  real<lower = 0> sigma_k;
  real<lower = 0> sigma_scst_k;
}
transformed parameters {
  vector[N] alpha;
  
  for (i in 1:N) {
    alpha[i] = b_0 + b_j[jj[i]]*sigma_j + b_k[kk[i]]*sigma_k + x_gs[i]*(b_sdo_d*x_sdo_d[i] + b_sdo_e*x_sdo_e[i]);
    if (kk[i] <= 4) {
      alpha[i] = alpha[i] + x_scst[i]*(b_scst_k[kk[i]]*sigma_scst_k);
    }
  }
}
model {
  b_0 ~ student_t(2.5, 0, 1);
  b_sdo_d ~ student_t(2.5, 0, 1);
  b_sdo_e ~ student_t(2.5, 0, 1);
  b_j ~ normal(0, 1);
  b_k ~ normal(0, 1);
  b_scst_k ~ normal(0, 1);
  sigma_j ~ cauchy(0, 1);
  sigma_k ~ cauchy(0, 1);
  sigma_scst_k ~ cauchy(0, 1);
  
  for (i in 1:N_t)
    y[ii_t[i]] ~ bernoulli_logit(alpha[ii_t[i]]);
}
generated quantities {
  vector[N_t] log_lik_t;
  vector[N_h] log_lik_h;
  for (i in 1:N_t)
    log_lik_t[i] = bernoulli_logit_lpmf(y[ii_t[i]] | alpha[ii_t[i]]);  
  for (i in 1:N_h)
    log_lik_h[i] = bernoulli_logit_lpmf(y[ii_h[i]] | alpha[ii_h[i]]);
}
