data {
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N];
  int<lower = 1, upper = K> kk[N];
  int<lower = 0, upper = 1> y[N];
  // k-fold cross-validation
  int<lower = 1> N_t;
  int<lower = 1> N_h;
  int<lower = 1, upper = N> ii_t[N_t];
  int<lower = 1, upper = N> ii_h[N_h];
}
parameters {
  real b_0;
  vector[J] b_j;
  vector[K] b_k;
  real<lower = 0> sigma_j;
  real<lower = 0> sigma_k;
}
transformed parameters {
  vector[N] alpha;
  
  for (i in 1:N)
    alpha[i] = b_0 + b_j[jj[i]]*sigma_j + b_k[kk[i]]*sigma_k;
}
model {
  b_0 ~ student_t(2.5, 0, 1);
  b_j ~ normal(0, 1);
  b_k ~ normal(0, 1);
  sigma_j ~ cauchy(0, 1);
  sigma_k ~ cauchy(0, 1);

  for (i in ii_t)
    y[i] ~ bernoulli_logit(alpha[i]);
}
generated quantities {
  vector[N_h] log_lik;
  for (i in 1:N_h)
    log_lik[i] = bernoulli_logit_lpmf(y[ii_h[i]] | alpha[ii_h[i]]);
}
