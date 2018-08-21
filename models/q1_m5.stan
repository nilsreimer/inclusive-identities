data {
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N];
  int<lower = 1, upper = K> kk[N];
  int<lower = 0, upper = 1> x_scst[N];
  int<lower = 0, upper = 1> x_og[N];
  vector[N] x_nc;
  vector[N] x_of;
  int<lower = 0, upper = 1> y[N];
}
parameters {
  real b_0;
  real b_nc;
  real b_of;
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
    alpha[i] = b_0 + b_j[jj[i]]*sigma_j + b_k[kk[i]]*sigma_k;
    if (kk[i] <= 4) {
      alpha[i] = alpha[i] + x_scst[i]*(b_scst_k[kk[i]]*sigma_scst_k) + x_og[i]*(b_nc*x_nc[i] + b_of*x_of[i]);
    }
  }
}
model {
  b_0 ~ student_t(2.5, 0, 1);
  b_nc ~ student_t(2.5, 0, 1);
  b_of ~ student_t(2.5, 0, 1);
  b_j ~ normal(0, 1);
  b_k ~ normal(0, 1);
  b_scst_k ~ normal(0, 1);
  sigma_j ~ cauchy(0, 1);
  sigma_k ~ cauchy(0, 1);
  sigma_scst_k ~ cauchy(0, 1);
  
  for (i in 1:N)
    y[i] ~ bernoulli_logit(alpha[i]);
}
generated quantities {
  vector[3] l_k[K];
  vector[N] log_lik;

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

  for (i in 1:N)
    log_lik[i] = bernoulli_logit_lpmf(y[i] | alpha[i]);
}
