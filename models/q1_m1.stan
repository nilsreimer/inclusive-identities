data {
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N];
  int<lower = 1, upper = K> kk[N];
  int<lower = 0, upper = 1> y[N];
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
  
  for (i in 1:N)
    y[i] ~ bernoulli_logit(alpha[i]);
}
generated quantities {
  vector<lower = 0, upper = 1>[K] p_k;
  vector[N] log_lik;
  
  for (k in 1:K)
    p_k[k] = inv_logit(b_0 + b_k[k]*sigma_k);
    
  for (i in 1:N)
    log_lik[i] = bernoulli_logit_lpmf(y[i] | alpha[i]);
}
