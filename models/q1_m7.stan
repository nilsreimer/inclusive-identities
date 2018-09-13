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
  
  for (i in 1:N)
    y[i] ~ bernoulli_logit(alpha[i]);
}
generated quantities {
  vector<lower = 0, upper = 1>[3] p_k[K];
  vector[N] log_lik;

  for (k in 1:K) {
    if (k <= 4) {
      p_k[k,1] = inv_logit(b_0 + b_k[k]*sigma_k);
      p_k[k,2] = inv_logit(b_0 + b_k[k]*sigma_k);
      p_k[k,3] = inv_logit(b_0 + b_k[k]*sigma_k + b_scst_k[k]*sigma_scst_k);      
    } else {
      p_k[k,1] = inv_logit(b_0 + b_k[k]*sigma_k);
      p_k[k,2] = inv_logit(b_0 + b_k[k]*sigma_k);
      p_k[k,3] = inv_logit(b_0 + b_k[k]*sigma_k);
    }
  }

  for (i in 1:N)
    log_lik[i] = bernoulli_logit_lpmf(y[i] | alpha[i]);
}
