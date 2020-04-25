data {
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N];
  int<lower = 1, upper = K> kk[N];
  int<lower = 0, upper = 1> y[N];
  vector[N] x_sim;
  vector[N] x_g1;
  vector[N] x_g2;}
parameters {
  real b_0;
  vector[J] b_j;
  vector[K] b_k;
  real<lower = 0> sigma_j;
  real<lower = 0> sigma_k;
  real b_x;
  real b_g1;
  real b_g2;
  vector[K] b_g1_k;
  vector[K] b_g2_k;
  real<lower = 0> sigma_g1_k;
  real<lower = 0> sigma_g2_k;
}
transformed parameters {
  vector[N] alpha;
  
  for (i in 1:N)
    alpha[i] = b_0 + b_x*x_sim[i] + b_j[jj[i]]*sigma_j + b_k[kk[i]]*sigma_k + x_g1[i]*(b_g1 + b_g1_k[kk[i]]*sigma_g1_k) + x_g2[i]*(b_g2 + b_g2_k[kk[i]]*sigma_g2_k);
}
model {
  // priors
  b_0 ~ student_t(2.5, 0, 1);
  b_j ~ normal(0, 1);
  b_k ~ normal(0, 1);
  sigma_j ~ cauchy(0, 1);
  sigma_k ~ cauchy(0, 1);
  b_x  ~ student_t(2.5, 0, 1);
  b_g1 ~ student_t(2.5, 0, 1);
  b_g2 ~ student_t(2.5, 0, 1);
  b_g1_k ~ normal(0, 1);
  b_g2_k ~ normal(0, 1);
  sigma_g1_k ~ cauchy(0, 1);
  sigma_g2_k ~ cauchy(0, 1);
  
  // likelihood
  for (i in 1:N)
    y[i] ~ bernoulli_logit(alpha[i]);
}
