data {
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N];
  int<lower = 1, upper = K> kk[N];
  int<lower = 0, upper = 1> y[N];
  // simulations
  int<lower = 1> N_sim;
  int<lower = 1> J_sim;
  int<lower = 1> K_sim;
  int<lower = 1, upper = J_sim> jj_sim[N_sim];
  int<lower = 1, upper = K_sim> kk_sim[N_sim];
  vector[N_sim] x_sim;
  vector[N_sim] x_g1;
  vector[N_sim] x_g2;
}
parameters {
  real b_0;
  vector[J] b_j;
  vector[K] b_k;
  real<lower = 0> sigma_j;
  real<lower = 0> sigma_k;
  // simulations
  real b_x;
  real b_g1;
  real b_g2;
  vector[J_sim] b_j_sim;
  vector[K_sim] b_g1_k;
  vector[K_sim] b_g2_k;
  real<lower = 0> sigma_g1_k;
  real<lower = 0> sigma_g2_k;
}
transformed parameters {
  vector[N] alpha;
  
  for (i in 1:N)
    alpha[i] = b_0 + b_j[jj[i]]*sigma_j + b_k[kk[i]]*sigma_k;
}
model {
  // priors
  b_0 ~ student_t(2.5, 0, 1);
  b_j ~ normal(0, 1);
  b_k ~ normal(0, 1);
  sigma_j ~ cauchy(0, 1);
  sigma_k ~ cauchy(0, 1);
  
  // simulations
  b_x  ~ normal(0, 1);
  b_g1 ~ normal(0, 1);
  b_g2 ~ normal(0, 1);
  b_j_sim ~ normal(0, 1);
  b_g1_k ~ normal(0, 1);
  b_g2_k ~ normal(0, 1);
  sigma_g1_k ~ student_t(10, 0, 0.2);
  sigma_g2_k ~ student_t(10, 0, 0.2);
  
  // likelihood
  for (i in 1:N)
    y[i] ~ bernoulli_logit(alpha[i]);
}
generated quantities {
  int<lower = 0, upper = 1> y_sim[N_sim];
  for (i in 1:N_sim)
    y_sim[i] = bernoulli_logit_rng(b_0 + b_x*x_sim[i] + b_j_sim[jj_sim[i]]*sigma_j + b_k[kk_sim[i]]*sigma_k + x_g1[i]*(b_g1 + b_g1_k[kk_sim[i]]*sigma_g1_k) + x_g2[i]*(b_g2 + b_g2_k[kk_sim[i]]*sigma_g2_k));
}
