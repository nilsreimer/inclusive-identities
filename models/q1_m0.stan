data {
  int<lower = 1> N;
  int<lower = 1> J;
  int<lower = 1> K;
  int<lower = 1, upper = J> jj[N];
  int<lower = 0, upper = 1> y[N];
}
parameters {
  real b_0;
  vector[J] b_j;
  real<lower = 0> sigma_j;
}
transformed parameters {
  vector[N] alpha;
  
  for (i in 1:N)
    alpha[i] = b_0 + b_j[jj[i]]*sigma_j;
}
model {
  b_0 ~ student_t(2.5, 0, 1);
  b_j ~ normal(0, 1);
  sigma_j ~ cauchy(0, 1);

  for (i in 1:N)
    y[i] ~ bernoulli_logit(alpha[i]);
}
generated quantities {
  vector<lower = 0, upper = 1>[3] p_k[K];
  for (k in 1:K) {
    for (l in 1:3) {
      p_k[k,l] = inv_logit(b_0);
    }
  }
}
