data {
  int<lower=1> N;                           // number of studies
  array[N] int<lower=0> x;                  // ABR cases
  array[N] int<lower=1> n;                  // number of patients 
}
// parameters {
//   real mu;                                  // population mean log-odds
//   real<lower=0> tau;                        // population SD
//   vector[N] theta;                          // study-level log-odds
// }
// model {
//   mu ~ normal(-2.5, 1);
//   tau ~ cauchy(0.0, 1);
//   theta ~ normal(mu, tau);
//   x ~ binomial_logit(n, theta);
// }
parameters {
  real mu;             // population mean log-odds
  real<lower=0> tau;   // population SD
  vector[N] z;         // standardized study effects
}
transformed parameters {
  vector[N] theta = mu + tau * z;
}
model {
  z ~ std_normal(); // implies theta ~ normal(mu, tau)
  // mu ~ normal(-2.5, 1);      // note, logit scale
  mu ~ normal(-2.5, 1);
  tau ~ cauchy(0, 1);
  x ~ binomial_logit(n, theta);
}
generated quantities {
  vector[N] p = inv_logit(theta);
  real p_mean_unweighted = mean(p);
  real p_mean_weighted = dot_product(to_vector(n), p) / sum(n);
}
