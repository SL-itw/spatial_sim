data {
  int<lower=1> N;                 // number of areas (299)
  array[N] int<lower=0> y;              // observed counts
  matrix[N, N] Q_star;            // INLA‐scaled CAR precision
}
parameters {
  real               intercept;
  vector[N]          eps_std;     // iid “epsilon” (var=1)
  vector[N]          u_std;       // CAR “u”   (var=1 under Q_star)
  real<lower=0>      tau;         // overall precision
  real<lower=0,upper=1> phi;      // mixing weight
}
transformed parameters {
  vector[N] b;                    // total spatial effect
  // mix them on the log‐rate scale
  b = (sqrt(1 - phi) * eps_std + sqrt(phi) * u_std) / sqrt(tau);
}
model {
  // priors
  intercept ~ normal(0, 5);
  eps_std   ~ normal(0, 1);
  u_std     ~ multi_normal_prec(rep_vector(0, N), Q_star);
  tau       ~ gamma(1, 0.01);     // weakly informative
  phi       ~ beta(2, 2);         // roughly uniform on [0,1]

  // likelihood
  y ~ poisson_log(intercept + b);
}
