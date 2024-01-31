data {
  int<lower=1> n;
  vector[n] M;
  real lmin;
  real lmax;
  matrix<lower=0, upper=1>[n, n] W;
  vector[n] y;
}
parameters {
  real<lower=0, upper=1> rho;
}

transformed parameters {
  vector[n] mu;
  vector[n] sigma;
  sigma = rep_vector(1/dot_self(M),n);
  mu = rho * (W * y) ./M;


}
model {
  rho ~ uniform(1/lmin,1/lmax);
  y ~ multi_normal(mu, diag_matrix(sigma));

}
