
functions {
  // sum-to-zero parameterisation
  matrix sum2zero_generating_matrix(int K) {
    matrix[K, K] A = diag_matrix(rep_vector(1, K));
    for (i in 1:K - 1) A[K, i] = -1;
    A[K, K] = 0;
    return qr_Q(A)[ , 1:(K - 1)];
  }
}
data{
  int<lower=1> N; // number of observations
  int<lower=0> P; // number of countries
  matrix[N,P] X; //predictors
  array[N] int<lower=0> y; // outcomes
}
transformed data{
  // sum-to-zero matrix
  real s2z_sd_b = inv(sqrt(1. - inv(P)));
  matrix[P, P-1] s2z_Q_b = sum2zero_generating_matrix(P);
}
parameters{
  real beta_0;
  real<lower=0> beta_sd;
  vector[P-1] beta_s2z_m1_rnde;
}
transformed parameters{
  vector[P] beta_nc_rnde = s2z_Q_b * beta_s2z_m1_rnde;
  vector[P] beta = beta_sd * beta_nc_rnde;
  vector[N] log_lambda;
  
  // linear predictor
    log_lambda = beta_0 + X * beta;
}
model{
  // likelihood
  y ~ poisson_log(log_lambda);
  
  // priors
  beta_0 ~ normal(0, 2);
  beta_s2z_m1_rnde ~ normal(0, s2z_sd_b);
  beta_sd ~ cauchy(0, 1);
}

