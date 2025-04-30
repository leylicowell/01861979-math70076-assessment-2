
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
  int<lower=0> M; // number of months
  array[N] int<lower=0> y; // outcomes
  array[N] int<lower=1, upper=M> month_id;
  
  // prediction data
  int<lower=1> N_star;
  array[N_star] int<lower=1, upper=M> month_id_star;
}
transformed data{
  int N_all;
  
  N_all= N + N_star;
  
  // sum-to-zero matrix
  real s2z_sd_b = inv(sqrt(1. - inv(M)));
  matrix[M, M-1] s2z_Q_b = sum2zero_generating_matrix(M);
}
parameters{
  real<lower=0> beta_sd;
  vector[M-1] beta_s2z_m1_rnde;
}
transformed parameters{
  vector[M] beta_nc_rnde = s2z_Q_b * beta_s2z_m1_rnde;
  vector[M] beta = beta_sd * beta_nc_rnde;
  vector[N] log_lambda;
  
  for (i in 1:N){
    log_lambda[i] = beta[month_id[i]];
  }
}
model{
  // likelihood
  y ~ poisson_log(log_lambda);
  
  // priors
  beta_s2z_m1_rnde ~ normal(0, s2z_sd_b);
  beta_sd ~ cauchy(0, 1);
}
generated quantities{
  vector[N_star] log_lambda_star;
  array[N_all] real y_all_pred;
  array [N] real log_lik;
  
  for(i in 1:N)
  {
    log_lik[i] = poisson_log_lpmf(y[i] | beta[month_id[i]]);
  }
  
  for (i in 1:N_star){
    log_lambda_star[i] = beta[month_id_star[i]];
  }
  y_all_pred = poisson_log_rng(append_row(log_lambda, log_lambda_star));
}

