
  functions {
    vector diagSPD_EQ(real alpha, real rho, real L, int M){
      return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * 
      linspaced_vector(M, 1, M)^2);}
    matrix PHI(int N, int M, real L, vector x){
      return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), 
      linspaced_vector(M, 1, M)))/sqrt(L);}
    // sum-to-zero parameterisation
    matrix sum2zero_generating_matrix(int K) {
      matrix[K, K] A = diag_matrix(rep_vector(1, K));
      for (i in 1:K - 1) A[K, i] = -1;
      A[K, K] = 0;
      return qr_Q(A)[ , 1:(K - 1)];}
  }
  data{
      int<lower=1> N; // number of observations
      int<lower=0> M; // number of months
      array[N] int<lower=0> y; // outputs
      array[N] real<lower=0> offs; // log offsets
      array[N] int<lower=1, upper=M> month_id;
      
      // data for GP
      vector[N] inputs_standardised_fit; // unique inputs for GP for fitting
      // additional distinct inputs for prediction
      int<lower=1> N_star;           // number of unique inputs for prediction
      vector[N_star] inputs_standardised_star; // unique inputs for GP prediction
      array[N_star] int<lower=1, upper=M> month_id_star;
      // HSGP arguments
      real<lower=0> hsgp_c;//factor c to determine the boundary value for the HSGP
      int<lower=1> hsgp_M;    // number of basis functions for the HSGP
  }
  transformed data{
    int N_all; 
    real hsgp_L;
    matrix[N, hsgp_M] hsgp_PHI;
    real hsgp_L_star;
    matrix[N_star, hsgp_M] hsgp_PHI_star;
    N_all = N + N_star;
    
    // precompute HSGP basis functions at inputs to fit
    hsgp_L = hsgp_c*max(inputs_standardised_fit);
    hsgp_PHI = PHI(N, hsgp_M, hsgp_L, inputs_standardised_fit);
    // precompute HSGP basis functions at inputs to predict
    hsgp_L_star = hsgp_c*max(inputs_standardised_star);
    hsgp_PHI_star = PHI(N_star, hsgp_M, hsgp_L_star, inputs_standardised_star);
    // beta sum-to-zero
    real s2z_sd_b; 
    matrix[R,R-1] s2z_Q_b;
    s2z_sd_b = inv(sqrt(1. - inv(M)));
    s2z_Q_b = sum2zero_generating_matrix(M);
  }
  parameters{
    // intercept and noise std
    real beta_0;
    real<lower=0> beta_sd;
    vector[M-1] beta_s2z_m1_rnde;
    real<lower=0> gp_lengthscale;
    real<lower=0> gp_sigma;
    vector[hsgp_M] z; 
  }
  transformed parameters{
    vector[N] f;
    vector[hsgp_M] hsgp_sqrt_spd;
    vector[M] beta_nc_rnde;
    vector[M] beta;
    vector[N] log_lambda;
    
    // sum to zero random effects
    beta_nc_rnde = s2z_Q_b * beta_s2z_m1_rnde;
    // non-centered random effect parameterisation 
    beta =  beta_sd * beta_nc_rnde;
    // square root of spectral densities
    hsgp_sqrt_spd = diagSPD_EQ( gp_sigma, gp_lengthscale, hsgp_L, hsgp_M);
    // construct HSGP at inputs for fitting
    f = hsgp_PHI * (hsgp_sqrt_spd .* z);
    // linear predictor with offset
    for (i in 1:N){
      log_lambda[i] = beta_0 + beta[month_id[i]] + f[i];
    }
  }
  model{
    // likelihood
    y ~ poisson_log( log_lambda );
    // priors
    beta_0 ~ normal( 0 , 10 );
    beta_s2z_m1_rnde ~ normal( 0, s2z_sd_b );
    beta_sd ~ cauchy(0,1);
    // priors for GP
    gp_lengthscale ~ inv_gamma( 5, 1 );
    gp_sigma ~ cauchy( 0 , 1 );
    z ~ std_normal();
  }
  generated quantities{
    vector[N_star] log_lambda_star;
    array[N_all] real y_all_pred;
    vector[N_star] hsgp_pred;{
      // sample GP at prediction inputs
      hsgp_pred = hsgp_PHI_star * (hsgp_sqrt_spd .* z);
      for (i in 1:N_star){
        log_lambda_star[i] = beta_0 + beta[month_id_star[i]] + hsgp_pred[i];}
      // predict
      y_all_pred = poisson_log_rng( append_row(log_lambda, log_lambda_star));
    }
}
