
// for horseshoe prior
functions {
    // square root of a vector (elementwise)
    vector sqrt_vec(vector x) {
        vector[dims(x)[1]] res;
        for (m in 1:dims(x)[1])
        {
            res[m] = sqrt(x[m]);}
        return res;
    }
}
data {
  int<lower=0> N; // number of observations
  int<lower=0> P; // number of predictors
  vector<lower=0>[N] y; // total damage
  matrix[N, P] X; //predictors (wind, duration, total bearing change and interaction effects)
  real<lower=1> nu; // degrees of freedom for the half t-priors
}
parameters {
  real beta_0; // baseline parameter
  real<lower=0> sigma;  
  // auxiliary variables for the variance parameters
  vector[P] z;
  // now the global scale parameter
  // it is computionally much better to write the Cauchy/Student
  // as infinite normal-inverse-gamma mixtures
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  vector<lower=0>[P] r1_local;
  vector<lower=0>[P] r2_local;
}

transformed parameters {
    // global and local variance parameters, and the input weights
    real<lower=0> tau;
    vector<lower=0>[P] lambda;
    vector[P] beta;

    // the normal-inverse-gamma mixture implementation
    tau = r1_global * sqrt(r2_global); 
    // the normal-inverse-gamma mixture implementation
    lambda = r1_local .* sqrt_vec(r2_local); 
    beta = z .* lambda * tau; // non-centered parameterisation
}
model {
    // likelihood
    y ~ lognormal(beta_0 + X*beta, sigma);
    // weakly informative prior for the intercept
    beta_0 ~ normal(0,5);   
    // local scale parameters, nu = 1 corresponds to horseshoe
    z ~ normal(0, 1);
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5*nu, 0.5*nu);
    // global scale parameter
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);
    sigma ~ cauchy(0,1);
}
