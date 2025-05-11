# our posterior predictive checks for our HSGP model displayed
# small yearly changes so in this script we design a simpler model, investigating
# monthly effects only, from count data stratified by month and year

#==============================================================================
# load packages and functions
#==============================================================================

library(data.table)  # for data mangling
library(ggplot2)  # for plotting
library(ggsci)  # for plotting colors
library(hexbin)  # for plotting pair plots
library(bayesplot)  # for plotting Stan outputs
library(kableExtra)  # for tables
library(cmdstanr)  # for Stan
library(webshot2) # for saving kbl tables
library(here)
library(dplyr)
library(magrittr)
library(tidyr)

source(here("src", "check-model-diagnostics.R"))
source(here("src", "check-posterior-predictions.R"))

# set colour scheme for Stan
bayesplot::color_scheme_set("brewer-RdYlBu")

#==============================================================================
# we load our data set
#==============================================================================

landfalls <- as.data.table(read.csv(here("data", 
                                         "derived", 
                                         "stan-landfall-monthly-freq.csv")))

#==============================================================================
# Run poisson regression model without year effect on data stratified by
# year and month
#==============================================================================

# define poisson regression model 
logpoi_no_year_effect_model_text <- "
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
"

#-------------------------------------------------------------------------------
# compile model
#-------------------------------------------------------------------------------

logpoi_no_year_effect_model_filename <- cmdstanr::write_stan_file(
  gsub('\t',' ', logpoi_no_year_effect_model_text),
  dir = here("outputs", "stan-models"),
  basename = NULL,
  force_overwrite = FALSE,
  hash_salt = "")

# compile Stan model
logpoi_no_year_effect_model_compiled <- cmdstanr::cmdstan_model(logpoi_no_year_effect_model_filename)

# define our maps from 1 to 60 race & ethnicity and years
stan_data <- list()
stan_data$N <- nrow( subset(landfalls, !is.na(OBS_ID)) )
stan_data$y <- landfalls[!is.na(OBS_ID), LANDFALL_COUNT]
stan_data$N_star <- nrow(subset(landfalls, !is.na(PRED_ID)))
stan_data$month_id <- subset(landfalls, !is.na(OBS_ID))$MONTH_ID
stan_data$month_id <- as.vector(unname(as.matrix(stan_data$month_id)))
stan_data$month_id_star <- subset(landfalls, !is.na(PRED_ID))$MONTH_ID
stan_data$month_id_star <- as.vector(unname(as.matrix(stan_data$month_id_star)))
stan_data$M <- 12

# sample from joint posterior
logpoi_no_year_effect_model_fit <- logpoi_no_year_effect_model_compiled$sample(
  data = stan_data,
  seed = 42,
  chains = 3,
  parallel_chains = 3,
  iter_warmup = 500,
  iter_sampling = 3500,
  refresh = 500,
  save_warmup = TRUE)

# load CmdStan output files into the fitted model object using qs to reduce file size
logpoi_no_year_effect_model_fit$draws() # load posterior draws into the object.
try(logpoi_no_year_effect_model_fit$sampler_diagnostics(), silent = TRUE) # load sampler diagnostics.
try(logpoi_no_year_effect_model_fit$init(), silent = TRUE) # load user-defined initial values.
try(logpoi_no_year_effect_model_fit$profiles(), silent = TRUE) # load profiling samples.

# save the object
qs::qsave(x = logpoi_no_year_effect_model_fit, 
          file = file.path(here("outputs", 
                                "stan-models", 
                                "logpoi_no_year_effect_model_fit.qs")))


#logpoi_no_year_effect_model_fit$save_object(file = file.path(
 # here("outputs", 
       #"stan-models", 
      # "no-year-effect-model-cmdstanr.rds")))


#-------------------------------------------------------------------------------
# check mixing and convergence
#-------------------------------------------------------------------------------

# prints table of all parameter convergence diagnostics, a trace plot of worst 
# performing parameter and pairwise plot of 4 worst performing parameters, then
# saves table of 5 worst parameters (in terms of their convergence diagnostics)
model_diagnostics <- check_model_diagnostics(logpoi_no_year_effect_model_fit,
                                             c("beta", "beta_sd"))

model_diagnostics

save_kable(model_diagnostics, file = here("outputs", 
                           "bayesian-analysis-monthly-freq",
                           "model-no-year-effect", 
                           "no-year-model-diagnostics.html"))

# use webshot to capture the html table as a pdf
webshot(here("outputs",
             "bayesian-analysis-monthly-freq",
             "model-no-year-effect",  
             "no-year-model-diagnostics.html"), 
        here("outputs", 
             "bayesian-analysis-monthly-freq",
             "model-no-year-effect", 
             "no-year-model-diagnostics.png"),
        selector = "table",
        zoom = 2)

#-------------------------------------------------------------------------------
# posterior predictive checks
#-------------------------------------------------------------------------------

post_checks <- check_posterior_predictions(logpoi_no_year_effect_model_fit, 
                                           landfalls[1:stan_data$N],
                                           "ALL_ID",
                                           "LANDFALL_COUNT", 
                                           'MONTH', 
                                           'YEAR')

post_checks

post_checks <- ggsave(here("outputs", 
                  "bayesian-analysis-monthly-freq",
                  "model-no-year-effect", 
                  "no-year-post-pred-checks.pdf"), 
             plot = post_checks, 
             height = 10, 
             width = 10, 
             dpi=600)




