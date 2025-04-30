# in this script, we implement a bayesian regression model using our cyclone
# data from the HURDAT2 data set
# we are interested in the number of landfalls per year and month, 
# investigating both monthly and yearly effects

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

# set colour scheme for Stan
bayesplot::color_scheme_set("brewer-RdYlBu")

source(here("src", "check-model-diagnostics.R"))
source(here("src", "check-posterior-predictions.R"))

#==============================================================================
# we load our data set
#==============================================================================

hurricane_data <- read.csv(here("data", "derived", "hurricane-data.csv"))
landfall_data <- read.csv(here("data", "derived", "landfall-data.csv"))

#===============================================================================
# Bayesian poisson random effects model for number of landfalls per year
#===============================================================================

# we use vectorized calculations, sum to zero constraints and 
# non-centered parameterisations to optimise our algorithm's efficiency
# we model the number of landfalls for each month from count data 
# stratified by month and year with zero-mean HSGP prior for time effect

# define poisson regression model 
logpoi_hsgp_model_text <- "
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
    matrix[M,M-1] s2z_Q_b;
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
    beta_0 ~ normal( 0 , 2 );
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
    array [N] real log_lik;
  
    for(i in 1:N){
      log_lik[i] = poisson_log_lpmf(y[i] | beta_0 + beta[month_id[i]]);
    }
    
    vector[N_star] hsgp_pred;{
      // sample GP at prediction inputs
      hsgp_pred = hsgp_PHI_star * (hsgp_sqrt_spd .* z);
      
      for (i in 1:N_star){
        log_lambda_star[i] = beta_0 + beta[month_id_star[i]] + hsgp_pred[i];}
      // predict
      y_all_pred = poisson_log_rng( append_row(log_lambda, log_lambda_star));
    }
}" 

#-------------------------------------------------------------------------------
# define data in format needed for Stan model
#-------------------------------------------------------------------------------

landfalls <- landfall_data %>%
  group_by(YEAR, MONTH) %>%
  summarise(LANDFALL_COUNT = n(), .groups = "drop")

all_combinations <- CJ(YEAR = 2000:2024, MONTH = 1:12)

# we keep all possible combinations and we fill in data set with landfalls data
landfalls <- merge(all_combinations, landfalls, by = c('MONTH','YEAR'), all.x = TRUE)
# we replace NA values (if any) with 0 
set(landfalls, landfalls[, which(is.na(LANDFALL_COUNT))], 'LANDFALL_COUNT', 0)
# replace month numbers with month names for readability
landfalls <- landfalls %>%
  mutate(MONTH = month.name[MONTH])

# first define integer IDs
landfalls$OBS_ID = 1:nrow(landfalls)

months <- unique(subset(landfalls, select = MONTH))
months$MONTH_ID = 1:nrow(months)
landfalls <- merge(landfalls, months, by = 'MONTH')

# we create a dataset that will contain both observed and forecast attributes
forecasts <- CJ(YEAR = 2000:2029, MONTH = unique(landfalls$MONTH))

forecasts <- left_join(forecasts,
                       unique(subset(landfalls, select = c(MONTH, MONTH_ID))),
                       by = 'MONTH')

landfalls <- left_join(forecasts, landfalls, by = c('YEAR','MONTH', 'MONTH_ID'))

forecasts_id <- subset(landfalls, is.na(OBS_ID))
forecasts_id$PRED_ID = 1:nrow(forecasts_id)
landfalls <- left_join(landfalls, forecasts_id, by = colnames(landfalls))

# define standardised inputs, so off-the-shelf GP priors can be used for fitting
# and for predictions
landfalls$DATE = (landfalls$YEAR - min(landfalls$YEAR)) / (max(landfalls$YEAR) - min(landfalls$YEAR))
set(landfalls, NULL, 'DATE', landfalls[, DATE*2 - 1])
landfalls<- landfalls[order(OBS_ID),]
landfalls[, ALL_ID := 1:nrow(landfalls)]

# we save this data set as we will use it later on, in our model without the year effect
write.csv(landfalls, 
          file = file.path(here("data", "derived"), "stan-landfall-freq-data.csv"), 
          row.names = FALSE)

#-------------------------------------------------------------------------------
# compile model
#-------------------------------------------------------------------------------

logpoi_hsgp_model_filename <- cmdstanr::write_stan_file(
  gsub('\t',' ', logpoi_hsgp_model_text),
  dir = here("outputs", "stan-models"),
  basename = NULL,
  force_overwrite = FALSE,
  hash_salt = "")

# compile Stan model
logpoi_hsgp_model_compiled <- cmdstanr::cmdstan_model(logpoi_hsgp_model_filename)

# define our maps from 1 to 60 race & ethnicity and years
stan_data <- list()
stan_data$N <- nrow( subset(landfalls, !is.na(OBS_ID)) )
stan_data$y <- landfalls[!is.na(OBS_ID), LANDFALL_COUNT]
stan_data$inputs_standardised_fit <- landfalls[!is.na(OBS_ID), DATE]
stan_data$N_star <- nrow(subset(landfalls, !is.na(PRED_ID)))
stan_data$inputs_standardised_star <- landfalls[!is.na(PRED_ID), DATE]
stan_data$hsgp_c <- 1.2
stan_data$hsgp_M <- 30
stan_data$month_id <- subset(landfalls, !is.na(OBS_ID))$MONTH_ID
stan_data$month_id <- as.vector(unname(as.matrix(stan_data$month_id)))
stan_data$month_id_star <- subset(landfalls, !is.na(PRED_ID))$MONTH_ID
stan_data$month_id_star <- as.vector(unname(as.matrix(stan_data$month_id_star)))
stan_data$M <- 12

# sample from joint posterior
logpoi_hsgp_model_fit <- logpoi_hsgp_model_compiled$sample(
  data = stan_data,
  seed = 42,
  chains = 3,
  parallel_chains = 3,
  iter_warmup = 500,
  iter_sampling = 4e3,
  refresh = 500,
  save_warmup = TRUE)

# save output to RDS
logpoi_hsgp_model_fit$save_object(file = file.path(here("outputs", 
                                                        "stan-models", 
                                                        "hsgp-model-cmdstanr.rds")))

#-------------------------------------------------------------------------------
# check mixing and convergence
#-------------------------------------------------------------------------------

# prints table of all parameter convergence diagnostics, a trace plot of worst 
# performing parameter and pairwise plot of 4 worst performing parameters, then
# saves table of 5 worst parameters (in terms of their convergence diagnostics)
model_diagnostics <- check_model_diagnostics(logpoi_hsgp_model_fit,
                                             c("beta_0", 
                                               "beta", 
                                               "beta_sd",
                                               "z",
                                               "gp_lengthscale", 
                                               "gp_sigma"))


model_diagnostics

save_kable(model_diagnostics, file = here("outputs", 
                           "bayesian-analysis-landfall-freq",
                           "model-HSGP", 
                           "HSGP-model-diagnostics.html"))

# use webshot to capture the html table as a pdf
webshot(here("outputs", 
             "bayesian-analysis-landfall-freq",
             "model-HSGP", 
             "HSGP-model-diagnostics.html"), 
        here("outputs",
             "bayesian-analysis-landfall-freq",
             "model-HSGP", 
             "HSGP-model-diagnostics.pdf"))

#-------------------------------------------------------------------------------
# posterior predictive checks
#-------------------------------------------------------------------------------

post_checks <- check_posterior_predictions(logpoi_hsgp_model_fit, 
                                           'log_lambda', 
                                           landfalls[1:300],
                                           "ALL_ID",
                                           "LANDFALL_COUNT", 
                                           'MONTH', 
                                           'YEAR')

post_checks 

post_checks <- ggsave(here("outputs", 
                  "bayesian-analysis-landfall-freq",
                  "model-HSGP",
                  "HSGP-post-pred-checks.pdf"), 
            plot = post_checks, 
            height = 10, 
            width = 10, 
            dpi=300)
                                  