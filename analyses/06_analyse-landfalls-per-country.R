# In this script, we analyse the number of landfalls per country and year in 
# our hurricane data.
# we use a Bayesian Poisson regression model for the top 10 countries with the most
# landfalls since 2000

#==============================================================================
# load packages and functions
#==============================================================================

library(data.table) 
library(ggplot2)  
library(ggsci)  
library(hexbin)  
library(bayesplot) 
library(kableExtra)  
library(cmdstanr) 
library(webshot2) 
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

landfall_data <- read.csv(here("data", "derived", "landfall-data.csv"))

#==============================================================================
# Run poisson regression model on data stratified by
# country and year since 2000
#==============================================================================

# define poisson regression model 
logpoi_country_model_text <- "
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
"
#-------------------------------------------------------------------------------
# define data in format needed for Stan model
#-------------------------------------------------------------------------------

# we count the number of landfalls per country to find countries with the
# most landfalls since 2000

landfalls_per_location <- landfall_data %>%
  group_by(LOCATION) %>%
  summarise(LANDFALL_NBR = n(), .groups = "drop")

# find 10 countries with the most landfalls and save these in a dataset
top_ten_landfall_locations <- landfalls_per_location %>%
  arrange(desc(LANDFALL_NBR))

top_ten_landfall_locations <- top_ten_landfall_locations[1:10,]


landfall_freq <- landfall_data %>%
  filter(LOCATION %in% top_ten_landfall_locations$LOCATION) %>%
  group_by(YEAR, LOCATION) %>%
  summarise(LANDFALL_NBR = n(), .groups = "drop")

all_combinations <- CJ(YEAR = 2000:2024, LOCATION = top_ten_landfall_locations$LOCATION)

# we keep all possible combinations and we fill in data set with landfalls data
landfall_freq <- merge(all_combinations, landfall_freq, by = c('LOCATION','YEAR'), all.x = TRUE)
# we replace NA values (if any) with 0 
set(landfall_freq, landfall_freq[, which(is.na(LANDFALL_NBR))], 'LANDFALL_NBR', 0)

# convert our data sets to data.table objects
landfall_freq <- as.data.table(landfall_freq)
top_ten_landfall_locations <- as.data.table(top_ten_landfall_locations)
# first define integer IDs
landfall_freq$OBS_ID = 1:nrow(landfall_freq)
top_ten_landfall_locations <- top_ten_landfall_locations[order(LOCATION),]
top_ten_landfall_locations[, LOC_ID := 1:nrow(top_ten_landfall_locations)]
landfall_freq <- merge(landfall_freq, 
                       subset(top_ten_landfall_locations, select = c(LOCATION, LOC_ID)), 
                              by = 'LOCATION')

# now we make one hot design matrix for predictors
predictors <- data.table::dcast(landfall_freq, 
                                OBS_ID ~ LOC_ID, 
                                value.var = 'LANDFALL_NBR', 
                                fun.aggregate = length)
landfall_freq <- merge(landfall_freq, predictors, by = 'OBS_ID')

#-------------------------------------------------------------------------------
# compile model
#-------------------------------------------------------------------------------

logpoi_country_model_filename <- cmdstanr::write_stan_file(
  gsub('\t',' ', logpoi_country_model_text),
  dir = here("outputs", "stan-models"),
  basename = NULL,
  force_overwrite = FALSE,
  hash_salt = "")

# compile Stan model
logpoi_country_model_compiled <- cmdstanr::cmdstan_model(logpoi_country_model_filename)

# define our Stan data
stan_data <- list()
stan_data$N <- nrow(landfall_freq)
stan_data$P <- 10
stan_data$X <- subset(landfall_freq, select = -c(OBS_ID, 
                                                 LOCATION, 
                                                 YEAR, 
                                                 LANDFALL_NBR,
                                                 LOC_ID))
stan_data$X <- unname(as.matrix(stan_data$X))
stan_data$y <- landfall_freq$LANDFALL_NBR


# sample from joint posterior
logpoi_country_model_fit <- logpoi_country_model_compiled$sample(
  data = stan_data,
  seed = 42,
  chains = 3,
  parallel_chains = 3,
  iter_warmup = 500,
  iter_sampling = 6e3,
  refresh = 500,
  save_warmup = TRUE)

# save output to RDS
logpoi_country_model_fit$save_object(file = file.path(here("outputs",
                                                           "stan-models", 
                                                           "country-model-cmdstanr.rds")))


#-------------------------------------------------------------------------------
# check mixing and convergence of model parameters
#-------------------------------------------------------------------------------

# prints table of all parameter convergence diagnostics, a trace plot of worst 
# performing parameter and pairwise plot of 4 worst performing parameters, then
# saves table of 5 worst parameters (in terms of their convergence diagnostics)
model_diagnostics <- check_model_diagnostics(logpoi_country_model_fit,
                                       c("beta_0", "beta", "beta_sd"))

model_diagnostics

save_kable(model_diagnostics, file = here("outputs", 
                           "bayesian-analysis-country-freq",
                           "country-model-diagnostics.html"))

# use webshot to capture the html table as a pdf
webshot(here("outputs", 
             "bayesian-analysis-country-freq",
             "country-model-diagnostics.html"), 
        here("outputs", 
             "bayesian-analysis-country-freq",
             "country-model-diagnostics.pdf"))

#-------------------------------------------------------------------------------
# posterior predictive checks
#-------------------------------------------------------------------------------

# plot posterior predictive check for each location and each year
post_checks <-  check_posterior_predictions(logpoi_country_model_fit, 
                                   'log_lambda', 
                                   landfall_freq,
                                   'OBS_ID', 
                                   "LANDFALL_NBR", 
                                   'LOCATION', 
                                   'YEAR')
post_checks

post_checks <- ggsave(here("outputs", 
                  "bayesian-analysis-country-freq",
                  "country-post-pred-checks.pdf"), 
             plot = post_checks, 
             height = 10, 
             width = 10, 
             dpi=300)
