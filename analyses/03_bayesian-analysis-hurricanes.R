# in this script, we implement a bayesian regression model using our cyclone
# data from the HURDAT2 data set
# we are interested in the number of landfalls per year

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
library(readr)


# set colour scheme for Stan
bayesplot::color_scheme_set("brewer-RdYlBu")

#==============================================================================
# we load our data set
#==============================================================================

hurricane_data <- read.csv(here("data", "derived", "hurricane-data.csv"))

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
}" 

#-------------------------------------------------------------------------------
# define data in format needed for Stan model
#-------------------------------------------------------------------------------


landfalls <- hurricane_data %>%
  filter(RECORD_ID == "L") %>%
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

#-------------------------------------------------------------------------------
# compile model
#-------------------------------------------------------------------------------

logpoi_hsgp_model_filename <- cmdstanr::write_stan_file(
  gsub('\t',' ', logpoi_hsgp_model_text),
  dir = here("analyses", "stan-models"),
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
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 4e3,
  refresh = 500,
  save_warmup = TRUE)

# save output to RDS
logpoi_hsgp_model_fit$save_object(file = file.path(here("analyses", 
                                                        "stan-models", 
                                                        "hsgp_model_cmdstanr.rds")))


#-------------------------------------------------------------------------------
# check mixing and convergence
#-------------------------------------------------------------------------------

model_summary <- logpoi_hsgp_model_fit$summary(
  variables = c("beta_0", "beta","z", "beta_sd", "gp_lengthscale", "gp_sigma"),
  posterior::default_summary_measures(),
  posterior::default_convergence_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975)))

# sort by smallest ess_bulk
model_summary  <- as.data.table(model_summary )
model_summary  <- model_summary[order(ess_bulk),]

# plot table
kbl(model_summary, caption = 'Model diagnostics', longtable = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                font_size = 12) 

# table of 5 worst performing parameters
t1 <- kbl(subset(model_summary[1:5,], 
    select = c(variable, rhat, ess_bulk, ess_tail)),
    caption = 'Model diagnostics for 5 worst performing parameters', 
    longtable = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                font_size = 12) 

t1

save_kable(t1, file = here("outputs", 
                           "bayesian-analysis-hurricanes", 
                           "model-diagnostics-landfalls-nbr.html"))

# use webshot to capture the html table as a pdf
webshot(here("outputs", 
             "bayesian-analysis-hurricanes", 
             "model-diagnostics-landfalls-nbr.html"), 
        here("outputs", 
             "bayesian-analysis-hurricanes", 
             "model-diagnostics-landfalls-nbr.pdf"))

# plot traces of parameter with smallest ess_bulk
# extract samples
model_draws <- logpoi_hsgp_model_fit$draws(
  variables = c('lp__',
                "beta_0", 
                "beta", "z", 
                "beta_sd", 
                "gp_lengthscale", 
                "gp_sigma"),
  inc_warmup = TRUE, format = "draws_array")

# make trace plots of worst performing parameter
model1_worst_var <- model_summary$variable[ which.min(model_summary$ess_bulk)]

bayesplot:::mcmc_trace(model_draws,  
                            pars = model1_worst_var, 
                            n_warmup = 500, 
                            facet_args = list(nrow = 1))+ 
  theme_bw() + 
  coord_cartesian(ylim = c(-0.5, 10))

#-------------------------------------------------------------------------------
# pairwise posterior geometry of 4 worst performing parameters
#-------------------------------------------------------------------------------

worst_param <- model_summary[1:4,]$variable
model_draws <- logpoi_hsgp_model_fit$draws(
  variables = worst_param,inc_warmup = FALSE,
  format = "draws_array")

bayesplot::color_scheme_set('viridisC')
bayesplot::mcmc_pairs(model_draws, 
                      pars = worst_param, 
                      diag_fun = "dens", 
                      off_diag_fun = "hex") 

# no tight hyperplanes and identifiability issues

#-------------------------------------------------------------------------------
# posterior predictive checks
#-------------------------------------------------------------------------------

# extract Monte Carlo samples of log_lambda vector
model <- list()
model$log_lambda <- logpoi_hsgp_model_fit$draws(variables = "log_lambda", 
                                             inc_warmup = FALSE,
                                             format = "draws_df")

# MELT log_lambda into long format, and link to all_id in the data
model$log_lambda <- data.table::melt(as.data.table(model$log_lambda), 
                                     id.vars = c('.chain','.iteration','.draw'))
set(model$log_lambda, 
    NULL, 
    'ALL_ID', 
    gsub('log_lambda\\[([0-9]+)\\]','\\1',as.character(model$log_lambda$variable)))
set(model$log_lambda, NULL, 'ALL_ID', as.integer(model$log_lambda$ALL_ID))
setnames(model$log_lambda, c('value'), c('log_lambda'))
set(model$log_lambda, NULL, 'variable', NULL)
logpoi_hsgp_model_lambda <- model$log_lambda

# clean up so we don't have millions of data points stored multiple times
model <- NULL
gc()

# make posterior predictions, conditional on the joint posterior for log_lambda
set.seed(42L)
logpoi_hsgp_model_lambda[ , post_pred := rpois( nrow(logpoi_hsgp_model_lambda), 
                                                exp(log_lambda))]

# create median and 95\% credible intervals for variables 
logpoi_hsgp_model_summary <- logpoi_hsgp_model_lambda[,list(
  summary_value = quantile(post_pred, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
  summary_name = c('q_lower','iqr_lower', 'median','iqr_upper', 'q_upper')),
  by = 'ALL_ID']

logpoi_hsgp_model_summary<- data.table::dcast(logpoi_hsgp_model_summary,
                                                 ALL_ID ~ summary_name, 
                                              value.var = 'summary_value')

month_effects <- subset(landfalls[1:300], 
                        select = c('MONTH', 'LANDFALL_COUNT', 'YEAR', 'ALL_ID'))

model_summary <- merge(month_effects,logpoi_hsgp_model_summary, by = 'ALL_ID')


# make posterior predictive check
model_summary[, IN_PPI := LANDFALL_COUNT >= q_lower & LANDFALL_COUNT <= q_upper]
cat(round(model_summary[, mean( as.numeric( IN_PPI ) )]*100, 1),
    "% of all 300 observed values are within our model's 95% posterior predictive intervals.")

# reorder months in chronological order
model_summary$MONTH <- factor(model_summary$MONTH, 
                              levels = month.name)  

# plot posterior predictive check for each year and each race & ethnicity
p1 <- ggplot(model_summary, aes(x = as.factor(YEAR))) + 
  geom_boxplot( aes( group = YEAR, 
                     ymin = q_lower, 
                     lower = iqr_lower,
                     middle = median, 
                     upper = iqr_upper, 
                     ymax = q_upper), 
                stat = 'identity') +
  geom_point( aes(y = LANDFALL_COUNT, colour = IN_PPI ) ) +
  scale_y_continuous() + 
  ggsci::scale_color_npg() +
  labs(x = 'Year', 
       y = 'Number of landfalls', 
       colour = 'within\n95% posterior\nprediction\ninterval') +
  facet_grid(MONTH ~ ., labeller = label_wrap_gen(20)) + 
  theme_bw()

p1

p1 <- ggsave(here("outputs", "bayesian-analysis-hurricanes","post-pred-checks.pdf"), 
            plot = p1, 
            height = 10, 
            width = 10, 
            dpi=300)
                                  

#-------------------------------------------------------------------------------
# forecasts
#-------------------------------------------------------------------------------

model <- list()
model$log_lambda_star <- logpoi_hsgp_model_fit$draws(variables = "log_lambda_star", 
                                                     inc_warmup = FALSE, 
                                                     format = "draws_df")

# MELT log_lambda_star into long format, and link to pred_idx in the data
model$log_lambda_star <- data.table::melt(as.data.table(model$log_lambda_star), 
                                          id.vars = c('.chain','.iteration','.draw'))
set(model$log_lambda_star, 
    NULL, 
    'PRED_ID',
    gsub('log_lambda_star\\[([0-9]+)\\]','\\1',as.character(model$log_lambda_star$variable)))
set(model$log_lambda_star, NULL, 'PRED_ID', as.integer(model$log_lambda_star$PRED_ID))
setnames(model$log_lambda_star, c('value'), c('log_lambda_star'))
set(model$log_lambda_star, NULL, 'variable', NULL)
logpoi_hsgp_model_lambda_star <- model$log_lambda_star

# clean up so we don't have millions of data points stored multiple times
model <- NULL
gc()

# make posterior predictions for each data point, conditional on the 
# joint posterior for log_lambda
set.seed(42L)
logpoi_hsgp_model_lambda_star[,post_pred := rpois( nrow(logpoi_hsgp_model_lambda_star), 
                                                   exp(log_lambda_star))]

# create median and 95\% credible intervals for variables 
# using AGGREGATION over .draws
logpoi_hsgp_model_star_summary <- logpoi_hsgp_model_lambda_star[,list( 
  summary_value = quantile(post_pred, prob =c(0.025,  0.5, 0.975)),
  summary_name = c('2.5% quantile','median','97.5% quantile')),
  by = 'PRED_ID']

# reshape using CASTING
logpoi_hsgp_model_star_summary<- data.table::dcast(
  logpoi_hsgp_model_star_summary, 
  PRED_ID ~ summary_name, 
  value.var = 'summary_value')

# merge data using INNER JOIN
forecast_summary <- subset(landfalls[301:360], select = c('MONTH', 'YEAR', 'PRED_ID'))
forecast_summary <- merge(forecast_summary,logpoi_hsgp_model_star_summary, by = 'PRED_ID')
set(forecast_summary, NULL, 'PRED_ID', NULL)

# reorder months in chronological order
forecast_summary$MONTH <- factor(forecast_summary$MONTH, levels = month.name)  
forecast_summary <- forecast_summary[order(forecast_summary$YEAR, forecast_summary$MONTH), ]

# plot table
kbl(forecast_summary, 
    caption = 'Forecasted number of landfalls per month in the next 5 years', 
    longtable = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12)


# create median and 95\% credible intervals for variables 
# using AGGREGATION over .draws for total monthly forecasts over 2025-2029
total_monthly_draws <- left_join(logpoi_hsgp_model_lambda_star,
                                   subset(landfalls[301:360], select = c(PRED_ID,MONTH)), 
                                   by = "PRED_ID" )

monthly_sums <- total_monthly_draws %>%
  group_by(.draw, MONTH) %>%
  summarise(month_total = sum(post_pred), .groups = "drop")

monthly_sums <- as.data.table(monthly_sums)

total_monthly_summary <- monthly_sums[,list( 
  summary_value = quantile(month_total, prob =c(0.025, 0.25,  0.5, 0.75, 0.975)),
  summary_name = c('2.5% quantile', 'IQR_lower','median', 'IQR_upper', '97.5% quantile')),
  by = 'MONTH']

total_monthly_summary <- data.table::dcast(
  total_monthly_summary, 
  MONTH ~ summary_name, 
  value.var = 'summary_value')

# reorder months in chronological order
total_monthly_summary$MONTH <- factor(total_monthly_summary$MONTH, levels = month.name)  
total_monthly_summary <- total_monthly_summary[order(total_monthly_summary$MONTH), ]

t2 <- kbl(total_monthly_summary, 
    caption = 'Total forecasted number of landfalls per month in the next 5 years', 
    longtable = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12)

t2

save_kable(t2, file = here("outputs", 
                           "bayesian-analysis-hurricanes", 
                           "forecasts-landfalls.html"))

# use webshot to capture the html table as a pdf
webshot(here("outputs", 
             "bayesian-analysis-hurricanes", 
             "forecasts-landfalls.html"), 
        here("outputs", 
             "bayesian-analysis-hurricanes", 
             "forecasts-landfalls.pdf"))


#==============================================================================
# posterior predictive density for number of landfalls 2025
#==============================================================================

forecast_2025 <- subset(landfalls[301:360], select = c('MONTH', 'YEAR', 'PRED_ID')) %>%
  filter(YEAR == 2025)


forecast_2025 <- merge(logpoi_hsgp_model_lambda_star,
                                    forecast_2025,
                                    by = "PRED_ID")

# reorder months in chronological order
forecast_2025$MONTH <- factor(forecast_2025$MONTH, levels = month.name)  
forecast_2025 <- forecast_2025[order(forecast_2025$MONTH), ]

p2 <- ggplot(forecast_2025, aes(x = post_pred)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "steelblue", alpha = 0.7, color = "black") +
  facet_wrap(~ MONTH, ncol = 3) +  # Create a separate plot for each month
  theme_bw() +
  labs(title = "Posterior Density Plots for 2025 Landfall Predictions", 
       x = "Number of Landfalls, L", y = "P(L)") +
  theme(legend.position = "none")

p2

ggsave(here("outputs", "bayesian-analysis-hurricanes", "landfall-density-plots.pdf"), 
       plot =p2,
       height = 10, 
       width = 10, 
       dpi=300)
