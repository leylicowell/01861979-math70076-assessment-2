# In this script, we analyse the number of landfalls per country and year in 
# our hurricane data.
# we use a Bayesian Poisson regression model for the top 10 countries with the most
# landfalls since 2000

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

# now make one hot design matrix for predictors
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

# define our maps from 1 to 60 race & ethnicity and years
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
# check mixing and convergence
#-------------------------------------------------------------------------------

model_summary <- logpoi_country_model_fit$summary(
  variables = c("beta_0", "beta", "beta_sd"),
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
                           "bayesian-analysis-landfall-freq",
                           "no-year-effect-model", 
                           "no-year-model-diagnostics.html"))

# use webshot to capture the html table as a pdf
webshot(here("outputs",
             "bayesian-analysis-landfall-freq",
             "no-year-effect-model",  
             "no-year-model-diagnostics.html"), 
        here("outputs", 
             "bayesian-analysis-landfall-freq",
             "no-year-effect-model", 
             "no-year-model-diagnostics.pdf"))



# plot traces of parameter with smallest ess_bulk
# extract samples
model_draws <- logpoi_country_model_fit$draws(
  variables = c('lp__',
                "beta_0", 
                "beta", 
                "beta_sd"),
  inc_warmup = TRUE, format = "draws_array")

# make trace plots of worst performing parameter
model1_worst_var <- model_summary$variable[ which.min(model_summary$ess_bulk)]

bayesplot:::mcmc_trace(model_draws,  
                       pars = model1_worst_var, 
                       n_warmup = 500, 
                       facet_args = list(nrow = 1))+ 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 4))

#-------------------------------------------------------------------------------
# pairwise posterior geometry of 4 worst performing parameters
#-------------------------------------------------------------------------------

worst_param <- model_summary[1:4,]$variable
model_draws <- logpoi_country_model_fit$draws(
  variables = worst_param,inc_warmup = FALSE,
  format = "draws_array")

bayesplot::color_scheme_set('viridisC')
bayesplot::mcmc_pairs(model_draws, 
                      pars = worst_param, 
                      diag_fun = "dens", 
                      off_diag_fun = "hex") 


#-------------------------------------------------------------------------------
# posterior predictive checks
#-------------------------------------------------------------------------------

# extract Monte Carlo samples of log_lambda vector
model <- list()
model$log_lambda <- logpoi_country_model_fit$draws(variables = "log_lambda", 
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
logpoi_no_year_effect_model_lambda <- model$log_lambda

# clean up so we don't have millions of data points stored multiple times
model <- NULL
gc()

# make posterior predictions, conditional on the joint posterior for log_lambda
set.seed(42L)
logpoi_no_year_effect_model_lambda[ , post_pred := rpois( nrow(logpoi_no_year_effect_model_lambda), 
                                                          exp(log_lambda))]

# create median and 95\% credible intervals for variables 
logpoi_no_year_effect_model_summary <- logpoi_no_year_effect_model_lambda[,list(
  summary_value = quantile(post_pred, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
  summary_name = c('q_lower','iqr_lower', 'median','iqr_upper', 'q_upper')),
  by = 'ALL_ID']

logpoi_no_year_effect_model_summary<- data.table::dcast(logpoi_no_year_effect_model_summary,
                                                        ALL_ID ~ summary_name, 
                                                        value.var = 'summary_value')

month_effects <- subset(landfalls[1:300], 
                        select = c('MONTH', 'LANDFALL_COUNT', 'YEAR', 'ALL_ID'))

model_summary <- merge(month_effects,logpoi_no_year_effect_model_summary, by = 'ALL_ID')


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

p1 <- ggsave(here("outputs", 
                  "bayesian-analysis-landfall-freq",
                  "no-year-effect-model", 
                  "no-year-post-pred-checks.pdf"), 
             plot = p1, 
             height = 10, 
             width = 10, 
             dpi=300)

