# our posterior predictive checks for our HSGP model displayed
# small yearly changes so in this script we design a simpler model, investigating
# monthly effects only, from count data stratified by month and year

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

landfalls <- as.data.table(read.csv(here("data", "derived", "stan-landfall-freq-data.csv")))

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
  real beta_0;
  real<lower=0> beta_sd;
  vector[M-1] beta_s2z_m1_rnde;
}
transformed parameters{
  vector[M] beta_nc_rnde = s2z_Q_b * beta_s2z_m1_rnde;
  vector[M] beta = beta_sd * beta_nc_rnde;
  vector[N] log_lambda;
  
  for (i in 1:N){
    log_lambda[i] = beta_0 + beta[month_id[i]];
  }
}
model{
  // likelihood
  y ~ poisson_log(log_lambda);
  
  // priors
  beta_0 ~ normal(0, 2);
  beta_s2z_m1_rnde ~ normal(0, s2z_sd_b);
  beta_sd ~ cauchy(0, 1);
}
generated quantities{
  vector[N_star] log_lambda_star;
  array[N_all] real y_all_pred;
  array [N] real log_lik;
  
  for(i in 1:N)
  {
    log_lik[i] = poisson_log_lpmf(y[i] | beta_0 + beta[month_id[i]]);
  }
  
  for (i in 1:N_star){
    log_lambda_star[i] = beta_0 + beta[month_id_star[i]];
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
  iter_sampling = 6e3,
  refresh = 500,
  save_warmup = TRUE)

# save output to RDS
logpoi_no_year_effect_model_fit$save_object(file = file.path(here("outputs", 
                                                                  "stan-models", 
                                                                  "no-year-effect-model-cmdstanr.rds")))


#-------------------------------------------------------------------------------
# check mixing and convergence
#-------------------------------------------------------------------------------

model_summary <- logpoi_no_year_effect_model_fit$summary(
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
model_draws <- logpoi_no_year_effect_model_fit$draws(
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
  coord_cartesian(ylim = c(-0.5, 10))

#-------------------------------------------------------------------------------
# pairwise posterior geometry of 4 worst performing parameters
#-------------------------------------------------------------------------------

worst_param <- model_summary[1:4,]$variable
model_draws <- logpoi_no_year_effect_model_fit$draws(
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
model$log_lambda <- logpoi_no_year_effect_model_fit$draws(variables = "log_lambda", 
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




