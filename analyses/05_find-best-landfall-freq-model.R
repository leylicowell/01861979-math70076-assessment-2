# in this script, we check the predictive significance of our yearly effect 
# by using approximate leave one out cross-validation to quantify the
# prediction accuracy of both our no-year and HSGP models
# we then use the best model to forecast landfall frequency over the next 5 years

library(data.table)  
library(ggplot2) 
library(ggsci) 
library(kableExtra)  
library(cmdstanr)  
library(webshot2)
library(here)
library(dplyr)
library(magrittr)
library(tidyr)

#-------------------------------------------------------------------------------
# load model fits and data
#-------------------------------------------------------------------------------

logpoi_hsgp_model_fit <- readRDS(file = file.path(here("outputs", 
                                                       "stan-models", 
                                                       "hsgp-model-cmdstanr.rds")))
logpoi_no_year_effect_model_fit <- readRDS(file = file.path(here("outputs", 
                                                          "stan-models", 
                                                          "no-year-effect-model-cmdstanr.rds")))

landfalls <- as.data.table(read.csv(here("data", "derived", "stan-landfall-freq-data.csv")))

#-------------------------------------------------------------------------------
# compare prediction accuracy of both models
#-------------------------------------------------------------------------------

# calculate approximate LOO cross validation for model without year effect
logpoi_no_year_effect_model_fit_loo <- logpoi_no_year_effect_model_fit$loo()
logpoi_hsgp_model_fit_loo <- logpoi_hsgp_model_fit$loo()

# compare ELPD across models
comp <- loo::loo_compare(logpoi_no_year_effect_model_fit_loo, 
                         logpoi_hsgp_model_fit_loo
)
print(comp, simplify = FALSE)

# we see that the model without yearly effects results in better expected log 
# posterior predictive accuracy than the HSGP model


#-------------------------------------------------------------------------------
# we use the model without yearly effects to forecast 
# number of landfalls over the next 5 years
#-------------------------------------------------------------------------------

model <- list()
model$log_lambda_star <- logpoi_no_year_effect_model_fit$draws(variables = "log_lambda_star", 
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
logpoi_no_year_effect_model_lambda_star <- model$log_lambda_star

# clean up so we don't have millions of data points stored multiple times
model <- NULL
gc()

# make posterior predictions for each data point, conditional on the 
# joint posterior for log_lambda
set.seed(42L)
logpoi_no_year_effect_model_lambda_star[,post_pred := rpois( nrow(logpoi_no_year_effect_model_lambda_star), 
                                                   exp(log_lambda_star))]

# create median and 95\% credible intervals for variables 
# using AGGREGATION over .draws
logpoi_no_year_effect_model_star_summary <- logpoi_no_year_effect_model_lambda_star[,list( 
  summary_value = quantile(post_pred, prob =c(0.025,  0.5, 0.975)),
  summary_name = c('2.5% quantile','median','97.5% quantile')),
  by = 'PRED_ID']

# reshape using CASTING
logpoi_no_year_effect_model_star_summary<- data.table::dcast(
  logpoi_no_year_effect_model_star_summary, 
  PRED_ID ~ summary_name, 
  value.var = 'summary_value')

# merge data using INNER JOIN
forecast_summary <- subset(landfalls[301:360], select = c('MONTH', 'YEAR', 'PRED_ID'))
forecast_summary <- merge(forecast_summary,logpoi_no_year_effect_model_star_summary, by = 'PRED_ID')
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
total_monthly_draws <- left_join(logpoi_no_year_effect_model_lambda_star,
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

t1 <- kbl(total_monthly_summary, 
          caption = 'Total forecasted number of landfalls per month in the next 5 years', 
          longtable = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12)

t1

save_kable(t2, file = here("outputs", 
                           "bayesian-analysis-landfall-freq", 
                           "landfalls-forecasts.html"))

# use webshot to capture the html table as a pdf
webshot(here("outputs", 
             "bayesian-analysis-landfall-freq", 
             "landfalls-forecasts.html"), 
        here("outputs", 
             "bayesian-analysis-landfall-freq", 
             "landfalls-forecasts.pdf"))


#-------------------------------------------------------------------------------
# posterior predictive density for number of landfalls in 2025
#-------------------------------------------------------------------------------

forecast_2025 <- subset(landfalls[301:360], select = c('MONTH', 'YEAR', 'PRED_ID')) %>%
  filter(YEAR == 2025)


forecast_2025 <- merge(logpoi_no_year_effect_model_lambda_star,
                       forecast_2025,
                       by = "PRED_ID")

# reorder months in chronological order
forecast_2025$MONTH <- factor(forecast_2025$MONTH, levels = month.name)  
forecast_2025 <- forecast_2025[order(forecast_2025$MONTH), ]

p1 <- ggplot(forecast_2025, aes(x = post_pred)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 1, 
                 fill = "steelblue",
                 alpha = 0.7, 
                 color = "black") +
  facet_wrap(~ MONTH, ncol = 3) +  # Create a separate plot for each month
  theme_bw() +
  labs(title = "Posterior Density Plots for 2025 Landfall Predictions", 
       x = "Number of Landfalls, L", 
       y = "P(L)") +
  theme(legend.position = "none")

p1

ggsave(here("outputs", 
            "bayesian-analysis-landfall-freq", 
            "landfall-density-plots.pdf"), 
       plot =p1,
       height = 10, 
       width = 10, 
       dpi=300)

