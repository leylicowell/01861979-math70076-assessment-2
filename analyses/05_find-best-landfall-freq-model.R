# in this script, we check the predictive significance of our yearly effect 
# by using approximate leave one out cross-validation to quantify the
# prediction accuracy of both our no-year and HSGP models
# we then use the best model to forecast landfall frequency over the next 5 years

#===============================================================================
# load packages and functions
#===============================================================================

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

source(here("src", "forecast-landfalls.R"))

#===============================================================================
# load model fits and data
#===============================================================================

logpoi_hsgp_model_fit <- readRDS(file = file.path(here("outputs", 
                                                       "stan-models", 
                                                       "hsgp-model-cmdstanr.rds")))

logpoi_no_year_effect_model_fit <- readRDS(file = file.path(here("outputs", 
                                                          "stan-models", 
                                                          "no-year-effect-model-cmdstanr.rds")))

landfalls <- as.data.table(read.csv(here("data", 
                                         "derived", 
                                         "stan-landfall-freq-data.csv")))

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
# and plot posterior predictive density for number of landfalls in 2025
#-------------------------------------------------------------------------------

forecast_data <- forecast_landfalls(logpoi_no_year_effect_model_fit, landfalls, "MONTH")
forecast_table <- forecast_data[[1]]
forecast_25_plot <- forecast_data[[2]]

forecast_table
forecast_25_plot

save_kable(forecast_table, file = here("outputs", 
                           "bayesian-analysis-landfall-freq", 
                           "landfalls-forecasts.html"))

# use webshot to capture the html table as a pdf
webshot(here("outputs", 
             "bayesian-analysis-landfall-freq", 
             "landfalls-forecasts.html"), 
        here("outputs", 
             "bayesian-analysis-landfall-freq", 
             "landfalls-forecasts.pdf"))


ggsave(here("outputs", 
            "bayesian-analysis-landfall-freq", 
            "landfall-density-plots.pdf"), 
       plot = forecast_25_plot,
       height = 10, 
       width = 10, 
       dpi=300)

