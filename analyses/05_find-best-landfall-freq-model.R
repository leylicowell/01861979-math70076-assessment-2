# in this script, we check the predictive significance of our yearly effect 
# by using approximate leave one out cross-validation to quantify the
# prediction accuracy of both our no-year and HSGP models

#===============================================================================
# load packages
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

#===============================================================================
# load model fits and data
#===============================================================================

logpoi_hsgp_model_fit <- readRDS(file = file.path(here("outputs", 
                                                       "stan-models", 
                                                       "hsgp-model-cmdstanr.rds")))

logpoi_no_year_effect_model_fit <- readRDS(file = file.path(here("outputs", 
                                                          "stan-models", 
                                                          "no-year-effect-model-cmdstanr.rds")))

landfalls_monthly <- as.data.table(read.csv(here("data", 
                                         "derived", 
                                         "stan-landfall-monthly-freq.csv")))

#-------------------------------------------------------------------------------
# compare prediction accuracy of both models
#-------------------------------------------------------------------------------

# calculate approximate LOO cross validation for model without year effect
logpoi_no_year_effect_model_fit_loo <- logpoi_no_year_effect_model_fit$loo()
logpoi_hsgp_model_fit_loo <- logpoi_hsgp_model_fit$loo()

# compare ELPD across models
comp <- loo::loo_compare(list(no_year_effect = logpoi_no_year_effect_model_fit_loo, 
                              hsgp = logpoi_hsgp_model_fit_loo))

print(comp, simplify = FALSE)

# we see that the model without yearly effects results in better expected log 
# posterior predictive accuracy than the HSGP model

# we now want to save this as a table for our report
comp_data <- as.data.frame(comp)

# create table
comp_table <- kbl(comp_data, caption = "LOO Model Comparison") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

print(comp_table)

# save table
save_kable(comp_table, file = here("outputs", 
                                       "bayesian-analysis-landfall-freq", 
                                       "best-model-comp.html"))

# use webshot to capture the html table as a pdf
webshot(here("outputs", 
             "bayesian-analysis-landfall-freq", 
             "best-model-comp.html"), 
        here("outputs", 
             "bayesian-analysis-landfall-freq", 
             "best-model-comp.pdf"))



