# In this script, we analyse the number of landfalls per country and year in our hurricane data
# we use a Bayesian Poisson regression model for the top 5 countries with the most
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

hurricane_data <- read.csv(here("data", "derived", "hurricane-data.csv"))

landfalls <- hurricane_data %>%
  filter(RECORD_ID == "L")

head(hurricane_data)
tail(hurricane_data)
