# in this script, we implement a bayesian regression model for economical damages
# using our merged data set

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
library(readr)
library(stringr)
library(geosphere)
library(corrr)
library(forcats)


# set colour scheme for Stan
bayesplot::color_scheme_set("brewer-RdYlBu")

#==============================================================================
# we load our data set
#==============================================================================

merged_data <- read.csv(here("data", "derived", "merged-data.csv"))

#==============================================================================
# we make a pairwise correlation plot of wind at time of landfall, 
# overall curvature of hurricane path, duration of hurricane and adjusted total 
# damages
#==============================================================================

# reformat date and time to calculate duration
merged_data <- merged_data %>%
  mutate(TIME = sprintf("%04s", as.integer(TIME)),
         DATE_TIME = as.POSIXct(paste(YEAR, 
                                      MONTH, 
                                      DAY, 
                                      str_sub(TIME, 1, 2), 
                                      str_sub(TIME, 3, 4)),
                                format = "%Y %m %d %H %M"))


storm_data <- merged_data %>%
  arrange(STORM_ID, DATE_TIME) %>%
  group_by(STORM_ID, NAME, YEAR)%>%
  mutate(
    NEXT_LAT = lead(LAT),
    NEXT_LON = lead(LON),
    BEARING = bearing(cbind(LON, LAT), cbind(NEXT_LON, NEXT_LAT)),
    # change bearings from -180 to 180 to 0 to 360 degrees
    BEARING = BEARING + 180,
    DIFF_BEARING = abs(lead(BEARING)-BEARING)) %>%
  summarise(
    total_bearing_change = sum(DIFF_BEARING, na.rm = TRUE),
    avg_landfall_wind = mean(WIND[RECORD_ID =="L"]),
    total_damage = mean(ADJ_TOTAL_DAMAGE), # damage is the same for all rows for each storm
    duration_h = as.numeric(difftime(max(DATE_TIME), min(DATE_TIME), units = "hours")),
    .groups = "drop")

# extract relevant columns for pairwise correlation plot
corr_matrix <- storm_data %>%
  select("avg_landfall_wind", 
         "total_bearing_change", 
         "duration_h", 
         "total_damage")%>%
  correlate(diagonal = 1) %>%
  shave(upper = FALSE)

corr_matrix<- corr_matrix %>%
  pivot_longer(cols = -term,
               names_to = "colname",
               values_to = "corr") %>%
  mutate(rowname = fct_inorder(term),
         colname = fct_inorder(colname),
         label = if_else(is.na(corr), "", sprintf("%1.2f", corr)))

p1 <- ggplot(corr_matrix, aes(rowname, fct_rev(colname),
                 fill = corr)) +
  geom_tile() +
  geom_text(aes(
    label = label,
    color = abs(corr) < .75
  )) +
  coord_fixed(expand = FALSE) +
  scale_color_manual(
    values = c("white", "black"),
    guide = "none"
  ) +
  scale_fill_distiller(
    palette = "PuOr", na.value = "white",
    direction = 1, limits = c(-1, 1),
    name = "Pearson\nCorrelation:"
  ) +
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(color = NA, fill = NA),
        legend.position = c(.85, .8))

ggsave(here("outputs", "bayesian-analysis-damages", "correlation-plot.pdf"), 
       plot =p1,
       height = 7, 
       width = 7, 
       dpi=600)

#===============================================================================
# Bayesian regression analysis of cyclone economical damage
#===============================================================================

log_norm_horseshoe_txt <-"
// for horseshoe prior
functions {
    // square root of a vector (elementwise)
    vector sqrt_vec(vector x) {
        vector[dims(x)[1]] res;
        for (m in 1:dims(x)[1])
        {
            res[m] = sqrt(x[m]);}
        return res;
    }
}
data {
  int<lower=0> N; // number of observations
  int<lower=0> P; // number of predictors
  vector<lower=0>[N] y; // total damage
  matrix[N, P] X; //predictors (wind, duration, total bearing change and interaction effects)
  real<lower=1> nu; // degrees of freedom for the half t-priors
}
parameters {
  real beta_0; // baseline parameter
  real<lower=0> sigma;  
  // auxiliary variables for the variance parameters
  vector[P] z;
  // now the global scale parameter
  // it is computionally much better to write the Cauchy/Student
  // as infinite normal-inverse-gamma mixtures
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  vector<lower=0>[P] r1_local;
  vector<lower=0>[P] r2_local;
}

transformed parameters {
    // global and local variance parameters, and the input weights
    real<lower=0> tau;
    vector<lower=0>[P] lambda;
    vector[P] beta;

    // the normal-inverse-gamma mixture implementation
    tau = r1_global * sqrt(r2_global); 
    // the normal-inverse-gamma mixture implementation
    lambda = r1_local .* sqrt_vec(r2_local); 
    beta = z .* lambda * tau; // non-centered parameterisation
}
model {
    // likelihood
    y ~ lognormal(beta_0 + X*beta, sigma);
    // weakly informative prior for the intercept
    beta_0 ~ normal(0,5);   
    // local scale parameters, nu = 1 corresponds to horseshoe
    z ~ normal(0, 1);
    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5*nu, 0.5*nu);
    // global scale parameter
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);
    sigma ~ cauchy(0,1);
}"


#-------------------------------------------------------------------------------
# define data in format needed for Stan model
#-------------------------------------------------------------------------------

#storm_data$wind_duration_effect = storm_data$avg_landfall_wind * storm_data$duration_h
#storm_data$wind_bearings_effect = storm_data$avg_landfall_wind * storm_data$total_bearing_change

#-------------------------------------------------------------------------------
# compile model
#-------------------------------------------------------------------------------

log_norm_horseshoe_filename <- cmdstanr::write_stan_file(
  gsub('\t',' ', log_norm_horseshoe_txt),
  dir = here("analyses", "stan-models"),
  basename = NULL,
  force_overwrite = FALSE,
  hash_salt = "")

# compile Stan model
log_norm_horseshoe_compiled <- cmdstanr::cmdstan_model(log_norm_horseshoe_filename)

# define our data for stan
stan_data <- list()
stan_data$N <- nrow(storm_data)
stan_data$X <- subset(storm_data, select = c(-STORM_ID, -NAME, -YEAR, -total_damage))
stan_data$X <- unname(as.matrix(stan_data$X))
stan_data$P <- 3
stan_data$y <- storm_data$total_damage
stan_data$nu <- 2.0

# sample
log_norm_horseshoe_fit <- log_norm_horseshoe_compiled$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1e3,
  iter_sampling = 10e3,
  refresh = 500, 
  save_warmup = TRUE,
  adapt_delta = 0.8)

# save output to RDS
log_norm_horseshoe_fit$save_object(file = file.path(here("analyses", 
                                                         "stan-models", 
                                                         "lognorm_model_cmdstanr.rds")))
model_summary <- log_norm_horseshoe_fit$summary(
  variables = c("beta_0", "beta"),
  posterior::default_summary_measures(),
  posterior::default_convergence_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975)))

# sort by smallest ess_bulk
model_summary  <- as.data.table(model_summary )
model_summary  <- model_summary[order(ess_bulk),]





    