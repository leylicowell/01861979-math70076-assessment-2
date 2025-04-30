library(cmdstanr)
library(data.table)
library(ggplot2)

#check_posterior_predictions runs and returns a ggplot
data = data.table(OBS_ID =1:5, 
                  YEAR = 2001:2005, 
                  COUNT = c(1,2,4,2,5), 
                  COUNTRY = c("UK", "FRANCE", "ALBANIA", "PORTUGAL", "RUSSIA"))

model_text <- '
  data{
  int<lower=1> N;
  array[N] int<lower=0> y;
  }
  parameters{
  vector[N] log_lambda;
  }
  model {
  y ~ poisson_log(log_lambda);
  }'

model <- write_stan_file( model_text,
                          dir = getOption("cmdstanr_write_stan_file_dir", tempdir()),
                          basename = NULL,
                          force_overwrite = FALSE,
                          hash_salt = "")

compiled <- cmdstan_model(model)
fit <- compiled$sample(data = list(N = 5, y = data$COUNT),
                       seed = 123,
                       chains = 2,
                       iter_warmup = 500, 
                       iter_sampling = 1e3, 
                       refresh = 0)

plot <- check_posterior_predictions(model_fit = fit,
                                    parameter = "log_lambda",
                                    data = data,
                                    observation_col = "OBS_ID",
                                    count_col = "COUNT",
                                    facet_col = "COUNTRY",
                                    x_col = "YEAR")

testthat::expect_s3_class(plot, "ggplot")

# visual check that our plot looks correct
plot 
