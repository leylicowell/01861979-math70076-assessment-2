library(cmdstanr)
library(data.table)
library(kableExtra)

# check_model_diagnostics example
# should print trace plot without warmup and no pairs plot
cmdstanr::cmdstanr_example("logistic")
testthat::expect_s3_class(check_model_diagnostics(fit, 
                                                  "alpha", 
                                                  warmup = FALSE, 
                                                  warmup_nbr = 0), 
                          "knitr_kable")

# check_model_diagnostics with invalid parameters
fit <- cmdstanr::cmdstanr_example("logistic")
testthat::expect_error(check_model_diagnostics(fit, "not_a_real_param"))


# check_model_diagnostics with warmup and more than one parameter
fit <- cmdstanr::cmdstanr_example("logistic", save_warmup = TRUE, iter_warmup = 500)
testthat::expect_s3_class(check_model_diagnostics(fit, c("alpha", "beta")),
                                 "knitr_kable")


