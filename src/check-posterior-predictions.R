#' Plot posterior predictive checks
#'
#' @param model_fit CmdStanMCMC fitted Stan model object
#' @param parameter logarithm of poisson model parameter 
#' @param data Stan dataset used in fitted model
#' @param observation_col column name of data ID
#' @param count_col column name of count data
#' @param facet_col column name of categories we want to facet in our plot
#' @param x_col column name of x-axis data in our plot
#'
#' @returns plot of posterior predictive checks
#' @import cmdstanr
#' @import data.table
#' @import ggplot2
#' @export
#'
#' @examples
#' data = data.table(OBS_ID =1:5, 
#' YEAR = 2001:2005, 
#' COUNT = c(1,2,4,2,5), 
#' COUNTRY = c("UK", "FRANCE", "ALBANIA", "PORTUGAL", "RUSSIA"))
#' 
#' model_text <- '
#' data{
#' int<lower=1> N;
#' array[N] int<lower=0> y;
#' }
#' parameters {
#' vector[N] log_lambda;
#' }
#' model {
#' y ~ poisson_log(log_lambda);
#' }'
#' 
#' model <- write_stan_file(
#' model_text,
#' dir = getOption("cmdstanr_write_stan_file_dir", tempdir()),
#' basename = NULL,
#' force_overwrite = FALSE,
#' hash_salt = ""
#' )
#' 
#' compiled <- cmdstan_model(model)
#' fit <- compiled$sample(
#' data = list(N = 5, y = data$COUNT),
#' seed = 123, 
#' chains = 2,
#' iter_warmup = 500, 
#' iter_sampling = 1e3, 
#' refresh = 0
#' )
#' 
#' check_posterior_predictions(
#' model_fit = fit,
#' data = data,
#' observation_col = "OBS_ID",
#' count_col = "COUNT",
#' facet_col = "COUNTRY",
#' x_col = "YEAR"
#' )


check_posterior_predictions <- function(model_fit, 
                                        data,
                                        observation_col,
                                        count_col,
                                        facet_col,
                                        x_col,
                                        parameter = "log_lambda"){
  
  # extract Monte Carlo samples of vector log_lambda
  # each dimension in one column, each draw in one row
  sample <- list()
  sample[[parameter]] <- model_fit$draws(variables = parameter, 
                                                      inc_warmup = FALSE,format = "draws_df")
  
  # MELT parameter into long format, and link to observation ID column in the data
  sample[[parameter]] <- data.table::melt(as.data.table(sample[[parameter]]), 
                                        id.vars = c('.chain','.iteration','.draw')
  )
  
  set(sample$log_lambda, 
      NULL, 
      observation_col, 
      as.integer(gsub(paste0(parameter,'\\[([0-9]+)\\]'),
                      '\\1',
                      as.character(sample[[parameter]]$variable))))
  
  setnames(sample[[parameter]], c('value'), c(parameter))
  
  set(sample[[parameter]], NULL, 'variable', NULL)
  
  model_sample <- sample[[parameter]]
  
  # clean up so we don't have millions of data points stored multiple times
  sample <- NULL
  gc()
  
  # make posterior predictions for each data point, conditional on the 
  # joint posterior for log_lambda
  set.seed(42L)
  model_sample[ , post_pred := rpois( nrow(model_sample), exp(get(parameter)))]
  
  # create median and 95% credible intervals for variables 
  # using AGGREGATION over .draws
  model_sample_summary <- model_sample[,list(
    summary_value = quantile(post_pred, prob = c(0.025, 0.25, 0.5, 0.75, 0.975)),
    summary_name = c('q_lower','iqr_lower','median','iqr_upper', 'q_upper')),
    by = observation_col]
  
  # reshape using CASTING
  model_sample_summary <- data.table::dcast(model_sample_summary,
                                            get(observation_col) ~ summary_name, 
                                            value.var = 'summary_value')
  
  setnames(model_sample_summary, old = "observation_col", new = observation_col)
  
  # merge data using INNER JOIN
  model_sample_summary <- merge(model_sample_summary, 
                                subset(data, select = c(observation_col,
                                                        count_col,
                                                        facet_col,
                                                        x_col)), 
                                by = observation_col)
  
  # if we have different plots for each month, ensure these are chronological
  if(facet_col == 'MONTH'){
    # reorder months in chronological order
    model_sample_summary[[facet_col]] <- factor(model_sample_summary[[facet_col]], 
                                                levels = month.name)  
  }
  
  # make posterior predictive check
  model_sample_summary <- as.data.table(model_sample_summary)
  model_sample_summary[,IN_CI := get(count_col) >= q_lower & get(count_col) <= q_upper]
  cat(round(model_sample_summary[, mean( as.numeric( IN_CI ) )]*100, 1),
      "% of all observed values are within our model's posterior predictive intervals.")
  

  # plot posterior predictive check
  post_plot <- ggplot(model_sample_summary, aes(x = as.factor(get(x_col)))) + 
    geom_boxplot( aes( group = get(x_col), 
                       ymin = q_lower, 
                       lower = iqr_lower,
                       middle = median, 
                       upper = iqr_upper, 
                       ymax = q_upper), 
                  stat = 'identity') +
    geom_point( aes(y = get(count_col), colour = IN_CI ) ) +
    scale_y_continuous() + 
    ggsci::scale_color_npg() +
    labs(x = 'Year', 
         y = 'Number of landfalls', 
         colour = 'within\n95% posterior\nprediction\ninterval') +
    facet_grid(get(facet_col) ~ ., labeller = label_wrap_gen(10)) + 
    theme_bw()+
    theme(axis.title = element_text(size = 18), 
          axis.text.x = element_text(size = 16, angle = 45,vjust = 1,hjust = 1),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12.5))
  
  return(post_plot)
}





