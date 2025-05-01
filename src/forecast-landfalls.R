#' Forecast landfall frequency from fitted poisson regression Stan model
#'
#' @param model_fit CmdStanMCMC fitted Stan model object
#' @param data Stan data used to fit model
#' @param category_name column name for group categories
#'
#' @returns table of forecasted landfall numbers per category over the next 5 years
#' @import KableExtra
#' @import data.table
#' @import dplyr
#' @export
#'

forecast_landfalls <- function(model_fit, data, category_name){
  model <- list()
  model$log_lambda_star <- model_fit$draws(
    variables = "log_lambda_star", 
    inc_warmup = FALSE, 
    format = "draws_df")
  
  # MELT log_lambda_star into long format, and link to pred_idx in the data
  model$log_lambda_star <- data.table::melt(as.data.table(model$log_lambda_star), 
                                            id.vars = c('.chain','.iteration','.draw'))
  set(model$log_lambda_star, 
      NULL, 
      'PRED_ID',
      gsub('log_lambda_star\\[([0-9]+)\\]','\\1',
           as.character(model$log_lambda_star$variable)))
  
  set(model$log_lambda_star, NULL, 'PRED_ID', as.integer(model$log_lambda_star$PRED_ID))
  setnames(model$log_lambda_star, c('value'), c('log_lambda_star'))
  set(model$log_lambda_star, NULL, 'variable', NULL)
  model_lambda_star <- model$log_lambda_star
  
  # clean up so we don't have millions of data points stored multiple times
  model <- NULL
  gc()
  
  # make posterior predictions for each data point, conditional on the 
  # joint posterior for log_lambda
  set.seed(42L)
  model_lambda_star[,post_pred := rpois(
    nrow(model_lambda_star), 
    exp(log_lambda_star))]
  
  # create median and 95% credible intervals for variables 
  # using AGGREGATION over .draws
  model_star_summary <- model_lambda_star[,list( 
    summary_value = quantile(post_pred, prob =c(0.025,  0.5, 0.975)),
    summary_name = c('2.5% quantile','median','97.5% quantile')),
    by = 'PRED_ID']
  
  # reshape using CASTING
  model_star_summary<- data.table::dcast(
    model_star_summary, 
    PRED_ID ~ summary_name, 
    value.var = 'summary_value')
  
  # merge data using INNER JOIN
  data
  forecast_summary <- subset(data[!is.na(PRED_ID)], select = c(category_name, 'YEAR', 'PRED_ID'))
  forecast_summary <- merge(forecast_summary, model_star_summary, by = 'PRED_ID')
  set(forecast_summary, NULL, 'PRED_ID', NULL)
  
  # plot table
  print(kbl(forecast_summary, 
      caption = paste0('Forecasted number of landfalls per ', 
                       tolower(category_name), 
                       ' in the next 5 years'), 
      longtable = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12))
  
  
  # create median and 95% credible intervals for variables 
  # using AGGREGATION over .draws for total monthly forecasts over 2025-2029
  total_group_draws <- left_join(model_lambda_star,
                                    subset(data[!is.na(PRED_ID)], 
                                           select = c("PRED_ID", category_name)), 
                                    by = "PRED_ID" )
  
  group_sums <- total_group_draws[, .(
    group_total = sum(post_pred)
  ), by = c(".draw", category_name)]
  
  group_sums <- as.data.table(group_sums)
  
  total_group_summary <- group_sums[,list( 
    summary_value = quantile(group_total, prob =c(0.025, 0.25,  0.5, 0.75, 0.975)),
    summary_name = c('2.5% quantile', 'IQR_lower','median', 'IQR_upper', '97.5% quantile')),
    by = category_name]
  
  total_group_summary <- data.table::dcast(
    total_group_summary, 
    as.formula(paste0(category_name, " ~ summary_name")), 
    value.var = 'summary_value')
  
  if(category_name == 'MONTH'){
    # reorder months in chronological order
    total_group_summary[[category_name]]<- factor(total_group_summary[[category_name]], 
                                                  levels = month.name)  
    total_group_summary <- total_group_summary[order(total_group_summary[[category_name]]), ]
  }
  
  
  forecast_table <- kbl(total_group_summary, 
                caption = paste0('Total forecasted number of landfalls per ', 
                tolower(category_name), ' in the next 5 years'), 
                longtable = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 12)
  
  #  plot posterior predictive density for number of landfalls in 2025
  forecasts_2025 <- subset(data[!is.na(PRED_ID)], 
                          select = c(category_name, 'YEAR', 'PRED_ID')) %>%
    filter(YEAR == 2025)
  
  
  forecasts_2025 <- merge(model_lambda_star,
                         forecasts_2025,
                         by = "PRED_ID")
  
  if(category_name == 'MONTH'){
    # reorder months in chronological order
    forecasts_2025[[category_name]]<- factor(forecasts_2025[[category_name]], 
                                                  levels = month.name)  
    forecasts_2025 <- forecasts_2025[order(forecasts_2025[[category_name]]), ]
  }
  
  forecasts_25_plot <- ggplot(forecasts_2025, aes(x = post_pred)) +
    geom_histogram(aes(y = after_stat(density)), 
                   binwidth = 1, 
                   fill = "steelblue",
                   alpha = 0.7, 
                   color = "black") +
    facet_wrap(~ get(category_name), ncol = 3) +  # Create a separate plot for each month
    theme_bw() +
    labs(title = "Posterior Density Plots for 2025 Landfall Predictions", 
         x = "Number of Landfalls, L", 
         y = "P(L)") +
    theme(legend.position = "none")
  
  
  return(list(forecast_table, forecasts_25_plot))
  
}



