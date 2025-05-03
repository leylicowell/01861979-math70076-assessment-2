#' Check model convergence diagnostics for Stan models
#'
#' @param model_fit CmdStanMCMC fitted Stan model object
#' @param var variables/parameters from our Stan model
#' @param warmup TRUE if warmup was included in stan model or FALSE otherwise
#' @param warmup_nbr number of iterations in model which were part of the warmup
#'
#' @returns table of model convergence diagnostics for 5 worst chosen parameters
#' 
#' @import cmdstanr
#' @importFrom data.table as.data.table
#' @importFrom kableExtra kable_styling, kbl
#' @import ggplot2
#' @export 
#'
#' @examples 
#' fit <- cmdstanr_example("logistic")
#' check_model_diagnostics(fit, c("alpha", "beta"), warmup = FALSE, warmup_nbr = 0)
#' 
check_model_diagnostics <- function(model_fit, var, warmup = TRUE, warmup_nbr = 500){
  # summarise model parameter measures for diagnostics table
  model_summary <- model_fit$summary(
    variables = var,
    posterior::default_summary_measures(),
    posterior::default_convergence_measures())
  
  # extract samples for trace plot and pairwise plot
  model_draws <- model_fit$draws(
    variables = var,
    inc_warmup = warmup, format = "draws_array")
  
  # sort by smallest ess_bulk
  model_summary  <- as.data.table(model_summary )
  model_summary  <- model_summary[order(ess_bulk),]
  
  # print table of parameter convergence diagnostics
  print(kbl(model_summary, caption = 'Model diagnostics', longtable = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                  font_size = 12))
  
  # save table of 5 (or all parameters if there are less than 5) worst performing parameters
  t1 <- kbl(subset(model_summary[1:min(5, nrow(model_summary)),], 
                   select = c(variable, rhat, ess_bulk, ess_tail)), 
            longtable = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), 
                  font_size = 16) 
  
  
  # make trace plot of worst performing parameter
  model1_worst_var <- model_summary$variable[ which.min(model_summary$ess_bulk)]
  
  print(bayesplot:::mcmc_trace(model_draws,  
                         pars = model1_worst_var, 
                         n_warmup = warmup_nbr, 
                         facet_args = list(nrow = 1))+ 
    theme_bw())
  
  # only print pairwise plot if we have more than one parameter
  if(nrow(model_summary) > 1){
    # make pairwise plot of 4 (or all, if less than 4) model parameters
    worst_param <- model_summary[1:min(4, nrow(model_summary)),]$variable
    model_draws <- model_fit$draws(
      variables = worst_param,inc_warmup = FALSE,
      format = "draws_array")
    
    bayesplot::color_scheme_set('viridisC')
    print(bayesplot::mcmc_pairs(model_draws, 
                          pars = worst_param, 
                          diag_fun = "dens", 
                          off_diag_fun = "hex"))
  }
    
  # return diagnostics table for report
  return(t1)
}





