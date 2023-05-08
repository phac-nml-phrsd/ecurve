#' @title Print the parameters of an ESC model.
#'
#' @param model A model of class \code{esc}.
#'
#' @return Print the model parameters.
#'
#' @export
#'
#' @examples
#' esc_data = data.frame(
#'      concentrations = c(1,1,10, 10, 100, 500, 500),
#'      cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
#' mod = esc_mle(esc_data)
#' print_esc_model(mod)
#'
print_esc_model <- function(model) {

  msgdata = 'No data is attached to the model.\n\n'
  if(!is.null(model$data)){
    msgdata = paste0('Number of data points: ', nrow(model$data),'\n\n')
  }

  cat(paste0(
    '---- ESC model ----\n\n',
    'intercept = ', model$intercept, '\n',
    'slope     = ', model$slope, '\n',
    'efficacy  = ', model$eff, '\n',
    'sigma     = ', model$sigma, '\n',
    msgdata
  ))
}
