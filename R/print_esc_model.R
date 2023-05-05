#' @title Print the parameters of an ESC model.
#'
#' @param model
#' @param show.data
#'
#' @return Print the model parameters.
#' @export
#'
#' @examples
#'
print_esc_model <- function(model) {
  cat(paste0(
    '---- ESC model ----\n\n',
    'intercept = ', model$intercept, '\n',
    'slope     = ', model$slope, '\n',
    'efficacy  = ', model$eff, '\n',
    'sigma     = ', model$sigma, '\n',
    'Number of data points: ', nrow(model$data),'\n\n'
  ))
}
