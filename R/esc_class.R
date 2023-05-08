#' @title Create a New ESC Object
#'
#' @description Creates new object of class \code{esc},
#' representing a ESC model as a list containing
#' intercept, slope and sigma parameters, as well as the calculated efficiency.
#'
#' @param intercept intercept parameter of ESC model.
#' @param slope slope parameter of ESC model.
#' @param sigma sigma parameter of ESC model.
#' @param data data frame containing concentration and Cq values used to fit the
#' model. Optional parameter, defaults to \code{NULL} if not provided.
#'
#' @return Object of class \code{esc}.
#'
new_esc <- function(intercept, slope, sigma, data = NULL) {
  res <- list(intercept = intercept,
              slope     = slope,
              eff       = exp(-log(10)/slope) - 1,
              sigma     = sigma,
              data      = data)

  # Check inputs
  if(!is.numeric(intercept)) stop('`intercept` must be a number.')
  if(!is.numeric(slope)) stop('`slope` must be a number.')
  if(!is.numeric(sigma)) stop('`sigma` must be a number.')
  if(intercept < 0) stop('`intercept` must be a positive number.')
  if(sigma < 0) stop('`sigma` must be a positive number.')

  # Return class object
  structure(res, class = "esc")
}

#' Create ESC Objects With Specified Parameters
#'
#' Creates \code{esc} object, given parameters intercept, eff, and sigma
#'
#' @param intercept Numeric. Intercept parameter of ESC model.
#' @param eff Numeric. qPCR efficiency.
#' @param sigma Numeric. Sigma parameter of ESC model.
#'
#' @return \code{esc} object
#' @export
#'
#' @examples
#'
#' m = gen_esc(intercept = 39, eff = 0.98, sigma = 0.2)
#' print_esc_model(m)
#'
gen_esc <- function(intercept, eff, sigma) {
  # Check input
  # (leave `intercept` and `sigma` to the nested function `new_esc()`)
  if(eff <= 0) stop('Efficiency `eff` must be positive.')

  return(new_esc(intercept, slope = -1 / log10(1 + eff), sigma))
}
