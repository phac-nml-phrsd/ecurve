#' Create New esc Object
#'
#' Creates new object of class esc, representing ESC model as a list containing
#' intercept, slope and sigma parameters, as well as the calculated efficiency
#'
#' @param intercept intercept parameter of ESC model
#' @param slope slope parameter of ESC model
#' @param sigma sigma parameter of ESC model
#' @param data data frame containing concentration and cq values used to fit the
#' model. Optional parameter, defaults to NULL if not provided.
#'
#' @return Object of class esc
#'
#' @examples
#'
new_esc <- function(intercept, slope, sigma, data = NULL) {
  res <- list(intercept = intercept, slope = slope, eff = exp(-1/slope) - 1,
              sigma = sigma, data = data)
  structure(res, class = "esc")
}
