#' Create New esc Object
#'
#' Creates new object of class esc, representing ESC model as a list containing
#' intercept, slope and sigma parameters, as well as the calculated efficiency
#'
#' @param intercept intercept parameter of ESC model
#' @param slope slope parameter of ESC model
#' @param sigma sigma parameter of ESC model
#'
#' @return Object of class esc
#'
#' @examples
#'
new_esc <- function(intercept, slope, sigma) {
  res <- list(intercept = intercept, slope = slope, eff = exp(-1/slope) - 1,
              sigma = sigma)
  structure(res, class = "esc")
}
