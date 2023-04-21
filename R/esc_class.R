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
#' @examples
#'
new_esc <- function(intercept, slope, sigma, data = NULL) {
  res <- list(intercept = intercept,
              slope = slope,
              eff = exp(-log(10)/slope) - 1,
              sigma = sigma,
              data = data)
  structure(res, class = "esc")
}
