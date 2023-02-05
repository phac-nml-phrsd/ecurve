#' Simulate Cq values based on ESC Model
#'
#' Generate simulated Cq values given vector of concentrations and parameters
#' for the ESC Model. Non-detects are coded as NaN
#'
#' @param concentrations Vector of input concentrations
#' @param eff PCR Efficiency
#' @param cq1 Cq value for initial copy number of 1
#' @param sigma Standard Deviation of Gaussian error in Cq values
#'
#' @return Numeric vector of Cq Values
#' @export
#'
#' @examples
sim_cqs <- function(concentrations, eff, cq1, sigma){
  N0s <- rpois(length(concentrations), concentrations)
  detects <- which(N0s > 0)
  cqs <- rep(NaN, length(N0s))
  cqs[detects] <- rnorm(length(detects), mean = cq1 - log(N0s[detects])/log(1 + eff),
                        sd = sigma)
  cqs
}
