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
  if(!all(is.numeric(concentrations))) {stop("concentrations must be numeric")}
  if(!all(concentrations >= 0 & is.finite(concentrations))) {
    stop("concentrations must be non-negative real numbers")
  }
  if(!is.numeric(eff)) {stop("eff must be numeric")}
  if(!is.numeric(cq1)) {stop("cq1 must be numeric")}
  if(!is.numeric(sigma)) {stop("sigma must be numeric")}
  if(!is.finite(eff) | eff < 0) {stop("invalid efficiency ", eff, " provided")}
  if(!is.finite(cq1)) {stop("invalid intercept ", cq1, " provided")}
  if(!is.finite(sigma) | sigma < 0) {stop("invalid sigma ", sigma, " provided")}
  if(eff > 1) {warning("efficiency greater than 1 provided")}
  if(cq1 < 0) {warning("negative intercept provided")}
  N0s <- rpois(length(concentrations), concentrations)
  detects <- which(N0s > 0)
  cqs <- rep(NaN, length(N0s))
  cqs[detects] <- rnorm(length(detects), mean = cq1 - log(N0s[detects])/log(1 + eff),
                        sd = sigma)
  cqs
}
