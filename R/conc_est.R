conc_likelihood_factory <- function(cqs, intercept, slope, sigma) {
  detects <- which(!is.nan(cqs))
  num_non_detects <- length(cqs) - length(detects)
  cqs <- cqs[detects]
  N0max <- N0min <- round(exp(mean(cqs) - intercept)/slope)
  gen_norm_densities <- function(N0) {
    dnorm(cqs, mean = intercept + slope * log(N0), sd = sigma)
  }
  norm_densities <- gen_norm_densities(N0max)
  function(concentration) {
    if(concentration < 0) {return(0)}
    N0start <- max(qpois(1E-15, concentration), 1)
    N0end <- qpois(1E-15, concentration, lower.tail = FALSE)
    if (N0start < N0min) {
      norm_densities <<- cbind(sapply(N0start:(N0min - 1), gen_norm_densities),
                               norm_densities)
      N0min <<- N0start
    }
    if(N0end > N0max) {
      norm_densities <<- cbind(norm_densities,
                               sapply((N0max + 1):N0end, gen_norm_densities))
      N0max <<- N0end
    }
    N0s <- N0start:N0end
    prod(norm_densities[,N0s - N0min + 1] %*% dpois(N0s, concentration),
         exp(-num_non_detects * concentration))
  }
}

#' Estimate concentration from Cq values
#'
#' Given list of Cq values from a set of technical replicates and fitted ESC
#' model, generates maximum likelihood estimate of concentration
#'
#' @param cqs Numeric vector of Cq values from sample replicates, non-detects
#' coded as NaN
#' @param model esc object representing fitted model to use for estimation
#'
#' @return MLE of concentration
#' @export
#'
#' @examples
conc_mle <- function(cqs, model) {
  if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(class(model) != "esc") {stop("model is not an esc object")}
  if(all(is.nan(cqs))) {return(0)}
  conc_likelihood <- conc_likelihood_factory(cqs, model$intercept, model$slope,
                                             model$sigma)
  nlm(function(conc) {-log(conc_likelihood(conc))},
      exp((mean(cqs, na.rm = TRUE) - model$intercept)/model$slope))$estimate
}

#' Generate credible interval for concentration analytically
#'
#' Given list of Cq values from a set of technical replicates and fitted ESC
#' model, generates credible interval of desired level through numerical
#' integration of the posterior distribution
#'
#' @param cqs Numeric vector of Cq values from sample replicates, non-detects
#' coded as NaN
#' @param model esc object representing fitted model to use for estimation
#' @param level Desired credibily level, defaults to 0.95
#'
#' @return Numeric vector specifying the lower and upper bounds of desired
#' credible interval
#' @export
#'
#' @examples
conc_interval <- function(cqs, model, level = 0.95) {
  if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(class(model) != "esc") {stop("model is not an esc object")}
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}
  if(all(is.nan(cqs))) {
    return(c(0, -log(1 - level)/length(cqs)))
  }
  conc_likelihood <- conc_likelihood_factory(cqs, model$intercept, model$slope,
                                             model$sigma)
  mle <- nlm(function(conc) {-log(conc_likelihood(conc))},
             exp((mean(cqs, na.rm = TRUE) - model$intercept)/model$slope))$estimate
  threshold <- conc_likelihood(mle) / 1E9
  lb <- uniroot(function(conc) {conc_likelihood(conc) - threshold}, lower = 0,
                upper = mle)$root
  ub <- uniroot(function(conc) {conc_likelihood(conc) - threshold}, lower = mle,
                upper = 2 * mle, extendInt = "downX")$root
  lb <- max(lb, ub/2001)
  grid <- seq(lb, ub, length.out = 1001)
  pdf <- sapply(grid, conc_likelihood)
  cdf <- cumsum(pdf)
  limits <- c((1 - level)/2, 1 - (1 - level)/2) * cdf[1001]
  bounds <- findInterval(limits, cdf) + 1
  grid[bounds] + (0.5 - (cdf[bounds] - limits)/pdf[bounds]) * (ub - lb) / 1000
}
