#' Factory Function to Create Concentration Log Likelihood Function
#'
#' Given ESC Model parameters and observed Cq data, returns function implementing
#' negative log likelihood calculation for concentration estimation
#'
#' @param cqs Numeric vector of observed Cq values, with non-detects coded as NaN
#' @param model `esc` object representing fitted model to use for estimation.
#'
#' @return Function that takes single input concentration and returns negative of
#' corresponding log likelihhod
conc_log_likelihood_factory <- function(cqs, model) {

  # --- unpack
  intercept = model$intercept
  slope     = model$slope
  sigma     = model$sigma

  # --- Identify non-detects (if any)
  detects <- which(!is.nan(cqs))
  num_non_detects <- length(cqs) - length(detects)
  cqs <- cqs[detects]

  N0max <- N0min <- max(round(exp((mean(cqs) - intercept)/slope)), 1)

  gen_norm_densities <- function(N0) {
    dnorm(cqs, mean = intercept + slope * log(N0), sd = sigma)
  }

  norm_densities <- gen_norm_densities(N0max)

  function(concentration) {

    if(concentration < 0) {return(Inf)}

    bounds <- get_bounds(concentration, cqs, intercept, slope, sigma)
    N0start <- bounds[1]
    N0end   <- bounds[2]

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

    # --- Likelihood components
    lk.norm = norm_densities[,N0s - N0min + 1]
    lk.pois = dpois(N0s, concentration)
    cond.detect = num_non_detects * concentration

    # --- Final likelihood
    res = sum(-log(lk.norm %*% lk.pois), cond.detect)
    return(res)
    # prod(norm_densities[,N0s - N0min + 1] %*% dpois(N0s, concentration),
    #      exp(-num_non_detects * concentration))
  }
}

#' Estimate concentration from Cq values
#'
#' Given list of Cq values from a set of technical replicates and fitted ESC
#' model, generates maximum likelihood estimate of concentration
#'
#' @param cqs Numeric vector of Cq values from sample replicates, non-detects
#' coded as NaN
#' @param model `esc` object representing fitted model to use for estimation.
#'
#' @return MLE of concentration
#' @export
#'
#' @examples
conc_mle <- function(cqs, model) {

  # --- Input checks

  if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(class(model) != "esc") {stop("model is not an esc object")}
  if(all(is.nan(cqs))) {return(0)}

  # --- Compute MLE

  conc_log_like <- conc_log_likelihood_factory(cqs,model)

  res = nlm(
    f = conc_log_like,
    p = exp((mean(cqs, na.rm = TRUE) - model$intercept)/model$slope) )

  return(res$estimate)
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
#' @return conc_int object, containing a list interval specifying the lower and
#' upper bounds of desired credible interval, as well as the maximum likelihood
#' concentration estimate, and a data frame distribution, containing values of
#' the pdf and cdf of the posterior distribution evaluated at specified points
#' @export
#'
#' @examples
conc_interval <- function(cqs, model, level = 0.95) {

  # --- Input checks
  if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(class(model) != "esc") {stop("model is not an esc object")}
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}

  # --- Deal with case of all non-detects seperately
  if(all(is.nan(cqs))) {
    n <- length(cqs)
    grid <- seq(from = 0, to = log(10) * 9 / n, length.out = 1001)
    return(new_conc_int(mle = 0, interval = c(0, -log(1 - level)/n), grid = grid,
                        pdf = sapply(grid, function(x) {n * exp(-x * n)}),
                        cdf = sapply(grid, function(x) {1 - exp(-x * n)})))
  }

  # --- Compute MLE
  conc_log_like<- conc_log_likelihood_factory(cqs, model)

  mle <- nlm(conc_log_like,
             exp((mean(cqs, na.rm = TRUE) - model$intercept)/model$slope))$estimate

  # Compute bounds for numerical intergration
  threshold <- conc_log_like(mle) + 9 * log(10)
  root_est_factor <- exp(4 * model$sigma / abs(model$slope))
  lb <- uniroot(function(conc) {conc_log_like(conc) - threshold}, upper = mle,
                lower = mle / root_est_factor, extendInt = "downX")$root
  ub <- uniroot(function(conc) {conc_log_like(conc) - threshold}, lower = mle,
                upper = mle * root_est_factor, extendInt = "upX")$root
  lb <- max(lb, ub/2001)

  # --- Perform numerical integration
  grid <- seq(lb, ub, length.out = 1001)
  pdf <- exp(conc_log_like(mle) - sapply(grid, conc_log_like))
  cdf <- cumsum(pdf) * (ub - lb) / 1000

  # --- Construct and return interval
  limits <- c((1 - level)/2, 1 - (1 - level)/2) * cdf[1001]
  interval <- approx(x = cdf, y = grid, xout = limits, rule = 2, ties = "ordered")$y

  res = new_conc_int(mle, interval, grid, pdf, cdf)

  return(res)
}
