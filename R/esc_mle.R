get_bounds <- function(conc, cqs, intercept, slope, sigma) {
  density_term <- function(N0) {
    -sum(dnorm(cqs, mean = intercept + slope * log(N0), sd = sigma, log = TRUE),
      dpois(N0, conc, log = TRUE))
  }
  density_interp <- function(N0) {
    weight <- N0 - floor(N0)
    density_term(floor(N0)) * (1 - weight) + density_term(ceiling(N0)) * weight
  }
  cq_conc_ests <- exp((cqs - intercept)/slope)
  lb <- max(1, min(conc, cq_conc_ests))
  ub <- max(conc, cq_conc_ests, lb + 1)
  optN0 <- optimize(density_interp, lower = lb, upper = ub, tol = 0.01)
  threshold <- optN0$objective + 6 * log(10)
  if (density_term(1) < threshold) N0min <- 1
  else {
    N0min <- floor(uniroot(function(conc) {density_interp(conc) - threshold},
                           lower = 1, upper = optN0$minimum, tol = 0.01)$root)
  }
  N0max <- ceiling(uniroot(function(conc) {density_interp(conc) - threshold},
                           lower = optN0$minimum, upper = 2 * optN0$minimum,
                           extendInt = "upX", tol = 0.01)$root)
  return(c(N0min, max(N0min + 1, N0max)))
}

likelihood_est <- function(conc, cqs, intercept, slope, sigma) {
  dens_est <- function(N0) {
    prod(dnorm(cqs, mean = intercept + slope * log(N0), sd = sigma),
         dnorm(N0, mean = conc, sd = sqrt(conc)))
  }
  bounds <- get_bounds(conc, cqs, intercept, slope, sigma)
  integrate(Vectorize(dens_est), lower = bounds[1], upper = bounds[2], abs.tol = 1e-100)$value
}

#' Calculate cq PDFs under ESC Model
#'
#' Given ESC model parameters and vector of specified concentrations, calculates
#' ESC PDF for given Cq values (one for each concentration)
#'
#' @param concentrations numeric vector of specified concentrations
#' @param cqs numeric vector containing Cq value at which to compute pdf for
#' each specificed concentration, with non-detects coded as NaN
#' @param intercept intercept parameter of ESC model
#' @param slope slope parameter of ESC model
#' @param sigma sigma parameter of ESC model
#'
#' @return numeric vector of calculated pdf values
#'
esc_probdens <- function(concentrations, cqs, intercept, slope, sigma) {

  # DC: we may want to merge `intercept`, `slope` and `sigma`
  # in one single _named_ vector `params`

  bounds <- mapply(get_bounds, conc = concentrations, cqs = cqs, MoreArgs =
                     list(intercept = intercept, slope = slope, sigma = sigma))
  N0s <- unique(unlist(mapply(FUN = seq, sort(bounds[1,]), sort(bounds[2,]))))
  mean_cqs <- intercept + slope * log(N0s)
  N0starts <- match(bounds[1,], N0s)
  N0ends <- match(bounds[2,], N0s)

  probdens <- function(conc, cq, N0start, N0end) {

    tmp1 <- dpois(x      = N0s[N0start:N0end],
                 lambda = conc)

    tmp2 <- dnorm(x    = cq,
                 mean = mean_cqs[N0start:N0end],
                 sd   = sigma)

    return( sum( tmp1 * tmp2 ) )
  }
  res <- mapply(probdens, conc = concentrations, cq = cqs, N0start = N0starts,
                N0end = N0ends)
  return(res)
}

#' Log Likelihood function for ESC Modle fitting
#'
#' @param params  numeric vector containing parameters of the model: first
#' element is intercept, second is slope and third is sigma
#' @param concentrations numeric vector of concentrations that model is being
#' fit to
#' @param cqs numeric vector containing Cq value for each specified concentration,
#' with non-detects coded as NaN
#'
#' @return negative log likelihood evalueated at the given parameters
#'
esc_log_likelihood <- function(params, concentrations, cqs) {
  if(params[3] < 0) {return(Inf)}

  tmp <- esc_probdens(concentrations, cqs,
                     intercept = params[1],
                     slope = params[2],
                     sigma = params[3])
  res <- -sum(log(tmp))
  return(res)
}


#' Calculate Theoretical Quantiles for Cq Values
#'
#' Given ESC model parameters, quantile level alpha, and vector of concentrations
#' evaluates that alpha-th quantile of the theoretical distribution of cq values
#' at each of the specified concentrations
#'
#' @param concentrations numeric vector of concentrations at which to evaluate
#' specified quantile
#' @param alpha quantile level
#' @param intercept intercept parameter of ESC model
#' @param slope slope parameter of ESC model
#' @param sigma sigma parameter of ESC model
#'
#' @return Numeric vector containing the calculated quantiles
cq_quantile <- function(alpha,
                        concentrations,
                        intercept,
                        slope,
                        sigma) {
  N0mins <- pmax(qpois(1E-9, concentrations), 1)
  N0maxes <- pmax(qpois(1E-9, concentrations, lower.tail = FALSE), N0mins + 1)
  N0s <- unique(unlist(mapply(FUN = seq, sort(N0mins), sort(N0maxes))))
  means <- intercept + slope * log(N0s)
  N0starts <- match(N0mins, N0s)
  N0ends <- match(N0maxes, N0s)

  quant <- function(conc, N0start, N0end) {
    uniroot(function(cq) {
      sum(pnorm(cq, mean = means[N0start:N0end], sd = sigma) *
            dpois(N0s[N0start:N0end], conc)) - alpha * (1 - exp(-conc))},
      lower = qnorm(alpha, mean = means[N0end], sd = sigma),
      upper = qnorm(alpha, mean = means[N0start], sd = sigma))$root
  }
  res = mapply(quant, conc = concentrations, N0start = N0starts, N0end = N0ends)
  return(res)
}



#' Fit ESC Model Using MLE
#'
#' @param esc_data Data frame containing data to be used to fit the model. Must
#' contain a column named "concentrations" with known sample concentrations, and
#' a column named "cqs" with corresponding Cq values. Non-detects should be
#' encoded by a Cq value of NaN
#' @param CI Numeric. Width of the confidence interval.
#'
#' @return esc object representing fitted model
#' @export
#'
#' @example
#'
esc_mle <- function(esc_data, CI = 0.95) {

  # --- Inputs Checks

  if(!is.data.frame(esc_data)) {stop("esc_data must be a data frame")}
  if(!"concentrations" %in% names(esc_data)) {
    stop("esc_data must contain concentrations column")
  }
  if(!"cqs" %in% names(esc_data)) {
    stop("esc_data must contain cqs column")
  }

  concentrations <- esc_data[,"concentrations"]
  cqs <- esc_data[,"cqs"]
  if(!all(is.numeric(concentrations))) {stop("concentrations must be numeric")}
  if(!all(concentrations >= 0 & is.finite(concentrations))) {
    stop("concentrations must be non-negative real numbers")
  }
  if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(length(concentrations) != length(cqs)) {
    stop("concentrations and cqs must be the same length")
  }

  if(!is.numeric(CI)) {stop("CI must be numeric")}
  if(CI > 1 | CI < 0) {stop("CI must be between 0 and 1")}

  # --- Filter out the non-detects
  detects <- !is.nan(cqs)
  concentrations <- concentrations[detects]
  cqs <- cqs[detects]

  # --- Naive standard curve fitting, i.e., linear regression.
  # Will also be used as initial parameter values
  # for the ESC model fit (below).
  naive_sc <- lm(cqs ~ log(concentrations))

  # --- ESC model fitting

  init <- c(unname(coef(naive_sc)), sigma(naive_sc))

  res <- nlm(f = esc_log_likelihood,
             p = init,
             concentrations = concentrations,
             cqs = cqs)

  intercept = res$estimate[1]
  slope     = res$estimate[2]
  sigma     = res$estimate[3]

  # --- Quantiles for Cq values

  alphas = c(0.5 - CI/2, 0.5, 0.5 + CI/2)

  cq.quants = lapply(X = alphas, FUN = cq_quantile,
                     concentrations = concentrations,
                     intercept = intercept,
                     slope = slope,
                     sigma = sigma)

  names(cq.quants) = c('low', 'median', 'high')


  m = new_esc(intercept = res$estimate[1],
              slope = res$estimate[2],
              sigma = res$estimate[3],
              cq_quantiles = as.data.frame(cq.quants),
              cq_quantiles_CI = CI,
              data = data.frame(concentrations = concentrations,
                                cqs = cqs))
  return(m)
}
