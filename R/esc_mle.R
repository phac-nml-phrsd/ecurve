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

  N0starts <- pmax(qpois(1E-15, concentrations), 1)
  N0ends   <- qpois(1E-15, concentrations, lower.tail = FALSE)
  N0min    <- min(N0starts)
  N0max    <- max(N0ends)
  N0s      <- N0min:N0max
  mean_cqs <- intercept + slope * log(N0s)
  probdens <- function(conc, cq, N0start, N0end) {

    tmp1 <- dpois(x      = N0s[N0start:N0end],
                 lambda = conc)

    tmp2 <- dnorm(x    = cq,
                 mean = mean_cqs[N0start:N0end],
                 sd   = sigma)

    return( sum( tmp1 * tmp2 ) )
  }
  res <- mapply(probdens, conc = concentrations, cq = cqs,
                N0start = N0starts - N0min + 1, N0end = N0ends - N0min + 1)
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


#' Fit ESC Model Using MLE
#'
#' @param concentrations Numeric vector of known sample concentrations
#' @param cqs Numeric vector of corresponding Cq values, with non-detects coded
#' as NaN
#'
#' @return esc object representing fitted model
#' @export
#'
#' @example
#'
esc_mle <- function(concentrations, cqs) {

  # --- Inputs Checks

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

  # DC: why not `optim()`??
  res <- nlm(f = esc_log_likelihood,
             p = init,
             concentrations = concentrations,
             cqs = cqs)

  new_esc(intercept = res$estimate[1],
          slope = res$estimate[2],
          sigma = res$estimate[3])
}
