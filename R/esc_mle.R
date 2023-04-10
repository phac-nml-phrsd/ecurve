#' Faster Estimation of log likelihood values for ESC model
#'
#' Estimates likelihood of given set of observed Cq values given a concentration
#' and ESC model parameters. For efficiency, ensures that a maximum of 100 N0
#' values are considered when the computation is performed. N0 values used for
#' computation are equally spaced out between the computed lower and upper bounds,
#' and are used to estimate contribution of nearby N0 values not included in the
#' computation.
#'
#' @param conc numeric concentration
#' @param cqs numeric vector containing observed Cq values to estimate the
#' likelihood of, with non-detects encoded as NaN
#' @param intercept intercept parameter of ESC model
#' @param slope slope parameter of ESC model
#' @param sigma sigma parameter of ESC model
#'
#' @return estimated log likelihood values
#'
log_likelihood_est <- function(conc, cqs, intercept, slope, sigma) {
  if(conc < 0) return(Inf)

  # --- Identify non-detects (if any)
  detects <- which(!is.nan(cqs))
  num_non_detects <- length(cqs) - length(detects)
  cqs <- cqs[detects]

  #determine N0 values at which to perform evaluation
  N0start <- max(qpois(1e-15, conc), 1)
  N0end <- max(qpois(1e-15, conc, lower.tail = FALSE), N0start + 1)
  granularity <- ceiling((N0end - N0start)/100)
  N0s <- seq(N0start, N0end, by = granularity)

  #calculate intermediate results
  norm_dens <- sapply(N0s, function(N0){
    dnorm(cqs, mean = intercept + slope * log(N0), sd = sigma)
  })
  pois_dens <- dpois(N0s, conc)

  #calculate final log likelyhood
  sum(-log((norm_dens %*% pois_dens) * granularity), num_non_detects * conc)
}

#' Calculate cq PDFs under ESC Model
#'
#' Given ESC model parameters and vector of specified concentrations, calculates
#' ESC PDF for given Cq values (one for each concentration)
#'
#' @param concentrations numeric vector of specified concentrations
#' @param cqs numeric vector containing Cq value at which to compute pdf for
#' each specified concentration, with non-detects coded as NaN
#' @param intercept intercept parameter of ESC model
#' @param slope slope parameter of ESC model
#' @param sigma sigma parameter of ESC model
#'
#' @return numeric vector of calculated pdf values
#'
esc_probdens <- function(concentrations, cqs, intercept, slope, sigma) {

  # DC: we may want to merge `intercept`, `slope` and `sigma`
  # in one single _named_ vector `params`

  N0mins <- pmax(qpois(1E-15, concentrations), 1)
  N0maxes <- pmax(qpois(1E-15, concentrations, lower.tail = FALSE), N0mins)
  N0s <- unique(unlist(mapply(FUN = seq, sort(N0mins), sort(N0maxes))))
  mean_cqs <- intercept + slope * log(N0s)
  N0starts <- match(N0mins, N0s)
  N0ends <- match(N0maxes, N0s)

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
#' @param approximate logical. If TRUE (the default), a faster but potentially
#' less accurate approximation for the likelihood function will be used at high
#' concentrations
#'
#' @return negative log likelihood evaluated at the given parameters
#'
esc_log_likelihood <- function(params, concentrations, cqs, approximate = TRUE) {
  if(params[3] < 0) {return(Inf)}
  else if (approximate) {
    res <- sum(mapply(log_likelihood_est, conc = concentrations, cqs = cqs,
                      MoreArgs = list(intercept = params[1], slope = params[2],
                                      sigma = params[3])))
  }
  else {
    tmp <- esc_probdens(concentrations, cqs,
                        intercept = params[1],
                        slope = params[2],
                        sigma = params[3])
    res <- -sum(log(tmp))
  }
  return(res)
}



#' Fit ESC Model Using MLE
#'
#' @param esc_data Data frame containing data to be used to fit the model. Must
#' contain a column named "concentrations" with known sample concentrations, and
#' a column named "cqs" with corresponding Cq values. Non-detects should be
#' encoded by a Cq value of NaN
#' @param approximate logical. If TRUE (the default), a faster but potentially
#' less accurate approximation for the likelihood function will be used at high
#' concentrations
#'
#' @return esc object representing fitted model
#' @export
#'
#' @example
#'
esc_mle <- function(esc_data, approximate = TRUE) {

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
  if(!is.logical(approximate)) {stop("approximate must be logical")}

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
             cqs = cqs,
             approximate = approximate)

  intercept = res$estimate[1]
  slope     = res$estimate[2]
  sigma     = res$estimate[3]

  m = new_esc(intercept = res$estimate[1],
              slope = res$estimate[2] * log(10),
              sigma = res$estimate[3],
              data = data.frame(concentrations = concentrations,
                                cqs = cqs))
  return(m)
}

#' MCMC Estimation of ESC Model Parameters
#'
#' Performs Baysean estimation of the parameters of the ESC model by using MCMC
#' sampling to approximate the posterior distribution of those samples. For each
#' parameter, returns the mean and median of the posterior samples, as well as
#' the bounds of a credible interval which level specified by the optional level
#' parameter. Requires installation of JAGS and the rjags and runjags packages
#'
#' @param esc_data Data frame containing data to be used to fit the model. Must
#' contain a column named "concentrations" with known sample concentrations, and
#' a column named "cqs" with corresponding Cq values. Non-detects should be
#' encoded by a Cq value of NaN
#' @param level Desired credible level for interval estimates of parameters.
#' Defaults to 0.95
#'
#' @return A list of 5 elements. The first 4, named intercept, slope, eff, and
#' sigma, are lists summarizing the credible interval endpoints and posterior
#' mean and median for the giver parameter, while the 5th, named mcmc_samples,
#' is the object produced by the runjags package representing the generated
#' samples
#' @export
#'
#' @examples
esc_mcmc <- function(esc_data, level = 0.95) {
  #software checks
  if(!requireNamespace("runjags", quietly = TRUE)) {
    stop("package 'runjags' must be installed to use mcmc functions")
  }
  test <- runjags::testjags(silent = TRUE)
  if(!test$rjags.found) {
    stop("package 'rjags' must be installed to use mcmc functions")
  }
  if(!test$JAGS.found) {
    stop("JAGS software must be installed to use mcmc functions")
  }

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
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}

  # --- Filter out the non-detects
  detects <- !is.nan(cqs)
  concentrations <- concentrations[detects]
  cqs <- cqs[detects]

  #run mcmc sampling with jags
  results <- runjags::run.jags(system.file("jags-models", "esc-model.txt",
                                           package = "ecurve"),
                               data = list(n = length(cqs), cq = cqs,
                                           conc = concentrations),
                               monitor = c("alpha", "beta", "eff", "sigma"))

  #process and return results
  results <- runjags::add.summary(results, confidence = c(level))
  extract_int <- function(param) {
    interval <- as.list(results$summaries[param,c(1, 2, 4, 3)])
    names(interval) <- c("lower", "median", "mean", "upper")
    return(interval)
  }
  return(list(intercept = extract_int("alpha"), slope = extract_int("beta"),
              eff = extract_int("eff"), sigma = extract_int("sigma"),
              mcmc_samples = results))
}
