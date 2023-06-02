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

  # Determine N0 values at which to perform evaluation
  N0start <- max(stats::qpois(1e-15, conc), 1)
  N0end   <- max(stats::qpois(1e-15, conc, lower.tail = FALSE), N0start + 1)
  granularity <- ceiling((N0end - N0start)/100)
  N0s <- seq(N0start, N0end, by = granularity)

  # Calculate intermediate results
  norm_dens <- sapply(N0s, function(N0){
    stats::dnorm(x    = cqs,
                 mean = intercept + slope * log(N0),
                 sd   = sigma)
  })
  pois_dens <- dpois(N0s, conc)

  # Calculate final log likelihood
  tmp = norm_dens %*% pois_dens
  tmp = tmp[tmp>0]

  res = sum(-log(tmp * granularity), num_non_detects * conc)

  return(res)
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

  N0mins <- pmax(stats::qpois(1E-15, concentrations), 1)
  N0maxes <- pmax(stats::qpois(1E-15, concentrations, lower.tail = FALSE), N0mins)
  N0s <- unique(unlist(mapply(FUN = seq, sort(N0mins), sort(N0maxes))))
  mean_cqs <- intercept + slope * log(N0s)
  N0starts <- match(N0mins, N0s)
  N0ends <- match(N0maxes, N0s)

  probdens <- function(conc, cq, N0start, N0end) {

    tmp1 <- stats::dpois(x      = N0s[N0start:N0end], lambda = conc)

    tmp2 <- stats::dnorm(
      x    = cq,
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
esc_log_likelihood <- function(params,
                               concentrations,
                               cqs,
                               approximate = TRUE) {

  if(params[3] < 0) {return(Inf)}

  # Parameterization to constraint
  # the Efficiency between 0 and 1:
  slope = slope_from_theta(theta = params[2]) / log(10)


  if (approximate) {
    res <- sum( mapply(
      FUN      = log_likelihood_est,
      conc     = concentrations,
      cqs      = cqs,
      MoreArgs = list(intercept = params[1],
                      slope     = slope,
                      sigma     = params[3])
      ))
  }
  else {
    tmp <- esc_probdens(
      concentrations = concentrations,
      cqs            = cqs,
      intercept      = params[1],
      slope          = slope,
      sigma          = params[3])

    res <- -sum(log(tmp))
  }
  return(res)
}

#' @title Check the input of the function `esc_mle()`.
#'
#' @param esc_data Dataframe. See \code{esc_mle()}.
#' @param approximate Logical. See \code{esc_mle()}.
#' @param assumeND Logical. See \code{esc_mle()}.
#'
#' @return A list containing formated \code{concentrations} and \code{cqs}.
#'
check_input_esc_mle <- function(esc_data, approximate, assumeND) {

  # Handle tibbles:
  q = class(esc_data)
  if(q[1] == 'tbl_df') esc_data = as.data.frame(esc_data)

  if(!is.data.frame(esc_data)) {stop("esc_data must be a data frame")}
  if(!"concentrations" %in% names(esc_data)) {
    stop("esc_data must contain concentrations column")
  }
  if(!"cqs" %in% names(esc_data)) {
    stop("esc_data must contain cqs column")
  }

  concentrations <- esc_data[,"concentrations"]
  cqs            <- esc_data[,"cqs"]

  if(!all(is.numeric(concentrations))) {stop("concentrations must be numeric")}
  if(!all(concentrations >= 0 & is.finite(concentrations))) {
    stop("concentrations must be non-negative real numbers")
  }

  # Translate non-numeric Cqs into NaN which
  # is the coding for non-detect in this library:
  if(!is.logical(assumeND)) {stop("`assumeND` must be logical.")}
  if(assumeND){
    tmp = suppressWarnings(as.numeric(cqs))
    n.na = sum(is.na(tmp))
    if(n.na > 0){
      warning('`esc_mle()`: ',n.na,
              ' Cq values were assumed non-detects because they were inputed as non-numeric.')
    }
    tmp[is.na(tmp)] <- NaN
    cqs <- tmp
  }
  if(!assumeND){
    if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
    if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
      stop("cqs must be non-negative real numbers or NaN.")
    }
  }

  if(!all(cqs[!is.nan(cqs)]>0)) stop("cqs must be non-negative")

  if(!is.logical(approximate)) {stop("`approximate` must be logical.")}

  return(list(concentrations = concentrations, cqs=cqs))
}

#' @title Fit ESC Model Using Maximum Likelihood Estimation.
#'
#' @description Fit the parameters of an ESC model to qPCR standard curve data. The ESC model is
#' described in
#'  \href{https://doi.org/10.3389/fmicb.2023.1048661}{Schmidt et al. 2022}.
#'
#' @param esc_data Data frame containing data to be used to fit the model. Must
#' contain a column named \code{concentrations} with known sample concentrations, and
#' a column named \code{cqs} with corresponding Cq values. Non-detects should be
#' encoded by a Cq value of \code{NaN}.
#' @param approximate Logical. If \code{TRUE} (the default), a faster but potentially
#' less accurate approximation for the likelihood function will be used at high
#' concentrations.
#' @param assumeND Logical. If \code{TRUE} (the default), any non numerical
#' Cq values in the dataframe \code{esc_data} is assumed a non-detect (ND).
#'
#' @return An \code{esc} object representing fitted model.
#'
#' @export
#'
#' @examples
#' esc_data = data.frame(
#'      concentrations = c(1,1,10, 10, 100, 500, 500),
#'      cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
#' mod = esc_mle(esc_data)
#'
esc_mle <- function(esc_data,
                    approximate = TRUE,
                    assumeND = TRUE) {

  # --- Inputs Checks
  chk            = check_input_esc_mle(esc_data, approximate, assumeND)
  cqs            = chk$cqs
  concentrations = chk$concentrations

  # --- Filter out the non-detects
  detects <- !is.nan(cqs)
  concentrations <- concentrations[detects]
  cqs <- cqs[detects]

  # --- Naive standard curve fitting, i.e., linear regression.
  # Will also be used as initial parameter values
  # for the ESC model fit (below).
  naive_sc <- stats::lm(cqs ~ log(concentrations))

  # --- ESC model fitting

  # Note:
  # coef[1] = intercept
  # coef[2] = slope
  k = stats::coef(naive_sc)

  # Parameterization to constraint
  # the Efficiency between 0 and 1:
  theta.init = theta_from_slope(min(-1.1/log10(2), k[2]))

  # initial guess for optimization
  # Note that: init = c(intercept, theta, sigma)
  init <- c(k[1],
            theta.init,
            stats::sigma(naive_sc))

  res <- suppressWarnings(
    stats::nlm(
      f              = esc_log_likelihood,
      p              = init,
      concentrations = concentrations,
      cqs            = cqs,
      approximate    = approximate,
      # The argument `typsize` specifies the
      # expected "size" (order of magnitude)
      # of the solution (ie the parameter values
      # that minimize the likelihood).
      # In the context of qPCR, it is expected
      # that the intercept is in the 30s/40s,
      # that `sigma` should be around or smaller
      # than 1. Expected value for `theta` is more
      # uncertain, so set to "1" for lack of a
      # better guess. The precise values of `typsize`
      # should not matter tremendously, however
      # setting some reasonable values should
      # improve the robustness of the optimization.
      typsize        = c(35, 1, 0.1)
    )
  )

  # Extract the optimal slope from the optimal theta
  # (theta is used in the optimization to constraint
  # the implied Efficiency between 0 and 1)

  slope = slope_from_theta(res$estimate[2])

  m <- new_esc(
    intercept = res$estimate[1],
    slope     = slope,
    sigma     = res$estimate[3],
    data      = data.frame(concentrations = concentrations,
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
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}

  # --- Filter out the non-detects
  detects <- !is.nan(cqs)
  concentrations <- concentrations[detects]
  cqs <- cqs[detects]



  # --- run mcmc sampling with jags


  init_for_jags <- function(){

    naive_sc <- stats::lm(cqs ~ log(concentrations))

    res = list(
      alpha     = stats::coef(naive_sc)[1],
      eff       = min(1, exp(-1/stats::coef(naive_sc)[2]) - 1),
      log2sigma = log(stats::sigma(naive_sc), base = 2)
    )
    return(res)
  }


  results <- runjags::run.jags(
    model = system.file("jags-models", "esc-model.txt", package = "ecurve"),
    data  = list(
      n = length(cqs),
      cq = cqs,
      conc = concentrations),
    monitor  = c("alpha", "beta", "eff", "sigma"),
    thin     = 4,
    n.chains = 2,
    inits = init_for_jags
  )

  # check effective sample sizes
  if(any(results$summaries[,9] < 100)) {
    warning(paste0("Low effective sample sizes (less than 100) for estimates of",
                   "parametes ",
                   paste0(rownames(results$summaries)[which(results$summaries[,9] < 100)],
                          collapse = ", "),
                   ". Results may be unreliable."))
  }

  # process and return results
  results <- runjags::add.summary(results, confidence = c(level))

  extract_int <- function(param) {
    interval <- as.list(results$summaries[param,c(1, 2, 4, 3)])
    names(interval) <- c("lower", "median", "mean", "upper")
    return(interval)
  }

  return(list(
    intercept    = extract_int("alpha"),
    slope        = extract_int("beta"),
    eff          = extract_int("eff"),
    sigma        = extract_int("sigma"),
    mcmc_samples = results))
}

