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
  slope     = model$slope / log(10)
  sigma     = model$sigma

  # --- Identify non-detects (if any)
  detects <- which(!is.nan(cqs))
  num_non_detects <- length(cqs) - length(detects)
  cqs <- cqs[detects]

  #bouds of N0 values for which normal densities have been stored
  N0max <- N0min <- max(round(exp((mean(cqs) - intercept)/slope)), 1)

  #function for generating normal densities
  gen_norm_densities <- function(N0) {
    dnorm(cqs, mean = intercept + slope * log(N0), sd = sigma)
  }

  #stored normal densities, for reuse between calls of generated function
  norm_densities <- gen_norm_densities(N0max)

  function(concentration) {

    #check for negative concentration
    if(concentration < 0) {return(Inf)}

    #calculate bounds on N0 for summation
    N0start <- max(qpois(1e-15, concentration), 1)
    N0end   <- max(qpois(1e-15, concentration, lower.tail = FALSE), N0start + 1)

    #calculate and store new normal densities if required, and update bounds
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

    #values to sum over
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
#' Given a list of Cq values from a set of technical replicates and
#' a fitted ESC model, generates a maximum likelihood estimate
#' of the concentrations.
#'
#' @param cqs Numeric vector of Cq values from sample replicates, non-detects
#' coded as \code{NaN}.
#' @param model \code{esc} object representing a fitted model to use for estimation.
#' @param approximate Logical. If \code{TRUE} (the default), a faster but potentially
#' less accurate approximation for the likelihood function will be used at high
#' concentrations.
#'
#' @return The maximum likelihood estimation of the concentrations.
#'
#' @export
#'
#' @examples
#'
conc_mle <- function(cqs, model, approximate = TRUE) {

  # --- Input checks

  if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(!inherits(model, 'esc')) {stop("model is not an esc object")}
  if(!is.logical(approximate)) {stop("approximate must be logical")}
  if(all(is.nan(cqs))) {return(0)}

  # --- Compute MLE
  if(approximate){
    conc_log_like <- function(conc) {
      log_likelihood_est(conc, cqs, model$intercept, model$slope / log(10),
                         model$sigma)
    }
  }
  else {
    conc_log_like <- conc_log_likelihood_factory(cqs, model)
  }

  res <- suppressWarnings(nlm(f = conc_log_like,
                              p = exp((mean(cqs, na.rm = TRUE) - model$intercept)
                                      * log(10)/model$slope)))

  return(res$estimate)
}

#' @title Generate credible interval for concentration analytically
#'
#' @description Given list of Cq values from a set of technical replicates
#' and a fitted ESC model, generates credible interval of desired level
#' through numerical integration of the posterior distribution.
#'
#' @param cqs Numeric vector of Cq values from sample replicates, non-detects
#' coded as \code{NaN}.
#' @param model A \code{esc} object representing a fitted model to use for estimation.
#' @param level Desired credibility level, defaults to 0.95.
#' @param approximate Logical. If \code{TRUE} (the default), a faster but potentially
#' less accurate approximation for the likelihood function will be used at high
#' concentrations.
#'
#' @return A \code{conc_int} object, containing a list interval specifying the lower and
#' upper bounds of desired credible interval, as well as the maximum likelihood
#' concentration estimate, and a data frame distribution, containing values of
#' the probability density function  and the cumulative density function
#' of the posterior distribution evaluated at specified points.
#'
#'
#' @export
#'
#' @examples
#'
conc_interval <- function(cqs,
                          model,
                          level = 0.95,
                          approximate = TRUE) {

  # --- Input checks
  if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(! inherits(model, 'esc')) {stop("model is not an esc object")}
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}
  if(!is.logical(approximate)) {stop("approximate must be logical")}

  # --- Deal with case of all non-detects separately
  if(all(is.nan(cqs))) {
    n <- length(cqs)
    grid <- seq(from = 0, to = log(10) * 9 / n, length.out = 1001)
    return(new_conc_int(mle = 0, interval = c(0, -log(1 - level)/n), grid = grid,
                        pdf = sapply(grid, function(x) {n * exp(-x * n)}),
                        cdf = sapply(grid, function(x) {1 - exp(-x * n)}),
                        level = level))
  }

  # --- Compute MLE
  if(approximate){
    conc_log_like <- function(conc) {
      log_likelihood_est(conc, cqs,
                         intercept = model$intercept,
                         slope = model$slope / log(10),
                         sigma = model$sigma)
    }
  }
  else {
    conc_log_like <- conc_log_likelihood_factory(cqs, model)
  }

  mle <- suppressWarnings(stats::nlm(conc_log_like,
                                     exp((mean(cqs, na.rm = TRUE) - model$intercept) *
                                     log(10)/model$slope)))$estimate

  # Compute bounds for numerical integration
  threshold <- conc_log_like(mle) + 9 * log(10)
  root_est_factor <- exp(9.21 * model$sigma / abs(model$slope)) #empirical factor
  #for getting reasonable estimates of roots before using uniroot. Helps improve efficiency
  lb <- stats::uniroot(function(conc) {conc_log_like(conc) - threshold}, upper = mle,
                lower = mle / root_est_factor, extendInt = "downX")$root
  ub <- stats::uniroot(function(conc) {conc_log_like(conc) - threshold}, lower = mle,
                upper = mle * root_est_factor, extendInt = "upX")$root
  lb <- max(lb, ub/2001)

  # --- Perform numerical integration
  grid <- seq(lb, ub, length.out = 1001)
  pdf <- exp(conc_log_like(mle) - sapply(grid, conc_log_like))
  cdf <- cumsum(pdf) * (ub - lb) / 1000

  # --- Construct and return interval
  limits <- c((1 - level)/2, 1 - (1 - level)/2) * cdf[1001]
  interval <- stats::approx(x = cdf, y = grid, xout = limits, rule = 2, ties = "ordered")$y

  res = new_conc_int(mle, interval, grid, pdf, cdf, level)

  return(res)
}

#' Compute Multiple Concentration Credible Intervals
#'
#' Computes Baysian credible intervals and maximum likelihood estimates for the
#' concentrations of multiple samples at once, given Cq data for each sample and
#' a single \code{esc} model object and specified credible level to use for all the
#' intervals.
#'
#' @param cq_data data frame with a \code{sample} column specifying the names of the
#' samples from which reactions were generated, and a \code{cqs} column containing the
#' corresponding Cq values. Cq values must be numeric, with non-detects encoded
#' as \code{NaN}.
#' @param model A \code{esc} object representing fitted model to use for estimation.
#' @param level Numeric. Desired credibility level, defaults to 0.95.
#' @param approximate Logical. If \code{TRUE} (the default), a faster but potentially
#' less accurate approximation for the likelihood function will be used at high
#' concentrations.
#'
#' @return Data frame with a \code{sample} column specifying the sample index,
#' and \code{lower},
#' \code{mle}, and \code{upper} columns specifying the interval
#' lower bound, maximum likelihood
#' estimate, and interval lower bound for the gene concentration in that sample.
#' One row is generated for each unique sample index in the original data frame.
#'
#' @export
#'
#' @examples
#'
#'
multi_interval <- function(cq_data, model, level = 0.95, approximate = TRUE) {
  # input checks
  if(!"sample" %in% names(cq_data)) {
    stop("cq_data must contain sample column")
  }
  if(!"cqs" %in% names(cq_data)) {
    stop("cq_data must contain cqs column")
  }
  if(!all(is.numeric(cq_data$cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cq_data$cqs) | (cq_data$cqs >= 0 & is.finite(cq_data$cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(!inherits(model, 'esc')) {stop("model is not an esc object")}
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}
  if(!is.logical(approximate)) {stop("approximate must be logical")}

  #compute intervals
  res <- aggregate(cqs ~ sample,
                   data = cq_data,
                   na.action = NULL,
                   FUN =
                     function(cqs) {
                       c(conc_interval(cqs, model, level, approximate)$interval)
                     })

  #reshape result
  res <- data.frame(sample = res[,1], res[,2])
  res[] <- lapply(res, unlist)
  return(res)
}

#' MCMC Estimation of Concentrations Using ESC Model
#'
#' Given list of Cq values from a set of technical replicates and fitted ESC
#' model, generates credible interval of desired level by approximating the
#' posterior distribution for concentration using MCMC sampling. Requires
#' installation of JAGS and the rjags and runjags packages
#'
#' @param cqs Numeric vector of Cq values from sample replicates, non-detects
#' coded as NaN
#' @param model esc object representing fitted model to use for estimation
#' @param level Desired credible level, defaults to 0.95
#'
#' @return A list of two elements. The first, named interval, is a list
#' containing the bounds of the credible interval, along with the mean and median
#' of the generated posterior samples. The second, named mcmc_samples, is the
#' object produced by the runjags package representing the generated samples
#' @export
#'
#' @examples
conc_mcmc <- function(cqs, model, level = 0.95){
  #software checks
  if(!requireNamespace("runjags", quietly = TRUE)) {
    stop("package 'runjags' must be installed to use mcmc functions")
  }
  test <- runjags::testjags(silent = TRUE)
  if(!test$JAGS.found) {
    stop("JAGS software must be installed to use mcmc functions")
  }

  # --- Input checks
  if(!all(is.numeric(cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cqs) | (cqs >= 0 & is.finite(cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(class(model) != "esc") {stop("model is not an esc object")}
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}

  #unpack model
  intercept <- model$intercept
  slope <- model$slope / log(10)
  sigma <- model$sigma

  #pre-process inputs
  n <- length(cqs)
  nds <- is.nan(cqs)
  cqs[which(nds)] <- NA

  #run mcmc sampling with jags
  results <- runjags::run.jags(system.file("jags-models", "conc-model.txt",
                                           package = "ecurve"),
                               data = list(n = n, cq = cqs, ND = as.numeric(nds),
                                           k = 1, index = rep(1, n),
                                           alpha = intercept, beta = slope,
                                           sigma = sigma),
                               monitor = c("conc"),
                               thin = 4,
                               n.chains = 2,
                               inits = function() {
                                 mod_cqs <- cqs
                                 mod_cqs[which(nds)] <- intercept
                                 list(conc = exp((mean(mod_cqs) - intercept)/slope))
                               })

  #check effective sample size
  if(results$summaries[1,9] < 100) {
    warning(paset0("Low effective sample size (less than 100). Results may be",
            "unreliable - consider using numerical integration instead."))
  }

  #process and return results
  results <- runjags::add.summary(results, confidence = c(level))
  interval <- as.list(results$summaries[,c(1, 2, 4, 3)])
  names(interval) <- c("lower", "median", "mean", "upper")
  return(list(interval = interval, mcmc_samples = results))
}

#' MCMC Estimation of Concentrations for Multiple Samples at Once
#'
#' Computes Baysian credible intervals for the concentrations of multiple samples
#' at once, given Cq data for each sample and a single esc model object and
#' specified credible level to use for all the intervals. Uses MCMC sampling to
#' approximate the posterior distribution. Requires installation of JAGS and the
#' rjags and runjags packages
#'
#' @param cq_data data frame with an sample column specifying the names of the
#' samples from which reactions were generated, and a cqs column containing the
#' corresponding Cq values. Cq values must be numeric, with non0detects encoded
#' as NaN.
#' @param model esc object representing fitted model to use for estimation
#' @param level Desired credible level, defaults to 0.95
#'
#' @return A list of two elements. The first, named interval, is a data frame
#' containing the bounds of the credible intervals, along with the means and
#' medians of the generated posteriorsamples. The second, named mcmc_samples,
#' is the object produced by the runjags package representing the generated
#' samples.
#' @export
#'
#' @examples
multi_conc_mcmc <- function(cq_data, model, level = 0.95) {
  #software checks
  if(!requireNamespace("runjags", quietly = TRUE)) {
    stop("package 'runjags' must be installed to use mcmc functions")
  }
  test <- runjags::testjags(silent = TRUE)
  if(!test$JAGS.found) {
    stop("JAGS software must be installed to use mcmc functions")
  }

  #input checks
  if(!"sample" %in% names(cq_data)) {
    stop("cq_data must contain sample column")
  }
  if(!"cqs" %in% names(cq_data)) {
    stop("cq_data must contain cqs column")
  }
  if(!all(is.numeric(cq_data$cqs))) {stop("cqs must be numeric")}
  if(!all(is.nan(cq_data$cqs) | (cq_data$cqs >= 0 & is.finite(cq_data$cqs)))) {
    stop("cqs must be non-negative real numbers or NaN")
  }
  if(class(model) != "esc") {stop("model is not an esc object")}
  if(!is.numeric(level)) {stop("level must be numeric")}
  if(level > 1 | level < 0) {stop("level must be between 0 and 1")}

  #unpack model
  intercept <- model$intercept
  slope <- model$slope / log(10)
  sigma <- model$sigma

  #pre-process inputs
  nds <- is.nan(cq_data$cqs)
  cq_data$cqs[which(nds)] <- NA
  samples <- unique(cq_data$sample)

  #run mcmc sampling with jags
  results <- runjags::run.jags(system.file("jags-models", "conc-model.txt",
                                           package = "ecurve"),
                               data = list(n = dim(cq_data)[1], cq = cq_data$cqs,
                                           ND = as.numeric(nds), k = length(samples),
                                           index = match(cq_data$sample, samples),
                                           alpha = intercept, beta = slope,
                                           sigma = sigma),
                               monitor = c("conc"),
                               thin = 4,
                               n.chains = 2,
                               inits = function() {
                                 mod_cq_data <- cq_data
                                 mod_cq_data$cqs[which(nds)] <- intercept
                                 mean_cqs <- aggregate(cqs ~ sample,
                                                       data = mod_cq_data,
                                                       FUN = mean)$cqs
                                 list(conc = exp((mean_cqs - intercept)/slope))
                               })

  #check effective sample sizes
  if(any(results$summaries[,9] < 100)) {
    warning(paste0("Low effective sample sizes (less than 100) for concentration",
                   " estimates for samples ",
                   paste0(samples[which(results$summaries[,9] < 100)],
                          collapse = ", "),
                   ". Results may be unreliable - consider using numerical",
                   " integration instead"))
  }

  #process and return results
  results <- runjags::add.summary(results, confidence = c(level))
  intervals <- data.frame(sample = samples,
                          lower = results$summaries[,1],
                          median = results$summaries[,2],
                          mean = results$summaries[,4],
                          upper = results$summaries[,3])
  rownames(intervals) <- NULL
  return(list(intervals = intervals, mcmc_samples = results))
}
