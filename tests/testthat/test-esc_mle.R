test_that("esc MLE works", {
  concs <- 2^rep(seq(from = -2, to = 10), each = 3)
  cq1   <- 38
  sigma <- 0.5
  eff   <- 0.97

  cqs <- sim_cqs(concs, cq1 = cq1, eff = eff, sigma = sigma)
  df <- data.frame(concentrations = concs, cqs = cqs)

  model <- esc_mle(df)
  expect_s3_class(model, "esc")
  expect_true(abs(model$intercept - cq1)/cq1 < 0.1)
  expect_true(abs(model$sigma - sigma)/sigma < 0.5)
  expect_true(abs(model$eff - eff)/eff < 0.1)
  expect_equal(model$slope, -1/log10(1 + model$eff))
})

test_that("esc MLE works when approximation is turned off", {
  concs <- 2^rep(seq(from = -2, to = 10), each = 3)
  cq1 <- 38
  sigma <- 0.5
  eff <- 0.97
  cqs <- sim_cqs(concs, cq1 = cq1, eff = eff, sigma = sigma)
  df <- data.frame(concentrations = concs, cqs = cqs)
  model <- esc_mle(df, approximate = FALSE)
  expect_s3_class(model, "esc")
  expect_true(abs(model$intercept - cq1)/cq1 < 0.1)
  expect_true(abs(model$sigma - sigma)/sigma < 0.5)
  expect_true(abs(model$eff - eff)/eff < 0.25)
})

test_that("esc_mle error checks work", {
  df <- data.frame(concentrations = c(1, 2), cps = c(1, 2))
  expect_error(esc_mle(df), "must contain cqs column")
  df <- data.frame(concentraions = c(1, 2), cqs = c(1, 2))
  expect_error(esc_mle(df), "must contain concentrations column")
  df <- data.frame(concentrations = c(1, "a"), cqs = c(1, 2))
  expect_error(esc_mle(df), "concentrations must be numeric")
  df <- data.frame(concentrations = c(1, 2), cqs = c(1, -1))
  expect_error(esc_mle(df), "cqs must be non-negative")
  df <- data.frame(concentrations = c(1, -1), cqs = c(1, 2))
  expect_error(esc_mle(df), "concentrations must be non-negative")
  df <- data.frame(concentrations = c(1, 2), cqs = c(1, 2))
  expect_error(esc_mle(df, approximate = "a"), "`approximate` must be logical.")
  expect_error(esc_mle(df, assumeND = "a"), "`assumeND` must be logical.")
})

test_that("esc mcmc estimation works", {
  skip_if_not_installed("runjags")
  concs <- 2^rep(seq(from = -2, to = 10), each = 3)
  cq1 <- 38
  sigma <- 0.5
  eff <- 0.97
  cqs <- sim_cqs(concs, cq1 = cq1, eff = eff, sigma = sigma)
  df <- data.frame(concentrations = concs, cqs = cqs)
  res <- esc_mcmc(esc_data = df)
  expect_s3_class(res$mcmc_samples, "runjags")
  expect_true(abs(res$intercept$median - cq1)/cq1 < 0.1)
  expect_true(abs(res$intercept$mean - cq1)/cq1 < 0.1)
  expect_true(abs(res$sigma$median - sigma)/sigma < 0.5)
  expect_true(abs(res$sigma$mean - sigma)/sigma < 0.5)
  expect_true(abs(res$eff$median - eff)/eff < 0.1)
  expect_true(abs(res$eff$mean - eff)/eff < 0.1)
  expect_equal(res$slope$median, -1/log10(1 + res$eff$median), tolerance = 1e-5)
  expect_true(res$intercept$lower < res$intercept$median)
  expect_true(res$intercept$upper > res$intercept$median)
  expect_true(res$sigma$lower < res$sigma$median)
  expect_true(res$sigma$upper > res$sigma$median)
  expect_true(res$eff$lower < res$eff$median)
  expect_true(res$eff$upper > res$eff$median)
  expect_true(res$slope$lower < res$slope$median)
  expect_true(res$slope$upper > res$slope$median)
})

test_that("esc_mcmc error checks work", {
  skip_if_not_installed("runjags")
  df <- data.frame(concentrations = c(1, 2), cps = c(1, 2))
  expect_error(esc_mcmc(df), "must contain cqs column")
  df <- data.frame(concentraions = c(1, 2), cqs = c(1, 2))
  expect_error(esc_mcmc(df), "must contain concentrations column")
  df <- data.frame(concentrations = c(1, 2), cqs = c(1, "a"))
  expect_error(esc_mcmc(df), "cqs must be numeric")
  df <- data.frame(concentrations = c(1, "a"), cqs = c(1, 2))
  expect_error(esc_mcmc(df), "concentrations must be numeric")
  df <- data.frame(concentrations = c(1, 2), cqs = c(1, -1))
  expect_error(esc_mcmc(df), "cqs must be non-negative")
  df <- data.frame(concentrations = c(1, -1), cqs = c(1, 2))
  expect_error(esc_mcmc(df), "concentrations must be non-negative")
  df <- data.frame(concentrations = c(1, 2), cqs = c(1, 2))
  expect_error(esc_mcmc(df, level = "a"), "level must be numeric")
  expect_error(esc_mcmc(df, level = 1.2), "level must be between 0 and 1")
})
