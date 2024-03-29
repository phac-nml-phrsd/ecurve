test_that("esc MLE works", {

  concs <- 2^rep(seq(from = -2, to = 10), each = 2)
  cq1   <- 39
  sigma <- 0.10
  eff   <- 0.97

  cqs <- sim_cqs(concs, cq1 = cq1, eff = eff, sigma = sigma)
  esc_data <- data.frame(concentrations = concs, cqs = cqs)

  model <- esc_mle(esc_data)
  # plot_esc_model(model)

  expect_s3_class(model, "esc")
  expect_equal(model$intercept, cq1, tolerance = 0.1)
  expect_equal(model$sigma, sigma, tolerance = 0.6)
  expect_equal(model$eff,eff, tolerance =  0.1)
  expect_equal(model$slope, -1/log10(1 + model$eff))
})

test_that("esc MLE works when approximation is turned off", {
  concs <- 2^rep(seq(from = -2, to = 10), each = 3)
  cq1   <- 38
  sigma <- 0.02
  eff   <- 0.97
  cqs   <- sim_cqs(concs, cq1 = cq1, eff = eff, sigma = sigma)
  df    <- data.frame(concentrations = concs, cqs = cqs)

  model <- esc_mle(df, approximate = FALSE, nlm.print.level = 0)
  # plot_esc_model(model)

  expect_s3_class(model, "esc")
  expect_equal(model$intercept, cq1, tolerance = 0.1)
  expect_equal(model$sigma, sigma, tolerance = 0.6)
  expect_equal(model$eff,eff, tolerance =  0.1)
  expect_equal(model$slope, -1/log10(1 + model$eff))
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

test_that("esc MLE works on a specific example", {
  # This is the example of the vignette
  df <- data.frame(
    concentrations = c(200, 200, 200, 100, 100, 100, 50, 50, 50, 25, 25, 25, 12.5,
                       12.5, 12.5, 6.25, 6.25, 6.25, 3.125, 3.125, 3.125, 1.5625, 1.5625,
                       1.5625, 1.5625, 1.5625, 1.5625, 0.7813, 0.7813, 0.3906, 0.3906,
                       0.1953),
    cqs = c(30.99, 30.84, 30.73, 31.98, 31.68, 31.69, 33.27, 32.83, 33.01,
            33.51, 33.15, 33.83, 35.36, 34.64, 35.09, 36.09, 36.05, 36, 37.34,
            36.23, 38.36, 38.22, 36.5, 37.63, 37.61, 37.17, 37.21, 36.2,
            38.25, 39.15, 38.13, 39.17))

  model <- esc_mle(df)
  # plot_esc_model(model)
  # print_esc_model(model)

  model.benchmark = gen_esc(
    intercept = 38.50225,
    eff       = 0.9951575,
    sigma     = 0.3197636)

  expect_s3_class(model, "esc")

  for(p in c('intercept', 'eff', 'sigma')){
    # print(p)
    expect_equal(model[[p]],
                 expected = model.benchmark[[p]],
                 tolerance = ifelse(p=='sigma',0.5,0.05))
  }

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
