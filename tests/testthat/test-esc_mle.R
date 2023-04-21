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
