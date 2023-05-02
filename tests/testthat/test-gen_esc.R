
test_that("gen_esc() works", {

  intercept = 39
  eff       = 0.93
  sigma     = 0.30

  x = gen_esc(intercept, eff, sigma)

  expect_equal(x$intercept, intercept)
  expect_equal(x$eff, eff)
  expect_equal(x$sigma, sigma)
  expect_equal(x$intercept, intercept)


  expect_error(
    gen_esc(intercept = -1, eff, sigma),
    "`intercept` must be a positive number."
  )

  expect_error(
    gen_esc(intercept, eff = -1, sigma),
    "Efficiency `eff` must be positive."
  )

  expect_error(
    gen_esc(intercept, eff, sigma = -1),
    "`sigma` must be a positive number."
  )

  expect_error(
    gen_esc(intercept = "A", eff, sigma),
    "`intercept` must be a number."
  )

  expect_error(
    gen_esc(intercept, eff, sigma = "A"),
    "`sigma` must be a number."
  )

})
