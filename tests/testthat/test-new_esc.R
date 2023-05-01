test_that("new() works", {

  intercept = 39
  slope     = -1.5
  sigma     = 0.30

  x = new_esc(intercept = intercept,
              slope = slope,
              sigma = sigma)

  expect_equal(x$intercept, intercept)
  expect_equal(x$slope, slope)
  expect_equal(x$sigma, sigma)


  expect_error(
    new_esc(intercept = -1, slope, sigma),
    "`intercept` must be a positive number."
  )

  expect_error(
    new_esc(intercept, slope, sigma = -1),
    "`sigma` must be a positive number."
  )

})
