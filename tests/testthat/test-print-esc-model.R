test_that("print_esc_model works", {

  expect_warning(print_esc_model('foo'))

  esc_data = data.frame(
     concentrations = c(1,1,10, 10, 100, 500, 500),
     cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
  mod = esc_mle(esc_data)
  expect_output(print_esc_model(mod), regexp = 'ESC model')

})
