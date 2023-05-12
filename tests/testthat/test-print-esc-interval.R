test_that("print_esc_model works", {

  expect_warning(print_esc_interval('foo'))

  esc_data = data.frame(
    concentrations = c(1,1,10, 10, 100, 500, 500),
    cqs = c(40.2, 39.3, 35.9, 36.4, 32.6, 30.0, 31.1))
  model = esc_mle(esc_data)
  cqs = new.cqs = c( 34, NaN, 36)

  int = conc_interval(cqs = new.cqs, model = model)

  expect_output( print_esc_interval(int),
                 regexp = 'ESC interval')

})
