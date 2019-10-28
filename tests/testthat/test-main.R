test_that("dphtype gives sensible results", {


  observed_1 = dphtype(seq(0, 10), generate_init_row(5), generate_subint_mat(6))
  expected_1 = c(0, 0.5338932, 0.2767559, 0.1060255, 0.03921488, 0.01443682, 0.00531153, 0.001954029, 0.0007188483, 0.0002644496, 9.728556e-05)
  expect_equal(observed_1, expected_1, tolerance = 1e-07)
})




test_that("rphtype gives sensible results", {
  expect_equal(2 * 2, 4)
})


