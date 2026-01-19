# Comprehensive tests for univariate moment-matching approximation

test_that("mm_tdiff_univariate mathematical correctness", {
  # Test case from paper example
  result <- mm_tdiff_univariate(
    mu1 = 0, sigma1 = 1, nu1 = 10,
    mu2 = 0, sigma2 = 1.5, nu2 = 15
  )

  # Verify structure
  expect_type(result, "list")
  expect_s3_class(result, "mm_tdiff_univariate")
  expect_equal(length(result), 5)
  expect_named(result, c("mu_diff", "sigma_star", "nu_star", "input_params", "method"))

  # Verify computations
  expect_equal(result$mu_diff, 0)

  # Verify nu_star is positive and reasonable
  expect_true(result$nu_star > 0)
  expect_true(result$nu_star > 4)
})

test_that("mm_tdiff_univariate validates input correctly", {
  # Test nu <= 4 rejection
  expect_error(
    mm_tdiff_univariate(0, 1, 3, 0, 1, 10),
    "must be greater than 4"
  )
  expect_error(
    mm_tdiff_univariate(0, 1, 10, 0, 1, 4),
    "must be greater than 4"
  )

  # Test negative sigma rejection
  expect_error(
    mm_tdiff_univariate(0, -1, 10, 0, 1, 10),
    "positive"
  )
  expect_error(
    mm_tdiff_univariate(0, 1, 10, 0, 0, 10),
    "positive"
  )

  # Test non-numeric input
  expect_error(
    mm_tdiff_univariate("a", 1, 10, 0, 1, 10),
    "numeric"
  )
})

test_that("asymptotic behavior is correct", {
  # As nu -> infinity, should approach normal distribution
  result_large_nu <- mm_tdiff_univariate(
    mu1 = 0, sigma1 = 1, nu1 = 10000,
    mu2 = 0, sigma2 = 1.5, nu2 = 10000
  )

  # Variance should approach sigma1^2 + sigma2^2
  expected_var <- 1^2 + 1.5^2  # 3.25
  actual_var <- result_large_nu$sigma_star^2 *
    result_large_nu$nu_star / (result_large_nu$nu_star - 2)
  expect_equal(actual_var, expected_var, tolerance = 0.01)

  # nu_star should be very large
  expect_true(result_large_nu$nu_star > 1000)
})

test_that("equal parameters special case works", {
  # When parameters are equal, use the simplified formula
  result <- mm_tdiff_univariate(
    mu1 = 5, sigma1 = 2, nu1 = 20,
    mu2 = 3, sigma2 = 2, nu2 = 20
  )

  # For equal parameters: nu_star = 2*nu - 4
  expect_equal(result$nu_star, 2 * 20 - 4)
  expect_equal(result$mu_diff, 5 - 3)
})

# Tests for distribution functions

test_that("dtdiff works correctly", {
  result <- mm_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  # Density at mean should be positive
  d_mean <- dtdiff(result$mu_diff, result)
  expect_true(d_mean > 0)

  # Density should integrate to approximately 1
  x_vals <- seq(result$mu_diff - 5 * result$sigma_star,
                result$mu_diff + 5 * result$sigma_star,
                length.out = 1000)
  densities <- dtdiff(x_vals, result)
  integral <- sum(densities) * diff(x_vals[1:2])
  expect_equal(integral, 1, tolerance = 0.01)
})

test_that("ptdiff works correctly", {
  result <- mm_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  # CDF at -Inf should be 0, at +Inf should be 1
  expect_equal(ptdiff(-Inf, result), 0)
  expect_equal(ptdiff(Inf, result), 1)

  # CDF at mean should be approximately 0.5
  expect_equal(ptdiff(result$mu_diff, result), 0.5, tolerance = 0.01)

  # CDF should be monotonically increasing
  x_vals <- seq(-5, 5, by = 0.5)
  p_vals <- ptdiff(x_vals, result)
  expect_true(all(diff(p_vals) >= 0))
})

test_that("qtdiff works correctly", {
  result <- mm_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  # Median should be close to mean
  q_median <- qtdiff(0.5, result)
  expect_equal(q_median, result$mu_diff, tolerance = 0.01)

  # 2.5% and 97.5% quantiles should give 95% interval
  q_lower <- qtdiff(0.025, result)
  q_upper <- qtdiff(0.975, result)
  expect_true(q_lower < result$mu_diff)
  expect_true(q_upper > result$mu_diff)

  # Inverse relationship with ptdiff
  x <- 1.5
  p <- ptdiff(x, result)
  q <- qtdiff(p, result)
  expect_equal(q, x, tolerance = 1e-6)

  # Test error for invalid probabilities
  expect_error(qtdiff(-0.1, result), "between 0 and 1")
  expect_error(qtdiff(1.1, result), "between 0 and 1")
})

test_that("rtdiff works correctly", {
  result <- mm_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  # Generate samples
  set.seed(123)
  samples <- rtdiff(1000, result)

  # Check basic properties
  expect_equal(length(samples), 1000)
  expect_true(is.numeric(samples))

  # Sample mean should be close to mu_diff
  expect_equal(mean(samples), result$mu_diff, tolerance = 0.2)

  # Sample variance should be close to theoretical variance
  theoretical_var <- result$sigma_star^2 * result$nu_star / (result$nu_star - 2)
  expect_equal(var(samples), theoretical_var, tolerance = 0.5)
})

test_that("distribution functions validate input", {
  result <- mm_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  # Test with wrong class
  expect_error(dtdiff(0, list(mu = 0)), "must be output from mm_tdiff_univariate")
  expect_error(ptdiff(0, list(mu = 0)), "must be output from mm_tdiff_univariate")
  expect_error(qtdiff(0.5, list(mu = 0)), "must be output from mm_tdiff_univariate")
  expect_error(rtdiff(10, list(mu = 0)), "must be output from mm_tdiff_univariate")
})

test_that("print method works", {
  result <- mm_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  expect_output(print(result), "Moment-Matching Approximation")
  expect_output(print(result), "Location")
  expect_output(print(result), "Scale")
  expect_output(print(result), "df")
})

test_that("nu_star is non-integer and used exactly", {
  result <- mm_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  # nu_star should be non-integer for these parameters
  expect_false(result$nu_star == round(result$nu_star))

  # All distribution functions should accept non-integer df
  expect_no_error(dtdiff(0, result))
  expect_no_error(ptdiff(0, result))
  expect_no_error(qtdiff(0.5, result))
  expect_no_error(rtdiff(10, result))
})
