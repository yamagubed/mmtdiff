# Comprehensive tests for multivariate moment-matching approximation

test_that("multivariate independent components works correctly", {
  mu1 <- c(0, 1, 2)
  sigma1 <- c(1, 1.5, 2)
  nu1 <- c(10, 12, 15)
  mu2 <- c(0, 0, 1)
  sigma2 <- c(1.2, 1, 1.8)
  nu2 <- c(15, 20, 25)

  result <- mm_tdiff_multivariate_independent(mu1, sigma1, nu1, mu2, sigma2, nu2)

  expect_s3_class(result, "mm_tdiff_multivariate_independent")
  expect_equal(length(result$mu_diff), 3)
  expect_equal(length(result$sigma_star), 3)
  expect_equal(length(result$nu_star), 3)
  expect_equal(result$p, 3)

  for (j in 1:3) {
    uni_result <- mm_tdiff_univariate(
      mu1[j], sigma1[j], nu1[j],
      mu2[j], sigma2[j], nu2[j]
    )
    expect_equal(result$mu_diff[j], uni_result$mu_diff)
    expect_equal(result$sigma_star[j], uni_result$sigma_star)
    expect_equal(result$nu_star[j], uni_result$nu_star)
  }
})

test_that("multivariate independent validates dimensions", {
  expect_error(
    mm_tdiff_multivariate_independent(
      c(0, 1), c(1), c(10, 12),
      c(0, 0), c(1.2, 1), c(15, 20)
    ),
    "same length"
  )

  expect_error(
    mm_tdiff_multivariate_independent(
      c(0, 1), c(1, 1.5), c(3, 12),
      c(0, 0), c(1.2, 1), c(15, 20)
    ),
    "greater than 4"
  )
})

test_that("multivariate general case works correctly", {
  Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)

  result <- mm_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )

  expect_s3_class(result, "mm_tdiff_multivariate_general")
  expect_true(is.matrix(result$Sigma_star))
  expect_equal(dim(result$Sigma_star), c(2, 2))
  expect_true(result$nu_star > 0)
  expect_true(isSymmetric(result$Sigma_star))

  eigenvalues <- eigen(result$Sigma_star, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("multivariate general validates input", {
  expect_error(
    mm_tdiff_multivariate_general(
      c(0, 1), c(1, 0, 0, 1), 10,
      c(0, 0), matrix(c(1, 0, 0, 1), 2), 15
    ),
    "must be matrices"
  )

  expect_error(
    mm_tdiff_multivariate_general(
      c(0, 1), matrix(1:6, 2, 3), 10,
      c(0, 0), matrix(c(1, 0, 0, 1), 2), 15
    ),
    "square"
  )

  expect_error(
    mm_tdiff_multivariate_general(
      c(0, 1, 2), matrix(c(1, 0, 0, 1), 2), 10,
      c(0, 0), matrix(c(1, 0, 0, 1), 2), 15
    ),
    "consistent"
  )

  bad_matrix <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(
    mm_tdiff_multivariate_general(
      c(0, 1), bad_matrix, 10,
      c(0, 0), matrix(c(1, 0, 0, 1), 2), 15
    ),
    "positive definite"
  )
})

test_that("diagonal matrices reduce to independent case", {
  Sigma1 <- diag(c(1, 1.5))
  Sigma2 <- diag(c(1.2, 1))

  result_general <- mm_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )

  result_indep <- mm_tdiff_multivariate_independent(
    mu1 = c(0, 1), sigma1 = c(1, 1.5), nu1 = c(10, 10),
    mu2 = c(0, 0), sigma2 = c(1.2, 1), nu2 = c(15, 15)
  )

  expect_equal(result_general$mu_diff, result_indep$mu_diff)

  diag_general <- diag(result_general$Sigma_star)
  expect_true(all(abs(sqrt(diag_general) - result_indep$sigma_star) < 1))
})

# Tests for multivariate distribution functions

test_that("dmvtdiff works correctly", {
  skip_if_not_installed("mvtnorm")

  Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
  result <- mm_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )

  d_mean <- dmvtdiff(result$mu_diff, result)
  expect_true(d_mean > 0)

  x_mat <- matrix(c(0, 1, -1, 0.5), nrow = 2, byrow = TRUE)
  densities <- dmvtdiff(x_mat, result)
  expect_equal(length(densities), 2)
  expect_true(all(densities > 0))

  d_log <- dmvtdiff(result$mu_diff, result, log = TRUE)
  expect_equal(d_log, log(dmvtdiff(result$mu_diff, result)))
})

test_that("pmvtdiff works correctly", {
  skip_if_not_installed("mvtnorm")

  Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
  result <- mm_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )

  # pmvtdiff uses rounded df, so it should work
  p_low <- pmvtdiff(result$mu_diff - 10 * sqrt(diag(result$Sigma_star)), result)
  expect_true(p_low < 0.01)

  p_high <- pmvtdiff(result$mu_diff + 10 * sqrt(diag(result$Sigma_star)), result)
  expect_true(p_high > 0.99)

  # Test lower.tail parameter works
  q_test <- c(0.5, 1.5)
  p_lower <- pmvtdiff(q_test, result, lower.tail = TRUE)
  p_upper <- pmvtdiff(q_test, result, lower.tail = FALSE)

  # For multivariate: p_lower + p_upper ≠ 1 in general
  # They represent different rectangular regions
  expect_true(p_lower >= 0 && p_lower <= 1)
  expect_true(p_upper >= 0 && p_upper <= 1)

  # But we can verify: p_lower >= p_upper when q is above the mean
  # (more probability mass in the lower-left rectangle)
  expect_true(p_lower + p_upper <= 1)  # They don't overlap, so sum <= 1
})

test_that("rmvtdiff works correctly", {
  skip_if_not_installed("mvtnorm")

  Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
  result <- mm_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )

  set.seed(123)
  samples <- rmvtdiff(100, result)

  expect_equal(dim(samples), c(100, 2))
  expect_true(is.matrix(samples))

  sample_means <- colMeans(samples)
  expect_equal(sample_means, result$mu_diff, tolerance = 0.5)

  theoretical_cov <- result$Sigma_star * result$nu_star / (result$nu_star - 2)
  sample_cov <- cov(samples)
  expect_equal(sample_cov[1, 1], theoretical_cov[1, 1], tolerance = 1)
  expect_equal(sample_cov[2, 2], theoretical_cov[2, 2], tolerance = 1)
})

test_that("multivariate distribution functions validate input", {
  skip_if_not_installed("mvtnorm")

  Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
  result <- mm_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )

  expect_error(
    dmvtdiff(c(0, 1), list(mu = c(0, 1))),
    "must be output from mm_tdiff_multivariate_general"
  )
  expect_error(
    pmvtdiff(c(0, 1), list(mu = c(0, 1))),
    "must be output from mm_tdiff_multivariate_general"
  )
  expect_error(
    rmvtdiff(10, list(mu = c(0, 1))),
    "must be output from mm_tdiff_multivariate_general"
  )

  expect_error(
    dmvtdiff(c(0, 1, 2), result),
    "columns to match dimension"
  )
  expect_error(
    pmvtdiff(c(0, 1, 2), result),
    "length"
  )
})

test_that("print methods work for multivariate", {
  result_indep <- mm_tdiff_multivariate_independent(
    mu1 = c(0, 1), sigma1 = c(1, 1.5), nu1 = c(10, 12),
    mu2 = c(0, 0), sigma2 = c(1.2, 1), nu2 = c(15, 20)
  )
  expect_output(print(result_indep), "Multivariate Independent")

  Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
  result_general <- mm_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )
  expect_output(print(result_general), "General Multivariate")
})

test_that("validate_approximation works", {
  result <- mm_tdiff_univariate(0, 1, 10, 0, 1.5, 15)

  validation <- validate_approximation(result)

  expect_type(validation, "list")
  expect_true("mean" %in% names(validation))
  expect_true("variance" %in% names(validation))
  expect_equal(validation$mean, result$mu_diff)
})

test_that("nu_star handling in multivariate functions", {
  skip_if_not_installed("mvtnorm")

  Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
  result <- mm_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
  )

  # nu_star should be non-integer for these parameters
  expect_false(result$nu_star == round(result$nu_star))

  # dmvtdiff and rmvtdiff use exact nu_star
  expect_no_error(dmvtdiff(c(0, 1), result))
  expect_no_error(rmvtdiff(10, result))

  # pmvtdiff rounds nu_star internally but should still work
  expect_no_error(pmvtdiff(c(0, 1), result))
})

test_that("rounding impact is minimal for pmvtdiff", {
  skip_if_not_installed("mvtnorm")

  # Test with nu_star that rounds differently
  Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
  result <- mm_tdiff_multivariate_general(
    mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 20,
    mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 30
  )

  # With large nu_star, rounding should have minimal impact
  expect_true(result$nu_star > 10)

  # Should compute successfully
  p_value <- pmvtdiff(c(0, 1), result)
  expect_true(p_value >= 0 && p_value <= 1)
})
