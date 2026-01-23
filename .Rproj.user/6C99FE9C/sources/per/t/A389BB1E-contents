#' @name mvtdiff_distributions
#' @title Distribution Functions for Multivariate Approximated t-Difference
#'
#' @param x Matrix of quantiles (n x p) or vector for single point
#' @param q Vector of quantiles (length p) for cumulative probability
#' @param n Number of observations
#' @param mm_result Result from mm_tdiff_multivariate_general()
#' @param log Logical; if TRUE, returns log density
#' @param lower.tail Logical; if TRUE (default), probabilities are P(X <= x)
#'
#' @return
#' For \code{dmvtdiff}: Numeric vector of density values.
#' For \code{pmvtdiff}: Numeric scalar of cumulative probability.
#' For \code{rmvtdiff}: Matrix of random samples (n x p).
#'
#' @details
#' These functions implement the distribution functions for the approximated
#' multivariate t-difference based on Theorem 3 from the paper.
#'
#' **Note on degrees of freedom:**
#' \itemize{
#'   \item \code{dmvtdiff} uses the exact (non-integer) nu_star from the paper
#'   \item \code{pmvtdiff} rounds nu_star to the nearest integer due to
#'         mvtnorm::pmvt requirements. This introduces minimal approximation
#'         error when nu_star > 10 (the recommended range).
#'   \item \code{rmvtdiff} uses the exact (non-integer) nu_star
#' }
#'
#' @examples
#' # Setup
#' Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
#' Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
#' result <- mm_tdiff_multivariate_general(
#'   mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
#'   mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
#' )
#'
#' # Density at a point
#' dmvtdiff(c(0, 1), result)
#'
#' # Density at multiple points
#' x_mat <- matrix(c(0, 1, -1, 0.5), nrow = 2, byrow = TRUE)
#' dmvtdiff(x_mat, result)
#'
#' # Cumulative probability
#' pmvtdiff(c(0, 1), result)
#'
#' # Random samples
#' samples <- rmvtdiff(100, result)
#' head(samples)
#'
#' @importFrom mvtnorm dmvt pmvt rmvt
NULL

#' @rdname mvtdiff_distributions
#' @export
dmvtdiff <- function(x, mm_result, log = FALSE) {
  if (!inherits(mm_result, "mm_tdiff_multivariate_general")) {
    stop("mm_result must be output from mm_tdiff_multivariate_general()")
  }

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  p <- length(mm_result$mu_diff)
  if (ncol(x) != p) {
    stop("x must have ", p, " columns to match dimension of mm_result")
  }

  # dmvt accepts non-integer df - use exact nu_star
  mvtnorm::dmvt(x,
                delta = mm_result$mu_diff,
                sigma = mm_result$Sigma_star,
                df = mm_result$nu_star,
                log = log)
}

#' @rdname mvtdiff_distributions
#' @export
pmvtdiff <- function(q, mm_result, lower.tail = TRUE) {
  if (!inherits(mm_result, "mm_tdiff_multivariate_general")) {
    stop("mm_result must be output from mm_tdiff_multivariate_general()")
  }

  if (is.matrix(q)) {
    if (nrow(q) != 1) {
      stop("q must be a vector or single-row matrix for pmvtdiff")
    }
    q <- as.vector(q)
  }

  p <- length(mm_result$mu_diff)
  if (length(q) != p) {
    stop("q must have length ", p, " to match dimension of mm_result")
  }

  if (lower.tail) {
    lower <- rep(-Inf, p)
    upper <- q - mm_result$mu_diff
  } else {
    lower <- q - mm_result$mu_diff
    upper <- rep(Inf, p)
  }

  # pmvt requires integer df - round nu_star
  # This introduces minimal error when nu_star > 10
  result <- mvtnorm::pmvt(lower = lower,
                          upper = upper,
                          delta = rep(0, p),
                          sigma = mm_result$Sigma_star,
                          df = round(mm_result$nu_star))

  as.numeric(result)
}

#' @rdname mvtdiff_distributions
#' @export
rmvtdiff <- function(n, mm_result) {
  if (!inherits(mm_result, "mm_tdiff_multivariate_general")) {
    stop("mm_result must be output from mm_tdiff_multivariate_general()")
  }

  # rmvt accepts non-integer df - use exact nu_star
  samples <- mvtnorm::rmvt(n = n,
                           sigma = mm_result$Sigma_star,
                           df = mm_result$nu_star,
                           delta = rep(0, length(mm_result$mu_diff)))

  sweep(samples, 2, mm_result$mu_diff, "+")
}
