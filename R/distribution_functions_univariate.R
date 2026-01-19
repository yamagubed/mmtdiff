#' @name tdiff_distributions
#' @title Distribution Functions for Approximated t-Difference
#'
#' @param x,q Vector of quantiles
#' @param p Vector of probabilities
#' @param n Number of observations
#' @param mm_result Result from mm_tdiff_univariate()
#'
#' @return
#' For \code{dtdiff}: Numeric vector of density values.
#' For \code{ptdiff}: Numeric vector of cumulative probabilities.
#' For \code{qtdiff}: Numeric vector of quantiles.
#' For \code{rtdiff}: Numeric vector of random samples from the approximated
#' t-difference distribution.
#'
#' @examples
#' result <- mm_tdiff_univariate(0, 1, 10, 0, 1.5, 15)
#' dtdiff(0, result)
#' ptdiff(0, result)
#' qtdiff(c(0.025, 0.975), result)
#' samples <- rtdiff(100, result)
NULL

#' @rdname tdiff_distributions
#' @export
dtdiff <- function(x, mm_result) {
  if (!inherits(mm_result, "mm_tdiff_univariate")) {
    stop("mm_result must be output from mm_tdiff_univariate()")
  }
  standardized <- (x - mm_result$mu_diff) / mm_result$sigma_star
  dt(standardized, df = mm_result$nu_star) / mm_result$sigma_star
}

#' @rdname tdiff_distributions
#' @export
ptdiff <- function(q, mm_result) {
  if (!inherits(mm_result, "mm_tdiff_univariate")) {
    stop("mm_result must be output from mm_tdiff_univariate()")
  }
  standardized <- (q - mm_result$mu_diff) / mm_result$sigma_star
  pt(standardized, df = mm_result$nu_star)
}

#' @rdname tdiff_distributions
#' @export
qtdiff <- function(p, mm_result) {
  if (!inherits(mm_result, "mm_tdiff_univariate")) {
    stop("mm_result must be output from mm_tdiff_univariate()")
  }
  if (any(p < 0) || any(p > 1)) {
    stop("Probabilities must be between 0 and 1")
  }
  mm_result$mu_diff + mm_result$sigma_star * qt(p, df = mm_result$nu_star)
}

#' @rdname tdiff_distributions
#' @export
rtdiff <- function(n, mm_result) {
  if (!inherits(mm_result, "mm_tdiff_univariate")) {
    stop("mm_result must be output from mm_tdiff_univariate()")
  }
  mm_result$mu_diff + mm_result$sigma_star * rt(n, df = mm_result$nu_star)
}
