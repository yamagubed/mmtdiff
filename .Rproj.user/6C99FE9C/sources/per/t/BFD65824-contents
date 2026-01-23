#' Moment-Matching Approximation for Univariate t-Differences
#'
#' @description
#' Approximates the distribution of the difference between two independent
#' non-standardized t-distributed random variables using the moment-matching
#' method.
#'
#' @details
#' For two independent non-standardized t-distributed random variables:
#' \itemize{
#'   \item X1 ~ t(mu1, sigma1^2, nu1)
#'   \item X2 ~ t(mu2, sigma2^2, nu2)
#' }
#'
#' The difference Z = X1 - X2 is approximated as:
#' Z ~ t(mu1 - mu2, sigma_star^2, nu_star)
#'
#' where the effective parameters are computed through moment matching:
#' \itemize{
#'   \item sigma_star is derived from the second moment matching
#'   \item nu_star is derived from the fourth moment matching
#' }
#'
#' The method requires nu1 > 4 and nu2 > 4 for the existence of fourth moments.
#' The approximation quality improves as degrees of freedom increase and
#' approaches exactness as nu -> infinity (normal limit).
#'
#' @param mu1 Location parameter of first distribution
#' @param sigma1 Scale parameter of first distribution (must be > 0)
#' @param nu1 Degrees of freedom of first distribution (must be > 4)
#' @param mu2 Location parameter of second distribution
#' @param sigma2 Scale parameter of second distribution (must be > 0)
#' @param nu2 Degrees of freedom of second distribution (must be > 4)
#'
#' @return An S3 object of class "mm_tdiff_univariate" containing:
#'   \item{mu_diff}{Location parameter of difference (mu1 - mu2)}
#'   \item{sigma_star}{Scale parameter}
#'   \item{nu_star}{Degrees of freedom}
#'   \item{input_params}{List of input parameters for reference}
#'   \item{method}{Character string "univariate"}
#'
#' @references
#' Yamaguchi, Y., Homma, G., Maruo, K., & Takeda, K.
#' Moment-Matching Approximation for Difference of Non-Standardized
#' t-Distributed Variables. (unpublished).
#'
#' @examples
#' # Example 1: Different scale parameters
#' result <- mm_tdiff_univariate(
#'   mu1 = 0, sigma1 = 1, nu1 = 10,
#'   mu2 = 0, sigma2 = 1.5, nu2 = 15
#' )
#' print(result)
#'
#' # Example 2: Equal parameters (special case)
#' result_equal <- mm_tdiff_univariate(
#'   mu1 = 5, sigma1 = 2, nu1 = 20,
#'   mu2 = 3, sigma2 = 2, nu2 = 20
#' )
#' print(result_equal)
#'
#' @seealso
#' \code{\link{dtdiff}}, \code{\link{ptdiff}}, \code{\link{qtdiff}}, \code{\link{rtdiff}}
#' for density, distribution function, quantile function, and random generation
#' respectively
#'
#' @export
#' @importFrom stats dt pt qt rt
mm_tdiff_univariate <- function(mu1, sigma1, nu1, mu2, sigma2, nu2) {
  # Input validation
  if (nu1 <= 4 || nu2 <= 4) {
    stop("Both nu1 and nu2 must be greater than 4 for fourth moment to exist")
  }
  if (sigma1 <= 0 || sigma2 <= 0) {
    stop("Scale parameters must be positive")
  }
  if (!is.numeric(c(mu1, sigma1, nu1, mu2, sigma2, nu2))) {
    stop("All parameters must be numeric")
  }

  # Compute variance components (Var(X) = sigma^2 * nu/(nu-2) for t-distribution)
  var1 <- sigma1^2 * nu1 / (nu1 - 2)
  var2 <- sigma2^2 * nu2 / (nu2 - 2)

  # Degrees of freedom
  Qu_star <- (var1 + var2)^2
  term1 <- (sigma1^2)^2 * nu1^2 / ((nu1 - 2) * (nu1 - 4))
  term2 <- (sigma2^2)^2 * nu2^2 / ((nu2 - 2) * (nu2 - 4))
  term3 <- 2 * var1 * var2
  Qu <- term1 + term2 + term3
  nu_star <- (2 * Qu_star - 4 * Qu) / (Qu_star - Qu)

  # Scale parameter
  sigma_star <- sqrt((var1 + var2) * (nu_star - 2) / nu_star)

  # Prepare result object
  result <- list(
    mu_diff = mu1 - mu2,
    sigma_star = sigma_star,
    nu_star = nu_star,
    input_params = list(
      mu1 = mu1, sigma1 = sigma1, nu1 = nu1,
      mu2 = mu2, sigma2 = sigma2, nu2 = nu2
    ),
    method = "univariate"
  )
  class(result) <- c("mm_tdiff_univariate", "list")
  return(result)
}
