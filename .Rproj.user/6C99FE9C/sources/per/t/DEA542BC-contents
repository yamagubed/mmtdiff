#' Moment-Matching Approximation for Multivariate t-Differences (Independent)
#'
#' @description
#' Approximates the distribution of differences between two independent
#' p-dimensional vectors with independent t-distributed components.
#'
#' @details
#' This function applies the univariate moment-matching approximation
#' component-wise when all components are mutually independent. Each
#' component difference Zj = X1j - X2j is approximated independently using
#' the univariate method.
#'
#' This approach is optimal for:
#' \itemize{
#'   \item Marginal inference on specific components
#'   \item Cases where components have different tail behaviors
#'   \item Maintaining computational efficiency in high dimensions
#' }
#'
#' @param mu1 Location vector of first distribution (length p)
#' @param sigma1 Scale vector of first distribution (length p, all > 0)
#' @param nu1 Degrees of freedom vector of first distribution (length p, all > 4)
#' @param mu2 Location vector of second distribution (length p)
#' @param sigma2 Scale vector of second distribution (length p, all > 0)
#' @param nu2 Degrees of freedom vector of second distribution (length p, all > 4)
#'
#' @return An S3 object of class "mm_tdiff_multivariate_independent" containing:
#'   \item{mu_diff}{Location vector of difference}
#'   \item{sigma_star}{Vector of scale parameters}
#'   \item{nu_star}{Vector of degrees of freedom}
#'   \item{p}{Dimension of the vectors}
#'   \item{method}{Character string "multivariate_independent"}
#'
#' @examples
#' result <- mm_tdiff_multivariate_independent(
#'   mu1 = c(0, 1), sigma1 = c(1, 1.5), nu1 = c(10, 12),
#'   mu2 = c(0, 0), sigma2 = c(1.2, 1), nu2 = c(15, 20)
#' )
#' print(result)
#'
#' @seealso \code{\link{mm_tdiff_multivariate_general}} for correlated components
#'
#' @export
mm_tdiff_multivariate_independent <- function(mu1, sigma1, nu1, mu2, sigma2, nu2) {
  p <- length(mu1)

  # Input validation
  if (length(sigma1) != p || length(nu1) != p ||
      length(mu2) != p || length(sigma2) != p || length(nu2) != p) {
    stop("All input vectors must have the same length")
  }
  if (any(nu1 <= 4) || any(nu2 <= 4)) {
    stop("All degrees of freedom must be greater than 4")
  }
  if (any(sigma1 <= 0) || any(sigma2 <= 0)) {
    stop("All scale parameters must be positive")
  }

  # Apply univariate approximation to each component
  sigma_star <- numeric(p)
  nu_star <- numeric(p)

  for (j in 1:p) {
    result <- mm_tdiff_univariate(mu1[j], sigma1[j], nu1[j],
                                  mu2[j], sigma2[j], nu2[j])
    sigma_star[j] <- result$sigma_star
    nu_star[j] <- result$nu_star
  }

  result <- list(
    mu_diff = mu1 - mu2,
    sigma_star = sigma_star,
    nu_star = nu_star,
    p = p,
    method = "multivariate_independent"
  )
  class(result) <- c("mm_tdiff_multivariate_independent", "list")
  return(result)
}

#' Moment-Matching Approximation for General Multivariate t-Differences
#'
#' @description
#' Approximates the distribution of differences between two independent
#' multivariate t-distributed random vectors with arbitrary covariance
#' structure.
#'
#' @details
#' This function handles the general case where components may be correlated
#' within each multivariate t-distribution. The approximation uses a single
#' scalar degrees of freedom parameter to capture the overall tail behavior.
#'
#' @param mu1 Location vector of first distribution (length p)
#' @param Sigma1 Scale matrix of first distribution (p x p, positive definite)
#' @param nu1 Degrees of freedom of first distribution (must be > 4)
#' @param mu2 Location vector of second distribution (length p)
#' @param Sigma2 Scale matrix of second distribution (p x p, positive definite)
#' @param nu2 Degrees of freedom of second distribution (must be > 4)
#'
#' @return An S3 object of class "mm_tdiff_multivariate_general" containing:
#'   \item{mu_diff}{Location vector of difference}
#'   \item{Sigma_star}{Scale matrix}
#'   \item{nu_star}{Degrees of freedom (scalar)}
#'   \item{method}{Character string "multivariate_general"}
#'
#' @examples
#' Sigma1 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
#' Sigma2 <- matrix(c(1.5, 0.5, 0.5, 1.2), 2, 2)
#' result <- mm_tdiff_multivariate_general(
#'   mu1 = c(0, 1), Sigma1 = Sigma1, nu1 = 10,
#'   mu2 = c(0, 0), Sigma2 = Sigma2, nu2 = 15
#' )
#' print(result)
#'
#' @export
mm_tdiff_multivariate_general <- function(mu1, Sigma1, nu1, mu2, Sigma2, nu2) {
  # Input validation
  if (nu1 <= 4 || nu2 <= 4) {
    stop("Both nu1 and nu2 must be greater than 4")
  }
  if (!is.matrix(Sigma1) || !is.matrix(Sigma2)) {
    stop("Sigma1 and Sigma2 must be matrices")
  }
  if (nrow(Sigma1) != ncol(Sigma1) || nrow(Sigma2) != ncol(Sigma2)) {
    stop("Sigma matrices must be square")
  }

  p <- length(mu1)
  if (nrow(Sigma1) != p || nrow(Sigma2) != p || length(mu2) != p) {
    stop("Dimensions of mu and Sigma must be consistent")
  }
  if (!isSymmetric(Sigma1,tol=1e-10) || !isSymmetric(Sigma2,tol=1e-10)) {
    warning("Sigma matrices should be symmetric; using symmetrized version")
    Sigma1 <- (Sigma1 + t(Sigma1)) / 2
    Sigma2 <- (Sigma2 + t(Sigma2)) / 2
  }

  # Check positive definiteness
  ev1 <- eigen(Sigma1, only.values = TRUE)$values
  ev2 <- eigen(Sigma2, only.values = TRUE)$values
  if (any(ev1 <= 0) || any(ev2 <= 0)) {
    stop("Sigma matrices must be positive definite")
  }

  # Constants
  beta1  <- nu1 / (nu1 - 2)
  beta2  <- nu2 / (nu2 - 2)
  alpha1 <- nu1^2 / ((nu1 - 2) * (nu1 - 4))
  alpha2 <- nu2^2 / ((nu2 - 2) * (nu2 - 4))
  gamma  <- 2 * nu1 * nu2 / ((nu1 - 2) * (nu2 - 2))

  V <- Sigma1 * beta1 + Sigma2 * beta2

  # Use a stable inverse via Cholesky when possible
  A <- tryCatch({
    R <- chol(V)
    chol2inv(R)
  }, error = function(e) {
    solve(V)
  })

  # Helper traces
  tr_AS1 <- sum(diag(A %*% Sigma1))
  tr_AS2 <- sum(diag(A %*% Sigma2))

  tr_AS1_AS1 <- sum(diag(A %*% Sigma1 %*% A %*% Sigma1))
  tr_AS2_AS2 <- sum(diag(A %*% Sigma2 %*% A %*% Sigma2))
  tr_AS1_AS2 <- sum(diag(A %*% Sigma1 %*% A %*% Sigma2))

  # Assemble Qm
  Qm_num <- alpha1 * (tr_AS1^2 + 2 * tr_AS1_AS1) +
            alpha2 * (tr_AS2^2 + 2 * tr_AS2_AS2) +
            gamma  * (tr_AS1 * tr_AS2 + 2 * tr_AS1_AS2)

  Qm <- Qm_num / (p * (p + 2))

  # Degrees of freedom
  nu_star <- (2 - 4 * Qm) / (1 - Qm)

  # Scale parameter
  Sigma_star <- V * (nu_star - 2) / nu_star

  result <- list(
    mu_diff = mu1 - mu2,
    Sigma_star = Sigma_star,
    nu_star = nu_star,
    method = "multivariate_general"
  )
  class(result) <- c("mm_tdiff_multivariate_general", "list")
  return(result)
}
