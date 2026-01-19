#' @export
print.mm_tdiff_univariate <- function(x, ...) {
  cat("Moment-Matching Approximation (Univariate)\n")
  cat("==============================================\n")
  cat(sprintf("Location (mu1 - mu2): %.4f\n", x$mu_diff))
  cat(sprintf("Scale (sigma*): %.4f\n", x$sigma_star))
  cat(sprintf("df (nu*): %.4f\n", x$nu_star))
  invisible(x)
}

#' @export
print.mm_tdiff_multivariate_independent <- function(x, ...) {
  cat("Moment-Matching Approximation (Multivariate Independent)\n")
  cat("=============================================================\n")
  cat("Location difference:\n"); print(x$mu_diff)
  cat("\nScale:\n"); print(x$sigma_star)
  cat("\ndf:\n"); print(x$nu_star)
  invisible(x)
}

#' @export
print.mm_tdiff_multivariate_general <- function(x, ...) {
  cat("Moment-Matching Approximation (General Multivariate)\n")
  cat("=========================================================\n")
  cat("Location difference:\n"); print(x$mu_diff)
  cat("\nScale matrix:\n"); print(x$Sigma_star)
  cat(sprintf("\ndf: %.4f\n", x$nu_star))
  invisible(x)
}
