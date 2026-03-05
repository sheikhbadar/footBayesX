#' Print method for stanFootMulti objects
#'
#' Thiss will provide a compact summary of the fitted footBayesX model.
#'
#' @param x A stanFootMulti object returned by stan_foot_multi().
#' @param ... Unused.
#'
#' @export
print.stanFootMulti <- function(x, ...) {
  
  cat("\nfootBayesX model fit\n")
  cat("--------------------\n")
  
  N_train <- x$stan_data$N
  N_pred  <- x$stan_data$N_prev
  K       <- x$stan_data$K
  teams   <- length(x$teams)
  
  cat("Teams            :", teams, "\n")
  cat("Training matches :", N_train, "\n")
  cat("Prediction set   :", N_pred, "\n")
  cat("Covariates (K)   :", K, "\n")
  
  if (!is.null(x$X_names)) {
    cat("Covariate names  :", paste(x$X_names, collapse = ", "), "\n")
  }
  
  cat("\nStan sampling:\n")
  
  draws <- x$fit$metadata()
  
  cat("Chains           :", draws$chains, "\n")
  cat("Iterations       :", draws$iter_sampling, "sampling /",
      draws$iter_warmup, "warmup\n")
  
  cat("\nUse:\n")
  cat("  x$fit$summary()\n")
  cat("  posterior::summarise_draws(x$fit$draws())\n")
  cat("  foot_abilities(x, data)\n")
  cat("  foot_rank(x, data)\n")
  
  invisible(x)
}