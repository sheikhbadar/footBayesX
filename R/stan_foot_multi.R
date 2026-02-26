#' Multi-Covariate Extension of footBayes Double Poisson
#'
#' Extends the footBayes double Poisson model to allow multiple covariates
#' through a linear predictor term X * beta.
#'
#' @param data Data frame with columns:
#' periods, home_team, away_team, home_goals, away_goals.
#' @param model Optional model label (kept for compatibility).
#' @param model_file Path to the Stan model file.
#' @param X Numeric matrix of covariates (N x K).
#' @param X_prev Optional matrix for out-of-sample prediction.
#' @param predict Number of matches to predict (last rows).
#' @param ranking Optional ranking matrix (currently not supported).
#' @param dynamic_type Placeholder for compatibility.
#' @param dynamic_weight Placeholder for compatibility.
#' @param dynamic_par Placeholder list for dynamic priors.
#' @param prior_par Placeholder list for prior settings.
#' @param home_effect Logical, include home advantage.
#' @param norm_method Placeholder normalization method.
#' @param ranking_map Placeholder ranking mapping.
#' @param method Estimation method (default "MCMC").
#' @param prior_dist_num Prior type for team effects.
#' @param prior_dist_sd_num Prior type for standard deviations.
#' @param hyper_df Degrees of freedom for Student-t priors.
#' @param hyper_location Location parameter for priors.
#' @param hyper_sd_df Degrees of freedom for sd priors.
#' @param hyper_sd_location Location for sd priors.
#' @param hyper_sd_scale Scale for sd priors.
#' @param mean_home Prior mean for home effect.
#' @param sd_home Prior sd for home effect.
#' @param ... Additional arguments passed to cmdstanr::sample().
#'
#' @return A list containing:
#' \itemize{
#'   \item fit: CmdStanMCMC object
#'   \item data: original data
#'   \item stan_data: list passed to Stan
#'   \item stan_code: Stan code used
#' }
#'
#' @export



stan_foot_multi <- function (
    data,
    model,
    model_file = NULL,
    X,
    X_prev = NULL,
    predict = 0,
    ranking,
    dynamic_type,
    dynamic_weight = FALSE,
    dynamic_par = list(
      common_sd = FALSE,
      spike = footBayes::normal(200, 0.1),
      slab  = footBayes::normal(0, 10),
      spike_prob = 0.2
    ),
    prior_par = list(
      ability    = footBayes::normal(0, NULL),
      ability_sd = footBayes::cauchy(0, 5),
      home       = footBayes::normal(0, 5)
    ),
    home_effect = TRUE,
    norm_method = "none",
    ranking_map = NULL,
    method = "MCMC",

    # -----------------------
    # NEW: Prior controls (these go into Stan data)
    # -----------------------
    prior_dist_num    = 2,   # 1=Normal, 2=Student-t, 3=Cauchy, 4=Laplace
    prior_dist_sd_num = 2,   # 1=Normal, 2=Student-t, 3=Cauchy, 4=Laplace

    hyper_df          = 4,
    hyper_location    = 0,

    hyper_sd_df       = 3,
    hyper_sd_location = 0,
    hyper_sd_scale    = 0.5,

    mean_home         = 0,
    sd_home           = 0.5,
    ...
){

  # -----------------------
  # Packages (safer than require())
  # -----------------------
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Install tidyr")
  if (!requireNamespace("cmdstanr", quietly = TRUE)) stop("Install cmdstanr")

  # -----------------------
  # BASIC CHECKS
  # -----------------------
  if (!is.data.frame(data))
    stop("data must be a data.frame")

  required_cols <- c("periods", "home_team", "away_team", "home_goals", "away_goals")
  if (!all(required_cols %in% names(data)))
    stop("Missing required columns: ", paste(setdiff(required_cols, names(data)), collapse = ", "))

  if (is.null(X))
    stop("X must be supplied")

  if (!is.matrix(X))
    X <- as.matrix(X)

  if (predict > 0 && is.null(X_prev))
    stop("X_prev required when predict > 0")

  if (!is.null(X_prev) && !is.matrix(X_prev))
    X_prev <- as.matrix(X_prev)

  # -----------------------
  # SPLIT TRAIN / TEST
  # -----------------------
  if (predict == 0){
    N <- nrow(data)
    N_prev <- 0
  } else {
    if (predict >= nrow(data)) stop("predict must be < nrow(data)")
    N <- nrow(data) - predict
    N_prev <- predict
  }

  # -----------------------
  # TEAMS (include home + away teams)
  # -----------------------
  teams <- sort(unique(c(data$home_team, data$away_team)))
  nteams <- length(teams)

  team_home <- match(data$home_team, teams)
  team_away <- match(data$away_team, teams)

  if (anyNA(team_home) || anyNA(team_away))
    stop("Some teams could not be matched to team index (NA found).")

  team1 <- team_home[1:N]
  team2 <- team_away[1:N]

  if (N_prev > 0){
    team1_prev <- team_home[(N+1):(N+N_prev)]
    team2_prev <- team_away[(N+1):(N+N_prev)]
  } else {
    team1_prev <- integer(0)
    team2_prev <- integer(0)
  }

  # -----------------------
  # GOALS
  # -----------------------
  y <- cbind(
    data$home_goals[1:N],
    data$away_goals[1:N]
  )

  # -----------------------
  # RANKING
  # -----------------------
  if (missing(ranking)){
    ntimes_rank <- 1
    ranking_matrix <- matrix(0, 1, nteams)
    instants_rank <- rep(1, N)
  } else {
    stop("Ranking not supported in this custom wrapper yet")
  }

  # -----------------------
  # HOME EFFECT
  # -----------------------
  ind_home <- as.integer(home_effect)

  # -----------------------
  # COVARIATES
  # -----------------------
  K <- ncol(X)

  X_train <- X[1:N, , drop = FALSE]

  if (N_prev > 0){
    if (is.null(X_prev)) stop("X_prev required when N_prev > 0")
    if (ncol(X_prev) != K) stop("X_prev must have same number of columns as X (K)")
  } else {
    X_prev <- matrix(0, 0, K)
  }

  # (Optional) ensure colnames propagate to X_prev for easier debugging
  if (!is.null(colnames(X_train)) && is.null(colnames(X_prev)) && nrow(X_prev) > 0){
    colnames(X_prev) <- colnames(X_train)
  }

  # -----------------------
  # STAN DATA
  # -----------------------
  data_stan <- list(
    N = N,
    N_prev = N_prev,
    y = y,

    nteams = nteams,
    team1 = team1,
    team2 = team2,
    team1_prev = team1_prev,
    team2_prev = team2_prev,

    ntimes_rank = ntimes_rank,
    instants_rank = instants_rank,
    ranking = ranking_matrix,

    ind_home = ind_home,
    mean_home = mean_home,
    sd_home = sd_home,

    prior_dist_num = prior_dist_num,
    prior_dist_sd_num = prior_dist_sd_num,
    hyper_df = hyper_df,
    hyper_location = hyper_location,
    hyper_sd_df = hyper_sd_df,
    hyper_sd_location = hyper_sd_location,
    hyper_sd_scale = hyper_sd_scale,

    K = K,
    X = X_train,
    X_prev = X_prev
  )

  # -----------------------
  # LOAD STAN MODEL
  # -----------------------
  stan_path <- system.file(
    "stan",
    "double_pois_multi.stan",
    package = "footBayesX"
  )

  if (is.null(model_file)) {
    model_file <- system.file("stan/double_pois_multi.stan",
                              package = "footBayesX")
  }

  model_stan <- cmdstanr::cmdstan_model(model_file)

  # -----------------------
  # FIT
  # -----------------------
  user_dots <- list(...)
  fit <- do.call(
    model_stan$sample,
    c(list(data = data_stan), user_dots)
  )

  out <- list(
    fit = fit,
    data = data,
    stan_data = data_stan,
    stan_code = fit$code()
  )

  class(out) <- "stanFoot"
  return(out)
}
