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
#' @param ind_home indicator for home advantage (1 = include home effect)
#' @param beta_prior_mean mean of Normal prior for regression coefficients beta
#' @param beta_prior_sd standard deviation of beta prior
#' @param chains number of MCMC chains
#' @param iter_sampling number of sampling iterations
#' @param iter_warmup number of warmup iterations
#' @param seed random seed
#' @param adapt_delta Stan NUTS adapt_delta
#' @param max_treedepth Stan max tree depth
#' @return A list containing:
#' \itemize{
#'   \item fit: CmdStanMCMC object
#'   \item data: original data
#'   \item stan_data: list passed to Stan
#'   \item stan_code: Stan code used
#' }
#'
#' @export


stan_foot_multi <- function(
    data,
    model_file,
    X,
    X_prev = NULL,
    predict = 0,

    # Home prior
    ind_home = 1,
    mean_home = 0,
    sd_home = 0.5,

    # Beta prior
    beta_prior_mean = 0,
    beta_prior_sd = 1,

    # Ability priors
    prior_dist_num = 1,
    prior_dist_sd_num = 1,
    hyper_df = 4,
    hyper_location = 0,
    hyper_sd_df = 3,
    hyper_sd_location = 0,
    hyper_sd_scale = 0.5,

    chains = 4,
    iter_sampling = 1000,
    iter_warmup = 1000,
    seed = 123,
    adapt_delta = 0.95,
    max_treedepth = 12
){

  if(!requireNamespace("cmdstanr", quietly = TRUE))
    stop("Package 'cmdstanr' is required.")

  # ----------------------------------------------------
  # Validate required columns
  # ----------------------------------------------------

  required_cols <- c("home_team","away_team","home_goals","away_goals")
  missing <- setdiff(required_cols, names(data))

  if(length(missing) > 0)
    stop("data is missing required columns: ",
         paste(missing, collapse=", "))

  # ----------------------------------------------------
  # Validate X
  # ----------------------------------------------------

  if(!is.matrix(X) || !is.numeric(X))
    stop("X must be a numeric matrix. Encode categorical variables first.")

  if(anyNA(X))
    stop("X contains NA values. Remove or impute them before fitting.")

  N <- nrow(data)

  if(nrow(X) != N)
    stop(paste0(
      "X must have ", N, " rows (one per match) but has ", nrow(X), "."
    ))

  # ----------------------------------------------------
  # Warn about X_prev behaviour
  # ----------------------------------------------------

  if(predict > 0 && !is.null(X_prev))
    warning("X_prev argument ignored when predict > 0: X_prev is sliced from X.")

  K <- ncol(X)

  # ----------------------------------------------------
  # Team indexing
  # ----------------------------------------------------

  teams <- sort(unique(c(data$home_team, data$away_team)))
  nteams <- length(teams)

  team_index <- setNames(seq_along(teams), teams)

  team1 <- team_index[data$home_team]
  team2 <- team_index[data$away_team]

  # ----------------------------------------------------
  # Outcomes
  # ----------------------------------------------------

  y <- cbind(data$home_goals, data$away_goals)

  # ----------------------------------------------------
  # Train / prediction split
  # ----------------------------------------------------

  if(predict > 0){

    N_prev <- predict
    N_train <- N - predict

    X_train <- X[1:N_train, , drop = FALSE]
    X_prev <- X[(N_train+1):N, , drop = FALSE]

    y <- y[1:N_train, , drop = FALSE]

    team1_prev <- team1[(N_train+1):N]
    team2_prev <- team2[(N_train+1):N]

    team1 <- team1[1:N_train]
    team2 <- team2[1:N_train]

    N <- N_train

  } else {

    N_prev <- 0
    team1_prev <- integer(0)
    team2_prev <- integer(0)

    X_train <- X
    X_prev <- matrix(0,0,K)

  }

  # ----------------------------------------------------
  # Ranking placeholder
  # ----------------------------------------------------

  instants_rank <- rep(1, N)
  ntimes_rank <- 1
  ranking <- matrix(0, ntimes_rank, nteams)

  # ----------------------------------------------------
  # Stan data
  # ----------------------------------------------------

  stan_data <- list(

    N = N,
    N_prev = N_prev,

    y = y,

    nteams = nteams,
    team1 = team1,
    team2 = team2,

    team1_prev = team1_prev,
    team2_prev = team2_prev,

    instants_rank = instants_rank,
    ntimes_rank = ntimes_rank,
    ranking = ranking,

    K = K,
    X = X_train,
    X_prev = X_prev,

    ind_home = ind_home,
    mean_home = mean_home,
    sd_home = sd_home,

    beta_prior_mean = beta_prior_mean,
    beta_prior_sd = beta_prior_sd,

    prior_dist_num = prior_dist_num,
    prior_dist_sd_num = prior_dist_sd_num,

    hyper_df = hyper_df,
    hyper_location = hyper_location,
    hyper_sd_df = hyper_sd_df,
    hyper_sd_location = hyper_sd_location,
    hyper_sd_scale = hyper_sd_scale
  )

  # ----------------------------------------------------
  # Compile and run Stan
  # ----------------------------------------------------

  mod <- cmdstanr::cmdstan_model(model_file)

  fit <- mod$sample(
    data = stan_data,
    chains = chains,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )

  # ----------------------------------------------------
  # Return object
  # ----------------------------------------------------

  out <- list(
    fit = fit,
    stan_data = stan_data,
    data = data,
    teams = teams
  )

  class(out) <- "stanFootMulti"

  return(out)
}

