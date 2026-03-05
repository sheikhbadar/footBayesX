#' Multi-Covariate Extension of footBayes Double Poisson
#'
#' Extends the footBayes double Poisson model to allow multiple match-level
#' covariates through a symmetric linear predictor term X * beta.
#'
#' @param data Data frame with columns:
#'   home_team, away_team, home_goals, away_goals.
#' @param model_file Path to the compiled Stan model file.
#'   Use \code{system.file("stan","double_pois_multi.stan", package="footBayesX")}.
#' @param X Numeric matrix of covariates (N x K). Should be standardised
#'   via \code{scale()} before passing. One row per match.
#' @param X_prev Optional. Ignored when \code{predict > 0} (X_prev is
#'   sliced automatically from X).
#' @param predict Integer. Number of last rows of \code{data} to use as
#'   out-of-sample prediction set. Default \code{0} (no prediction).
#' @param ind_home Integer 0 or 1. Include home advantage effect. Default \code{1}.
#' @param mean_home Prior mean for home effect. Default \code{0}.
#' @param sd_home Prior SD for home effect. Default \code{0.5}.
#' @param beta_prior_mean Prior mean for covariate coefficients beta. Default \code{0}.
#' @param beta_prior_sd Prior SD for covariate coefficients beta. Default \code{1}.
#'   Assumes X is standardised. Increase if X is not scaled.
#' @param prior_dist_num Integer 1-4. Prior family for team abilities.
#'   1=Normal, 2=Student-t, 3=Cauchy, 4=Laplace. Default \code{2}.
#' @param prior_dist_sd_num Integer 1-4. Prior family for ability SDs.
#'   1=Normal, 2=Student-t, 3=Cauchy, 4=Laplace. Default \code{2}.
#' @param hyper_df Degrees of freedom for Student-t ability prior.
#'   Only used when \code{prior_dist_num = 2}. Default \code{4}.
#' @param hyper_location Location hyperparameter for ability prior. Default \code{0}.
#' @param hyper_sd_df Degrees of freedom for Student-t SD prior.
#'   Only used when \code{prior_dist_sd_num = 2}. Default \code{3}.
#' @param hyper_sd_location Location hyperparameter for SD prior. Default \code{0}.
#' @param hyper_sd_scale Scale hyperparameter for SD prior. Default \code{0.5}.
#' @param chains Number of MCMC chains. Default \code{4}.
#' @param iter_sampling Number of post-warmup samples per chain. Default \code{1000}.
#' @param iter_warmup Number of warmup samples per chain. Default \code{1000}.
#' @param seed Random seed for reproducibility. Default \code{123}.
#' @param adapt_delta Stan NUTS target acceptance rate. Default \code{0.95}.
#'   Increase to \code{0.99} if divergences occur.
#' @param max_treedepth Stan NUTS maximum tree depth. Default \code{12}.
#' @param parallel_chains Number of chains to run in parallel. Default \code{1}.
#'   Cannot exceed \code{chains}. Set equal to \code{chains} for maximum speed.
#' @param init Initial values for Stan parameters. Must be \code{NULL} (Stan
#'   default), a named list, or a function returning a named list. Default \code{NULL}.
#'
#' @return A list of class \code{"stanFootMulti"} containing:
#' \itemize{
#'   \item \code{fit}: CmdStanMCMC object from cmdstanr
#'   \item \code{stan_data}: list passed to Stan
#'   \item \code{data}: original input data frame
#'   \item \code{teams}: sorted vector of team names
#' }
#'
#' @examples
#' \dontrun{
#' stan_file <- system.file("stan", "double_pois_multi.stan", package = "footBayesX")
#'
#' X <- scale(cbind(psych_diff = rnorm(100), tact_diff = rnorm(100)))
#'
#' fit <- stan_foot_multi(
#'   data              = my_data,
#'   model_file        = stan_file,
#'   X                 = X,
#'   predict           = 10,
#'   prior_dist_num    = 2,
#'   prior_dist_sd_num = 2,
#'   hyper_df          = 4,
#'   chains            = 4,
#'   parallel_chains   = 4,
#'   seed              = 42
#' )
#' }
#'
#' @export

stan_foot_multi <- function(
    data,
    model_file,
    X,
    X_prev            = NULL,
    predict           = 0,

    # ── Home effect #
    ind_home          = 1,
    mean_home         = 0,
    sd_home           = 0.5,

    # ── Beta (covariate) prior #
    beta_prior_mean   = 0,
    beta_prior_sd     = 1,

    # ── Team ability priors #
    prior_dist_num    = 2,
    prior_dist_sd_num = 2,
    hyper_df          = 4,
    hyper_location    = 0,
    hyper_sd_df       = 3,
    hyper_sd_location = 0,
    hyper_sd_scale    = 0.5,

    # ── Sampler settings #
    chains            = 4,
    iter_sampling     = 1000,
    iter_warmup       = 1000,
    seed              = 123,
    adapt_delta       = 0.95,
    max_treedepth     = 12,
    parallel_chains   = 1,
    init              = NULL
) {

  # ── 1. Check cmdstanr #
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop(
      "Package 'cmdstanr' is required.\n",
      "Install with: install.packages('cmdstanr')"
    )

  # ── 2. Validate model file #
  if (!file.exists(model_file))
    stop(
      "Stan model file not found: ", model_file, "\n",
      "Use system.file('stan', 'double_pois_multi.stan', package = 'footBayesX')"
    )

  # ── 3. Validate data #
  if (!is.data.frame(data))
    stop("'data' must be a data.frame.")

  required_cols <- c("home_team", "away_team", "home_goals", "away_goals")
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop("data is missing required columns: ",
         paste(missing_cols, collapse = ", "))

  if (!is.numeric(data$home_goals) || !is.numeric(data$away_goals))
    stop("home_goals and away_goals must be numeric.")

  # ── 4. Validate predict #
  if (!is.numeric(predict) || predict < 0 || predict %% 1 != 0)
    stop("'predict' must be a non-negative integer.")

  if (predict >= nrow(data))
    stop("'predict' (", predict, ") must be less than nrow(data) (",
         nrow(data), ").")

  # ── 5. Validate X #
  if (!is.matrix(X) || !is.numeric(X))
    stop("X must be a numeric matrix. Use cbind() then scale().")

  if (anyNA(X))
    stop("X contains NA values. Remove or impute before fitting.")

  if (!all(is.finite(X)))
    stop("X contains Inf or NaN values. Check your covariate construction.")

  if (nrow(X) != nrow(data))
    stop("X must have ", nrow(data), " rows (one per match in data) ",
         "but has ", nrow(X), " rows.")

  # ── 6. Validate prior arguments #
  if (!(prior_dist_num %in% 1:4))
    stop("prior_dist_num must be 1 (Normal), 2 (Student-t), ",
         "3 (Cauchy), or 4 (Laplace).")

  if (!(prior_dist_sd_num %in% 1:4))
    stop("prior_dist_sd_num must be 1 (Normal), 2 (Student-t), ",
         "3 (Cauchy), or 4 (Laplace).")

  if (beta_prior_sd <= 0)
    stop("beta_prior_sd must be strictly positive.")

  if (sd_home <= 0)
    stop("sd_home must be strictly positive.")

  if (hyper_sd_scale <= 0)
    stop("hyper_sd_scale must be strictly positive.")

  if (hyper_df <= 0)
    stop("hyper_df must be strictly positive.")

  if (hyper_sd_df <= 0)
    stop("hyper_sd_df must be strictly positive.")

  # ── 7. Validate ind_home (fixed precedence bug)#
  if (!(ind_home %in% c(0L, 1L)))
    stop("ind_home must be 0 or 1.")

  # ── 8. Validate sampler arguments #
  if (chains < 1)
    stop("chains must be >= 1.")

  if (iter_sampling < 1)
    stop("iter_sampling must be >= 1.")

  if (iter_warmup < 1)
    stop("iter_warmup must be >= 1.")

  if (parallel_chains < 1)
    stop("parallel_chains must be >= 1.")

  if (parallel_chains > chains)
    stop("parallel_chains (", parallel_chains,
         ") cannot exceed chains (", chains, ").")

  if (adapt_delta <= 0 || adapt_delta >= 1)
    stop("adapt_delta must be strictly between 0 and 1.")

  # ── 9. Validate init #
  if (!is.null(init) && !is.list(init) && !is.function(init))
    stop("init must be NULL, a named list, or a function returning a named list.")

  # ── 10. Warn about X_prev #
  if (predict > 0 && !is.null(X_prev))
    warning("X_prev is ignored when predict > 0. ",
            "X_prev is sliced automatically from X.")

  # ── 11. Covariate column names (fix NULL display) #
  cov_names <- colnames(X)
  if (is.null(cov_names))
    cov_names <- paste0("X", seq_len(ncol(X)))

  K <- ncol(X)

  # ── 12. Inform user of settings #
  prior_names <- c("Normal", "Student-t", "Cauchy", "Laplace")
  message(
    "\n[footBayesX] Running stan_foot_multi\n",
    "  Matches (train)   : ", nrow(data) - predict, "\n",
    "  Matches (predict) : ", predict, "\n",
    "  Covariates (K)    : ", K,
    " [", paste(cov_names, collapse = ", "), "]\n",
    "  Ability prior     : ", prior_names[prior_dist_num],
    if (prior_dist_num == 2) paste0(" (df=", hyper_df, ")") else "", "\n",
    "  SD prior          : ", prior_names[prior_dist_sd_num],
    if (prior_dist_sd_num == 2) paste0(" (df=", hyper_sd_df, ")") else "", "\n",
    "  beta prior        : Normal(", beta_prior_mean, ", ", beta_prior_sd, ")\n",
    "  Home effect       : ", ifelse(ind_home == 1, "Yes", "No"), "\n",
    "  Chains            : ", chains,
    " | Parallel: ", parallel_chains, "\n",
    "  Iterations        : warmup=", iter_warmup,
    " | sampling=", iter_sampling, "\n",
    "  adapt_delta       : ", adapt_delta,
    " | max_treedepth: ", max_treedepth, "\n"
  )

  # ── 13. Team indexing #
  teams      <- sort(unique(c(data$home_team, data$away_team)))
  nteams     <- length(teams)
  team_index <- setNames(seq_along(teams), teams)
  team1      <- unname(team_index[data$home_team])
  team2      <- unname(team_index[data$away_team])

  # ── 14. Outcome matrix #
  y <- cbind(as.integer(data$home_goals),
             as.integer(data$away_goals))

  # ── 15. Train / prediction split #
  N_total <- nrow(data)

  if (predict > 0) {

    N_train  <- N_total - predict
    N_prev   <- predict

    y_train      <- y[1:N_train,            , drop = FALSE]
    X_train      <- X[1:N_train,            , drop = FALSE]
    X_pred       <- X[(N_train + 1):N_total, , drop = FALSE]

    team1_prev   <- team1[(N_train + 1):N_total]
    team2_prev   <- team2[(N_train + 1):N_total]
    team1        <- team1[1:N_train]
    team2        <- team2[1:N_train]
    N            <- N_train

  } else {

    N            <- N_total
    N_prev       <- 0L
    y_train      <- y
    X_train      <- X
    X_pred       <- matrix(0.0, 0L, K)
    team1_prev   <- integer(0)
    team2_prev   <- integer(0)

  }

  # ── 16. Ranking placeholder #
  instants_rank <- rep(1L, N)
  ntimes_rank   <- 1L
  ranking       <- matrix(0.0, ntimes_rank, nteams)

  # ── 17. Build Stan data list #
  stan_data <- list(

    N            = N,
    N_prev       = N_prev,
    y            = y_train,

    nteams       = nteams,
    team1        = team1,
    team2        = team2,
    team1_prev   = team1_prev,
    team2_prev   = team2_prev,

    instants_rank = instants_rank,
    ntimes_rank   = ntimes_rank,
    ranking       = ranking,

    K            = K,
    X            = X_train,
    X_prev       = X_pred,

    ind_home     = ind_home,
    mean_home    = mean_home,
    sd_home      = sd_home,

    beta_prior_mean   = beta_prior_mean,
    beta_prior_sd     = beta_prior_sd,

    prior_dist_num    = prior_dist_num,
    prior_dist_sd_num = prior_dist_sd_num,

    hyper_df          = hyper_df,
    hyper_location    = hyper_location,
    hyper_sd_df       = hyper_sd_df,
    hyper_sd_location = hyper_sd_location,
    hyper_sd_scale    = hyper_sd_scale
  )

  # ── 18. Compiling the Stan model (quiet — no recompile if up to date) #
  mod <- cmdstanr::cmdstan_model(model_file, quiet = TRUE)

  # ── 19. Sample #
  fit <- mod$sample(
    data            = stan_data,
    chains          = chains,
    parallel_chains = parallel_chains,
    iter_sampling   = iter_sampling,
    iter_warmup     = iter_warmup,
    seed            = seed,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    init            = init
  )

  # ── 20. Return #
  out <- list(
    fit       = fit,
    stan_data = stan_data,
    data      = data,
    teams     = teams
  )

  class(out) <- "stanFootMulti"

  return(out)
}
