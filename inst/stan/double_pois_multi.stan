
data{

  int N;
  int<lower=0> N_prev;

  array[N,2] int y;

  int nteams;
  array[N] int team1;
  array[N] int team2;

  array[N_prev] int team1_prev;
  array[N_prev] int team2_prev;

  array[N] int instants_rank;
  int ntimes_rank;
  matrix[ntimes_rank,nteams] ranking;

  // ----------------------------------------------------
  // MATCH COVARIATES
  // ----------------------------------------------------
  // K = number of covariates
  // X  = covariates for training matches
  // X_prev = covariates for prediction matches
  // Covariates should typically be standardized in R
  // ----------------------------------------------------

  int<lower=1> K;
  matrix[N, K] X;
  matrix[N_prev, K] X_prev;

  // ----------------------------------------------------
  // HOME EFFECT CONTROL
  // ----------------------------------------------------

  int<lower=0, upper=1> ind_home;
  real mean_home;
  real<lower=1e-8> sd_home;

  // ----------------------------------------------------
  // BETA PRIOR (user controlled via wrapper)
  // ----------------------------------------------------
  // Default: Normal(0,1) assuming X is standardized
  // beta_prior_sd must be strictly positive
  // ----------------------------------------------------

  real beta_prior_mean;
  real<lower=1e-8> beta_prior_sd;

  // ----------------------------------------------------
  // TEAM EFFECT PRIORS
  // ----------------------------------------------------

  int<lower=1,upper=4> prior_dist_num;
  int<lower=1,upper=4> prior_dist_sd_num;

  real<lower=0> hyper_df;
  real hyper_location;

  real<lower=0> hyper_sd_df;
  real hyper_sd_location;
  real<lower=1e-8> hyper_sd_scale;
}

parameters{

  vector[nteams] att_raw;
  vector[nteams] def_raw;

  real<lower=1e-8> sigma_att;
  real<lower=1e-8> sigma_def;

  real home;
  real gamma;

  // ----------------------------------------------------
  // COVARIATE COEFFICIENTS
  // ----------------------------------------------------

  vector[K] beta;
}

transformed parameters{

  real adj_h_eff;

  vector[nteams] att;
  vector[nteams] def;

  array[N] vector[2] theta;

  // sum-to-zero constraint for identifiability
  for (t in 1:nteams){
    att[t] = att_raw[t] - mean(att_raw);
    def[t] = def_raw[t] - mean(def_raw);
  }

  adj_h_eff = home * ind_home;

  // ----------------------------------------------------
  // SCORING INTENSITIES
  // ----------------------------------------------------

  for (n in 1:N){

    theta[n,1] =
      exp(adj_h_eff +
          att[team1[n]] +
          def[team2[n]] +
          X[n] * beta +
          (gamma/2) *
          (ranking[instants_rank[n],team1[n]] -
           ranking[instants_rank[n],team2[n]]));

    theta[n,2] =
      exp(att[team2[n]] +
          def[team1[n]] -
          X[n] * beta -
          (gamma/2) *
          (ranking[instants_rank[n],team1[n]] -
           ranking[instants_rank[n],team2[n]]));
  }
}

model{

  // ----------------------------------------------------
  // TEAM ABILITY PRIORS
  // ----------------------------------------------------

  for (t in 1:nteams){

    if (prior_dist_num == 1){
      target += normal_lpdf(att_raw[t] | hyper_location, sigma_att);
      target += normal_lpdf(def_raw[t] | hyper_location, sigma_def);
    }
    else if (prior_dist_num == 2){
      target += student_t_lpdf(att_raw[t] | hyper_df,
                               hyper_location,
                               sigma_att);
      target += student_t_lpdf(def_raw[t] | hyper_df,
                               hyper_location,
                               sigma_def);
    }
    else if (prior_dist_num == 3){
      target += cauchy_lpdf(att_raw[t] | hyper_location, sigma_att);
      target += cauchy_lpdf(def_raw[t] | hyper_location, sigma_def);
    }
    else{
      target += double_exponential_lpdf(att_raw[t] |
                                        hyper_location,
                                        sigma_att);
      target += double_exponential_lpdf(def_raw[t] |
                                        hyper_location,
                                        sigma_def);
    }
  }

  // ----------------------------------------------------
  // SD PRIORS
  // ----------------------------------------------------

  if (prior_dist_sd_num == 1){
    target += normal_lpdf(sigma_att |
                          hyper_sd_location,
                          hyper_sd_scale);
    target += normal_lpdf(sigma_def |
                          hyper_sd_location,
                          hyper_sd_scale);
  }
  else if (prior_dist_sd_num == 2){
    target += student_t_lpdf(sigma_att |
                             hyper_sd_df,
                             hyper_sd_location,
                             hyper_sd_scale);
    target += student_t_lpdf(sigma_def |
                             hyper_sd_df,
                             hyper_sd_location,
                             hyper_sd_scale);
  }
  else if (prior_dist_sd_num == 3){
    target += cauchy_lpdf(sigma_att |
                          hyper_sd_location,
                          hyper_sd_scale);
    target += cauchy_lpdf(sigma_def |
                          hyper_sd_location,
                          hyper_sd_scale);
  }
  else{
    target += double_exponential_lpdf(sigma_att |
                                      hyper_sd_location,
                                      hyper_sd_scale);
    target += double_exponential_lpdf(sigma_def |
                                      hyper_sd_location,
                                      hyper_sd_scale);
  }

  // ----------------------------------------------------
  // FIXED EFFECTS
  // ----------------------------------------------------

  target += normal_lpdf(home | mean_home, sd_home);
  target += normal_lpdf(gamma | 0,1);

  // ----------------------------------------------------
  // COVARIATE PRIOR
  // ----------------------------------------------------

  beta ~ normal(beta_prior_mean, beta_prior_sd);

  // ----------------------------------------------------
  // LIKELIHOOD
  // ----------------------------------------------------

  target += poisson_lpmf(y[,1] | theta[,1]);
  target += poisson_lpmf(y[,2] | theta[,2]);
}

generated quantities{

  array[N,2] int y_rep;
  array[N_prev,2] int y_prev;

  array[N_prev] vector[2] theta_prev;

  vector[N] log_lik;
  vector[N] diff_y_rep;

  for (n in 1:N){

    y_rep[n,1] = poisson_rng(theta[n,1]);
    y_rep[n,2] = poisson_rng(theta[n,2]);

    log_lik[n] =
      poisson_lpmf(y[n,1] | theta[n,1]) +
      poisson_lpmf(y[n,2] | theta[n,2]);

    diff_y_rep[n] = y_rep[n,1] - y_rep[n,2];
  }

  if (N_prev > 0){

    for (n in 1:N_prev){

      theta_prev[n,1] =
        exp(adj_h_eff +
            att[team1_prev[n]] +
            def[team2_prev[n]] +
            X_prev[n] * beta +
            (gamma/2) *
            (ranking[instants_rank[N],team1_prev[n]] -
             ranking[instants_rank[N],team2_prev[n]]));

      theta_prev[n,2] =
        exp(att[team2_prev[n]] +
            def[team1_prev[n]] -
            X_prev[n] * beta -
            (gamma/2) *
            (ranking[instants_rank[N],team1_prev[n]] -
             ranking[instants_rank[N],team2_prev[n]]));

      y_prev[n,1] = poisson_rng(theta_prev[n,1]);
      y_prev[n,2] = poisson_rng(theta_prev[n,2]);
    }
  }
}
