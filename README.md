# footBayesX  
### Multi-Covariate Extension of the footBayes Double Poisson Model

footBayesX extends the `footBayes` package by allowing **multiple covariates** to be incorporated into the hierarchical Bayesian double Poisson football model.

---

## 📌 Motivation

The original `footBayes` double Poisson model defines goal intensity as:

λ_home = exp(h + α_i + δ_j)

Where:
- h = home advantage  
- α_i = attacking strength of team i  
- δ_j = defensive strength of team j  

### Limitation

The original implementation does **not support a general covariate matrix**.  
Users cannot directly incorporate multiple engineered features into the model.

---

## 🚀 What footBayesX Adds

We generalize the intensity to:

λ_home = exp(h + α_i + δ_j + X_n β)

Where:
- X_n = row of a covariate matrix (N × K)
- β = vector of coefficients
- K = number of covariates

This enables inclusion of:

- Match-level engineered features  
- Rest days  
- Travel distance  
- Squad rotation  
- Weather  
- Market value  
- Expected goals  
- Any custom-designed covariate  

---

# 🏗 Architecture

footBayesX consists of two components:

---

## 1️⃣ Modified Stan Model

File location:
```
inst/stan/double_pois_multi.stan
```

### Added to Stan:

```stan
int<lower=1> K;
matrix[N, K] X;
matrix[N_prev, K] X_prev;
vector[K] beta;
```

Linear predictor modification:

```stan
X[n] * beta
```

The hierarchical structure remains identical to footBayes:
- Team attack random effects
- Team defense random effects
- Home advantage
- Prior structure

Only the linear predictor is generalized.

---

## 2️⃣ Wrapper Function: stan_foot_multi()

File location:
```
R/stan_foot_multi.R
```

### What the wrapper does:

1. Validates input data
2. Converts team names to numeric indices
3. Handles train/test splitting
4. Constructs the Stan data list
5. Passes covariate matrix X
6. Calls cmdstanr
7. Returns fitted model object

---

## 🔧 Why a Wrapper Was Necessary

The original `footBayes::stan_foot()`:

- Internally builds Stan data
- Does not accept arbitrary X matrices
- Limits extensibility

The custom wrapper:

- Accepts any numeric matrix X
- Supports K ≥ 1 covariates
- Preserves hierarchical structure
- Maintains compatibility with footBayes framework

---

# 🧠 How To Use With Your Own Covariates

You can use **any engineered features**.

---

## Step 1 — Prepare Data

Data must contain:

- periods  
- home_team  
- away_team  
- home_goals  
- away_goals  

---

## Step 2 — Create Covariate Matrix

Example:

```r
X <- scale(cbind(
  xG_Diff            = home_xG - away_xG,
  ShotDistance_Diff  = home_avg_shot_distance - away_avg_shot_distance,
  PressingIntensity  = home_ppda - away_ppda,
  RestDaysDiff       = home_rest_days - away_rest_days
))
```

Rules:
- Must be numeric matrix
- Rows must match match data
- Columns = number of covariates
- Standardization recommended

---

## Step 3 — Fit Model

```r
library(footBayesX)

stan_file <- system.file(
  "stan",
  "double_pois_multi.stan",
  package = "footBayesX"
)

fit <- stan_foot_multi(
  data = match_data,
  model_file = stan_file,
  X = X,
  predict = 0,
  chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)
```

---

## Step 4 — Inspect Covariate Effects

```r
posterior::summarise_draws(
  fit$fit$draws(),
  variable = grep("beta", posterior::variables(fit$fit$draws()), value = TRUE)
)
```

Each beta[k] corresponds to column k of your X matrix.

---


# 📌 Important Notes

- Users construct X externally.
- Covariates must align row-wise with match data.
- Scaling strongly recommended.
- Interpretation of β is on log-goal intensity scale.
- Hierarchical attack/defense structure remains unchanged.

---

# 📦 Installation

```r
devtools::install_github("https://github.com/sheikhbadar/footBayesX")
```

---

# 📝 License

MIT

---

# 🎓 Research Use

This extension enables:

- Multi-feature football modeling
- Bayesian inference on match-level covariates
- LOO / WAIC model comparison
- Extension of footBayes without modifying original package

---

footBayesX preserves the structure of footBayes while removing its single-covariate limitation.
