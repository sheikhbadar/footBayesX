# footBayesX  
### Multi-Covariate Extension of the footBayes Double Poisson Model

footBayesX extends the `footBayes` package by allowing **multiple covariates** to be incorporated into the hierarchical Bayesian double Poisson football model.

---

##  Motivation

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

##  What footBayesX Adds

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

#  Architecture

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

##  Why a Wrapper Was Necessary

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

  How To Use With Your Own Covariates

---

# 3️⃣ Example: Using Multiple Covariates
The original `footBayes` double Poisson model explains match outcomes using:

- Team attacking strength  
- Team defensive strength  
- Home advantage  
- (Optional) ranking dynamics  

However, modern football performance is often influenced by additional contextual features such as:

- Form trends  
- Momentum  
- Tactical intensity  
- Match progression  
- Interaction effects between variables  

The `footBayesX` extension allows users to incorporate **multiple match-level covariates** directly into the log-intensity of the scoring model.

This example demonstrates how to:

Mathematically, the scoring intensity is extended as:

\[
\log(\theta_{home}) = \text{home effect} + \text{attack}_i + \text{defense}_j + X\beta
\]

\[
\log(\theta_{away}) = \text{attack}_j + \text{defense}_i - X\beta
\]

Where:

- `X` is an \( N \times K \) design matrix of covariates  
- `β` are regression coefficients estimated in Stan  
- Each covariate can positively or negatively influence scoring intensity  

This example shows how to:

• Construct custom match-level covariates  
• Build the design matrix `X`  
• Fit the extended multi-covariate model  
• Compare against the baseline `footBayes` model  
• Perform posterior predictive checks  
• Reconstruct league rankings  

---


## 1) Load Packages

```r
library(dplyr)
library(footBayes)
library(footBayesX)
library(posterior)
library(loo)
```

---

## 2) Load Example Data

```r
data("italy")

italy_df <- as_tibble(italy) %>%
  select(Season, home, visitor, hgoal, vgoal) %>%
  filter(Season %in% c("2000","2001","2002"))

colnames(italy_df) <- c(
  "periods",
  "home_team",
  "away_team",
  "home_goals",
  "away_goals"
)
```

---

## 3) Create Custom Covariates

Here we construct:

- Goal Difference Trend  
- Season Progress  
- Interaction Term  

```r
italy_df <- italy_df %>%
  arrange(periods)

italy_df$GD_diff <- italy_df$home_goals - italy_df$away_goals
italy_df$SeasonProgress <- seq_len(nrow(italy_df)) / nrow(italy_df)
italy_df$Interaction <- italy_df$GD_diff * italy_df$SeasonProgress
```

---

## 4) Build Design Matrix (N × K)

```r
X <- scale(cbind(
  GD_diff = italy_df$GD_diff,
  SeasonProgress = italy_df$SeasonProgress,
  Interaction = italy_df$Interaction
))
```

- Rows must match match data  
- Columns represent covariates  
- Numeric matrix required  

---

## 5) Fit Baseline Model

```r
fit_baseline <- stan_foot(
  data = italy_df,
  model = "double_pois",
  predict = 0,
  chains = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  seed = 123
)
```

---

## 6) Fit Multi-Covariate Model

```r
stan_file <- system.file(
  "stan",
  "double_pois_multi.stan",
  package = "footBayesX"
)

fit_multi <- stan_foot_multi(
  data = italy_df,
  model_file = stan_file,
  X = X,
  predict = 0,
  chains = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  adapt_delta = 0.95,
  seed = 123
)
We can visually analyze the attack and defense effects using foot_abilities function

foot_abilities(fit_multi, italy_df2)

## Attack and Defense Effects

<!--- [Attack and Defense Plot](figures/figure1.png) --->!
```

---

#  Important Notes

- Users construct X externally.
- Covariates must align row-wise with match data.
- Scaling strongly recommended.
- Interpretation of β is on log-goal intensity scale.
- Hierarchical attack/defense structure remains unchanged.

---

#  Installation

```r
devtools::install_github("https://github.com/sheikhbadar/footBayesX")
```

---

#  License

MIT

---

#  Research Use

This extension enables:

- Multi-covariate football modeling
- Bayesian inference on match-level covariates
- LOO / WAIC model comparison
- Extension of footBayes without modifying original package

---

footBayesX preserves the structure of footBayes while removing its single-covariate limitation.



