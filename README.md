# footBayesX  
### Multi-Covariate Extension of the footBayes Double Poisson Model
footBayesX extends the `footBayes` package by allowing **multiple covariates** to be incorporated into the hierarchical Bayesian double Poisson football model.
---
##  Motivation
The original `footBayes` double Poisson model defines goal scoring intensity as:
```
λ_home = exp(h + α_i + δ_j)
λ_away = exp(α_j + δ_i)
```
Where:
- `h` = home advantage
- `α_i` = attacking strength of team i (home)
- `δ_j` = defensive weakness of team j (away)
### Limitation
The original implementation does **not support a general covariate matrix**.
Users cannot directly incorporate multiple engineered match-level features into the model.
---
##  What footBayesX Adds
footBayesX generalises the scoring intensity to:
```
log(θ_home) = h + α_i + δ_j + X[n] β
log(θ_away) =     α_j + δ_i - X[n] β
```
Where:
- `X[n]` = row n of an N × K covariate matrix (one row per match)
- `β` = vector of K coefficients estimated by Stan
- `K` = number of user-supplied features
### The Symmetric Sign Design
The covariate enters with a **plus sign** in the home intensity and a **minus sign** in the away intensity.
This is not arbitrary. It follows the exact same logic already present in footBayes for the ranking term:
```
footBayes ranking:      + (γ/2) × (rank_home − rank_away)   in θ_home
                        − (γ/2) × (rank_home − rank_away)   in θ_away
footBayesX covariates:  + X[n] × β                          in θ_home
                        − X[n] × β                          in θ_away
```
footBayesX simply **replaces the ranking difference with user-supplied features** — same principle, generalised.
### The Key Requirement
> **X must be constructed as (home value − away value) difference features.**
When `X[n, k] > 0`, the home team has an advantage on feature k:
- `+Xβ` increases home scoring intensity 
- `−Xβ` decreases away scoring intensity 
When `X[n, k] < 0`, the away team has an advantage:
- `+Xβ` decreases home scoring intensity 
- `−Xβ` increases away scoring intensity 
---
##  The X Matrix — What It Is and How to Build It
### Structure
`X` is a numeric matrix of size **N × K**, where:
- **N** = number of matches (rows must align exactly with match data)
- **K** = number of features (columns)
Each row `X[n]` contains the feature values for match `n`.
### What β Means
Each coefficient `β[k]` describes how feature `k` affects the scoreline:
| β[k] value | Meaning |
|---|---|
| `β[k] > 0` | Home team having more of feature k → more home goals, fewer away goals |
| `β[k] < 0` | Home team having more of feature k → fewer home goals, more away goals |
| `β[k] ≈ 0` | Feature k has no meaningful effect on the scoreline |
**Example:** If `β[1] = 0.25` and feature 1 is a form difference:
- A home team with +1 SD better recent form scores `exp(0.25) ≈ 1.28×` more goals
- The same away team scores `exp(−0.25) ≈ 0.78×` fewer goals
### Feature Categories Users Can Build
All features must be constructed as **(home value − away value)** before passing to Stan.
#### Category 1: Form and Performance Differences
| Feature | How to construct in R |
|---|---|
| Recent points difference | `home_last5_points − away_last5_points` |
| Goals scored difference | `home_avg_goals − away_avg_goals` |
| Goal difference trend | `home_avg_GD − away_avg_GD` |
| Win streak difference | `home_win_streak − away_win_streak` |
| Shot accuracy difference | `home_shots_on_target − away_shots_on_target` |
#### Category 2: Tactical and Physical Differences
| Feature | How to construct in R |
|---|---|
| Possession difference | `home_avg_possession − away_avg_possession` |
| Expected goals difference | `home_xG_last5 − away_xG_last5` |
| Pressing intensity difference | `home_PPDA − away_PPDA` |
| Corners difference | `home_avg_corners − away_avg_corners` |
#### Category 3: Squad and Fitness Differences
| Feature | How to construct in R |
|---|---|
| Injuries difference | `home_injured − away_injured` |
| Days rest difference | `home_days_rest − away_days_rest` |
| Squad value difference | `home_squad_value − away_squad_value` |
| Suspensions difference | `home_suspended − away_suspended` |
#### Category 4: Contextual and Market Differences
| Feature | How to construct in R |
|---|---|
| Elo rating difference | `home_elo − away_elo` |
| League position difference | `home_position − away_position` |
| Points in table difference | `home_points − away_points` |
| Market value difference | `home_market_value − away_market_value` |
### Rules for Building X
1. **Always use lagged features** — compute rolling statistics using only past matches to avoid data leakage
2. **Always standardise using training data** — `scale(X_train)` then apply the same centre and scale to `X_test`
3. **Features must be differences** — absolute values (e.g. home stadium capacity alone) are not compatible with the `+Xβ / −Xβ` design
4. **Handle missing values** — early season rows may have NA due to insufficient history; replace with 0 (neutral) or league average
5. **Check collinearity** — if two features are highly correlated (|r| > 0.7), consider dropping one
```r
# Check collinearity before fitting
cor(X_train)
```
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
All covariates must be:
- Constructed as **(home − away) differences**
- Built using **only past match information** (lagged) to avoid data leakage
- **Standardised** using training data statistics only
Here we construct two difference features:
- **Form Difference**: rolling 5-match points average (home − away)
- **Goal Difference Trend**: rolling 5-match average goal difference (home − away)
```r
library(dplyr)
library(zoo)
# Step 1 — Build team-level panel (one row per team per match)
team_long <- italy_df %>%
  pivot_longer(
    c(home_team, away_team),
    names_to  = "side",
    values_to = "team"
  ) %>%
  mutate(
    is_home = side == "home_team",
    GF      = ifelse(is_home, home_goals, away_goals),
    GA      = ifelse(is_home, away_goals, home_goals),
    Points  = case_when(GF > GA ~ 3, GF == GA ~ 1, TRUE ~ 0),
    GD      = GF - GA
  ) %>%
  arrange(team, periods)
# Step 2 — Compute lagged rolling features (lag prevents data leakage)
team_long <- team_long %>%
  group_by(team) %>%
  mutate(
    form_raw = lag(rollmean(Points, 5, fill = NA, align = "right")),
    gd_raw   = lag(rollapply(GD, 5, mean, fill = NA, align = "right"))
  ) %>%
  ungroup()
# Step 3 — Join back to match level as (home − away) differences
home_feat <- team_long %>%
  filter(is_home) %>%
  select(periods, home_team = team, FormHome = form_raw, GDHome = gd_raw)
away_feat <- team_long %>%
  filter(!is_home) %>%
  select(periods, away_team = team, FormAway = form_raw, GDAway = gd_raw)
italy_df <- italy_df %>%
  left_join(home_feat, by = c("periods", "home_team")) %>%
  left_join(away_feat, by = c("periods", "away_team")) %>%
  mutate(
    FormDiff = coalesce(FormHome - FormAway, 0),  # 0 = neutral when history unavailable
    GDDiff   = coalesce(GDHome   - GDAway,   0)
  )
```
---
## 4) Build Design Matrix (N × K)
Split into training and test, then standardise using **training statistics only**:
```r
train_df <- italy_df %>% filter(periods < 2002)
test_df  <- italy_df %>% filter(periods == 2002)
# Standardise on training data
form_scale <- scale(train_df$FormDiff)
gd_scale   <- scale(train_df$GDDiff)
X_train <- cbind(
  form = form_scale[, 1],
  gd   = gd_scale[, 1]
)
# Apply SAME centre and scale to test data (no leakage)
X_test <- cbind(
  form = scale(test_df$FormDiff,
               center = attr(form_scale, "scaled:center"),
               scale  = attr(form_scale, "scaled:scale")),
  gd   = scale(test_df$GDDiff,
               center = attr(gd_scale, "scaled:center"),
               scale  = attr(gd_scale, "scaled:scale"))
)
# Combine in chronological order
all_df <- bind_rows(train_df, test_df) %>% arrange(periods)
X_all  <- rbind(X_train, X_test)
# Sanity checks
stopifnot(nrow(all_df) == nrow(X_all))
stopifnot(ncol(X_all) == 2)
anyNA(X_all)  # should be FALSE
```
**Requirements:**
- `X` must be a numeric matrix
- Rows must align exactly with `data` row order
- Columns represent features (K ≥ 1)
- Standardisation is strongly recommended when `beta_prior_sd = 1`
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
### The (home − away) Requirement
`X` must contain **(home − away) difference features**. The model uses `+Xβ` for home intensity and `−Xβ` for away intensity. This means a positive value in `X[n, k]` is interpreted as a home-team advantage on feature k. Absolute features (e.g. home stadium capacity without subtracting away capacity) are not compatible with this design.

### Collinearity
If two features are highly correlated (|r| &gt; 0.7), consider dropping one. High collinearity inflates posterior uncertainty on `β`. Check with `cor(X_train)` before fitting.
### Handling Missing Values
Early-season matches may have NA values due to insufficient rolling history. Replace these with `0` (neutral difference) using `coalesce(..., 0)` in R, or impute with the training league average.
### Interpreting β
`β` coefficients are on the **log-goal intensity scale**. A coefficient of `β[k] = 0.2` means a 1-unit increase in feature k multiplies expected home goals by `exp(0.2) ≈ 1.22` and divides expected away goals by the same factor.
### Hierarchical structure is unchanged
`att`, `def`, `home`, and `gamma` remain identical to the original footBayes model. `Xβ` is an additive term on top of the existing structure — it does not replace team abilities.
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
