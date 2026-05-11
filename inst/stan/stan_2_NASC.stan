// Step 2: Network-Aware Synthetic Control (NASC) Model
// The NASC penalty is always active in this file. If the user does not want
// the penalty, model1.stan is used instead (routed by Base.R on nasc_penalty).
//
// The bias-correction toggle remains: penalty-on, bias-off is a meaningful
// configuration (NASC weights that account for network contamination via the
// penalty, but with tau left on the untransformed scale).
//
// Generated-quantity names (y_sim_pre, y_counterfactual) match model1.stan so
// a single R post-processing path handles both.
//
// COVARIATE MATCHING:
//   This file mirrors model1.stan's optional augmented matching: when
//   covariates are supplied, the simplex weights w must reconcile not only
//   the pre-treatment outcome trajectory but also the pre-period covariate
//   means of the treated unit against a w-weighted mix of donor covariate
//   means. Each covariate is standardized by its donor-side mean/SD (the
//   treated value is standardized by the same donor-derived statistics) so
//   it enters the likelihood on the same dimensionless scale as the outcome
//   features. K_cov = 0 disables the augmented matching and reproduces the
//   pre-existing behavior exactly.
//
// PREDICTOR IMPORTANCE WEIGHTS (v_cov):
//   Optional length-K_cov non-negative vector that re-weights each covariate
//   row's contribution to the matching loss, mirroring the V-matrix in the
//   classical Abadie-Diamond-Hainmueller SCM and bsynth's `vs` argument.
//   The k-th covariate row is treated as Normal(., sigma_sc / sqrt(v_cov[k]))
//   so v_cov[k] = 1 (default) recovers equal weighting and v_cov[k] = 0
//   removes covariate k from the matching loss without dropping it from the
//   data. Outcome rows always have implicit weight 1.
//
// LAMBDA SCALING (unit-information reparameterization):
//   The NASC penalty term -lambda * <w, |s|> is on a different natural scale
//   than the Gaussian likelihood, which scales linearly in the effective
//   number of likelihood contributions
//
//       n_eff = T0 + sum(v_cov)
//
//   (T0 outcome rows each contributing weight 1, plus K_cov covariate rows
//   contributing weight v_cov[k] to the Fisher information in w). Sampling
//   lambda directly on the raw scale puts the prior on a parameter whose
//   data-informed region is O(n_eff), causing severe prior sensitivity when
//   the prior (e.g. gamma(2,1)) places negligible mass there.
//
//   The fix is to sample lambda_tilde on a unit-information scale and
//   recover lambda = lambda_tilde * n_eff. With this rescaling,
//   lambda_tilde ~ O(1) means "penalty competes with the likelihood at
//   O(1)", which is what a unit-scale prior like gamma(2,1) intends. The
//   raw lambda is still exposed in generated quantities so downstream R
//   code that extracts "lambda" continues to work unchanged.
//
//   Note: we do NOT also divide by max|s| (or any other summary of |s|).
//   The contamination vector s = rho * (I - rho W)^{-1} w_{J1} already
//   shrinks toward zero as rho -> 0, which is the correct statement of
//   "no contamination problem when there is no spatial spillover".
//   Dividing the penalty by max|s| would cancel exactly this rho-dependence
//   and force the penalty to act with the same strength regardless of how
//   much contamination is actually present. The current scaling preserves
//   that: when rho is small, the penalty is correctly weak (and lambda is
//   accordingly weakly identified, with its posterior reflecting the prior);
//   when rho is large, the penalty has real bite and lambda becomes
//   data-informed. The pooled MI posterior across rho draws is then a
//   honest mixture that reflects rho's role in identification.

data {
  // NASC Specific Data
  int<lower=0> J;
  int<lower=1> T0;
  matrix[J + 1, T0] Y_panel;     // Pre-treatment panel: rows 1..J donors, row J+1 treated

  matrix[J, J] W_J;              // Donor-donor block of W
  vector[J] w_J1;                // Donor-to-treated column of W

  int<lower=1> T_post;
  matrix[T_post, J] Y0_post;     // Post-treatment donor outcomes (raw scale)
  vector[T_post] Y1_post;        // Post-treatment treated outcome (raw scale)

  // rho is provided as data (one draw from Stage 1, or an exogenous value)
  real<lower=-1, upper=1> rho;

  // Bias-correction toggle. Penalty toggle has been removed: this file is
  // only used when the penalty is on.
  int<lower=0, upper=1> use_bias_correction;

  // ---- Optional covariate-matching inputs ----
  // K_cov = 0 disables the augmented matching (zero-sized arrays).
  int<lower=0> K_cov;
  matrix[K_cov, J] X_cov0;       // Donor covariate means (pre-treatment)
  vector[K_cov] X_cov1;          // Treated covariate means (pre-treatment)
  vector<lower=0>[K_cov] v_cov;  // Per-covariate importance weights (V matrix diagonal)
}

transformed data {
  matrix[J, J] I_J = diag_matrix(rep_vector(1.0, J));

  // -------------------------------------------------------------
  // PREPARE OUTCOMES FOR SC MATCHING (pre-treatment, standardized)
  // -------------------------------------------------------------
  vector[T0] y_pre = to_vector(Y_panel[J + 1, ]);
  real mean_y = mean(y_pre);
  real sd_y   = sd(y_pre);
  vector[T0] y_pre_std = (y_pre - mean_y) / sd_y;

  // Donor standardization logic aligned with model1.stan
  matrix[T0, J] X_pre = Y_panel[1:J, ]';
  matrix[T0, J] X_pre_std;
  matrix[T_post, J] Y0_post_std;
  vector[J] mean_x;
  vector[J] sd_x;

  for (j in 1:J) {
    mean_x[j] = mean(X_pre[, j]);
    sd_x[j]   = sd(X_pre[, j]);
    X_pre_std[, j] = (X_pre[, j] - mean_x[j]) / sd_x[j];
    Y0_post_std[, j] = (Y0_post[, j] - mean_x[j]) / sd_x[j];
  }

  // -------------------------------------------------------------
  // PREPARE COVARIATE-MATCHING ROWS (standardized)
  // For each covariate k, compute donor-side mean/SD across the J
  // donor values, then standardize both the donor row (X_cov0) and
  // the treated scalar (X_cov1) by the same statistics. Guard SDs
  // near zero so a constant covariate produces a zero residual
  // rather than a NaN.
  // -------------------------------------------------------------
  matrix[K_cov, J] X_cov0_std;
  vector[K_cov] X_cov1_std;
  for (k in 1:K_cov) {
    real m_k = mean(X_cov0[k, ]);
    real s_k = sd(X_cov0[k, ]);
    if (s_k > 1e-12) {
      X_cov0_std[k, ] = (X_cov0[k, ] - m_k) / s_k;
      X_cov1_std[k]   = (X_cov1[k]    - m_k) / s_k;
    } else {
      X_cov0_std[k, ] = rep_row_vector(0.0, J);
      X_cov1_std[k]   = 0.0;
    }
  }

  // -------------------------------------------------------------
  // CONTAMINATION VECTOR s (always needed: penalty uses |s|)
  // -------------------------------------------------------------
  vector[J] s = rho * mdivide_left(I_J - rho * W_J, w_J1);
  vector[J] s_abs = fabs(s);

  // -------------------------------------------------------------
  // EFFECTIVE LIKELIHOOD SIZE FOR LAMBDA RESCALING
  //   n_eff = T0 + sum(v_cov)
  // T0 outcome rows each contribute Fisher information of weight 1 in w;
  // each covariate row k contributes v_cov[k] (since its precision is
  // v_cov[k] / sigma_sc^2). Rows with v_cov[k] = 0 contribute nothing,
  // consistent with their omission from the likelihood loop below.
  // -------------------------------------------------------------
  real n_eff = T0 + sum(v_cov);
}

parameters {
  simplex[J] w;                  // NASC donor weights
  real<lower=0> sigma_sc;        // Residual SD on standardized scale
  real<lower=0> lambda_tilde;    // NASC penalty strength on unit-information scale
}

transformed parameters {
  // Recover the raw penalty strength so downstream R code that extracts
  // "lambda" from the fit continues to work unchanged. Both lambda and
  // lambda_tilde are exposed in generated quantities for diagnostics.
  real<lower=0> lambda = lambda_tilde * n_eff;
}

model {
  // Priors
  sigma_sc     ~ normal(0, 1);
  lambda_tilde ~ cauchy(0, 1);     // Unit-information scale: O(1) is meaningful.

  // NASC Penalty (always on in this file). The penalty contribution is
  // -lambda_tilde * n_eff * <w, |s|>, so for lambda_tilde ~ O(1) the
  // penalty competes with the likelihood at O(1) per likelihood unit.
  target += -lambda * dot_product(w, s_abs);

  // Likelihood, augmented by covariate-matching rows when K_cov > 0.
  // Outcome rows: y_pre_std[t] ~ Normal(X_pre_std[t, ] * w, sigma_sc)
  target += normal_lpdf(y_pre_std | X_pre_std * w, sigma_sc);
  // Covariate rows: X_cov1_std[k] ~ Normal(X_cov0_std[k, ] * w,
  //                                        sigma_sc / sqrt(v_cov[k]))
  // Higher v_cov[k] => tighter implied SD => stronger matching pressure
  // on covariate k. v_cov[k] = 0 => infinite SD => row contributes 0
  // (covariate effectively excluded). The K_cov = 0 case skips the loop.
  for (k in 1:K_cov) {
    if (v_cov[k] > 0)
      target += normal_lpdf(X_cov1_std[k] |
                            dot_product(X_cov0_std[k, ], w),
                            sigma_sc / sqrt(v_cov[k]));
  }
}

generated quantities {
  vector[T0] y_sim_pre;             // Pre-treatment fitted values (raw scale)
  vector[T_post] y_counterfactual;  // Post-treatment counterfactual (raw scale)

  real wts_dot = dot_product(w, s);
  // bias_correction is 1.0 when toggled off, so downstream code applies it
  // unconditionally without changing tau.
  real bias_correction = use_bias_correction == 1 ? 1.0 / (1.0 - wts_dot) : 1.0;

  for (t in 1:T0) {
    y_sim_pre[t] = normal_rng(dot_product(X_pre_std[t, ], w), sigma_sc) * sd_y + mean_y;
  }

  for (t in 1:T_post) {
    y_counterfactual[t] = normal_rng(dot_product(Y0_post_std[t, ], w), sigma_sc) * sd_y + mean_y;
  }
}
