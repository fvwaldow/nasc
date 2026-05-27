// Step 2: Network-Aware Synthetic Control (NASC) — unified model.

data {
  // ---- Core panel dimensions ----
  int<lower=0> J;                    // number of donor units
  int<lower=1> T0;                   // number of pre-treatment periods

  // Pre-treatment panel: rows 1..J are donors, row J+1 is the treated unit.
  matrix[J + 1, T0] Y_panel;

  // Spatial weight sub-matrices (always required; used for bias correction
  // and/or the NASC penalty).
  matrix[J, J] W_J;                  // donor-donor block of W
  vector[J] w_J1;                    // donor-to-treated column of W

  // Post-treatment data (used only in generated quantities).
  int<lower=1> T_post;
  matrix[T_post, J] Y0_post;         // post-treatment donor outcomes (raw scale)
  vector[T_post] Y1_post;            // post-treatment treated outcome (kept for R compatibility)

  // rho: one draw from Stage 1, or an exogenous value.
  real<lower=-1, upper=1> rho;

  // ---- Model configuration flags ----
  int<lower=0, upper=1> use_bias_correction;  // compute bias_correction in GQ
  int<lower=0, upper=1> use_penalty;          // add NASC penalty to target

  // ---- Optional covariate-matching inputs (K_cov = 0 disables) ----
  int<lower=0> K_cov;
  matrix[K_cov, J] X_cov0;           // donor covariate means (pre-treatment)
  vector[K_cov] X_cov1;              // treated covariate means (pre-treatment)
  vector<lower=0>[K_cov] v_cov;      // per-covariate importance weights (V-matrix diagonal)
}

transformed data {
  matrix[J, J] I_J = diag_matrix(rep_vector(1.0, J));

  // ------------------------------------------------------------------
  // OUTCOME STANDARDIZATION
  // Treated unit: standardize by the treated unit's own pre-treatment
  // mean/SD (T0 observations only, no covariate rows mixed in).
  // ------------------------------------------------------------------
  vector[T0] y_pre    = to_vector(Y_panel[J + 1, ]);
  real        mean_y  = mean(y_pre);
  real        sd_y    = sd(y_pre);
  vector[T0]  y_pre_std = (y_pre - mean_y) / sd_y;

  // Donors: each donor column standardized by its own pre-treatment
  // outcome mean/SD (T0 observations only).
  matrix[T0, J] X_pre     = Y_panel[1:J, ]';   // T0 x J
  matrix[T0, J] X_pre_std;
  matrix[T_post, J] Y0_post_std;
  vector[J] mean_x;
  vector[J] sd_x;

  for (j in 1:J) {
    mean_x[j] = mean(X_pre[, j]);
    sd_x[j]   = sd(X_pre[, j]);
    X_pre_std[, j]   = (X_pre[, j]   - mean_x[j]) / sd_x[j];
    Y0_post_std[, j] = (Y0_post[, j] - mean_x[j]) / sd_x[j];
  }

  // ------------------------------------------------------------------
  // COVARIATE STANDARDIZATION
  // Each covariate row is standardized by the donor-side mean/SD of
  // that covariate (across the J donors).  The treated value is
  // standardized by the same donor-derived statistics so that the
  // residual X_cov1_std[k] - dot(X_cov0_std[k,], w) is dimensionless
  // and comparable to the outcome residuals.  A near-zero SD guard
  // sets the row to zero (contributes nothing to the matching loss).
  // ------------------------------------------------------------------
  matrix[K_cov, J] X_cov0_std;
  vector[K_cov]    X_cov1_std;
  for (k in 1:K_cov) {
    real m_k = mean(X_cov0[k, ]);
    real s_k = sd(X_cov0[k, ]);
    if (s_k > 1e-12) {
      X_cov0_std[k, ] = (X_cov0[k, ] - m_k) / s_k;
      X_cov1_std[k]   = (X_cov1[k]   - m_k) / s_k;
    } else {
      X_cov0_std[k, ] = rep_row_vector(0.0, J);
      X_cov1_std[k]   = 0.0;
    }
  }

  // ------------------------------------------------------------------
  // CONTAMINATION VECTOR s (needed for bias correction and/or penalty).
  // s = rho * (I_J - rho * W_J)^{-1} * w_{J1}.
  // When rho = 0, s = 0 and the penalty has no effect, which is
  // correct: there is no contamination without spatial spillover.
  // ------------------------------------------------------------------
  vector[J] s     = rho * mdivide_left(I_J - rho * W_J, w_J1);
  vector[J] s_abs = fabs(s);

  // ------------------------------------------------------------------
  // EFFECTIVE LIKELIHOOD SIZE FOR LAMBDA RESCALING
  //   n_eff = T0 + sum(v_cov)
  // T0 outcome rows each contribute Fisher information weight 1 in w;
  // each covariate row k contributes v_cov[k] (precision is
  // v_cov[k] / sigma_sc^2).  Rows with v_cov[k] = 0 contribute nothing.
  // ------------------------------------------------------------------
  real n_eff = T0 + sum(v_cov);
}

parameters {
  simplex[J] w;                  // SC donor weights
  real<lower=0> sigma_sc;        // residual SD on the standardized scale

  // Penalty strength on the unit-information scale.
  // When use_penalty = 0 this parameter is sampled from its prior but
  // does not influence w or sigma_sc; downstream code may ignore it.
  real<lower=0> lambda_tilde;
}

transformed parameters {
  // Recover raw penalty strength so R code that extracts "lambda"
  // continues to work unchanged for the penalty case.
  real<lower=0> lambda = lambda_tilde * n_eff;
}

model {
  // ---- Priors ----
  sigma_sc     ~ normal(0, 1);
  lambda_tilde ~ gamma(2, 1);    // unit-information scale: O(1) is meaningful

  // ---- NASC penalty (active only when use_penalty = 1) ----
  // -lambda * <w, |s|> penalizes donor weights that lean heavily on
  // units whose outcomes are contaminated by spatial spillover from
  // the treated unit.  When use_penalty = 0 this block is skipped and
  // lambda_tilde is irrelevant to the posterior on w.
  if (use_penalty)
    target += -lambda * dot_product(w, s_abs);

  // ---- Outcome likelihood ----
  // y_pre_std[t] ~ Normal(X_pre_std[t, ] * w, sigma_sc)
  target += normal_lpdf(y_pre_std | X_pre_std * w, sigma_sc);

  // ---- Covariate-matching rows (K_cov = 0 skips the loop) ----
  // X_cov1_std[k] ~ Normal(X_cov0_std[k, ] * w, sigma_sc / sqrt(v_cov[k]))
  // Higher v_cov[k] => tighter implied SD => stronger matching pressure.
  // v_cov[k] = 0 => row excluded from the likelihood.
  for (k in 1:K_cov) {
    if (v_cov[k] > 0)
      target += normal_lpdf(X_cov1_std[k] |
                            dot_product(X_cov0_std[k, ], w),
                            sigma_sc / sqrt(v_cov[k]));
  }
}

generated quantities {
  // Pre-treatment fitted values and post-treatment counterfactual,
  // both back-transformed to the original scale of Y using the treated
  // unit's pre-treatment mean/SD (not contaminated by covariate rows).
  vector[T0]     y_sim_pre;
  vector[T_post] y_counterfactual;

  // Bias correction factor: 1 / (1 - <w, s>) when active, else 1.
  // Applied unconditionally by downstream R so the formula for tau is
  // the same regardless of whether bias correction is on.
  real wts_dot         = dot_product(w, s);
  real bias_correction = use_bias_correction == 1
                         ? 1.0 / (1.0 - wts_dot)
                         : 1.0;

  for (t in 1:T0)
    y_sim_pre[t] = normal_rng(dot_product(X_pre_std[t, ], w), sigma_sc)
                   * sd_y + mean_y;

  for (t in 1:T_post)
    y_counterfactual[t] = normal_rng(dot_product(Y0_post_std[t, ], w), sigma_sc)
                          * sd_y + mean_y;

  // Penalty diagnostics (always present; meaningful only when use_penalty = 1).
  real lambda_out       = lambda;
  real lambda_tilde_out = lambda_tilde;
  real n_eff_out        = n_eff;
}
