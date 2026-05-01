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
  // CONTAMINATION VECTOR s (always needed: penalty uses |s|)
  // -------------------------------------------------------------
  vector[J] s = rho * mdivide_left(I_J - rho * W_J, w_J1);
  vector[J] s_abs = fabs(s);
}

parameters {
  simplex[J] w;                  // NASC donor weights
  real<lower=0> sigma_sc;        // Residual SD on standardized scale
  real<lower=0> lambda;          // NASC penalty strength (always active here)
}

model {
  // Priors
  sigma_sc ~ normal(0, 1);
  lambda   ~ gamma(1, 0.5);

  // NASC Penalty (always on in this file)
  target += -lambda * dot_product(w, s_abs);

  // Likelihood: standardized treated outcome ~ Normal(standardized donor mix, sigma_sc)
  target += normal_lpdf(y_pre_std | X_pre_std * w, sigma_sc);
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
