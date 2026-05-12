// Panel Spatial Durbin Model (SDM) -- Explicit hierarchical fixed
// effects (Option B).
//
// Replaces the two-way within-demean specification, which suffered from
// Lee-Yu (2010) incidental-parameter bias on rho. Instead of projecting
// the unit/time fixed effects out (which produces a singular error
// covariance and an inconsistent QMLE), we sample alpha_i and gamma_t
// explicitly with hierarchical N(0, sigma_alpha) / N(0, sigma_gamma)
// priors. Partial pooling regularizes the high-dimensional FE problem
// for small T0, and the likelihood + Jacobian now live on the same
// (untransformed) scale so rho is unbiased.
//
// Structural form (on the standardized scale):
//   Y_std = rho * W Y_std + sum_k beta_k X_std_k + sum_k theta_k WX_std_k
//           + alpha_i + gamma_t + epsilon_{it}
// with epsilon_{it} ~ iid N(0, sigma_sdm^2).
//
// Standardization is retained: Y is divided by sd_y and each X_k by
// sd_x_k. WX_std_k is built from the already-standardized X, so it
// inherits the 1/sd_x_k factor, and beta_k / theta_k share the same
// back-transform. The N(0,1) priors on coefficients and SDs are then
// scale-free and HMC adaptation is well-behaved.
//
// Identifiability of the fixed effects:
//   alpha_i + gamma_t is only identified up to an additive constant.
//   With independent N(0, sigma_alpha) / N(0, sigma_gamma) priors this
//   would create a posterior ridge. We pin it down with soft sum-to-
//   zero anchors that are tight enough to eliminate the ridge but
//   orders of magnitude weaker than the hierarchical priors, so they
//   do not distort rho, beta, theta, or the FE SDs.
//
// Time-invariant covariates (SDM-specific note):
//   A column of X that is constant in time within every unit is
//   perfectly collinear with the unit fixed effects. Because W acts
//   only on the unit dimension, WX of a time-invariant X is *also*
//   time-invariant within units, so theta_k is unidentified for the
//   same reason beta_k is. We detect these columns via the within-
//   unit variance of X_std and zero out both X_std[k] and WX_std[k]
//   simultaneously. beta_k and theta_k then revert cleanly to their
//   priors instead of wandering an unidentified ridge.
//
// Note on SDM identifiability of beta vs theta:
//   Even when both X_k and WX_k are time-varying, the data may not
//   distinguish beta_k from theta_k well if W is close to a uniform-
//   averaging matrix or if X_k varies little across neighbors. This
//   is a feature of SDM, not a bug; the N(0,1) priors regularize
//   the weakly-identified direction and the posterior correctly
//   reflects the resulting uncertainty.
//
// Edge case K_pred = 0: rho is still identified from Y, WY, and the
// fixed effects alone.

data {
  int<lower=0> K_pred;
  int<lower=0> J;
  int<lower=1> T0;

  array[T0] vector[K_pred] X1;
  array[T0] matrix[K_pred, J] X0;

  matrix[J + 1, T0] Y_panel;
  matrix[J + 1, J + 1] W;

  // Eigenvalues of W, split into real and imaginary parts. The Jacobian
  // log|det(I - rho W)| uses |1 - rho * lambda_i| which for complex
  // lambda becomes sqrt((1 - rho Re)^2 + (rho Im)^2). For symmetric W
  // the imaginary parts are zero and this collapses to the real-only
  // sum log1m(rho * lambda).
  vector[J + 1] lambda_W_re;
  vector[J + 1] lambda_W_im;
}

transformed data {
  int N = J + 1;

  // ---- Standardize Y by grand pre-treatment mean / SD ----
  real mean_y = mean(to_vector(Y_panel));
  real sd_y   = sd(to_vector(Y_panel));
  matrix[N, T0] Y_std  = (Y_panel - mean_y) / sd_y;
  matrix[N, T0] WY_std = W * Y_std;

  // ---- Assemble + standardize the covariate tensor X_std ----
  vector[K_pred] sd_x;
  array[K_pred] matrix[N, T0] X_std;
  // Identifiability flag per covariate. A covariate is time-varying
  // (and so identifiable in the presence of unit FE) iff its within-
  // unit variance is non-negligible. When the flag is 0, both X_std[k]
  // and WX_std[k] are zeroed: beta[k] and theta[k] then revert cleanly
  // to their priors instead of wandering the alpha+beta+theta ridge.
  vector[K_pred] x_identified;
  for (k in 1:K_pred) {
    matrix[N, T0] Xk;
    for (t in 1:T0) {
      for (j in 1:J)
        Xk[j, t] = X0[t, k, j];
      Xk[N, t] = X1[t, k];
    }
    real m = mean(to_vector(Xk));
    real s = sd(to_vector(Xk));
    sd_x[k] = s;

    if (s > 1e-12) {
      matrix[N, T0] Xk_std = (Xk - m) / s;

      // Within-unit variance: mean over units of per-unit time
      // variance. Time-invariant covariates have within-unit variance
      // exactly zero (they are constant along the time dimension for
      // each unit). The threshold is on the dimensionless standardized
      // scale so it is independent of the user's units.
      real within_var = 0;
      for (i in 1:N) {
        real row_mean = mean(Xk_std[i, ]);
        within_var += variance(Xk_std[i, ] - row_mean);
      }
      within_var /= N;

      if (within_var > 1e-10) {
        X_std[k]        = Xk_std;
        x_identified[k] = 1;
      } else {
        X_std[k]        = rep_matrix(0.0, N, T0);
        x_identified[k] = 0;
      }
    } else {
      X_std[k]        = rep_matrix(0.0, N, T0);
      x_identified[k] = 0;
    }
  }

  // ---- Compute spatial lag WX in the standardized covariate space ----
  // Because W acts only on the unit (row) dimension, applying it to a
  // time-invariant standardized X column yields another time-invariant
  // column -- both beta[k] and theta[k] are then collinear with the
  // unit FE. The x_identified guard above already zeroed X_std[k] for
  // such k, so WX_std[k] = W * 0 = 0 and theta[k] reverts to its prior
  // alongside beta[k]. No separate WX-side guard is needed.
  array[K_pred] matrix[N, T0] WX_std;
  for (k in 1:K_pred)
    WX_std[k] = W * X_std[k];
}

parameters {
  real<lower=-1, upper=1> rho;
  vector[K_pred] beta;          // Direct effects of own X
  vector[K_pred] theta;         // Indirect (spillover) effects of WX
  vector[N] alpha_raw;          // unit fixed effects (non-centered)
  vector[T0] gamma_raw;         // time fixed effects (non-centered)
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_sdm;
}

transformed parameters {
  // Non-centered parameterization for the fixed effects. HMC mixes
  // much faster on the raw-scale parameters when sigma_alpha or
  // sigma_gamma get small (Neal's-funnel problem); multiplying by the
  // scale here restores the original interpretation.
  vector[N]  alpha = sigma_alpha * alpha_raw;
  vector[T0] gamma = sigma_gamma * gamma_raw;
}

model {
  // ---- Priors ----
  // Coefficients and residual SD on the standardized scale: N(0,1)
  // means "an SD-unit change in x_k shifts y by an SD-unit of y", a
  // weakly-informative scale-free prior. beta and theta share this
  // calibration since WX_std inherits the 1/sd_x_k factor from X_std.
  beta        ~ normal(0, 1);
  theta       ~ normal(0, 1);
  sigma_sdm   ~ normal(0, 1);

  // Half-normal hyperpriors on the FE SDs. Weak but proper; keep
  // sigma_alpha / sigma_gamma away from infinity.
  sigma_alpha ~ normal(0, 2);
  sigma_gamma ~ normal(0, 2);

  // Non-centered draws: alpha = sigma_alpha * alpha_raw is implied by
  // alpha_raw ~ N(0,1).
  alpha_raw ~ std_normal();
  gamma_raw ~ std_normal();

  // Soft sum-to-zero anchors. Without these, alpha and gamma are
  // jointly identified only up to a constant (alpha_i + gamma_t is
  // identified, but the split is not), creating a flat ridge through
  // the posterior and producing divergences. The anchors are tight
  // enough to remove the ridge but several orders of magnitude weaker
  // than the hierarchical priors so they do not bias the posterior on
  // rho, beta, theta, or the FE SDs.
  sum(alpha) ~ normal(0, 0.001 * N);
  sum(gamma) ~ normal(0, 0.001 * T0);

  // ---- Jacobian for Y -> (I - rho W) Y, applied each period ----
  // Now correct: the likelihood is on the (standardized) raw Y_std,
  // not a demeaned projection of it, so the full Jacobian over all
  // N eigenvalues and all T0 periods applies.
  real log_det_A = 0.5 * sum(log(
      square(1 - rho * lambda_W_re) + square(rho * lambda_W_im)
  ));
  target += T0 * log_det_A;

  // ---- SDM likelihood with explicit fixed effects ----
  // Y_std = rho * WY_std + sum_k beta_k X_std_k + sum_k theta_k WX_std_k
  //       + alpha_i + gamma_t + e
  matrix[N, T0] mu = rho * WY_std;
  for (k in 1:K_pred) {
    mu += beta[k]  * X_std[k];
    mu += theta[k] * WX_std[k];
  }
  // Broadcast alpha across columns and gamma across rows.
  for (i in 1:N)
    for (t in 1:T0)
      mu[i, t] += alpha[i] + gamma[t];

  target += normal_lpdf(to_vector(Y_std - mu) | 0, sigma_sdm);
}

generated quantities {
  real rho_out = rho;

  // ---- Back-transform beta and theta to the original scale ----
  // The model fits on Y_std = (Y - mean_y) / sd_y and X_std_k built
  // from (X_k - m_k) / sd_x_k (when identified), with WX_std_k built
  // from the already-standardized X (so WX_std_k also carries the
  // 1/sd_x_k factor). Both coefficients therefore satisfy
  //   param_stan_k = param_true_k * sd_x_k / sd_y,
  // and the back-transform is the same for beta and theta. Unidentified
  // covariates (constant globally or time-invariant within units) are
  // returned as not-a-number so downstream code can flag them.
  vector[K_pred] beta_orig;
  vector[K_pred] theta_orig;
  for (k in 1:K_pred) {
    if (x_identified[k] > 0.5) {
      beta_orig[k]  = beta[k]  * sd_y / sd_x[k];
      theta_orig[k] = theta[k] * sd_y / sd_x[k];
    } else {
      beta_orig[k]  = not_a_number();
      theta_orig[k] = not_a_number();
    }
  }

  // ---- Back-transform FE SDs and residual SD to the scale of Y ----
  // alpha and gamma are on the standardized-Y scale; multiplying by
  // sd_y restores their interpretation as deviations on Y's natural
  // scale.
  real sigma_alpha_orig = sigma_alpha * sd_y;
  real sigma_gamma_orig = sigma_gamma * sd_y;
  real sigma_sdm_orig   = sigma_sdm   * sd_y;
}
