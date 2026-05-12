// Panel Spatial Autoregressive (SAR) Model -- Explicit hierarchical fixed
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
// Standardization is retained: Y is divided by sd_y and each X_k by
// sd_x_k. This keeps the N(0,1) priors on beta/sigma_sar/sigma_alpha/
// sigma_gamma meaningfully calibrated and tightens HMC step-size
// adaptation. beta and the FE SDs are back-transformed to the original
// scale in generated quantities.
//
// Identifiability of the fixed effects:
//   alpha_i + gamma_t is only identified up to an additive constant
//   (alpha_i -> alpha_i + c, gamma_t -> gamma_t - c leaves the
//   likelihood unchanged). With independent N(0, sigma_alpha) and
//   N(0, sigma_gamma) priors this would create a ridge through the
//   posterior. We pin it down with soft sum-to-zero constraints on
//   alpha and gamma. The constraints are weak enough that they do
//   not distort the posterior on rho, beta, or the FE SDs, but tight
//   enough to eliminate the ridge and divergences it would cause.
//
// Time-invariant covariates:
//   A column of X that is constant in time within every unit is
//   perfectly collinear with the unit fixed effects: any shift in
//   beta[k] can be absorbed by a compensating shift in alpha. We
//   detect these columns inside transformed data via the within-unit
//   variance and zero them out, so beta[k] is sampled purely from its
//   prior rather than wandering an unidentified ridge. R-side code is
//   responsible for flagging the affected beta's as NA post-hoc; this
//   guard exists so the sampler does not see pathology even if the R
//   flagging is bypassed.
//
// Edge case K_pred = 0: no covariates supplied, rho is still
// identified from Y, WY, and the fixed effects alone.

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

  // ---- Assemble + standardize the covariate tensor ----
  vector[K_pred] sd_x;
  array[K_pred] matrix[N, T0] X_std;
  // Identifiability flag per covariate. A covariate is time-varying
  // (and so identifiable in the presence of unit FE) iff its within-
  // unit variance is non-negligible. Columns that fail this test are
  // zeroed: beta[k] then reverts cleanly to its prior instead of
  // wandering the alpha+beta ridge. We compute the flag after
  // standardization so the threshold is on a dimensionless scale.
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

      // Within-unit variance: mean over units of the per-unit time
      // variance. Time-invariant covariates have within-unit variance
      // exactly zero (they are constant along the time dimension for
      // each unit).
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
        // Time-invariant in standardized form: collinear with alpha.
        // Zero the column so beta[k] is sampled from its prior only.
        X_std[k]        = rep_matrix(0.0, N, T0);
        x_identified[k] = 0;
      }
    } else {
      // Globally constant covariate: cannot be standardized.
      X_std[k]        = rep_matrix(0.0, N, T0);
      x_identified[k] = 0;
    }
  }
}

parameters {
  real<lower=-1, upper=1> rho;
  vector[K_pred] beta;
  vector[N] alpha_raw;        // unit fixed effects (non-centered)
  vector[T0] gamma_raw;       // time fixed effects (non-centered)
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_sar;
}

transformed parameters {
  // Non-centered parameterization for the fixed effects. Stan's HMC
  // mixes much faster on the raw-scale parameters when sigma_alpha or
  // sigma_gamma get small (the funnel problem); multiplying by the
  // scale here restores the original interpretation.
  vector[N]  alpha = sigma_alpha * alpha_raw;
  vector[T0] gamma = sigma_gamma * gamma_raw;
}

model {
  // ---- Priors ----
  // Coefficients and residual SD on the standardized scale: N(0,1)
  // means "an SD-unit change in x_k shifts y by an SD-unit of y", a
  // genuinely weakly-informative scale-free prior.
  beta        ~ normal(0, 1);
  sigma_sar   ~ normal(0, 1);

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
  // identified, but the split is not), which creates a flat ridge
  // through the posterior and produces divergences. The anchors are
  // tight enough to remove the ridge but several orders of magnitude
  // weaker than the hierarchical priors so they do not bias the
  // posterior on the parameters of interest.
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

  // ---- SAR likelihood with explicit fixed effects ----
  // Y_std = rho * WY_std + sum_k beta_k X_std_k + alpha_i + gamma_t + e
  matrix[N, T0] mu = rho * WY_std;
  for (k in 1:K_pred)
    mu += beta[k] * X_std[k];
  // Broadcast alpha across columns and gamma across rows.
  for (i in 1:N)
    for (t in 1:T0)
      mu[i, t] += alpha[i] + gamma[t];

  target += normal_lpdf(to_vector(Y_std - mu) | 0, sigma_sar);
}

generated quantities {
  real rho_out = rho;

  // ---- Back-transform beta to the original (un-standardized) scale ----
  // The model fits on Y_std = (Y - mean_y) / sd_y and X_std_k built
  // from (X_k - m_k) / sd_x_k (when identified), so
  //   beta_stan_k = beta_true_k * sd_x_k / sd_y.
  // Unidentified covariates (constant globally or time-invariant
  // within units) are returned as not-a-number so downstream code can
  // flag them.
  vector[K_pred] beta_orig;
  for (k in 1:K_pred) {
    if (x_identified[k] > 0.5)
      beta_orig[k] = beta[k] * sd_y / sd_x[k];
    else
      beta_orig[k] = not_a_number();
  }

  // ---- Back-transform FE SDs to the original scale of Y ----
  // alpha and gamma are on the standardized-Y scale; multiplying by
  // sd_y restores their interpretation as deviations on Y's natural
  // scale, which is what users will want to report.
  real sigma_alpha_orig = sigma_alpha * sd_y;
  real sigma_gamma_orig = sigma_gamma * sd_y;
  real sigma_sar_orig   = sigma_sar   * sd_y;
}
