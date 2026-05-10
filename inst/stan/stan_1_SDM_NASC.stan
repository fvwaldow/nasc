// Panel Spatial Durbin Model (SDM) -- Within (two-way demeaned).
//
// Single regime: always applies the two-way within transform to Y, WY,
// X, and WX. Unit + time fixed effects are absorbed by the projection.
//
// Behavior on time-invariant covariates:
//   The within transform maps any time-invariant column of X to zero.
//   Since W does not introduce time variation, the corresponding column
//   of WX is also time-invariant and likewise zeroed. Both beta[k] and
//   theta[k] therefore drop out of the likelihood and their posteriors
//   revert to their priors. This is harmless for rho and for time-varying
//   coefficients; the calling R code is responsible for flagging the
//   affected beta[k] / theta[k] as not-identified (NA) post-hoc.
//
// Edge case K_pred = 0: no covariates supplied, model still identifies
// rho from the within-transformed Y alone. beta and theta are length-0.

data {
  int<lower=0> K_pred;
  int<lower=0> J;
  int<lower=1> T0;

  array[T0] vector[K_pred] X1;
  array[T0] matrix[K_pred, J] X0;

  matrix[J + 1, T0] Y_panel;
  matrix[J + 1, J + 1] W;
  vector[J + 1] lambda_W;
}

transformed data {
  int N = J + 1;

  // ---- Standardize Y by grand pre-treatment mean / SD ----
  // sd_y / sd_x are HOISTED out of any inner block so that `generated
  // quantities` can read them and back-transform beta and theta to the
  // original scale.
  real mean_y = mean(to_vector(Y_panel));
  real sd_y   = sd(to_vector(Y_panel));
  matrix[N, T0] Y_std  = (Y_panel - mean_y) / sd_y;
  matrix[N, T0] WY_std = W * Y_std;

  // ---- Assemble + standardize the covariate tensor X_full_std ----
  vector[K_pred] sd_x;                        // per-covariate grand SD, kept for back-transform
  array[K_pred] matrix[N, T0] X_full_std;
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
    X_full_std[k] = s > 1e-12 ? (Xk - m) / s : rep_matrix(0.0, N, T0);
  }

  // ---- Compute spatial lag WX in the standardized covariate space ----
  // (W applied along the unit dimension, period-by-period)
  array[K_pred] matrix[N, T0] WX_full_std;
  for (k in 1:K_pred)
    WX_full_std[k] = W * X_full_std[k];

  // ---- Two-way within transform on Y, WY, each X_k, each WX_k ----
  matrix[N, T0] Y_within;
  matrix[N, T0] WY_within;
  array[K_pred] matrix[N, T0] X_within;
  array[K_pred] matrix[N, T0] WX_within;
  {
    // Y
    vector[N]      row_means_Y  = Y_std * rep_vector(1.0 / T0, T0);
    row_vector[T0] col_means_Y  = rep_row_vector(1.0 / N, N) * Y_std;
    real           grand_Y      = mean(to_vector(Y_std));
    for (i in 1:N)
      for (t in 1:T0)
        Y_within[i, t] = Y_std[i, t] - row_means_Y[i] - col_means_Y[t] + grand_Y;

    // WY
    vector[N]      row_means_WY = WY_std * rep_vector(1.0 / T0, T0);
    row_vector[T0] col_means_WY = rep_row_vector(1.0 / N, N) * WY_std;
    real           grand_WY     = mean(to_vector(WY_std));
    for (i in 1:N)
      for (t in 1:T0)
        WY_within[i, t] = WY_std[i, t] - row_means_WY[i] - col_means_WY[t] + grand_WY;

    // X_k and WX_k
    for (k in 1:K_pred) {
      // X
      vector[N]      row_means_X = X_full_std[k] * rep_vector(1.0 / T0, T0);
      row_vector[T0] col_means_X = rep_row_vector(1.0 / N, N) * X_full_std[k];
      real           grand_X     = mean(to_vector(X_full_std[k]));
      for (i in 1:N)
        for (t in 1:T0)
          X_within[k, i, t] = X_full_std[k, i, t]
                              - row_means_X[i] - col_means_X[t] + grand_X;

      // WX
      vector[N]      row_means_WX = WX_full_std[k] * rep_vector(1.0 / T0, T0);
      row_vector[T0] col_means_WX = rep_row_vector(1.0 / N, N) * WX_full_std[k];
      real           grand_WX     = mean(to_vector(WX_full_std[k]));
      for (i in 1:N)
        for (t in 1:T0)
          WX_within[k, i, t] = WX_full_std[k, i, t]
                               - row_means_WX[i] - col_means_WX[t] + grand_WX;
    }
  }
}

parameters {
  real<lower=-1, upper=1> rho;
  vector[K_pred] beta;          // Direct effects of own X
  vector[K_pred] theta;         // Indirect (spillover) effects of WX
  real<lower=0> sigma_sdm;
}

model {
  // ---- Priors ----
  sigma_sdm ~ normal(0, 1);
  beta      ~ normal(0, 1);
  theta     ~ normal(0, 1);

  // ---- Jacobian for Y -> (I - rho W) Y, applied each period ----
  real log_det_A = sum(log1m(rho * lambda_W));
  target += T0 * log_det_A;

  // ---- Within-SDM likelihood ----
  // Q Y = rho * Q (W Y) + Q X beta + Q (W X) theta + Q e
  matrix[N, T0] mu_lp = rho * WY_within;
  for (k in 1:K_pred) {
    mu_lp += beta[k]  * X_within[k];
    mu_lp += theta[k] * WX_within[k];
  }

  target += normal_lpdf(to_vector(Y_within - mu_lp) | 0, sigma_sdm);
}

generated quantities {
  real rho_out = rho;

  // ---- Back-transform beta and theta to the original (un-standardized) scale ----
  // The model fits on Y_std = (Y - mean_y) / sd_y and X_std_k = (X_k - m_k) / sd_x_k,
  // with WX_full_std built from the already-standardized X (so WX_std also
  // carries the 1/sd_x_k factor). Both coefficients therefore satisfy
  //   param_stan_k = param_true_k * sd_x_k / sd_y,
  // and the back-transform is the same for beta and theta. Constant covariates
  // (sd_x = 0) cannot be identified; they are returned as not-a-number so
  // downstream code can flag them.
  vector[K_pred] beta_orig;
  vector[K_pred] theta_orig;
  for (k in 1:K_pred) {
    if (sd_x[k] > 1e-12) {
      beta_orig[k]  = beta[k]  * sd_y / sd_x[k];
      theta_orig[k] = theta[k] * sd_y / sd_x[k];
    } else {
      beta_orig[k]  = not_a_number();
      theta_orig[k] = not_a_number();
    }
  }
}
