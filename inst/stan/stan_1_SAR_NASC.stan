// Panel Spatial Autoregressive (SAR) Model -- Within (two-way demeaned).
//
// Single regime: always applies the two-way within transform to Y, WY,
// and X. Unit + time fixed effects are absorbed by the projection.
//
// Behavior on time-invariant covariates:
//   The within transform maps any time-invariant column of X to zero.
//   The corresponding beta[k] therefore drops out of the likelihood and
//   its posterior reverts to the prior. This is harmless for rho and
//   for time-varying beta's; the calling R code is responsible for
//   flagging the affected beta's as not-identified (NA) post-hoc.
//
// Edge case K_pred = 0: no covariates supplied, model still identifies
// rho from the within-transformed Y alone. beta is a length-0 vector.

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
  real mean_y = mean(to_vector(Y_panel));
  real sd_y   = sd(to_vector(Y_panel));
  matrix[N, T0] Y_std  = (Y_panel - mean_y) / sd_y;
  matrix[N, T0] WY_std = W * Y_std;

  // ---- Assemble + standardize the covariate tensor ----
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
    // Guard sd: a fully-constant covariate (all units, all periods)
    // would give s = 0 and divide-by-zero. Keep the column at zero
    // in that case; the within transform would zero it out anyway.
    X_full_std[k] = s > 1e-12 ? (Xk - m) / s : rep_matrix(0.0, N, T0);
  }

  // ---- Two-way within transform applied to Y, WY, each X_k ----
  matrix[N, T0] Y_within;
  matrix[N, T0] WY_within;
  array[K_pred] matrix[N, T0] X_within;
  {
    vector[N]      row_means_Y  = Y_std * rep_vector(1.0 / T0, T0);
    row_vector[T0] col_means_Y  = rep_row_vector(1.0 / N, N) * Y_std;
    real           grand_Y      = mean(to_vector(Y_std));
    for (i in 1:N)
      for (t in 1:T0)
        Y_within[i, t] = Y_std[i, t] - row_means_Y[i] - col_means_Y[t] + grand_Y;

    vector[N]      row_means_WY = WY_std * rep_vector(1.0 / T0, T0);
    row_vector[T0] col_means_WY = rep_row_vector(1.0 / N, N) * WY_std;
    real           grand_WY     = mean(to_vector(WY_std));
    for (i in 1:N)
      for (t in 1:T0)
        WY_within[i, t] = WY_std[i, t] - row_means_WY[i] - col_means_WY[t] + grand_WY;

    for (k in 1:K_pred) {
      vector[N]      row_means_X = X_full_std[k] * rep_vector(1.0 / T0, T0);
      row_vector[T0] col_means_X = rep_row_vector(1.0 / N, N) * X_full_std[k];
      real           grand_X     = mean(to_vector(X_full_std[k]));
      for (i in 1:N)
        for (t in 1:T0)
          X_within[k, i, t] = X_full_std[k, i, t]
                              - row_means_X[i] - col_means_X[t] + grand_X;
    }
  }
}

parameters {
  real<lower=-1, upper=1> rho;
  vector[K_pred] beta;
  real<lower=0> sigma_sar;
}

model {
  // ---- Priors ----
  sigma_sar ~ normal(0, 1);
  beta      ~ normal(0, 1);

  // ---- Jacobian for Y -> (I - rho W) Y, applied each period ----
  real log_det_A = sum(log1m(rho * lambda_W));
  target += T0 * log_det_A;

  // ---- Within-SAR likelihood ----
  matrix[N, T0] mu_lp = rho * WY_within;
  for (k in 1:K_pred)
    mu_lp += beta[k] * X_within[k];

  target += normal_lpdf(to_vector(Y_within - mu_lp) | 0, sigma_sar);
}

generated quantities {
  real rho_out = rho;
}
