// Step 1: Pure Spatial Autoregressive (SAR) Panel Model
// Standardized internally using ALL units in the pre-treatment period

data {
  int<lower=1> K_pred;
  vector[K_pred] X1;            // Treated unit covariates
  int<lower=0> J;
  matrix[K_pred, J] X0;         // Donor unit covariates

  int<lower=1> T0;
  matrix[J + 1, T0] Y_panel;    // Pre-treatment outcomes for ALL units
  matrix[J + 1, J + 1] W;
  vector[J + 1] lambda_W;
}

transformed data {
  // 1. Combine raw covariates into a single matrix
  matrix[J + 1, K_pred] X_full_raw;
  for (k in 1:K_pred) {
    X_full_raw[1:J, k] = to_vector(X0[k, ]);
    X_full_raw[J + 1, k] = X1[k];
  }

  // 2. Standardize Covariates (using grand mean and SD across all units)
  matrix[J + 1, K_pred] X_full_std;
  for (k in 1:K_pred) {
    real mean_k = mean(X_full_raw[, k]);
    real sd_k = sd(X_full_raw[, k]);
    X_full_std[, k] = (X_full_raw[, k] - mean_k) / sd_k;
  }

  // 3. Standardize Outcomes (using grand pre-treatment mean and SD)
  matrix[J + 1, T0] Y_panel_std;

  // Since Y_panel only contains T0, we can use the whole matrix!
  real mean_y = mean(to_vector(Y_panel));
  real sd_y = sd(to_vector(Y_panel));

  for (j in 1:(J + 1)) {
    for (t in 1:T0) {
      Y_panel_std[j, t] = (Y_panel[j, t] - mean_y) / sd_y;
    }
  }

  // 4. Precompute the spatial lag on the standardized Y
  matrix[J + 1, T0] WY_panel = W * Y_panel_std;
}

parameters {
  real<lower=-1, upper=1> rho;
  vector[K_pred] beta;
  real<lower=0> sigma_sar;
  // Intercept omitted because the data is globally centered
}

model {
  // Priors
  sigma_sar ~ normal(0, 1);
  beta ~ normal(0, 1);

  // SAR Likelihood
  real log_det_A = sum(log(1 - rho * lambda_W));
  target += T0 * log_det_A;

  matrix[J + 1, T0] AY = Y_panel_std - rho * WY_panel;
  vector[J + 1] sar_mean = X_full_std * beta;

  target += normal_lpdf(to_vector(AY) | to_vector(rep_matrix(sar_mean, T0)), sigma_sar);
}
