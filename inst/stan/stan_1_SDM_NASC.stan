// Step 1: Pure Spatial Durbin Model (SDM) Panel Estimation
// Fully standardized internally using the global pre-treatment period
// Used to get an unbiased posterior distribution of 'rho'

data {
  int<lower=1> K_pred;
  vector[K_pred] X1;            // Treated unit covariates
  int<lower=0> J;
  matrix[K_pred, J] X0;         // Donor unit covariates

  int<lower=1> T0;
  matrix[J + 1, T0] Y_panel;    // Pre-treatment outcomes for all units
  matrix[J + 1, J + 1] W;       // Spatial weight matrix
  vector[J + 1] lambda_W;       // Eigenvalues of W
}

transformed data {
  // --- 1. Combine and Standardize Covariates ---
  matrix[J + 1, K_pred] X_full_raw;
  matrix[J + 1, K_pred] X_full_std;

  for (k in 1:K_pred) {
    X_full_raw[1:J, k] = to_vector(X0[k, ]);
    X_full_raw[J + 1, k] = X1[k];

    // Calculate grand mean and SD across all units
    real mean_k = mean(X_full_raw[, k]);
    real sd_k = sd(X_full_raw[, k]);

    // Mean-center and scale
    X_full_std[, k] = (X_full_raw[, k] - mean_k) / sd_k;
  }

  // --- 2. Standardize Outcomes ---
  matrix[J + 1, T0] Y_panel_std;
  real mean_y = mean(to_vector(Y_panel));
  real sd_y = sd(to_vector(Y_panel));

  for (j in 1:(J + 1)) {
    for (t in 1:T0) {
      Y_panel_std[j, t] = (Y_panel[j, t] - mean_y) / sd_y;
    }
  }

  // --- 3. Precompute Spatial Lags (Highly Optimized) ---
  // Apply W to the standardized Y and standardized X
  matrix[J + 1, T0] WY_panel = W * Y_panel_std;
  matrix[J + 1, K_pred] WX_full_std = W * X_full_std;
}

parameters {
  real<lower=-1, upper=1> rho;  // Global spatial autocorrelation
  vector[K_pred] beta;          // Direct effects
  vector[K_pred] theta;         // Indirect (spillover) effects
  real<lower=0> sigma_sdm;      // Error standard deviation

  // NOTE: 'alpha' is completely removed because all data has a mean of 0.
}

model {
  // --- Priors ---
  // Adjusted to normal(0, 1) because all inputs are strictly standardized
  sigma_sdm ~ normal(0, 1);
  beta ~ normal(0, 1);
  theta ~ normal(0, 1);

  // --- OPTIMIZED SDM Likelihood ---
  // Jacobian penalty to guarantee an unbiased rho
  real log_det_A = sum(log(1 - rho * lambda_W));
  target += T0 * log_det_A;

  // Isolate the errors (AY) and calculate the expected mean
  matrix[J + 1, T0] AY = Y_panel_std - rho * WY_panel;

  // Expected mean based on own traits (beta) and neighbors' traits (theta)
  vector[J + 1] sdm_mean = X_full_std * beta + WX_full_std * theta;

  // Evaluate likelihood
  target += normal_lpdf(to_vector(AY) | to_vector(rep_matrix(sdm_mean, T0)), sigma_sdm);
}
