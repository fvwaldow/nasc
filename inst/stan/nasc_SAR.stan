// Bayesian Network-Aware Synthetic Control (NASC) Model
// Spatial Autoregressive (SAR) Panel Model

data {
  // 1. Predictor Match Data
  int<lower=1> K_pred;           // Number of predictors
  vector[K_pred] X1;             // Treated unit predictors
  int<lower=0> J;                // Number of donors
  matrix[K_pred, J] X0;          // Donor units predictors
  vector<lower=0>[K_pred] vs;    // External importance for predictors

  // 2. SAR Panel Data (to estimate rho, beta)
  int<lower=1> T0;               // Pre-treatment time periods
  matrix[J + 1, T0] Y_panel;     // Full outcome matrix (Treated + Donors)
  matrix[J + 1, J + 1] W;        // Full row-standardized weight matrix
  vector[J + 1] lambda_W; // NEW: Pre-calculated eigenvalues passed from R

  // 3. Submatrices for NASC penalty (derived from W)
  matrix[J, J] W_J;              // Donor-to-donor submatrix
  vector[J] w_J1;                // Treated-to-donor exposure vector

  // 4. Post-intervention Data
  int<lower=1> T_post;
  matrix[T_post, J] Y0_post;
  vector[T_post] Y1_post;
}

transformed data {
  matrix[J, J] I_J = diag_matrix(rep_vector(1.0, J));
  matrix[J + 1, J + 1] I_full = diag_matrix(rep_vector(1.0, J + 1));

  // Standardize predictors for SC Matching
  matrix[K_pred, J] X0_std;
  vector[K_pred] X1_std;
  for (k in 1:K_pred) {
    real sd_k = sd(X0[k, ]);
    X0_std[k, ] = X0[k, ] / sd_k;
    X1_std[k] = X1[k] / sd_k;
  }

  // --- Prepare Covariates for SAR ---
  // Combine treated and donor covariates into one (J+1 x K) matrix
  matrix[J + 1, K_pred] X_full_std;
  for (k in 1:K_pred) {
    X_full_std[1:J, k] = to_vector(X0_std[k, ]); // Donors 1 to J
    X_full_std[J + 1, k] = X1_std[k];            // Treated is J+1
  }

  // --- NEW OPTIMIZATION: Precompute the spatial lag of Y ---
  matrix[J + 1, T0] WY_panel = W * Y_panel;

  real epsilon = 1e-8;
}

parameters {
  // --- SAR Parameters ---
  real<lower=-1, upper=1> rho;   // Spatial lag of Y
  vector[K_pred] beta;           // Covariate coefficients
  vector[J + 1] mu;              // Unit fixed effects
  real<lower=0> sigma_sar;       // Noise in spatial process

  // --- SC Parameters ---
  simplex[J] w;                  // NASC weights
  simplex[K_pred] gamma_imp;     // Predictor importance
  real<lower=0> sigma_sc;        // Noise in predictor match
  real<lower=0> lambda;          // NASC penalty strength
}

transformed parameters {
  // s depends on the ESTIMATED rho
  vector[J] s = rho * mdivide_left(I_J - rho * W_J, w_J1);

  vector[K_pred] Omega;
  for (k in 1:K_pred) {
    Omega[k] = sigma_sc / gamma_imp[k];
  }
}

model {
  // --- Priors ---
  rho ~ normal(0, 0.5);
  mu ~ normal(0, 10);
  sigma_sar ~ normal(0, 5);
  beta ~ normal(0, 5);

  sigma_sc ~ normal(0, 5);
  lambda ~ normal(0, 10);
  gamma_imp ~ dirichlet(vs);

  // --- 1. OPTIMIZED SAR Likelihood ---
  {
    // Fast Log-Determinant using Eigenvalues (No matrix inversion!)
    real log_det_A = sum(log(1 - rho * lambda_W));
    target += T0 * log_det_A;

    // Evaluate (I - rho*W)Y efficiently using the precomputed WY_panel
    matrix[J + 1, T0] AY = Y_panel - rho * WY_panel;

    // Calculate the expected mean for all units
    vector[J + 1] sar_mean = mu + X_full_std * beta;

    // Fully Vectorized Likelihood:
    // rep_matrix copies the mean vector T0 times so it matches the dimensions of AY.
    // to_vector() flattens matrices into 1D vectors for lightning-fast computation.
    target += normal_lpdf(to_vector(AY) | to_vector(rep_matrix(sar_mean, T0)), sigma_sar);
  }

  // --- 2. NASC Regularization ---
  target += -lambda * dot_product(w, sqrt(square(s) + epsilon));

  // --- 3. SC Likelihood (Predictor Match) ---
  target += normal_lpdf(X1_std | X0_std * w, Omega);
}

generated quantities {
  vector[T_post] tau_nasc;
  real bias_correction = 1.0 / (1.0 - dot_product(w, s));

  for (t in 1:T_post) {
    real raw_effect = Y1_post[t] - dot_product(Y0_post[t, ], w);
    tau_nasc[t] = raw_effect * bias_correction;
  }
}
