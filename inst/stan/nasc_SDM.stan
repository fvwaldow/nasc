// Bayesian Network-Aware Synthetic Control (NASC) Model
// Spatial Durbin Model (SDM) Estimation

data {
  // 1. Predictor Match Data
  int<lower=1> K_pred;           // Number of predictors
  vector[K_pred] X1;             // Treated unit predictors
  int<lower=0> J;                // Number of donors
  matrix[K_pred, J] X0;          // Donor units predictors
  vector<lower=0>[K_pred] vs;    // External importance for predictors

  // 2. SDM / Panel Data (to estimate rho, beta, theta)
  int<lower=1> T0;               // Pre-treatment time periods
  matrix[J + 1, T0] Y_panel;     // Full outcome matrix (Treated + Donors)
  matrix[J + 1, J + 1] W;        // Full row-standardized weight matrix

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

  // --- NEW: Prepare Covariates for SDM ---
  // Combine treated and donor covariates into one (J+1 x K) matrix
  matrix[J + 1, K_pred] X_full_std;
  for (k in 1:K_pred) {
    X_full_std[1:J, k] = to_vector(X0_std[k, ]); // Donors 1 to J
    X_full_std[J + 1, k] = X1_std[k];            // Treated is J+1
  }

  // Pre-calculate the spatial lag of covariates (W * X)
  matrix[J + 1, K_pred] WX_full_std = W * X_full_std;
}

parameters {
  // --- SDM Parameters ---
  real<lower=-1, upper=1> rho;   // Spatial lag of Y
  vector[K_pred] beta;           // Covariate coefficients (own)
  vector[K_pred] theta;          // Spatial lag of X coefficients (neighbors)
  vector[J + 1] mu;              // Unit fixed effects
  real<lower=0> sigma_sdm;       // Noise in spatial process

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
  rho ~ uniform(-1, 1);
  mu ~ normal(0, 10);
  sigma_sdm ~ cauchy(0, 5);
  beta ~ normal(0, 5);           // Prior for X coefficients
  theta ~ normal(0, 5);          // Prior for WX coefficients

  sigma_sc ~ cauchy(0, 5);
  lambda ~ cauchy(0, 10);
  gamma_imp ~ dirichlet(vs);

  // --- 1. SDM Likelihood (Estimates rho, beta, theta) ---
  {
    matrix[J + 1, J + 1] A = I_full - rho * W;
    real log_det_A = log_determinant(A);

    // Calculate the SDM expected mean: mu + X*beta + WX*theta
    vector[J + 1] sdm_mean = mu + X_full_std * beta + WX_full_std * theta;

    for (t in 1:T0) {
      target += log_det_A; // Jacobian adjustment
      // SDM: (I - rho*W)Y_t = mu + X*beta + WX*theta + epsilon
      target += normal_lpdf(A * Y_panel[, t] | sdm_mean, sigma_sdm);
    }
  }

  // --- 2. NASC Regularization ---
  target += -lambda * dot_product(w, fabs(s));

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
