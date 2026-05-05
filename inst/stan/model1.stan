// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Bayesian Synthetic control model with optional spatial bias correction.
//
// The likelihood is the standard SC simplex-weight model on standardized
// pre-treatment outcomes. The optional bias-correction block computes
//   s = rho * (I_J - rho * W_J)^{-1} * w_{J1}
//   bias_correction = 1 / (1 - <w, s>)
// in generated quantities, conditional on use_bias_correction. When that flag
// is 0, no spatial inputs are required (J_bc is set to 0) and bias_correction
// is the constant 1.0.
//
// Generated-quantity names (y_sim_pre, y_counterfactual) are harmonized with
// stan_2_NASC.stan so that one R post-processing path handles both files.
//
// PREDICTOR IMPORTANCE WEIGHTS (v_aug):
//   When R augments X with appended covariate rows (use_covariate_rows path),
//   the bottom (N - N_outcome) rows of X are pre-treatment unit means of
//   covariates, not outcomes. v_aug provides per-row importance weights for
//   those appended rows, mirroring the V matrix in classical Abadie-Diamond-
//   Hainmueller SCM. Each augmented row k is treated as
//     Normal(., sigma / sqrt(v_aug[k]))
//   so v_aug[k] = 1 (default) recovers equal weighting and v_aug[k] = 0
//   removes that covariate from the matching loss. Outcome rows always have
//   implicit weight 1. When N_outcome == N (no covariate augmentation), v_aug
//   is length 0 and has no effect.

data {
   int<lower=1> N;                    // number of pre-intervention periods (or N_outcome + N_cov rows)
   vector[N] y;                       // outcome
   int<lower=0> K;                    // number of donors
   matrix[N,K] X;                     // pre-intervention outcome for donors
   int<lower=1> N_pred;               // number of post-intervention periods
   matrix[N_pred,K] X_pred;           // post-intervention outcome for donors

   // ---- Optional spatial bias-correction inputs ----
   int<lower=0, upper=1> use_bias_correction;
   int<lower=0> J_bc;                 // = K when active, 0 when inactive
   matrix[J_bc, J_bc] W_J;            // donor-donor block of W (zero-sized if off)
   vector[J_bc] w_J1;                 // donor-to-treated column of W
   real<lower=-1, upper=1> rho_bc;    // rho used for bias correction

   // ---- Optional predictor-importance weights for augmented matching ----
   // N_outcome = N when no covariate rows are appended (default).
   // When N_outcome < N, rows N_outcome+1..N are appended covariate-mean
   // rows and v_aug provides their per-row V-matrix entries.
   int<lower=1, upper=N> N_outcome;
   vector<lower=0>[N - N_outcome] v_aug;
}

transformed data{ // normalize using pre-treatment values
   matrix[N, K] X_std;
   matrix[N_pred, K] X_pred_std;
   vector[K] mean_X;
   vector[K] sd_X;
   real mean_y = mean(y);
   real sd_y = sd(y);
   vector[N] y_std = (y - mean_y) / sd_y;

   for (k in 1:K) {
      mean_X[k] = mean(X[,k]);
      sd_X[k] = sd(X[,k]);
      X_std[,k] = (X[,k] - mean_X[k]) / sd_X[k];
      X_pred_std[,k] = (X_pred[,k] - mean_X[k]) / sd_X[k];
   }

   // Compute s only when bias correction is requested. When off, J_bc = 0 and
   // s_bc is a zero-length vector (no matrix solve performed).
   vector[J_bc] s_bc;
   if (use_bias_correction == 1) {
      matrix[J_bc, J_bc] I_J = diag_matrix(rep_vector(1.0, J_bc));
      s_bc = rho_bc * mdivide_left(I_J - rho_bc * W_J, w_J1);
   }
}

parameters {
   real<lower=0> sigma;
   simplex[K] w;
}

model {
   // Priors.
   sigma ~ normal(0,1);

   // Outcome rows: equal weight (implicit V = 1).
   target += normal_lpdf(y_std[1:N_outcome] | X_std[1:N_outcome, ] * w, sigma);

   // Augmented covariate-matching rows: per-row V-weighted precision.
   // sigma_k = sigma / sqrt(v_aug[k]); v_aug[k] = 0 -> row dropped from loss.
   for (k in 1:(N - N_outcome)) {
      if (v_aug[k] > 0) {
         int row_idx = N_outcome + k;
         target += normal_lpdf(y_std[row_idx] |
                               dot_product(X_std[row_idx, ], w),
                               sigma / sqrt(v_aug[k]));
      }
   }
}

generated quantities {
   // Names harmonized with stan_2_NASC.stan: y_sim_pre and y_counterfactual.
   vector[N] y_sim_pre;
   vector[N_pred] y_counterfactual;

   // bias_correction is 1.0 when toggled off, so downstream code can apply it
   // unconditionally without changing tau.
   real bias_correction = use_bias_correction == 1
                          ? 1.0 / (1.0 - dot_product(w, s_bc))
                          : 1.0;

   for (i in 1:N) {
      y_sim_pre[i] = normal_rng(X_std[i,]*w, sigma) * sd_y + mean_y;
   }
   for (j in 1:N_pred) {
      y_counterfactual[j] = normal_rng(X_pred_std[j,]*w, sigma) * sd_y + mean_y;
   }
}
