// Model 1

data {
   int<lower=1> N;
   vector[N] y;
   int<lower=0> K;
   matrix[N,K] X;
   int<lower=1> N_pred;
   matrix[N_pred,K] X_pred;

   // optional spatial bias-correction
   int<lower=0, upper=1> use_bias_correction;
   int<lower=0> J_bc;
   matrix[J_bc, J_bc] W_J;
   vector[J_bc] w_J1;
   real<lower=-1, upper=1> rho_bc;

   // optional predictor-importance weights
   int<lower=1, upper=N> N_outcome;
   vector<lower=0>[N - N_outcome] v_aug;
}

transformed data{
   matrix[N, K] X_std;
   matrix[N_pred, K] X_pred_std;
   vector[N] y_std;

   // 1. Outcome standardization (bsynth individual-wise)
   real mean_y = mean(y[1:N_outcome]);
   real sd_y   = sd(y[1:N_outcome]);

   for (i in 1:N_outcome) {
      y_std[i] = (y[i] - mean_y) / sd_y;
   }

   for (k in 1:K) {
      real mean_X_k = mean(X[1:N_outcome, k]);
      real sd_X_k   = sd(X[1:N_outcome, k]);

      if (sd_X_k > 1e-12) {
         X_std[1:N_outcome, k] = (X[1:N_outcome, k] - mean_X_k) / sd_X_k;
         X_pred_std[, k]       = (X_pred[, k] - mean_X_k) / sd_X_k;
      } else {
         X_std[1:N_outcome, k] = rep_vector(0.0, N_outcome);
         X_pred_std[, k]       = rep_vector(0.0, N_pred);
      }
   }

   // 2. Covariate standardization (row-wise across donors, like Model 2)
   for (aug_idx in (N_outcome + 1):N) {
      real m_aug = mean(X[aug_idx, ]);
      real s_aug = sd(X[aug_idx, ]);

      if (s_aug > 1e-12) {
         X_std[aug_idx, ] = (X[aug_idx, ] - m_aug) / s_aug;
         y_std[aug_idx]   = (y[aug_idx] - m_aug) / s_aug;
      } else {
         X_std[aug_idx, ] = rep_row_vector(0.0, K);
         y_std[aug_idx]   = 0.0;
      }
   }

   // 3. Bias correction vector
   // Initialized to 0 to prevent uninitialized memory warnings
   vector[J_bc] s_bc = rep_vector(0.0, J_bc);
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
   sigma ~ normal(0,1);
   target += normal_lpdf(y_std[1:N_outcome] | X_std[1:N_outcome, ] * w, sigma);

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
   // FIX: Only loop up to N_outcome so covariates are not back-transformed
   vector[N_outcome] y_sim_pre;
   vector[N_pred] y_counterfactual;

   real bias_correction = use_bias_correction == 1
                          ? 1.0 / (1.0 - dot_product(w, s_bc))
                          : 1.0;

   for (i in 1:N_outcome) {
      y_sim_pre[i] = normal_rng(X_std[i,]*w, sigma) * sd_y + mean_y;
   }
   for (j in 1:N_pred) {
      y_counterfactual[j] = normal_rng(X_pred_std[j,]*w, sigma) * sd_y + mean_y;
   }
}
