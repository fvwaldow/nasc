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
   // standardize using pre-treatment values
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
   vector[N] y_sim_pre;
   vector[N_pred] y_counterfactual;

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
