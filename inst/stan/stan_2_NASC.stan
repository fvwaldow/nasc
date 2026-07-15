// Module 2
data {
  int<lower=0> J;
  int<lower=1> T0;
  matrix[J + 1, T0] Y_panel;
  matrix[J, J] W_J;
  vector[J] w_J1;
  int<lower=1> T_post;
  matrix[T_post, J] Y0_post;
  vector[T_post] Y1_post;
  int<lower=0> K_cov;
  matrix[K_cov, J] X_cov0;
  vector[K_cov] X_cov1;
  vector<lower=0>[K_cov] v_cov;
  real<lower=-1, upper=1> rho;
  int<lower=0, upper=1> use_bias_correction;
  int<lower=0, upper=1> use_penalty;
  real<lower=0> lambda;      // CV-selected (or user-fixed) penalty strength
  real<lower=0> sigma_ref;   // reference noise scale (data): calibrates the penalty
}
transformed data {
  matrix[J, J] I_J = diag_matrix(rep_vector(1.0, J));
  vector[T0] y_pre     = to_vector(Y_panel[J + 1, ]);
  real       mean_y    = mean(y_pre);
  real       sd_y      = sd(y_pre);
  vector[T0] y_pre_std = (y_pre - mean_y) / sd_y;

  matrix[T0, J]     X_pre = Y_panel[1:J, ]';   // T0 x J
  matrix[T0, J]     X_pre_std;
  matrix[T_post, J] Y0_post_std;
  for (j in 1:J) {
    X_pre_std[, j]   = (X_pre[, j]   - mean_y) / sd_y;
    Y0_post_std[, j] = (Y0_post[, j] - mean_y) / sd_y;
  }

  matrix[K_cov, J] X_cov0_std;
  vector[K_cov]    X_cov1_std;
  for (k in 1:K_cov) {
    real m_k = mean(X_cov0[k, ]);
    real s_k = sd(X_cov0[k, ]);
    if (s_k > 1e-12) {
      X_cov0_std[k, ] = (X_cov0[k, ] - m_k) / s_k;
      X_cov1_std[k]   = (X_cov1[k]   - m_k) / s_k;
    } else {
      X_cov0_std[k, ] = rep_row_vector(0.0, J);
      X_cov1_std[k]   = 0.0;
    }
  }

  // Contamination vector and rescaled-penalty normalizer.
  vector[J] s     = rho * mdivide_left(I_J - rho * W_J, w_J1);
  vector[J] s_abs = fabs(s);
  real n_eff = T0 + sum(v_cov);

  real sigma_floor = 0.05;

  // Penalty strength, calibrated at the REFERENCE noise scale sigma_ref (data,
  // = the residual scale of the unpenalized simplex fit) rather than at the
  // sampled sigma_sc. lambda is fixed by CV upstream, so lambda_tilde depends on
  // DATA ONLY: it lives here in transformed data, is computed once as a plain
  // double, and never enters the autodiff tape or the per-draw output. (In
  // transformed parameters it would be promoted to a var with zero derivatives
  // and written out every draw.) With lambda fixed, no ETDir normalizing
  // constant is needed: the model is a valid penalized (Gibbs) posterior.
  //
  // Why not square(sigma_sc): that makes the penalty term sigma-dependent, so it
  // enters sigma's own conditional and inflates it --
  //   sigma^2 = (SSR + 2*lambda*n_eff*(w's)) / n_eff
  // -- which corrupts sigma as an estimate of pre-fit noise and widens every
  // predictive interval. With sigma_ref the penalty drops out of sigma's
  // conditional entirely (sigma^2 = SSR / n_eff, as in the unpenalized model)
  // while remaining calibrated in likelihood-precision units at the scale the
  // data actually exhibit. The CV-selected lambda is then exactly the lambda
  // deployed, since CV calibrates at the same reference scale.
  real lambda_tilde = (lambda * n_eff) / square(sigma_ref);
}
parameters {
  simplex[J]    w;             // SC donor weights
  real<lower=0> sigma_raw;     // free part of the likelihood SD
}
transformed parameters {
  // Soft floor
  real<lower=0> sigma_sc = sqrt(square(sigma_floor) + square(sigma_raw));
}
model {
  sigma_raw    ~ normal(0, 1);       // half-normal (implicit <lower=0>)

  if (use_penalty)
    target += -lambda_tilde * dot_product(w, s_abs); // nasc penalty

  target += normal_lpdf(y_pre_std | X_pre_std * w, sigma_sc); // outcome likelihood

  for (k in 1:K_cov) {
    if (v_cov[k] > 0)
      target += normal_lpdf(X_cov1_std[k] |
                            dot_product(X_cov0_std[k, ], w),
                            sigma_sc / sqrt(v_cov[k]));
  }
}
generated quantities {
  vector[T0]     y_sim_pre;
  vector[T_post] y_counterfactual;
  real wts_dot         = dot_product(w, s);
  real bias_correction = use_bias_correction == 1
                         ? 1.0 / (1.0 - wts_dot)
                         : 1.0;
  for (t in 1:T0)
    y_sim_pre[t] = normal_rng(dot_product(X_pre_std[t, ], w), sigma_sc)
                   * sd_y + mean_y;
  for (t in 1:T_post)
    y_counterfactual[t] = normal_rng(dot_product(Y0_post_std[t, ], w), sigma_sc)
                          * sd_y + mean_y;
  real lambda_tilde_out       = lambda_tilde;
  real lambda_out = lambda;
  real n_eff_out        = n_eff;
  real sigma_sc_out     = sigma_sc;   // monitor: should sit near sigma_floor
  real sigma_ref_out    = sigma_ref;
}
