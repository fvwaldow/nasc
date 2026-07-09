// Module 1 - SAR

data {
  int<lower=0> K_pred;
  int<lower=0> J;
  int<lower=1> T0;

  array[T0] vector[K_pred] X1;
  array[T0] matrix[K_pred, J] X0;

  matrix[J + 1, T0] Y_panel;
  matrix[J + 1, J + 1] W;

  vector[J + 1] lambda_W_re;
  vector[J + 1] lambda_W_im;
}

transformed data {
  int N = J + 1;

  // standardize using pre-treatment values
  real mean_y = mean(to_vector(Y_panel));
  real sd_y   = sd(to_vector(Y_panel));
  matrix[N, T0] Y_std  = (Y_panel - mean_y) / sd_y;
  matrix[N, T0] WY_std = W * Y_std;

  vector[K_pred] sd_x;
  array[K_pred] matrix[N, T0] X_std;

  vector[K_pred] x_identified;
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

    if (s > 1e-12) {
      matrix[N, T0] Xk_std = (Xk - m) / s;

      real within_var = 0;
      for (i in 1:N) {
        real row_mean = mean(Xk_std[i, ]);
        within_var += variance(Xk_std[i, ] - row_mean);
      }
      within_var /= N;

      if (within_var > 1e-10) {
        X_std[k]        = Xk_std;
        x_identified[k] = 1;
      } else {
        X_std[k]        = rep_matrix(0.0, N, T0);
        x_identified[k] = 0;
      }
    } else {
      X_std[k]        = rep_matrix(0.0, N, T0);
      x_identified[k] = 0;
    }
  }
}

parameters {
  real<lower=-1, upper=1> rho;
  vector[K_pred] beta;
  vector[N] alpha_raw;
  vector[T0] gamma_raw;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_gamma;
  real<lower=0> sigma_sar;
}

transformed parameters {
  vector[N]  alpha = sigma_alpha * alpha_raw;
  vector[T0] gamma = sigma_gamma * gamma_raw;
}

model {
  beta        ~ normal(0, 1);
  sigma_sar   ~ normal(0, 1);
  sigma_alpha ~ normal(0, 2);
  sigma_gamma ~ normal(0, 2);
  alpha_raw ~ std_normal();
  gamma_raw ~ std_normal();

  sum(alpha) ~ normal(0, 0.001 * N);
  sum(gamma) ~ normal(0, 0.001 * T0);

  real log_det_A = 0.5 * sum(log(
      square(1 - rho * lambda_W_re) + square(rho * lambda_W_im)
  ));
  target += T0 * log_det_A;

  matrix[N, T0] mu = rho * WY_std;
  for (k in 1:K_pred)
    mu += beta[k] * X_std[k];

  for (i in 1:N)
    for (t in 1:T0)
      mu[i, t] += alpha[i] + gamma[t];

  target += normal_lpdf(to_vector(Y_std - mu) | 0, sigma_sar);
}

generated quantities {
  real rho_out = rho;

  // re-standardize coefficients
  vector[K_pred] beta_orig;
  for (k in 1:K_pred) {
    if (x_identified[k] > 0.5)
      beta_orig[k] = beta[k] * sd_y / sd_x[k];
    else
      beta_orig[k] = not_a_number();
  }

  real sigma_alpha_orig = sigma_alpha * sd_y;
  real sigma_gamma_orig = sigma_gamma * sd_y;
  real sigma_sar_orig   = sigma_sar   * sd_y;
}
