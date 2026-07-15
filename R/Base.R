# Bayesian NASC estimator
# BASE code

nascSynth <- R6::R6Class(
  classname = "nascSynth",
  private = list(
    data = NULL,
    covariates = NULL,
    W = NULL,
    spatial_model = NULL,
    rho_exogenous = NULL,
    bias_correction = NULL,
    nasc_penalty = NULL,
    time = NULL,
    id = NULL,
    treated = NULL,
    outcome = NULL,
    ci_width = NULL,
    intervention = NULL,
    fitted = NULL,
    plot_data = NULL,
    stan_model = NULL,
    y_synth_draws = NULL,
    treated_ids = NULL,
    donor_ids   = NULL,
    uses_rho   = NULL,
    uses_step1 = NULL,
    beta_identified = NULL,
    cov_names = NULL,
    predictor_weights = NULL,
    predictors_op = NULL,
    special_predictors = NULL,
    time_predictors_prior = NULL,
    predictor_labels = NULL,
    rho_covariates = NULL,
    rho_time_window = NULL,
    lambda = NULL,        # user-fixed penalty strength (NULL => CV)
    lambda_cv_grid = NULL,      # candidate lambda values for CV
    lambda_train_frac = NULL,   # share of the pre-period used for training
    lambda_used = NULL,   # value actually used in Step 2
    lambda_cv_table = NULL      # CV results (grid x fold RMSE)
  ),
  active = list(
    plotData = function() { return(private$plot_data) },
    interventionTime = function() { return(private$intervention) }
  ),
  public = list(


    # nascSynth object.
    initialize = function(data, time, id, treated, outcome, ci_width = 0.75,
                          covariates = NULL,
                          W = NULL, spatial_model = "none",
                          rho = NULL,
                          bias_correction = NULL,
                          nasc_penalty = NULL,
                          predictor_weights = NULL,
                          predictors.op = "mean",
                          special.predictors = NULL,
                          time.predictors.prior = NULL,
                          rho.covariates = NULL,
                          rho.time.window = NULL,
                          lambda = NULL,
                          lambda_cv_grid = seq(from=0,to=100,by=1),
                          lambda_train_frac = 0.8) {

      stopifnot(ci_width > 0 & ci_width < 1)

      if (!spatial_model %in% c("none", "SAR", "SDM", "exogenous")) {
        stop("spatial_model must be 'none', 'SAR', 'SDM', or 'exogenous'.")
      }

      if (is.null(bias_correction)) {
        bias_correction <- spatial_model %in% c("SAR", "SDM", "exogenous")
      }
      if (is.null(nasc_penalty)) {
        nasc_penalty <- spatial_model %in% c("SAR", "SDM", "exogenous")
      }

      if (!is.logical(bias_correction) || length(bias_correction) != 1L) {
        stop("'bias_correction' must be a single logical (TRUE/FALSE).")
      }
      if (!is.logical(nasc_penalty) || length(nasc_penalty) != 1L) {
        stop("'nasc_penalty' must be a single logical (TRUE/FALSE).")
      }

      if (!is.null(rho)) {
        if (!is.numeric(rho) || length(rho) != 1L || is.na(rho)) {
          stop("'rho' must be a single numeric value.")
        }
        if (rho <= -1 || rho >= 1) {
          stop("'rho' must lie strictly inside (-1, 1).")
        }
      }

      if (spatial_model == "exogenous" && is.null(rho)) {
        stop("spatial_model = 'exogenous' requires a user-supplied 'rho'.")
      }

      uses_rho <- bias_correction || nasc_penalty

      if (uses_rho && is.null(W)) {
        stop(
          "A spatial weights matrix 'W' is required whenever ",
          "bias_correction = TRUE or nasc_penalty = TRUE."
        )
      }

      if (uses_rho && spatial_model == "none" && is.null(rho)) {
        stop(
          "spatial_model = 'none' with bias_correction = TRUE or ",
          "nasc_penalty = TRUE requires an exogenous 'rho'."
        )
      }

      if (!is.null(W)) {
        W_mat <- as.matrix(W)
        rsums <- rowSums(W_mat, na.rm = TRUE)
        is_row_standardized <- all(abs(rsums - 1) < 1e-6 | abs(rsums) < 1e-6)
        if (!is_row_standardized) {
          stop("The spatial weights matrix 'W' must be row-standardized (row sums must equal 1, or 0 for isolated units).")
        }
      }

      private$time            <- rlang::enquo(time)
      private$id              <- rlang::enquo(id)
      private$treated         <- rlang::enquo(treated)
      private$outcome         <- rlang::enquo(outcome)
      private$ci_width        <- ci_width
      private$covariates      <- covariates

      # predictors.op / special.predictors / time.predictors.prior (module 2)
      if (!is.character(predictors.op) || length(predictors.op) != 1L ||
          is.na(predictors.op) || !nzchar(predictors.op)) {
        stop("'predictors.op' must be a single non-empty character string ",
             "(e.g. \"mean\", \"median\").")
      }
      tryCatch(match.fun(predictors.op),
               error = function(e) {
                 stop("predictors.op = '", predictors.op,
                      "' is not a callable function.")
               })

      if (!is.null(special.predictors)) {
        if (!is.list(special.predictors)) {
          stop("'special.predictors' must be a list of length-3 lists.")
        }
        for (i in seq_along(special.predictors)) {
          e_i <- special.predictors[[i]]
          if (!is.list(e_i) || length(e_i) != 3L) {
            stop(sprintf(
              "special.predictors[[%d]] must be a list of length 3: ",
              i),
              "list(<predictor>, <time periods>, <operator>).")
          }
        }
      }

      if (!is.null(time.predictors.prior)) {
        if (!is.numeric(time.predictors.prior) ||
            length(time.predictors.prior) < 1L ||
            anyNA(time.predictors.prior)) {
          stop("'time.predictors.prior' must be a non-empty numeric vector ",
               "without NAs.")
        }
      }

      private$predictors_op         <- predictors.op
      private$special_predictors    <- special.predictors
      private$time_predictors_prior <- time.predictors.prior

      # rho.covariates / rho.time.window (module 1 only)
      if (!is.null(rho.covariates)) {
        if (!is.character(rho.covariates) || anyNA(rho.covariates) ||
            !all(nzchar(rho.covariates))) {
          stop("'rho.covariates' must be a character vector of covariate ",
               "column names (no NA / empty strings).")
        }
        if (is.null(covariates)) {
          stop("'rho.covariates' was supplied but 'covariates' is NULL; ",
               "there is no covariate panel to select from.")
        }
        avail_cov <- setdiff(
          names(covariates),
          c(rlang::as_name(private$time), rlang::as_name(private$id))
        )
        missing_cov <- setdiff(rho.covariates, avail_cov)
        if (length(missing_cov) > 0L) {
          stop("'rho.covariates' names not found in the covariate panel: ",
               paste(missing_cov, collapse = ", "), ".")
        }
        if (anyDuplicated(rho.covariates)) {
          stop("'rho.covariates' contains duplicate names.")
        }
      }

      if (!is.null(rho.time.window)) {
        if (!is.numeric(rho.time.window) || length(rho.time.window) < 1L ||
            anyNA(rho.time.window)) {
          stop("'rho.time.window' must be a non-empty numeric vector ",
               "without NAs.")
        }
        rho.time.window <- sort(unique(rho.time.window))
      }

      if ((!is.null(rho.covariates) || !is.null(rho.time.window)) &&
          !(bias_correction || nasc_penalty) ||
          (!is.null(rho.covariates) || !is.null(rho.time.window)) &&
          (!is.null(rho) || !spatial_model %in% c("SAR", "SDM"))) {
        message("nascSynth: 'rho.covariates'/'rho.time.window' are only used ",
                "when Step 1 estimates rho (spatial_model 'SAR'/'SDM' with no ",
                "exogenous rho); they will be ignored for this configuration.")
      }

      private$rho_covariates  <- rho.covariates
      private$rho_time_window <- rho.time.window

      if (!is.null(predictor_weights)) {
        if (is.null(covariates) && is.null(special.predictors)) {
          stop("predictor_weights was supplied but neither 'covariates' nor ",
               "'special.predictors' is set.")
        }
        if (!is.numeric(predictor_weights) ||
            any(!is.finite(predictor_weights)) ||
            any(predictor_weights < 0)) {
          stop("predictor_weights must be a finite, non-negative numeric vector.")
        }
      }
      private$predictor_weights <- predictor_weights

      private$W               <- W
      private$spatial_model   <- spatial_model
      private$rho_exogenous   <- rho

      # lambda / CV settings (module 2 penalty only)
      if (!is.null(lambda)) {
        if (!is.numeric(lambda) || length(lambda) != 1L ||
            is.na(lambda) || lambda < 0) {
          stop("'lambda' must be a single non-negative numeric value.")
        }
      }
      if (!is.numeric(lambda_cv_grid) || length(lambda_cv_grid) < 2L ||
          anyNA(lambda_cv_grid) || any(lambda_cv_grid < 0)) {
        stop("'lambda_cv_grid' must be a numeric vector (length >= 2) of ",
             "non-negative candidate values.")
      }
      if (!is.numeric(lambda_train_frac) || length(lambda_train_frac) != 1L ||
          is.na(lambda_train_frac) || lambda_train_frac <= 0 ||
          lambda_train_frac >= 1) {
        stop("'lambda_train_frac' must be a single number strictly between ",
             "0 and 1 (share of the pre-period used for training).")
      }
      private$lambda    <- lambda
      private$lambda_cv_grid    <- sort(unique(as.numeric(lambda_cv_grid)))
      private$lambda_train_frac <- as.numeric(lambda_train_frac)

      private$bias_correction <- bias_correction
      private$nasc_penalty    <- nasc_penalty
      private$uses_rho        <- uses_rho

      private$uses_step1 <- uses_rho &&
        spatial_model %in% c("SAR", "SDM") &&
        is.null(rho)

      rho_src <- if (!uses_rho) {
        "NA"
      } else if (!is.null(rho)) {
        "exogenous"
      } else {
        "posterior distribution"
      }
      message(sprintf(
        "nascSynth: bias correction = %s, nasc penalty = %s, rho source = %s",
        bias_correction, nasc_penalty, rho_src
      ))

      trt_vec <- data |> dplyr::pull({{ treated }})
      if (is.logical(trt_vec)) {
        trt_num <- as.integer(trt_vec)
      } else if (is.factor(trt_vec)) {
        lv <- levels(trt_vec)
        if (!all(lv %in% c("0", "1"))) {
          stop("Treated identifier is a factor with levels other than ",
               "'0'/'1': ", paste(lv, collapse = ", "), ".")
        }
        trt_num <- as.integer(as.character(trt_vec))
      } else if (is.numeric(trt_vec)) {
        trt_num <- trt_vec
      } else {
        stop("Treated identifier must be numeric, integer, logical, or a ",
             "0/1 factor; got class '",
             paste(class(trt_vec), collapse = "/"), "'.")
      }
      uniq_non_na <- unique(trt_num[!is.na(trt_num)])
      if (!all(uniq_non_na %in% c(0, 1))) {
        stop("Treated identifier contains non-binary values: ",
             paste(setdiff(uniq_non_na, c(0, 1)), collapse = ", "),
             ". Expected only 0 and 1 (NA allowed).")
      }
      if (!all(c(0, 1) %in% uniq_non_na)) {
        stop("Treated identifier must contain both 0 and 1; ",
             "found only: ", paste(uniq_non_na, collapse = ", "), ".")
      }
      data <- data |>
        dplyr::mutate(!!rlang::quo_name(rlang::enquo(treated)) := trt_num)

      private$data <- data |>
        dplyr::mutate(
          status = dplyr::case_when(
            {{ treated }} == 1 ~ "Treated",
            {{ treated }} == 0 ~ "Untreated",
            is.na({{ treated }}) ~ "N/A"
          ),
          status = factor(status, levels = c("Treated", "Untreated", "N/A"))
        )

      if (!(data |> dplyr::pull({{ id }}) |> is.factor())) {
        private$data <- private$data |>
          dplyr::mutate(!!rlang::quo_name(private$id) := factor({{ id }}))
      }

      private$treated_ids <- private$data |>
        dplyr::filter(!!private$treated == 1) |>
        dplyr::select({{ id }}) |>
        dplyr::distinct() |>
        dplyr::pull({{ id }})

      if (length(private$treated_ids) != 1) {
        stop("Expected exactly one treated unit; found ",
             length(private$treated_ids), ".")
      }

      private$intervention <- private$data |>
        dplyr::filter(status == "Treated") |>
        dplyr::summarise(
          !!rlang::as_label(private$time) := min(!!private$time)
        ) |>
        dplyr::pull(!!private$time)

      private$stan_model <- list(
        step1 = if (private$uses_step1) {
          switch(
            spatial_model,
            SAR = stanmodels$stan_1_SAR_NASC,
            SDM = stanmodels$stan_1_SDM_NASC
          )
        } else {
          NULL
        },
        step2 = stanmodels$stan_2_NASC
      )
    },

    # Fit the Stan model via HMC
    fit = function(n_samples = 100,
                   n_samples_cap = 500L,
                   n_samples_min = 30L,
                   cores = parallel::detectCores() - 1,
                   worker_iter = 2000L, worker_warmup = 1000L, ...) {

      auto_n_samples <- FALSE
      if (is.character(n_samples) && length(n_samples) == 1L &&
          identical(n_samples, "auto")) {
        auto_n_samples <- TRUE
      } else if (is.numeric(n_samples) && length(n_samples) == 1L &&
                 !is.na(n_samples) && n_samples >= 1) {
        n_samples <- as.integer(n_samples)
      } else {
        stop("'n_samples' must be a positive integer or the string \"auto\".")
      }
      if (!is.numeric(n_samples_cap) || length(n_samples_cap) != 1L ||
          n_samples_cap < 1) {
        stop("'n_samples_cap' must be a positive integer.")
      }
      if (!is.numeric(n_samples_min) || length(n_samples_min) != 1L ||
          n_samples_min < 1) {
        stop("'n_samples_min' must be a positive integer.")
      }
      if (n_samples_min > n_samples_cap) {
        stop("'n_samples_min' cannot exceed 'n_samples_cap'.")
      }
      n_samples_cap <- as.integer(n_samples_cap)
      n_samples_min <- as.integer(n_samples_min)

      if (auto_n_samples && !isTRUE(private$uses_step1)) {
        message("'n_samples = \"auto\"' has no effect when Step 1 does not ",
                "run (exogenous rho or no rho); ignoring.")
        n_samples <- n_samples_min
        auto_n_samples <- FALSE
      }
      wide_df <- .makeWide(
        data      = private$data,
        id        = private$id,
        time      = private$time,
        outcome   = private$outcome,
        treatment = private$treated
      )

      pre_data  <- wide_df |> dplyr::filter(!!private$time <  private$intervention)
      post_data <- wide_df |> dplyr::filter(!!private$time >= private$intervention)

      X      <- pre_data  |> dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)
      X1     <- pre_data  |> dplyr::pull(!!private$outcome)
      X_pred <- post_data |> dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)

      n_pre_real <- nrow(pre_data)

      # `donor_ids` is needed by .build_predictor_matrix(); set it before the
      # augmentation step so both the NASC and non-NASC branches use the same
      # column ordering.
      donor_ids  <- colnames(X_pred)
      treated_id <- as.character(private$treated_ids)
      private$donor_ids <- donor_ids

      use_covariate_rows <- !private$nasc_penalty &&
        (!is.null(private$covariates) ||
           !is.null(private$special_predictors))

      pred_labels_step2 <- character(0)
      v_pred_step2      <- numeric(0)

      if (use_covariate_rows) {
        pred_mat <- .build_predictor_matrix(
          data               = private$data,
          covariates         = private$covariates,
          id                 = private$id,
          time               = private$time,
          outcome            = private$outcome,
          treated_id         = treated_id,
          donor_ids          = donor_ids,
          intervention       = private$intervention,
          predictors_op      = private$predictors_op,
          special_predictors = private$special_predictors,
          time_pred_prior    = private$time_predictors_prior
        )

        if (length(pred_mat$names) > 0L) {
          X  <- rbind(as.data.frame(X),
                      as.data.frame(pred_mat$X0,
                                    check.names = FALSE,
                                    row.names = NULL))
          X1 <- c(X1, pred_mat$X1)
          pred_labels_step2 <- pred_mat$names
          private$predictor_labels <- pred_labels_step2
          v_pred_step2 <- .resolve_predictor_weights(
            private$predictor_weights, pred_labels_step2
          )
        }
      }

      W_J      <- NULL
      w_J1     <- NULL
      W_full   <- NULL
      lambda_W_re <- NULL
      lambda_W_im <- NULL
      Y_panel  <- NULL
      X0_arr   <- NULL
      X1_arr   <- NULL
      K_pred   <- 0L

      if (private$uses_rho) {
        if (is.null(rownames(private$W)) || is.null(colnames(private$W))) {
          all_ids <- levels(private$data[[rlang::as_name(private$id)]])
          rownames(private$W) <- colnames(private$W) <- all_ids
        }
        W_full <- private$W[c(donor_ids, treated_id), c(donor_ids, treated_id)]
        J <- length(donor_ids)
        W_J  <- W_full[1:J, 1:J]
        w_J1 <- as.vector(W_full[1:J, J + 1])

        if (private$uses_step1) {
          ev <- eigen(W_full, only.values = TRUE)$values
          lambda_W_re <- Re(ev)
          lambda_W_im <- Im(ev)

          # module 1 covariate selection
          cov_names <- if (is.null(private$covariates)) {
            character(0)
          } else {
            all_cov <- setdiff(
              names(private$covariates),
              c(rlang::as_name(private$time), rlang::as_name(private$id))
            )
            if (is.null(private$rho_covariates)) {
              all_cov
            } else {
              private$rho_covariates[private$rho_covariates %in% all_cov]
            }
          }

          # module 1 time window
          pre_times <- sort(unique(pre_data |> dplyr::pull(!!private$time)))
          if (is.null(private$rho_time_window)) {
            step1_times <- pre_times
          } else {
            req <- private$rho_time_window
            post_req <- req[req >= private$intervention]
            if (length(post_req) > 0L) {
              warning("'rho.time.window' includes post-intervention period(s): ",
                      paste(post_req, collapse = ", "),
                      "; these are dropped from Step-1 rho estimation.")
            }
            not_found <- setdiff(req[req < private$intervention], pre_times)
            if (length(not_found) > 0L) {
              warning("'rho.time.window' includes period(s) absent from the ",
                      "pre-treatment panel: ",
                      paste(not_found, collapse = ", "), "; they are ignored.")
            }
            step1_times <- sort(intersect(req, pre_times))
            if (length(step1_times) == 0L) {
              stop("'rho.time.window' selects no pre-treatment periods for ",
                   "Step-1 rho estimation.")
            }
          }

          T0_n       <- length(step1_times)              # module 1 periods
          unit_order <- c(donor_ids, treated_id)
          N_units    <- length(unit_order)
          J_n        <- length(donor_ids)
          treated_i  <- N_units                          # treated is last entry in unit_order

          if (!is.null(private$rho_covariates) ||
              !is.null(private$rho_time_window)) {
            message(sprintf(
              "Step 1: estimating rho from %d covariate(s) (%s) over %d pre-period(s).",
              length(cov_names),
              if (length(cov_names) > 0L) paste(cov_names, collapse = ", ") else "none",
              T0_n
            ))
          }

          if (length(cov_names) > 0L) {
            cov_pre <- private$covariates |>
              dplyr::filter(!!private$time %in% step1_times)

            expected_rows <- N_units * T0_n
            if (nrow(cov_pre) != expected_rows) {
              stop(sprintf(
                "Covariate panel is unbalanced: expected %d rows (%d units x %d periods), got %d.",
                expected_rows, N_units, T0_n, nrow(cov_pre)
              ))
            }

            X_kit <- array(
              NA_real_,
              dim = c(length(cov_names), N_units, T0_n),
              dimnames = list(cov_names, unit_order, NULL)
            )
            for (k_idx in seq_along(cov_names)) {
              v <- cov_names[k_idx]
              wide_kt <- cov_pre |>
                dplyr::select(!!private$id, !!private$time, dplyr::all_of(v)) |>
                tidyr::pivot_wider(
                  names_from  = !!private$id,
                  values_from = dplyr::all_of(v)
                ) |>
                dplyr::arrange(!!private$time)
              mat_tn <- as.matrix(wide_kt[, unit_order, drop = FALSE])  # T0 x N
              X_kit[k_idx, , ] <- t(mat_tn)                              # N x T0
            }
            if (anyNA(X_kit)) {
              stop("Missing values in the covariate panel after reshape; please impute or drop.")
            }

            within_var <- vapply(seq_along(cov_names), function(k) {
              m  <- X_kit[k, , , drop = TRUE]                # N x T0
              rm <- rowMeans(m)
              mean(rowSums((m - rm)^2) / max(T0_n - 1L, 1L))
            }, numeric(1))
            beta_identified <- within_var > 1e-10
            names(beta_identified) <- cov_names
            if (any(!beta_identified)) {
              message(sprintf(
                "%d time-invariant %s covariate(s) detected; posteriors are not identified",
                sum(!beta_identified),
                paste(cov_names[!beta_identified], collapse = ", ")
              ))
            }
            private$beta_identified <- beta_identified
            private$cov_names       <- cov_names
          } else {
            private$beta_identified <- logical(0)
            private$cov_names       <- character(0)
          }
          K_pred <- length(cov_names)

          if (K_pred == 0L) {
            X1_arr <- matrix(numeric(0), nrow = T0_n, ncol = 0)
            X0_arr <- array(numeric(0), dim = c(T0_n, 0, J_n))
          } else {
            X1_arr <- matrix(NA_real_, nrow = T0_n, ncol = K_pred)
            X0_arr <- array(NA_real_, dim = c(T0_n, K_pred, J_n))
            for (t in seq_len(T0_n)) {
              X1_arr[t, ]   <- X_kit[, treated_i, t]
              X0_arr[t, , ] <- X_kit[, seq_len(J_n), t]
            }
          }

          # module 1 outcome panel
          Y_panel_step1 <- t(as.matrix(
            pre_data |>
              dplyr::filter(!!private$time %in% step1_times) |>
              dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
          ))

          # module 2 outcome panel
          Y_panel <- t(as.matrix(
            pre_data |> dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
          ))
        } else if (private$nasc_penalty) {
          Y_panel <- t(as.matrix(
            pre_data |> dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
          ))
        }
      }

      if (is.null(Y_panel)) {
        Y_panel <- t(as.matrix(
          pre_data |> dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
        ))
      }

      # Module 1
      sampled_rhos <- if (private$uses_rho) {
        if (!is.null(private$rho_exogenous)) {
          private$rho_exogenous
        } else if (private$uses_step1) {

          step1_data <- list(
            K_pred   = K_pred,
            J        = length(donor_ids),
            T0       = T0_n,
            X1       = X1_arr,         # T0 x K_pred
            X0       = X0_arr,         # T0 x K_pred x J
            Y_panel  = Y_panel_step1,  # (J+1) x T0
            W        = W_full,
            lambda_W_re = lambda_W_re,
            lambda_W_im = lambda_W_im
          )

          message(sprintf("Running Step 1: Estimating %s parameters...",
                          private$spatial_model))
          private$fitted <- rstan::sampling(
            private$stan_model$step1,
            data  = step1_data,
            cores = cores,
            ...
          )

          # Stratified sampling across chains
          rho_array <- rstan::extract(private$fitted, pars = "rho", permuted = FALSE)
          n_iters  <- dim(rho_array)[1]
          n_chains <- dim(rho_array)[2]
          total_rho_draws <- n_iters * n_chains

          if (auto_n_samples) {
            n_eff_rho <- tryCatch({
              s <- rstan::summary(private$fitted, pars = "rho")$summary
              as.numeric(s[, "n_eff"])
            }, error = function(e) NA_real_)

            if (!is.finite(n_eff_rho) || n_eff_rho <= 0) {
              warning("Could not read n_eff for rho from the Step-1 fit; ",
                      "falling back to n_samples = ", n_samples_min, ".")
              n_samples <- n_samples_min
            } else {
              proposed  <- max(n_samples_min,
                               min(n_samples_cap, as.integer(round(n_eff_rho))))
              n_samples <- min(proposed, total_rho_draws)
              message(sprintf(
                "Auto n_samples: n_eff(rho) = %.1f -> n_samples = %d (cap %d, min %d).",
                n_eff_rho, n_samples, n_samples_cap, n_samples_min
              ))
            }
          } else {
            if (n_samples > total_rho_draws) {
              warning(sprintf(
                "n_samples (%d) exceeds total Step-1 draws (%d); reducing to %d.",
                n_samples, total_rho_draws, total_rho_draws
              ))
              n_samples <- total_rho_draws
            }
          }

          samples_per_chain <- ceiling(n_samples / n_chains)
          rho_draws <- numeric(0)
          for (c_idx in seq_len(n_chains)) {
            chain_draws <- rho_array[, c_idx, 1]
            take_n <- min(samples_per_chain, n_iters)
            rho_draws <- c(rho_draws,
                           sample(chain_draws, size = take_n, replace = FALSE))
          }
          if (length(rho_draws) > n_samples) {
            rho_draws <- sample(rho_draws, size = n_samples, replace = FALSE)
          }
          rho_draws
        } else {
          stop("Internal: uses_rho but no exogenous rho and no Step 1.")
        }
      } else {
        NA_real_
      }

      # Module 2
      if (private$nasc_penalty) {
        if (!is.null(private$covariates) ||
            !is.null(private$special_predictors)) {
          pred_mat_step2 <- .build_predictor_matrix(
            data               = private$data,
            covariates         = private$covariates,
            id                 = private$id,
            time               = private$time,
            outcome            = private$outcome,
            treated_id         = treated_id,
            donor_ids          = donor_ids,
            intervention       = private$intervention,
            predictors_op      = private$predictors_op,
            special_predictors = private$special_predictors,
            time_pred_prior    = private$time_predictors_prior
          )
          K_cov_step2 <- length(pred_mat_step2$names)
          X_cov0_mat  <- pred_mat_step2$X0
          X_cov1_vec  <- pred_mat_step2$X1
          private$predictor_labels <- pred_mat_step2$names
        } else {
          K_cov_step2 <- 0L
          X_cov0_mat  <- matrix(numeric(0), nrow = 0, ncol = length(donor_ids))
          X_cov1_vec  <- numeric(0)
        }

        if (K_cov_step2 > 0L) {
          v_cov_vec <- .resolve_predictor_weights(
            private$predictor_weights, private$predictor_labels
          )
        } else {
          v_cov_vec <- numeric(0)
        }

        base_data <- list(
          J                   = length(donor_ids),
          T0                  = nrow(pre_data),
          Y_panel             = Y_panel,
          W_J                 = W_J,
          w_J1                = w_J1,
          T_post              = nrow(post_data),
          Y0_post             = as.matrix(X_pred[, donor_ids, drop = FALSE]),
          Y1_post             = as.array(as.numeric(post_data |> dplyr::pull(!!private$outcome))),
          use_bias_correction = as.integer(private$bias_correction),
          use_penalty         = 1L,
          K_cov               = K_cov_step2,
          X_cov0              = X_cov0_mat,
          X_cov1              = if (K_cov_step2 > 0L) as.array(X_cov1_vec) else numeric(0),
          v_cov               = if (K_cov_step2 > 0L) as.array(v_cov_vec) else numeric(0)
        )

        # Reference noise scale: calibrates the penalty and keeps it out of
        # sigma's conditional. Computed once from the unpenalized fit and used
        # both by CV and by the final sampler, so the selected lambda is exactly
        # the lambda deployed.
        base_data$sigma_ref <- .compute_sigma_ref(base_data)

        # --- lambda: user-fixed or CV-selected -------------------------
        lt_used <- private$lambda
        if (is.null(lt_used)) {
          rho_cv <- if (length(sampled_rhos) == 1L && !is.na(sampled_rhos[1])) {
            sampled_rhos[1]
          } else {
            mean(sampled_rhos, na.rm = TRUE)
          }
          message(sprintf(
            "Selecting lambda by hold-out CV (%.0f%% train / %.0f%% validation just before treatment, rho = %.3f)...",
            100 * private$lambda_train_frac,
            100 * (1 - private$lambda_train_frac), rho_cv))
          cv_res <- .cv_lambda(
            base_data  = base_data,
            rho        = rho_cv,
            grid       = private$lambda_cv_grid,
            train_frac = private$lambda_train_frac
          )
          lt_used <- cv_res$lambda
          private$lambda_cv_table <- cv_res$table
          message(sprintf(
            "CV-selected lambda = %.4g (validation RMSE = %.4g; train = periods 1-%d, validation = %d-%d)",
            lt_used, cv_res$rmse, max(cv_res$train),
            min(cv_res$validation), max(cv_res$validation)))
        }
        message(sprintf("Penalty calibrated at sigma_ref = %.4g", base_data$sigma_ref))
        private$lambda_used <- lt_used
        base_data$lambda    <- lt_used

        results <- .run_step2_loop(
          rhos          = sampled_rhos,
          base_data     = base_data,
          step2_mod     = private$stan_model$step2,
          rho_field     = "rho",
          cores         = cores,
          extra_args    = list(...),
          extract_pars  = c("y_counterfactual", "y_sim_pre", "w", "lambda_tilde",
                            "sigma_sc", "bias_correction"),
          worker_iter   = worker_iter,
          worker_warmup = worker_warmup
        )

        if (is.null(private$fitted)) private$fitted <- results$last_fit

        private$y_synth_draws <- list(
          y_counterfactual = results$y_counterfactual,
          y_sim_pre        = results$y_sim_pre,
          w                = results$w,
          lambda_tilde           = results$lambda_tilde,
          lambda     = lt_used,
          sigma_sc         = results$sigma_sc,
          bias_correction  = results$bias_correction,
          rhos_used        = results$rhos_used,
          worker_diagnostics = results$worker_diagnostics
        )

      } else {
        J_n       <- length(donor_ids)
        N_total   <- nrow(X)
        N_outcome <- n_pre_real
        N_aug     <- N_total - N_outcome

        if (N_aug > 0L) {
          if (length(v_pred_step2) != N_aug) {
            stop(sprintf(
              "Internal: predictor-weight vector length (%d) does not match ",
              length(v_pred_step2)),
              sprintf("number of augmented matching rows (%d).", N_aug))
          }
          X_cov0_mat_nop <- as.matrix(X[(N_outcome + 1L):N_total, donor_ids,
                                        drop = FALSE])  # K_cov x J
          X_cov1_vec_nop <- as.numeric(X1[(N_outcome + 1L):length(X1)])
          v_cov_vec_nop  <- v_pred_step2
        } else {
          X_cov0_mat_nop <- matrix(numeric(0), nrow = 0L, ncol = J_n)
          X_cov1_vec_nop <- numeric(0)
          v_cov_vec_nop  <- numeric(0)
        }

        if (private$bias_correction) {
          W_J_nop  <- W_J
          w_J1_nop <- w_J1
        } else {
          W_J_nop  <- diag(J_n)
          w_J1_nop <- rep(0, J_n)
        }

        base_data <- list(
          J                   = J_n,
          T0                  = N_outcome,
          Y_panel             = Y_panel,
          W_J                 = W_J_nop,
          w_J1                = w_J1_nop,
          T_post              = nrow(X_pred),
          Y0_post             = as.matrix(X_pred[, donor_ids, drop = FALSE]),
          Y1_post             = as.array(as.numeric(
            post_data |> dplyr::pull(!!private$outcome))),
          use_bias_correction = as.integer(private$bias_correction),
          use_penalty         = 0L,
          lambda        = 1.0,  # required by Stan data block; unused when use_penalty = 0
          sigma_ref     = 1.0,  # ditto
          K_cov               = N_aug,
          X_cov0              = X_cov0_mat_nop,
          X_cov1              = if (N_aug > 0L) as.array(X_cov1_vec_nop) else numeric(0),
          v_cov               = if (N_aug > 0L) as.array(v_cov_vec_nop) else numeric(0)
        )

        if (!private$bias_correction) {
          base_data$rho <- 0

          private$fitted <- rstan::sampling(
            private$stan_model$step2,
            data  = base_data,
            cores = cores,
            ...
          )

          draws <- rstan::extract(
            private$fitted,
            pars = c("y_counterfactual", "y_sim_pre", "w", "sigma_sc",
                     "bias_correction")
          )
          private$y_synth_draws <- list(
            y_counterfactual = draws$y_counterfactual,
            y_sim_pre        = draws$y_sim_pre,
            w                = draws$w,
            sigma_sc         = draws$sigma_sc,
            bias_correction  = draws$bias_correction,
            rhos_used        = NA_real_
          )

        } else {
          results <- .run_step2_loop(
            rhos          = sampled_rhos,
            base_data     = base_data,
            step2_mod     = private$stan_model$step2,
            rho_field     = "rho",
            cores         = cores,
            extra_args    = list(...),
            extract_pars  = c("y_counterfactual", "y_sim_pre", "w", "sigma_sc",
                              "bias_correction"),
            worker_iter   = worker_iter,
            worker_warmup = worker_warmup
          )

          if (is.null(private$fitted)) private$fitted <- results$last_fit

          private$y_synth_draws <- list(
            y_counterfactual   = results$y_counterfactual,
            y_sim_pre          = results$y_sim_pre,
            w                  = results$w,
            sigma_sc           = results$sigma_sc,
            bias_correction    = results$bias_correction,
            rhos_used          = results$rhos_used,
            worker_diagnostics = results$worker_diagnostics
          )
        }
      }

      private$plot_data <- .get_nasc_results(
        y_counterfactual_draws = private$y_synth_draws$y_counterfactual,
        bias_correction_draws  = private$y_synth_draws$bias_correction,
        y_sim_pre_draws        = private$y_synth_draws$y_sim_pre,
        pre_data               = pre_data,
        post_data              = post_data,
        time                   = private$time,
        outcome                = private$outcome,
        ci                     = private$ci_width
      )
    },

    # Update credible interval width
    updateWidth = function(ci_width = 0.75) {
      stopifnot(ci_width > 0, ci_width < 1)
      private$ci_width <- ci_width
      wide_df <- .makeWide(
        data      = private$data,
        id        = private$id,
        time      = private$time,
        outcome   = private$outcome,
        treatment = private$treated
      )
      pre_data  <- wide_df |> dplyr::filter(!!private$time <  private$intervention)
      post_data <- wide_df |> dplyr::filter(!!private$time >= private$intervention)

      private$plot_data <- .get_nasc_results(
        y_counterfactual_draws = private$y_synth_draws$y_counterfactual,
        bias_correction_draws  = private$y_synth_draws$bias_correction,
        y_sim_pre_draws        = private$y_synth_draws$y_sim_pre,
        pre_data               = pre_data,
        post_data              = post_data,
        time                   = private$time,
        outcome                = private$outcome,
        ci                     = private$ci_width
      )
    },

    # fit summary
    summary = function(ci_width = NULL, print = TRUE) {
      if (is.null(private$y_synth_draws)) {
        stop("Run $fit() before calling summary().")
      }
      ci <- if (is.null(ci_width)) private$ci_width else {
        stopifnot(is.numeric(ci_width), length(ci_width) == 1L,
                  ci_width > 0, ci_width < 1)
        ci_width
      }

      rho_source <- if (!private$uses_rho) {
        "n/a (no rho in use)"
      } else if (!is.null(private$rho_exogenous)) {
        sprintf("exogenous (rho = %.4f)", private$rho_exogenous)
      } else {
        sprintf("Step 1 %s posterior", private$spatial_model)
      }

      out <- .nasc_summary_stats(list(
        y_synth_draws   = private$y_synth_draws,
        plot_data       = private$plot_data,
        intervention    = private$intervention,
        time            = private$time,
        outcome         = private$outcome,
        ci_width        = ci,
        treated_ids     = private$treated_ids,
        donor_ids       = private$donor_ids,
        id_levels       = private$data[[rlang::as_name(private$id)]],
        spatial_model   = private$spatial_model,
        bias_correction = private$bias_correction,
        nasc_penalty    = private$nasc_penalty,
        uses_rho        = private$uses_rho,
        rho_source      = rho_source,
        fitted          = private$fitted,
        cov_names       = private$cov_names,
        beta_identified = private$beta_identified,
        model           = self
      ))

      if (isTRUE(print)) {
        print.summary.nascSynth(out)
      }
      invisible(out)
    },

    # summary of the indirect treatment effect
    indirectEffect = function() {
      if (is.null(private$y_synth_draws)) {
        stop("Run $fit() before calling indirectEffect().")
      }
      .nasc_indirect_draws(self)
    },

    # Plot of the observed and SC outcome
    syntheticPlot = function() {
      if (is.null(private$plot_data)) {
        stop("Run $fit() before calling syntheticPlot().")
      }

      df <- as.data.frame(private$plot_data)
      time_name    <- rlang::as_name(private$time)
      outcome_name <- rlang::as_name(private$outcome)

      x   <- df[[time_name]]
      obs <- df[[outcome_name]]
      syn <- df$y_synth
      lb  <- df$LB
      ub  <- df$UB

      ord <- order(x)
      x <- x[ord]; obs <- obs[ord]; syn <- syn[ord]
      lb <- lb[ord]; ub <- ub[ord]

      yrng <- range(c(obs, syn, lb, ub), na.rm = TRUE)

      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op))
      graphics::par(bty = "l", mar = c(5, 4, 2, 1))

      plot(x, obs, type = "n", ylim = yrng,
           xlab = time_name, ylab = outcome_name)
      graphics::grid(lty = "dotted", col = "gray80")
      graphics::polygon(c(x, rev(x)), c(lb, rev(ub)),
                        col = grDevices::adjustcolor("gray", alpha.f = 0.2),
                        border = NA)
      graphics::lines(x, obs, lty = 1, lwd = 2)
      graphics::lines(x, syn, lty = 2, lwd = 2)
      graphics::abline(v = private$intervention, lty = 3)
      graphics::legend("topleft",
                       legend = c("Observed", "Synthetic"),
                       lty    = c(1, 2), lwd = 2,
                       bg     = grDevices::adjustcolor("white", alpha.f = 0.85),
                       box.col = "gray70")
      invisible(NULL)
    },

    # Plot estimated direct treatment effect
    effectPlot = function(indirect = NULL) {
      if (is.null(private$plot_data)) {
        stop("Run $fit() before calling effectPlot().")
      }

      time_name <- rlang::as_name(private$time)

      .plot_tau(
        data       = private$plot_data,
        x          = time_name,
        y          = "tau",
        ymin       = "tau_LB",
        ymax       = "tau_UB",
        xintercept = private$intervention
      )

      invisible(NULL)
    },

    # Posterior density plots
    posteriorPlot = function() {
      if (is.null(private$fitted) && is.null(private$y_synth_draws)) {
        stop("Run $fit() before calling posteriorPlot().")
      }

      panels <- list()

      add_panel <- function(label, draws) {
        if (is.null(draws)) return(invisible(NULL))
        x <- as.numeric(draws)
        x <- x[is.finite(x)]
        if (length(x) < 2L || stats::sd(x) == 0) return(invisible(NULL))
        panels[[length(panels) + 1L]] <<- list(label = label, draws = x)
      }

      step1_draws <- if (!is.null(private$fitted)) {
        tryCatch(
          rstan::extract(private$fitted, permuted = TRUE),
          error = function(e) list()
        )
      } else {
        list()
      }

      if (!is.null(step1_draws$beta_orig)) {
        step1_draws$beta <- step1_draws$beta_orig
      }
      if (!is.null(step1_draws$theta_orig)) {
        step1_draws$theta <- step1_draws$theta_orig
      }

      if (!is.null(step1_draws$rho)) {
        add_panel("rho", step1_draws$rho)
      }

      if (!is.null(step1_draws$theta)) {
        th <- step1_draws$theta
        th_labels <- if (!is.null(private$cov_names) &&
                         is.matrix(th) &&
                         length(private$cov_names) == ncol(th)) {
          sprintf("theta[%s]", private$cov_names)
        } else if (is.matrix(th)) {
          sprintf("theta[%d]", seq_len(ncol(th)))
        } else {
          "theta"
        }
        if (is.matrix(th)) {
          for (k in seq_len(ncol(th))) {
            col_k <- th[, k]
            if (any(is.finite(col_k))) {
              add_panel(th_labels[k], col_k)
            }
          }
        } else {
          if (any(is.finite(th))) add_panel(th_labels, th)
        }
      }

      if (!is.null(step1_draws$beta)) {
        bt <- step1_draws$beta
        K_b <- if (is.matrix(bt)) ncol(bt) else 1L
        b_labels <- if (!is.null(private$cov_names) &&
                        length(private$cov_names) == K_b) {
          sprintf("beta[%s]", private$cov_names)
        } else if (K_b > 1L) {
          sprintf("beta[%d]", seq_len(K_b))
        } else {
          "beta"
        }
        b_keep <- if (!is.null(private$beta_identified) &&
                      length(private$beta_identified) == K_b) {
          private$beta_identified
        } else {
          rep(TRUE, K_b)
        }
        if (is.matrix(bt)) {
          for (k in seq_len(K_b)) {
            if (!isTRUE(b_keep[k])) next
            col_k <- bt[, k]
            if (any(is.finite(col_k))) {
              add_panel(b_labels[k], col_k)
            }
          }
        } else {
          if (isTRUE(b_keep[1]) && any(is.finite(bt))) {
            add_panel(b_labels, bt)
          }
        }
      }

      if (!is.null(step1_draws$sigma_sar)) {
        add_panel("sigma_step1", step1_draws$sigma_sar)
      } else if (!is.null(step1_draws$sigma_sdm)) {
        add_panel("sigma_step1", step1_draws$sigma_sdm)
      }

      if (!is.null(private$y_synth_draws)) {
        if (!is.null(private$y_synth_draws$lambda_tilde)) {
          add_panel("lambda_tilde", private$y_synth_draws$lambda_tilde)
        }
        if (!is.null(private$y_synth_draws$sigma_sc)) {
          add_panel("sigma_step2", private$y_synth_draws$sigma_sc)
        }
        if (!is.null(private$y_synth_draws$bias_correction)) {
          add_panel("bias_correction", private$y_synth_draws$bias_correction)
        }
      }

      if (length(panels) == 0L) {
        stop("No estimated scalar parameters available to plot.")
      }

      n <- length(panels)
      ncol_grid <- if (n <= 4L) 1L else 2L
      nrow_grid <- ceiling(n / ncol_grid)

      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op))
      graphics::par(mfrow = c(nrow_grid, ncol_grid),
                    mar = c(4, 5, 2, 1),
                    bty = "l")

      for (pn in panels) {
        d <- stats::density(pn$draws, na.rm = TRUE)
        plot(d, main = pn$label,
             xlab = pn$label, ylab = "density")
        graphics::polygon(d,
                          col = grDevices::adjustcolor("steelblue", alpha.f = 0.3),
                          border = "steelblue")
      }

      invisible(NULL)
    },

    # Posterior density plot of ATE
    attPlot = function(indirect = NULL) {
      if (is.null(private$y_synth_draws)) stop("Run $fit() before calling attPlot.")

      ycf <- private$y_synth_draws$y_counterfactual
      bc  <- private$y_synth_draws$bias_correction
      if (is.null(bc)) bc <- rep(1, ncol(ycf))

      wide_df <- .makeWide(
        data      = private$data,
        id        = private$id,
        time      = private$time,
        outcome   = private$outcome,
        treatment = private$treated
      )
      post_data <- wide_df |>
        dplyr::filter(!!private$time >= private$intervention)
      Y1_post <- post_data[[rlang::as_name(private$outcome)]]

      Y1_mat <- matrix(Y1_post, nrow = nrow(ycf), ncol = length(Y1_post),
                       byrow = TRUE)
      bc_mat <- matrix(bc, nrow = nrow(ycf), ncol = ncol(ycf), byrow = FALSE)
      tau_draws <- (Y1_mat - ycf) * bc_mat
      att_draws <- rowMeans(tau_draws)

      indirect_default <- isTRUE(private$uses_rho)
      if (is.null(indirect)) indirect <- indirect_default
      stopifnot(is.logical(indirect), length(indirect) == 1L)

      avg_indirect_draws <- NULL
      if (indirect) {
        ind <- tryCatch(.nasc_indirect_draws(self), error = function(e) NULL)
        if (!is.null(ind)) {

          J <- length(ind$donor_names)
          avg_indirect_draws <- ind$avg_total / J
        } else {
          if (indirect_default) {
            warning("Indirect-effect draws unavailable; drawing ATT only.")
          }
          indirect <- FALSE
        }
      }

      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op))
      if (indirect) {
        graphics::par(mfrow = c(1, 2), mar = c(4, 5, 2, 1), bty = "l")
      } else {
        graphics::par(mar = c(4, 5, 2, 1), bty = "l")
      }

      d_att <- stats::density(att_draws, na.rm = TRUE)
      plot(d_att, main = "ATT (direct)",
           xlab = "ATT", ylab = "density")
      graphics::polygon(d_att,
                        col = grDevices::adjustcolor("steelblue", alpha.f = 0.3),
                        border = "steelblue")
      graphics::abline(v = 0, lty = 3, col = "gray40")

      if (indirect) {
        d_ind <- stats::density(avg_indirect_draws, na.rm = TRUE)
        plot(d_ind, main = "Average indirect effect (donor-avg)",
             xlab = expression(bar(bar(delta))), ylab = "density")
        graphics::polygon(d_ind,
                          col = grDevices::adjustcolor("indianred", alpha.f = 0.3),
                          border = "indianred")
        graphics::abline(v = 0, lty = 3, col = "gray40")
      }

      invisible(NULL)
    },

    # Posterior density plots of TE
    tauPlot = function(indirect = NULL) {
      if (is.null(private$y_synth_draws)) stop("Run $fit() before calling tauPlot.")

      ycf <- private$y_synth_draws$y_counterfactual
      bc  <- private$y_synth_draws$bias_correction
      if (is.null(bc)) bc <- rep(1, ncol(ycf))

      wide_df <- .makeWide(
        data      = private$data,
        id        = private$id,
        time      = private$time,
        outcome   = private$outcome,
        treatment = private$treated
      )
      post_data <- wide_df |>
        dplyr::filter(!!private$time >= private$intervention)
      Y1_post <- post_data[[rlang::as_name(private$outcome)]]
      time_post <- post_data[[rlang::as_name(private$time)]]

      Y1_mat <- matrix(Y1_post, nrow = nrow(ycf), ncol = length(Y1_post),
                       byrow = TRUE)
      bc_mat <- matrix(bc, nrow = nrow(ycf), ncol = ncol(ycf), byrow = FALSE)
      tau_draws <- (Y1_mat - ycf) * bc_mat

      n_t <- ncol(tau_draws)

      indirect_default <- isTRUE(private$uses_rho)
      if (is.null(indirect)) indirect <- indirect_default
      stopifnot(is.logical(indirect), length(indirect) == 1L)

      delta_avg_draws <- NULL
      if (indirect) {
        ind <- tryCatch(.nasc_indirect_draws(self), error = function(e) NULL)
        if (!is.null(ind)) {
          J <- length(ind$donor_names)
          delta_avg_draws <- ind$delta_total / J         # n_draws x T_post
        } else if (indirect_default) {
          warning("Indirect-effect draws unavailable; drawing direct effect only.")
        }
        if (is.null(delta_avg_draws)) indirect <- FALSE
      }

      # Direct-effect densities
      dens_dir <- lapply(seq_len(n_t),
                         function(i) stats::density(tau_draws[, i], na.rm = TRUE))
      xr_dir <- range(sapply(dens_dir, function(d) d$x))

      # Indirect-effect densities
      dens_ind <- NULL
      xr_ind   <- NULL
      if (indirect) {
        dens_ind <- lapply(seq_len(n_t),
                           function(i) stats::density(delta_avg_draws[, i],
                                                      na.rm = TRUE))
        xr_ind <- range(sapply(dens_ind, function(d) d$x))
      }

      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op))

      ncol_layout <- if (indirect) 2L else 1L
      graphics::par(mfrow = c(n_t, ncol_layout),
                    mar = c(2, 5, 1.2, 1),
                    oma = c(3, 0, 2.4, 0),
                    bty = "l")

      for (i in seq_len(n_t)) {
        is_last <- i == n_t

        # Direct
        d_i <- dens_dir[[i]]
        plot(d_i,
             main = paste0("Period ", time_post[i]), font.main = 1, cex.main = 1,
             xlab = if (is_last) expression(tau ~ "(direct)") else "",
             ylab = "density",
             xlim = xr_dir,
             xaxt = if (is_last) "s" else "n")
        graphics::polygon(d_i,
                          col = grDevices::adjustcolor("steelblue", alpha.f = 0.3),
                          border = "steelblue")
        graphics::abline(v = 0, lty = 3, col = "gray40")

        if (indirect) {
          d_j <- dens_ind[[i]]
          plot(d_j,
               main = paste0("Period ", time_post[i]), font.main = 1, cex.main = 1,
               xlab = if (is_last) expression(bar(delta) ~ "(indirect, donor-avg)") else "",
               ylab = "density",
               xlim = xr_ind,
               xaxt = if (is_last) "s" else "n")
          graphics::polygon(d_j,
                            col = grDevices::adjustcolor("indianred", alpha.f = 0.3),
                            border = "indianred")
          graphics::abline(v = 0, lty = 3, col = "gray40")
        }
      }

      if (indirect) {
        graphics::mtext(expression(bold(tau) ~ "(direct)"), side = 3, outer = TRUE,
                        at = 0.25, line = 0.6, font = 2)
        graphics::mtext(expression(bold(bar(delta)) ~ "(indirect, donor-avg)"), side = 3, outer = TRUE,
                        at = 0.75, line = 0.6, font = 2)
      } else {
        graphics::mtext("",
                        side = 3, outer = TRUE, line = 0.5, font = 2)
      }

      invisible(NULL)
    },

    # Plot of the posterior weight density per donor
    # CV table and selected penalty strength (NULL until $fit() with
    # nasc_penalty = TRUE and lambda = NULL has been run).
    lambdaCV = function() {
      list(lambda = private$lambda_used,
           table        = private$lambda_cv_table)
    },

    weightDraws = function(overlap = 0.5, scale = 1.4, fill_alpha = 0.85) {
      if (is.null(private$fitted)) stop("Run $fit() before calling weightDraws().")

      stopifnot(
        is.numeric(overlap), length(overlap) == 1L,
        overlap >= 0, overlap < 1,
        is.numeric(scale), length(scale) == 1L, scale > 0,
        is.numeric(fill_alpha), length(fill_alpha) == 1L,
        fill_alpha > 0, fill_alpha <= 1
      )

      w_mat <- private$y_synth_draws$w

      treated_id <- as.character(private$treated_ids)
      donor_names <- private$donor_ids
      if (is.null(donor_names)) {
        all_ids     <- levels(private$data[[rlang::as_name(private$id)]])
        donor_names <- setdiff(all_ids, treated_id)
      }

      if (ncol(w_mat) != length(donor_names)) {
        stop("Mismatch between the number of weights and donor names.")
      }
      colnames(w_mat) <- donor_names

      n <- ncol(w_mat)

      dens <- lapply(seq_len(n), function(i) {
        stats::density(w_mat[, i], na.rm = TRUE)
      })

      x_range <- range(unlist(lapply(dens, function(d) d$x)))
      max_y   <- max(vapply(dens, function(d) max(d$y), numeric(1)))
      if (!is.finite(max_y) || max_y <= 0) max_y <- 1

      ridge_h <- scale * max_y
      step    <- ridge_h * (1 - overlap)
      ylim_top <- step * n + ridge_h * 1.05
      ylim_bot <- -ridge_h * 0.05

      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op))
      graphics::par(mar = c(4, 6, 2, 1))

      plot(NA,
           xlim = x_range, ylim = c(ylim_bot, ylim_top),
           xlab = "Donor Weight", ylab = "Donor Unit",
           yaxt = "n", bty = "l")
      graphics::axis(2, at = step * seq_len(n), labels = donor_names, las = 1)

      cols   <- rep(grDevices::adjustcolor("steelblue", alpha.f = fill_alpha * 0.4), n)
      border <- "steelblue"

      for (i in seq(n, 1L, by = -1L)) {
        d <- dens[[i]]
        baseline <- step * i
        y <- baseline + d$y * (ridge_h / max_y)
        graphics::polygon(
          x = c(d$x, rev(d$x)),
          y = c(y, rep(baseline, length(d$x))),
          col    = cols[i],
          border = border
        )
        graphics::segments(x_range[1], baseline, x_range[2], baseline,
                           col = "gray60", lwd = 0.5)
      }

      invisible(NULL)
    },

    # Plot of the correlations between weights across draws
    weightCorr = function() {
      if (is.null(private$fitted)) stop("Run $fit() before calling weightCorr().")

      w_mat <- private$y_synth_draws$w

      treated_id <- as.character(private$treated_ids)
      donor_names <- private$donor_ids
      if (is.null(donor_names)) {
        all_ids     <- levels(private$data[[rlang::as_name(private$id)]])
        donor_names <- setdiff(all_ids, treated_id)
      }

      colnames(w_mat) <- donor_names

      cormat <- round(cor(w_mat), 3)
      diag(cormat) <- NA

      n <- length(donor_names)

      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op))
      graphics::par(mar = c(6, 6, 2, 5))

      pal <- grDevices::colorRampPalette(c("red", "white", "green"))(101)
      mat_plot <- t(cormat)[, n:1, drop = FALSE]

      graphics::image(
        x = seq_len(n), y = seq_len(n), z = mat_plot,
        zlim = c(-1, 1), col = pal,
        xlab = "Donor Units", ylab = "Donor Units",
        axes = FALSE
      )
      graphics::axis(1, at = seq_len(n), labels = donor_names, las = 2)
      graphics::axis(2, at = seq_len(n), labels = rev(donor_names), las = 1)
      graphics::box()

      lg <- seq(-1, 1, length.out = length(pal))
      usr <- graphics::par("usr")
      xl <- usr[2] + (usr[2] - usr[1]) * 0.02
      xr <- usr[2] + (usr[2] - usr[1]) * 0.05
      yb <- seq(usr[3], usr[4], length.out = length(pal) + 1)
      graphics::rect(xl, yb[-length(yb)], xr, yb[-1],
                     col = pal, border = NA, xpd = TRUE)
      graphics::text(xr, c(usr[3], mean(usr[3:4]), usr[4]),
                     labels = c("-1", "0", "1"),
                     pos = 4, xpd = TRUE)
      graphics::text(mean(c(xl, xr)), usr[4], "Corr",
                     pos = 3, xpd = TRUE)
      invisible(NULL)
    }
  )
)
