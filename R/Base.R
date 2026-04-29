# Bayesian NASC estimator
# BASE code

#' Bayesian Network-Aware Synthetic Control
#'
#' @description
#' An R6 class implementing Bayesian Synthetic Control estimators with optional
#' network/spatial awareness.
#'
#' Three orthogonal options govern how the spatial structure enters the
#' estimator:
#'
#' \itemize{
#'   \item \code{spatial_model}: \code{"none"}, \code{"SAR"}, \code{"SDM"}, or
#'         \code{"exogenous"}. The first three behave as before; the last skips
#'         Step 1 and uses a user-supplied scalar \code{rho}.
#'   \item \code{rho}: optional numeric in \code{(-1, 1)}. When provided,
#'         Step 1 is skipped regardless of \code{spatial_model} and the supplied
#'         value is used directly.
#'   \item \code{bias_correction} (logical): whether the post-treatment effect
#'         is rescaled by \eqn{1/(1 - w's)}.
#'   \item \code{nasc_penalty} (logical): whether the NASC penalty
#'         \eqn{-\lambda \langle w, |s| \rangle} enters the likelihood.
#' }
#'
#' Model dispatch (single source of truth):
#' \itemize{
#'   \item \code{nasc_penalty = TRUE}  -> stan_2_NASC.stan (penalty always on
#'         in this file, bias_correction toggleable).
#'   \item \code{nasc_penalty = FALSE} -> model1.stan (no penalty; bias
#'         correction toggleable through the same generated-quantities mechanism).
#' }
#'
#' Step 1 (SAR or SDM rho estimation) runs only when (i) the chosen Step 2
#' configuration uses a non-trivial rho (penalty on, or bias correction on)
#' AND (ii) no exogenous rho was supplied.
#'
#' @export
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
    # Internal flags:
    #   uses_rho   - TRUE iff penalty or bias correction is active. Determines
    #                whether s (and therefore rho/W/w_J1) is needed at all.
    #   uses_step1 - TRUE iff Step 1 (SAR or SDM) must be fitted.
    uses_rho   = NULL,
    uses_step1 = NULL
  ),
  active = list(
    #' @field plotData Tibble with observed and counterfactual outcomes.
    plotData = function() { return(private$plot_data) },

    #' @field interventionTime The first treatment period.
    interventionTime = function() { return(private$intervention) },

    #' @field synthetic ggplot2 plot of observed vs synthetic outcomes.
    synthetic = function() {
      df_plot <- private$plot_data |>
        dplyr::rename(Observed = !!private$outcome, Synthetic = y_synth) |>
        dplyr::select(!!private$time, Observed, Synthetic, LB, UB) |>
        tidyr::pivot_longer(cols = c(Observed, Synthetic))

      ggplot2::ggplot(data = df_plot, ggplot2::aes(x = !!private$time)) +
        ggplot2::geom_line(ggplot2::aes(y = value, linetype = name)) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = LB, ymax = UB),
          color = "gray", alpha = 0.2
        ) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
          legend.title   = ggplot2::element_blank(),
          legend.position = c(0.9, .1),
          legend.background = ggplot2::element_rect(
            fill = ggplot2::alpha("white", 0)
          ),
          panel.border = ggplot2::element_blank(),
          axis.line    = ggplot2::element_line()
        ) +
        ggplot2::geom_vline(
          xintercept = private$intervention, linetype = "dashed"
        )
    }
  ),
  public = list(

    #' @description
    #' Create a new nascSynth object.
    #'
    #' @param spatial_model One of \code{"none"}, \code{"SAR"}, \code{"SDM"},
    #'   or \code{"exogenous"}.
    #' @param rho Optional scalar in \code{(-1, 1)}. When supplied, Step 1 is
    #'   skipped and this value is used in Step 2 regardless of
    #'   \code{spatial_model}. Required when \code{spatial_model = "exogenous"}.
    #' @param bias_correction Logical. If \code{TRUE}, the post-treatment
    #'   counterfactual is rescaled by \eqn{1/(1 - w's)}. Default is \code{TRUE}
    #'   for SAR/SDM/exogenous and \code{FALSE} for \code{spatial_model = "none"}.
    #' @param nasc_penalty Logical. If \code{TRUE}, the NASC penalty
    #'   \eqn{-\lambda \langle w, |s| \rangle} enters the likelihood. Same
    #'   model-aware default as \code{bias_correction}.
    initialize = function(data, time, id, treated, outcome, ci_width = 0.75,
                          covariates = NULL,
                          W = NULL, spatial_model = "none",
                          rho = NULL,
                          bias_correction = NULL,
                          nasc_penalty = NULL) {

      stopifnot(ci_width > 0 & ci_width < 1)

      if (!spatial_model %in% c("none", "SAR", "SDM", "exogenous")) {
        stop("spatial_model must be 'none', 'SAR', 'SDM', or 'exogenous'.")
      }

      # Model-aware defaults preserve legacy behaviour for plain `none` calls.
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

      # Validate exogenous rho when supplied
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

      # Does any spatial element of the estimator activate?
      uses_rho <- bias_correction || nasc_penalty

      # Need W whenever rho is in use. Don't need W otherwise.
      if (uses_rho && is.null(W)) {
        stop(
          "A spatial weights matrix 'W' is required whenever ",
          "bias_correction = TRUE or nasc_penalty = TRUE."
        )
      }

      # When rho is in use under spatial_model = "none", the user must supply
      # an exogenous rho (no Step 1 model to fit).
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
      private$W               <- W
      private$spatial_model   <- spatial_model
      private$rho_exogenous   <- rho
      private$bias_correction <- bias_correction
      private$nasc_penalty    <- nasc_penalty
      private$uses_rho        <- uses_rho

      # Step 1 runs only if rho is needed AND the user did not supply one.
      private$uses_step1 <- uses_rho &&
        spatial_model %in% c("SAR", "SDM") &&
        is.null(rho)

      # ---- Status message ----
      engine <- if (nasc_penalty) "stan_2_NASC" else "model1"
      rho_src <- if (!uses_rho) {
        "n/a (no rho in use)"
      } else if (!is.null(rho)) {
        sprintf("exogenous (rho = %.4f)", rho)
      } else {
        sprintf("Step 1 %s posterior", spatial_model)
      }
      message(sprintf(
        "nascSynth: engine = %s | nasc_penalty = %s | bias_correction = %s | rho source = %s",
        engine, nasc_penalty, bias_correction, rho_src
      ))

      if (!setequal(
        data |>
        dplyr::select(!!private$treated) |>
        dplyr::distinct() |>
        dplyr::pull({{ treated }}),
        c(1, 0)
      )) {
        stop("Treated identifier is not binary (1/0).")
      }

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

      # Stan model dispatch
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
        step2 = if (nasc_penalty) stanmodels$stan_2_NASC else stanmodels$model1
      )
    },

    #' @description
    #' Fit the Stan model via MCMC.
    #' @param n_samples Number of rho draws to propagate to Step 2 when Step 1
    #'   runs (default 100). Ignored when an exogenous rho was supplied or when
    #'   no rho is in use.
    #' @param cores Number of CPU cores for parallel execution.
    #' @param ... Additional arguments forwarded to [rstan::sampling()].
    fit = function(n_samples = 100, cores = parallel::detectCores() - 1, ...) {
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

      # Number of *real* pre-treatment time periods. The covariate-as-extra-rows
      # feature below (model1 engine only) appends synthetic rows to X/X1 so the
      # Stan model also matches predictors. Stan then returns y_sim_pre with one
      # entry per row of X (i.e. n_pre_real + n_covariate_rows). We need to
      # remember n_pre_real here so we can slice y_sim_pre back down to real
      # time periods before post-processing.
      n_pre_real <- nrow(pre_data)

      # Covariate-as-extra-rows feature (Abadie-style predictor matching).
      # Only applied for the model1 engine, since stan_2_NASC's Y_panel has
      # spatial-aligned dimensions that cannot be extended with covariate rows.
      use_covariate_rows <- !private$nasc_penalty && !is.null(private$covariates)

      if (use_covariate_rows) {
        cov_names <- setdiff(
          names(private$covariates),
          c(rlang::as_name(private$time), rlang::as_name(private$id))
        )
        cov_wide <- private$covariates |>
          dplyr::filter(!!private$time < private$intervention) |>
          dplyr::group_by(!!private$id) |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(cov_names), mean),
            .groups = "drop"
          ) |>
          tidyr::pivot_longer(
            cols      = dplyr::all_of(cov_names),
            names_to  = ".covariate",
            values_to = ".value"
          ) |>
          tidyr::pivot_wider(
            names_from  = !!private$id,
            values_from = .value
          ) |>
          dplyr::select(-.covariate)

        treated_col <- as.character(private$treated_ids)
        donor_cols  <- colnames(X)

        X  <- rbind(X,  as.data.frame(cov_wide[, donor_cols,  drop = FALSE]))
        X1 <- c(X1, unlist(cov_wide[, treated_col, drop = TRUE]))
      }

      donor_ids  <- colnames(X_pred)  # post-treatment donor list (pre-cov stack-safe)
      treated_id <- as.character(private$treated_ids)

      # ----------------------------------------------------------------------
      # Build spatial inputs (W block, lambda_W, optional Step-1 covariates)
      # only if we actually need them.
      # ----------------------------------------------------------------------
      W_J      <- NULL
      w_J1     <- NULL
      W_full   <- NULL
      lambda_W <- NULL
      Y_panel  <- NULL
      X0_mat   <- NULL
      X1_vec   <- NULL
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
          if (max(abs(Im(ev))) > 1e-8) {
            warning("W has non-trivial complex eigenvalues; SAR/SDM ",
                    "identification may be unreliable.")
          }
          lambda_W <- Re(ev)

          # Build covariate matrices for Step 1
          if (is.null(private$covariates)) {
            stop("SAR/SDM Step 1 requires 'covariates'. Either supply ",
                 "covariates or pass an exogenous 'rho'.")
          }
          cov_names <- setdiff(
            names(private$covariates),
            c(rlang::as_name(private$time), rlang::as_name(private$id))
          )
          if (length(cov_names) == 0L) {
            stop("'covariates' contains no predictor columns (only time/id found).")
          }
          cov_wide <- private$covariates |>
            dplyr::filter(!!private$time < private$intervention) |>
            dplyr::group_by(!!private$id) |>
            dplyr::summarise(
              dplyr::across(dplyr::all_of(cov_names), mean),
              .groups = "drop"
            ) |>
            tidyr::pivot_longer(
              cols      = dplyr::all_of(cov_names),
              names_to  = ".covariate",
              values_to = ".value"
            ) |>
            tidyr::pivot_wider(
              names_from  = !!private$id,
              values_from = .value
            )
          pred_names <- cov_wide$.covariate
          cov_wide   <- cov_wide |> dplyr::select(-.covariate)

          X0_mat <- as.matrix(cov_wide[, donor_ids,  drop = FALSE])
          X1_vec <- as.numeric(cov_wide[[treated_id]])
          K_pred <- length(pred_names)

          if (length(X1_vec) != K_pred) {
            stop("Treated unit predictor vector length does not match K_pred.")
          }
          if (anyNA(X0_mat) || anyNA(X1_vec)) {
            stop("Missing values found in predictor matrix; please impute or drop.")
          }

          # Y_panel for Step 1 (pre-cov stacking) and stan_2_NASC
          Y_panel <- t(as.matrix(
            pre_data |> dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
          ))
        } else if (private$nasc_penalty) {
          # Need Y_panel for stan_2_NASC even without Step 1
          Y_panel <- t(as.matrix(
            pre_data |> dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
          ))
        }
      }

      # ----------------------------------------------------------------------
      # STEP 1 (only if needed)
      # ----------------------------------------------------------------------
      sampled_rhos <- if (private$uses_rho) {
        if (!is.null(private$rho_exogenous)) {
          private$rho_exogenous   # single scalar
        } else if (private$uses_step1) {

          step1_data <- list(
            K_pred   = K_pred,
            X1       = X1_vec,
            J        = length(donor_ids),
            X0       = X0_mat,
            T0       = nrow(pre_data),
            Y_panel  = Y_panel,
            W        = W_full,
            lambda_W = lambda_W
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
        NA_real_   # marker: no rho needed
      }

      # ----------------------------------------------------------------------
      # STEP 2 - dispatch to model1 or stan_2_NASC depending on nasc_penalty
      # ----------------------------------------------------------------------
      if (private$nasc_penalty) {

        # ---- stan_2_NASC engine ----
        # Always loops over rho draws (parallel when >1, single when scalar).
        base_data <- list(
          J                   = length(donor_ids),
          T0                  = nrow(pre_data),
          Y_panel             = Y_panel,
          W_J                 = W_J,
          w_J1                = w_J1,
          T_post              = nrow(post_data),
          Y0_post             = as.matrix(X_pred[, donor_ids, drop = FALSE]),
          Y1_post             = post_data |> dplyr::pull(!!private$outcome),
          use_bias_correction = as.integer(private$bias_correction)
        )

        results <- .run_step2_loop(
          rhos       = sampled_rhos,
          base_data  = base_data,
          step2_mod  = private$stan_model$step2,
          rho_field  = "rho",
          cores      = cores,
          extra_args = list(...),
          extract_pars = c("y_counterfactual", "y_sim_pre", "w", "lambda",
                           "sigma_sc", "bias_correction")
        )

        # When no Step 1 ran, expose Step 2 as the user-facing fit.
        if (is.null(private$fitted)) private$fitted <- results$last_fit

        private$y_synth_draws <- list(
          y_counterfactual = results$y_counterfactual,
          y_sim_pre        = results$y_sim_pre,
          w                = results$w,
          lambda           = results$lambda,
          sigma_sc         = results$sigma_sc,
          bias_correction  = results$bias_correction,
          rhos_used        = results$rhos_used
        )

      } else {

        # ---- model1 engine (with optional bias correction) ----
        # Build the bias-correction inputs only if active. When inactive,
        # J_bc = 0 and the spatial inputs are zero-sized placeholders Stan
        # tolerates without computing anything.
        if (private$bias_correction) {
          # In model1, the simplex `w` has length K = ncol(X). We need W_J
          # to also be K x K aligned to the donor columns of X. The donors
          # in X_pred (donor_ids) match the donor columns of X by construction
          # (they are derived from the same wide_df).
          J_bc   <- ncol(X)
          W_J_bc <- W_J
          w_J1_bc <- w_J1
          rho_bc <- if (length(sampled_rhos) == 1L) sampled_rhos else NA_real_
          # When sampled_rhos is a vector (Step 1 path), each Step 2 worker
          # plugs in its own rho. The .run_step2_loop helper handles that via
          # the `rho_field` argument.
        } else {
          J_bc    <- 0L
          W_J_bc  <- matrix(0, 0, 0)
          w_J1_bc <- numeric(0)
          rho_bc  <- 0  # arbitrary value within (-1, 1); unused when J_bc = 0
        }

        base_data <- list(
          N                   = nrow(X),
          y                   = X1,
          K                   = ncol(X),
          X                   = as.matrix(X),
          N_pred              = nrow(X_pred),
          X_pred              = as.matrix(X_pred),
          use_bias_correction = as.integer(private$bias_correction),
          J_bc                = J_bc,
          W_J                 = W_J_bc,
          w_J1                = w_J1_bc,
          rho_bc              = if (private$bias_correction) rho_bc else 0
        )

        if (!private$bias_correction) {

          # No spatial element at all -> single Stan fit, no rho loop.
          private$fitted <- rstan::sampling(
            private$stan_model$step2,
            data  = base_data,
            cores = cores,
            ...
          )

          draws <- rstan::extract(
            private$fitted,
            pars = c("y_counterfactual", "y_sim_pre", "w", "sigma",
                     "bias_correction")
          )
          # y_sim_pre has one column per row of X. When covariate rows were
          # stacked above, those trailing columns aren't real time periods --
          # drop them so the post-processing aligns with pre_data.
          y_sim_pre_trim <- draws$y_sim_pre[, seq_len(n_pre_real), drop = FALSE]
          private$y_synth_draws <- list(
            y_counterfactual = draws$y_counterfactual,
            y_sim_pre        = y_sim_pre_trim,
            w                = draws$w,
            sigma            = draws$sigma,
            bias_correction  = draws$bias_correction,  # all 1.0 in this branch
            rhos_used        = NA_real_
          )

        } else {

          # Bias correction on -> loop model1 across rho draws (same
          # propagation pattern as the stan_2_NASC path).
          results <- .run_step2_loop(
            rhos       = sampled_rhos,
            base_data  = base_data,
            step2_mod  = private$stan_model$step2,
            rho_field  = "rho_bc",
            cores      = cores,
            extra_args = list(...),
            extract_pars = c("y_counterfactual", "y_sim_pre", "w", "sigma",
                             "bias_correction")
          )

          if (is.null(private$fitted)) private$fitted <- results$last_fit

          # See the no-bias-correction branch above: drop covariate rows from
          # y_sim_pre so it aligns with pre_data.
          y_sim_pre_trim <- results$y_sim_pre[, seq_len(n_pre_real), drop = FALSE]
          private$y_synth_draws <- list(
            y_counterfactual = results$y_counterfactual,
            y_sim_pre        = y_sim_pre_trim,
            w                = results$w,
            sigma            = results$sigma,
            bias_correction  = results$bias_correction,
            rhos_used        = results$rhos_used
          )
        }
      }

      # ----------------------------------------------------------------------
      # Unified post-processing (single helper, both engines feed in)
      # ----------------------------------------------------------------------
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

    #' @description
    #' Update the credible interval width and recompute plot data.
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

    #' @description
    #' Posterior summary of the fit: ATT, per-period treatment effects,
    #' donor weights, model parameters and MCMC diagnostics. Prints a
    #' formatted summary and invisibly returns a \code{summary.nascSynth}
    #' list. Also dispatched from \code{summary(obj)} via the S3 method
    #' \code{summary.nascSynth}.
    #' @param ci_width Optional override for the credible interval width.
    #' @param print Logical; if \code{TRUE} (default), print the formatted
    #'   summary to the console. Set to \code{FALSE} to retrieve the
    #'   structured list silently for programmatic use.
    summary = function(ci_width = NULL, print = TRUE) {
      if (is.null(private$y_synth_draws)) {
        stop("Run $fit() before calling summary().")
      }
      ci <- if (is.null(ci_width)) private$ci_width else {
        stopifnot(is.numeric(ci_width), length(ci_width) == 1L,
                  ci_width > 0, ci_width < 1)
        ci_width
      }

      # Reconstruct the rho-source string the constructor printed, so the
      # summary stays self-describing without needing to capture it earlier.
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
        id_levels       = private$data[[rlang::as_name(private$id)]],
        spatial_model   = private$spatial_model,
        bias_correction = private$bias_correction,
        nasc_penalty    = private$nasc_penalty,
        uses_rho        = private$uses_rho,
        rho_source      = rho_source,
        fitted          = private$fitted
      ))

      # Print the formatted summary explicitly. This works even if the S3
      # method `print.summary.nascSynth` was not registered in NAMESPACE,
      # because we're calling the function by its full name from the
      # package namespace rather than going through S3 dispatch. Return
      # invisibly so the structured list is still available via
      # `s <- synth_sc$summary()` without spamming the console twice.
      if (isTRUE(print)) {
        print.summary.nascSynth(out)
      }
      invisible(out)
    },

    #' @description
    #' Plot the estimated treatment effect (tau) over time.
    effectPlot = function() {
      .plot_tau(
        data       = private$plot_data,
        x          = private$time,
        y          = tau,
        ymin       = tau_LB,
        ymax       = tau_UB,
        xintercept = private$intervention
      )
    },

    #' @description
    #' Plot implicit weight distribution across draws.
    weightDraws = function() {
      if (is.null(private$fitted)) stop("Run fit() first.")

      w_mat <- private$y_synth_draws$w

      treated_id <- as.character(private$treated_ids)
      all_ids <- levels(private$data[[rlang::as_name(private$id)]])
      donor_names <- setdiff(all_ids, treated_id)

      if (ncol(w_mat) != length(donor_names)) {
        stop("Mismatch between the number of weights and donor names.")
      }

      colnames(w_mat) <- donor_names

      melt_w <- as.data.frame(w_mat) |>
        tidyr::pivot_longer(
          cols      = dplyr::everything(),
          names_to  = "ID",
          values_to = "weight"
        )

      melt_w |>
        ggplot2::ggplot(ggplot2::aes(x = weight, y = ID, fill = ID)) +
        ggridges::geom_density_ridges(alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = "Donor Weight", y = "Donor Unit")
    },

    #' @description
    #' Plot correlations between weights across draws.
    weightCorr = function() {
      if (is.null(private$fitted)) stop("Run fit() first.")

      w_mat <- private$y_synth_draws$w

      treated_id <- as.character(private$treated_ids)
      all_ids <- levels(private$data[[rlang::as_name(private$id)]])
      donor_names <- setdiff(all_ids, treated_id)

      colnames(w_mat) <- donor_names

      cormat <- round(cor(w_mat), 3)
      diag(cormat) <- NA

      melted_cormat <- as.data.frame(cormat) |>
        dplyr::mutate(X1 = factor(rownames(cormat), levels = donor_names)) |>
        tidyr::pivot_longer(
          cols      = -X1,
          names_to  = "X2",
          values_to = "value"
        ) |>
        dplyr::mutate(X2 = factor(X2, levels = donor_names))

      ggplot2::ggplot(data = melted_cormat, ggplot2::aes(x = X1, y = X2, fill = value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_fill_gradient2(
          low      = "red",
          mid      = "white",
          high     = "green",
          midpoint = 0,
          limit    = c(-1, 1),
          space    = "Lab",
          name     = "Corr",
          na.value = "transparent"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x      = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::labs(x = "Donor Units", y = "Donor Units")
    }
  )
)
