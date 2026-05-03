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
    donor_ids   = NULL,
    uses_rho   = NULL,
    uses_step1 = NULL
  ),
  active = list(
    #' @field plotData Tibble with observed and counterfactual outcomes
    plotData = function() { return(private$plot_data) },

    #' @field interventionTime The first treatment period
    interventionTime = function() { return(private$intervention) }
  ),
  public = list(

    #' @description
    #' Create a new nascSynth object.
    #'
    #' @param data A data frame in long format containing the panel data.
    #'   Must include columns for the unit identifier, time, outcome, and
    #'   a binary treatment indicator (1 = treated, 0 = untreated).
    #' @param time Unquoted column name for the time variable.
    #' @param id Unquoted column name for the unit identifier.
    #' @param treated Unquoted column name for the binary treatment indicator
    #'   (1 in treated periods for the treated unit, 0 otherwise).
    #' @param outcome Unquoted column name for the outcome variable.
    #' @param ci_width Numeric in \code{(0, 1)}. Width of the credible interval.
    #'   Default is \code{0.75}.
    #' @param covariates Optional data frame in long format with columns for
    #'   \code{id}, \code{time}, and one or more predictor variables. Required
    #'   when \code{spatial_model} is \code{"SAR"} or \code{"SDM"}.
    #' @param W Optional row-standardized spatial weights matrix. Row and column
    #'   names should match the unit identifiers in \code{data}. Required when
    #'   \code{bias_correction = TRUE} or \code{nasc_penalty = TRUE}.
    #' @param spatial_model One of \code{"none"}, \code{"SAR"}, \code{"SDM"},
    #'   or \code{"exogenous"}.
    #' @param rho Optional scalar in \code{(-1, 1)}. When supplied, Step 1 is
    #'   skipped and this value is used in Step 2 regardless of
    #'   \code{spatial_model}. Required when \code{spatial_model = "exogenous"}.
    #' @param bias_correction Logical. If \code{TRUE}, the post-treatment
    #'   counterfactual is rescaled by \eqn{1/(1 - w's)}. Default is \code{TRUE}
    #'   for SAR/SDM/exogenous and \code{FALSE} for
    #'   \code{spatial_model = "none"}.
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
      private$W               <- W
      private$spatial_model   <- spatial_model
      private$rho_exogenous   <- rho
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

      # ----------------------------------------------------------------
      # Robust binary check on the treated indicator.
      #
      # The previous setequal({0,1}, distinct(D)) check failed legitimate
      # inputs in three common cases:
      #   * NAs present: distinct returns {0, 1, NA}, setequal -> FALSE
      #   * factor with levels "0"/"1": setequal compares strings, OK by
      #     coincidence, but downstream `D == 1` may behave oddly
      #   * logical TRUE/FALSE: coerces, but fragile
      #
      # We accept numeric, integer, logical, and a "0"/"1" factor; we
      # require both 0 and 1 to actually appear (no degenerate columns);
      # NAs are tolerated and treated as "N/A" downstream.
      # ----------------------------------------------------------------
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
        step2 = if (nasc_penalty) stanmodels$stan_2_NASC else stanmodels$model1
      )
    },

    #' @description
    #' Fit the Stan model via MCMC.
    #' @param n_samples Number of rho draws to propagate to Step 2 when Step 1
    #'   runs. Either a positive integer (default 100), or the string
    #'   \code{"auto"}. When \code{"auto"}, the number is set to
    #'   \code{min(round(n_eff_rho), n_samples_cap)} after Step 1 completes,
    #'   floored at \code{n_samples_min}, so that the cut-posterior
    #'   approximation does not oversample a low-information rho posterior.
    #'   Ignored when an exogenous rho was supplied or when no rho is in use.
    #' @param n_samples_cap Upper bound on the auto-selected \code{n_samples}.
    #'   Default 200. Only used when \code{n_samples = "auto"}.
    #' @param n_samples_min Lower bound on the auto-selected \code{n_samples}.
    #'   Default 30. Only used when \code{n_samples = "auto"}.
    #' @param cores Number of CPU cores for parallel execution.
    #' @param ... Additional arguments forwarded to [rstan::sampling()].
    #' @param worker_iter Iterations per worker chain in the multi-rho parallel
    #'   loop. Default 2000. Increase if Stan reports low Bulk/Tail ESS.
    #' @param worker_warmup Warmup per worker chain in the multi-rho parallel
    #'   loop. Default 1000.
    fit = function(n_samples = 100,
                   n_samples_cap = 200L,
                   n_samples_min = 30L,
                   cores = parallel::detectCores() - 1,
                   worker_iter = 2000L, worker_warmup = 1000L, ...) {

      # Validate n_samples up front: integer >= 1 or the literal string "auto".
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

      # If auto was requested but Step 1 won't run, n_samples is unused
      # downstream anyway. Keep the variable numerically well-typed for
      # safety and let the user know we ignored the auto request.
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

      donor_ids  <- colnames(X_pred)
      treated_id <- as.character(private$treated_ids)

      # Store the canonical donor ordering so post-hoc helpers (summary,
      # contamination plots, etc.) can align posterior draws with W and
      # with donor names without re-deriving them and getting it wrong.
      private$donor_ids <- donor_ids

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

          if (is.null(private$covariates)) {
            stop("SAR/SDM Step 1 requires 'covariates'. Either supply ",
                 "covariates or pass an exogenous rho.")
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

          Y_panel <- t(as.matrix(
            pre_data |> dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
          ))
        } else if (private$nasc_penalty) {
          Y_panel <- t(as.matrix(
            pre_data |> dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
          ))
        }
      }

      # STEP 1 (only if needed)
      sampled_rhos <- if (private$uses_rho) {
        if (!is.null(private$rho_exogenous)) {
          private$rho_exogenous
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
          total_rho_draws <- n_iters * n_chains

          # ----------------------------------------------------------------
          # Auto-sizing of n_samples from the Step-1 ESS of rho.
          #
          # The cut posterior is approximated by a Monte Carlo average over
          # rho draws. There is no point integrating over more independent
          # rho values than the Step-1 chain actually produced, so we cap
          # the subsample size at round(n_eff_rho), then bound the result
          # below by n_samples_min and above by n_samples_cap. We also cap
          # at the total number of post-warmup draws so we never request
          # more than exists.
          # ----------------------------------------------------------------
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

      # STEP 2
      if (private$nasc_penalty) {
        base_data <- list(
          J                   = length(donor_ids),
          T0                  = nrow(pre_data),
          Y_panel             = Y_panel,
          W_J                 = W_J,
          w_J1                = w_J1,
          T_post              = nrow(post_data),
          Y0_post             = as.matrix(X_pred[, donor_ids, drop = FALSE]),
          # as.array() forces a length-1 numeric to serialize as a length-1
          # vector rather than a scalar (rstan would otherwise pass dims=()
          # and Stan would reject vector[1] declarations). No effect when
          # length > 1.
          Y1_post             = as.array(as.numeric(post_data |> dplyr::pull(!!private$outcome))),
          use_bias_correction = as.integer(private$bias_correction)
        )

        results <- .run_step2_loop(
          rhos          = sampled_rhos,
          base_data     = base_data,
          step2_mod     = private$stan_model$step2,
          rho_field     = "rho",
          cores         = cores,
          extra_args    = list(...),
          extract_pars  = c("y_counterfactual", "y_sim_pre", "w", "lambda",
                            "sigma_sc", "bias_correction"),
          worker_iter   = worker_iter,
          worker_warmup = worker_warmup
        )

        if (is.null(private$fitted)) private$fitted <- results$last_fit

        private$y_synth_draws <- list(
          y_counterfactual = results$y_counterfactual,
          y_sim_pre        = results$y_sim_pre,
          w                = results$w,
          lambda           = results$lambda,
          sigma_sc         = results$sigma_sc,
          bias_correction  = results$bias_correction,
          rhos_used        = results$rhos_used,
          worker_diagnostics = results$worker_diagnostics
        )

      } else {
        if (private$bias_correction) {
          J_bc   <- ncol(X)
          W_J_bc <- W_J
          w_J1_bc <- w_J1
          rho_bc <- if (length(sampled_rhos) == 1L) sampled_rhos else NA_real_
        } else {
          J_bc    <- 0L
          W_J_bc  <- matrix(0, 0, 0)
          w_J1_bc <- numeric(0)
          rho_bc  <- 0
        }

        base_data <- list(
          N                   = nrow(X),
          # as.array() forces length-1 numerics to serialize as length-1
          # vectors (Stan rejects scalars where vector[N] is declared).
          y                   = as.array(as.numeric(X1)),
          K                   = ncol(X),
          X                   = as.matrix(X),
          N_pred              = nrow(X_pred),
          X_pred              = as.matrix(X_pred),
          use_bias_correction = as.integer(private$bias_correction),
          J_bc                = J_bc,
          W_J                 = W_J_bc,
          w_J1                = if (J_bc > 0L) as.array(as.numeric(w_J1_bc)) else numeric(0),
          rho_bc              = if (private$bias_correction) rho_bc else 0
        )

        if (!private$bias_correction) {
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
          y_sim_pre_trim <- draws$y_sim_pre[, seq_len(n_pre_real), drop = FALSE]
          private$y_synth_draws <- list(
            y_counterfactual = draws$y_counterfactual,
            y_sim_pre        = y_sim_pre_trim,
            w                = draws$w,
            sigma            = draws$sigma,
            bias_correction  = draws$bias_correction,
            rhos_used        = NA_real_
          )

        } else {
          results <- .run_step2_loop(
            rhos          = sampled_rhos,
            base_data     = base_data,
            step2_mod     = private$stan_model$step2,
            rho_field     = "rho_bc",
            cores         = cores,
            extra_args    = list(...),
            extract_pars  = c("y_counterfactual", "y_sim_pre", "w", "sigma",
                              "bias_correction"),
            worker_iter   = worker_iter,
            worker_warmup = worker_warmup
          )

          if (is.null(private$fitted)) private$fitted <- results$last_fit

          y_sim_pre_trim <- results$y_sim_pre[, seq_len(n_pre_real), drop = FALSE]
          private$y_synth_draws <- list(
            y_counterfactual = results$y_counterfactual,
            y_sim_pre        = y_sim_pre_trim,
            w                = results$w,
            sigma            = results$sigma,
            bias_correction  = results$bias_correction,
            rhos_used        = results$rhos_used,
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
        # Pass `self` so .nasc_summary_stats can pull contamination /
        # indirect-effect draws via .nasc_indirect_draws() (which expects
        # a fitted nascSynth object).
        model           = self
      ))

      if (isTRUE(print)) {
        print.summary.nascSynth(out)
      }
      invisible(out)
    },

    #' @description
    #' Posterior draws and summaries of the indirect (spillover) treatment
    #' effect. By Proposition 6.2, the per-donor spillover at post-period
    #' \eqn{t} is \eqn{\delta_{i,t}^{NASC} = s_i \cdot \tau_{1t}^{NASC}},
    #' where \eqn{s = \rho (I_J - \rho W_J)^{-1} w_{J1}}.
    #'
    #' Returns \code{NULL} when the model is configured without a network
    #' (\code{uses_rho = FALSE}), since spillover is then identically zero
    #' by construction.
    #'
    #' @return Either \code{NULL} (no rho) or a named list with components:
    #'   \itemize{
    #'     \item \code{donor_names} -- character vector of donor IDs.
    #'     \item \code{time_post}  -- post-period time vector.
    #'     \item \code{delta_arr}  -- \code{[n_draws x T_post x J]} array
    #'           of per-draw, per-period, per-donor spillover.
    #'     \item \code{delta_total} -- \code{[n_draws x T_post]} matrix of
    #'           per-period total spillover (\eqn{\sum_i \delta_{i,t}}).
    #'     \item \code{avg_per_donor} -- \code{[n_draws x J]} matrix of
    #'           per-donor average spillover across post-periods.
    #'     \item \code{avg_total} -- \code{[n_draws]} vector of average
    #'           total spillover (the spillover analog of ATT).
    #'   }
    indirectEffect = function() {
      if (is.null(private$y_synth_draws)) {
        stop("Run $fit() before calling indirectEffect().")
      }
      .nasc_indirect_draws(self)
    },

    #' @description
    #' Base R plot of the observed and synthetic-control outcome trajectories,
    #' with a shaded credible-interval band around the synthetic series and a
    #' dotted vertical line at the intervention time.
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

    #' @description
    #' Plot the estimated treatment effect (tau) over time. When the model
    #' uses a network and \code{indirect = TRUE} (the default in that
    #' case), the per-period total indirect (spillover) effect is drawn
    #' in a second stacked panel below the direct effect. Set
    #' \code{indirect = FALSE} to suppress the indirect panel even when
    #' a network is in use.
    #'
    #' @param indirect Logical or \code{NULL}. When \code{NULL} (default),
    #'   the indirect panel is shown iff \code{uses_rho = TRUE}.
    effectPlot = function(indirect = NULL) {
      indirect_default <- isTRUE(private$uses_rho)
      if (is.null(indirect)) indirect <- indirect_default
      stopifnot(is.logical(indirect), length(indirect) == 1L)

      time_name <- rlang::as_name(private$time)

      # Direct effect: as before, from plot_data.
      if (indirect) {
        # Try to compute indirect-effect time series. If unavailable
        # (e.g. uses_rho = FALSE), fall back to single-panel layout.
        ind <- tryCatch(.nasc_indirect_draws(self), error = function(e) NULL)
        if (is.null(ind)) {
          if (indirect_default) {
            warning("Indirect-effect draws unavailable; drawing direct effect only.")
          }
          indirect <- FALSE
        }
      }

      if (!indirect) {
        # Original single-panel direct-effect plot.
        .plot_tau(
          data       = private$plot_data,
          x          = time_name,
          y          = "tau",
          ymin       = "tau_LB",
          ymax       = "tau_UB",
          xintercept = private$intervention
        )
        return(invisible(NULL))
      }

      # ---- Two-panel layout: direct on top, total indirect below ----
      ci    <- private$ci_width
      probs <- .ci_probs(ci)
      delta_total <- ind$delta_total       # [n_draws x T_post]
      ind_summary <- tibble::tibble(
        !!time_name := ind$time_post,
        delta    = apply(delta_total, 2, mean),
        delta_LB = apply(delta_total, 2, \(z) stats::quantile(z, probs[1], names = FALSE)),
        delta_UB = apply(delta_total, 2, \(z) stats::quantile(z, probs[2], names = FALSE))
      )

      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op))
      graphics::par(mfrow = c(2, 1),
                    mar = c(2.5, 5, 1.5, 1),
                    oma = c(3, 0, 0, 0),
                    bty = "l")

      # ---- Top: direct effect (replicates .plot_tau, top panel) ----
      pd <- as.data.frame(private$plot_data)
      pd <- pd[order(pd[[time_name]]), , drop = FALSE]
      xv  <- pd[[time_name]]
      yv  <- pd$tau
      lbv <- pd$tau_LB
      ubv <- pd$tau_UB
      yrng <- range(c(yv, lbv, ubv), na.rm = TRUE)
      plot(xv, yv, type = "n", ylim = yrng, xlab = "", ylab = "tau (direct)",
           xaxt = "n")
      graphics::grid(lty = "dotted", col = "gray80")
      graphics::polygon(c(xv, rev(xv)), c(lbv, rev(ubv)),
                        col = grDevices::adjustcolor("steelblue", alpha.f = 0.2),
                        border = NA)
      graphics::lines(xv, yv, lwd = 2, col = "steelblue4")
      graphics::abline(v = private$intervention, lty = 2)
      graphics::abline(h = 0, lty = 1, col = "black")

      # ---- Bottom: total indirect effect ----
      xv2  <- ind_summary[[time_name]]
      yv2  <- ind_summary$delta
      lb2  <- ind_summary$delta_LB
      ub2  <- ind_summary$delta_UB

      # Pre-treatment is identically zero by construction (no spillover
      # before intervention); we plot only the post-treatment series and
      # share the x-range with the direct panel for visual alignment.
      xrng_top <- range(xv, na.rm = TRUE)
      yrng2 <- range(c(yv2, lb2, ub2, 0), na.rm = TRUE)

      plot(xv2, yv2, type = "n",
           xlim = xrng_top, ylim = yrng2,
           xlab = time_name,
           ylab = expression(delta ~ "(indirect, total)"))
      graphics::grid(lty = "dotted", col = "gray80")
      graphics::polygon(c(xv2, rev(xv2)), c(lb2, rev(ub2)),
                        col = grDevices::adjustcolor("indianred", alpha.f = 0.2),
                        border = NA)
      graphics::lines(xv2, yv2, lwd = 2, col = "indianred4")
      graphics::abline(v = private$intervention, lty = 2)
      graphics::abline(h = 0, lty = 1, col = "black")

      invisible(NULL)
    },

    #' @description
    #' Plot the estimated indirect (spillover) treatment effect over time.
    #' By Proposition 6.2 the per-donor spillover at post-period t is
    #' \eqn{\delta_{i,t}^{NASC} = s_i \cdot \tau_{1t}^{NASC}}. This method
    #' supports two views:
    #'
    #' \itemize{
    #'   \item \code{which = "total"} (default): draws the per-period
    #'         total spillover \eqn{\sum_i \delta_{i,t}} as a single
    #'         line with a credible-interval band. This is the natural
    #'         population analog of the direct-effect time series.
    #'   \item \code{which = "by_donor"}: draws one line per donor for
    #'         the per-period donor-specific spillover
    #'         \eqn{\delta_{i,t}}. Useful for diagnosing which donors
    #'         absorb most of the spillover; a legend is omitted when
    #'         the donor pool is large.
    #' }
    #'
    #' Requires the model to be configured with a network
    #' (\code{uses_rho = TRUE}, i.e. either \code{bias_correction = TRUE}
    #' or \code{nasc_penalty = TRUE} -- otherwise spillover is identically
    #' zero by construction and the plot is uninformative).
    #'
    #' @param which One of \code{"total"} or \code{"by_donor"}.
    #' @param top_n Integer or \code{NULL}. Only used when
    #'   \code{which = "by_donor"}: if not \code{NULL}, restrict to the
    #'   \code{top_n} donors with largest \eqn{|s_i|} (i.e. those that
    #'   carry most of the spillover). Default \code{NULL} (all donors).
    indirectEffectPlot = function(which = c("total", "by_donor"),
                                  top_n = NULL) {
      which <- match.arg(which)
      ind <- .nasc_indirect_draws(self)
      if (is.null(ind)) {
        stop("Indirect-effect plotting requires a model with a non-trivial ",
             "rho. Re-fit with bias_correction = TRUE or nasc_penalty = TRUE.")
      }

      time_name <- rlang::as_name(private$time)
      ci <- private$ci_width
      probs <- .ci_probs(ci)

      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op))
      graphics::par(bty = "l", mar = c(4, 5, 2, 1))

      if (which == "total") {
        delta_total <- ind$delta_total
        xv  <- ind$time_post
        yv  <- apply(delta_total, 2, mean)
        lb  <- apply(delta_total, 2, \(z) stats::quantile(z, probs[1], names = FALSE))
        ub  <- apply(delta_total, 2, \(z) stats::quantile(z, probs[2], names = FALSE))

        ord <- order(xv)
        xv <- xv[ord]; yv <- yv[ord]; lb <- lb[ord]; ub <- ub[ord]

        yrng <- range(c(yv, lb, ub, 0), na.rm = TRUE)
        plot(xv, yv, type = "n", ylim = yrng,
             xlab = time_name,
             ylab = expression(delta ~ "(indirect, total)"))
        graphics::grid(lty = "dotted", col = "gray80")
        graphics::polygon(c(xv, rev(xv)), c(lb, rev(ub)),
                          col = grDevices::adjustcolor("indianred", alpha.f = 0.2),
                          border = NA)
        graphics::lines(xv, yv, lwd = 2, col = "indianred4")
        graphics::abline(v = private$intervention, lty = 2)
        graphics::abline(h = 0, lty = 1, col = "black")

      } else {
        # by_donor: one line per donor, no CI band (would be unreadable
        # for J > ~3). Lines are coloured from a qualitative palette.
        delta_arr  <- ind$delta_arr             # [n_draws x T_post x J]
        donor_nm   <- ind$donor_names
        xv         <- ind$time_post

        # Per-period mean per donor: collapse the draws dimension.
        per_donor_mean <- apply(delta_arr, c(2L, 3L), mean)
        # Returns a [T_post x J] matrix.
        if (!is.matrix(per_donor_mean)) {
          per_donor_mean <- matrix(per_donor_mean,
                                   nrow = length(xv),
                                   ncol = length(donor_nm))
        }
        colnames(per_donor_mean) <- donor_nm

        keep <- seq_len(ncol(per_donor_mean))
        if (!is.null(top_n)) {
          stopifnot(is.numeric(top_n), length(top_n) == 1L, top_n >= 1)
          # Rank by mean |delta| over time as a proxy for "absorbs
          # most spillover".
          rank_score <- apply(abs(per_donor_mean), 2, mean)
          keep <- order(rank_score, decreasing = TRUE)[
            seq_len(min(as.integer(top_n), length(rank_score)))
          ]
        }
        per_donor_mean <- per_donor_mean[, keep, drop = FALSE]
        donor_show     <- donor_nm[keep]
        n_show         <- length(donor_show)

        ord <- order(xv)
        xv <- xv[ord]
        per_donor_mean <- per_donor_mean[ord, , drop = FALSE]

        yrng <- range(c(as.numeric(per_donor_mean), 0), na.rm = TRUE)
        plot(xv, per_donor_mean[, 1], type = "n", ylim = yrng,
             xlab = time_name,
             ylab = expression(delta[i] ~ "(indirect, per donor)"))
        graphics::grid(lty = "dotted", col = "gray80")
        graphics::abline(v = private$intervention, lty = 2)
        graphics::abline(h = 0, lty = 1, col = "black")

        cols <- grDevices::hcl.colors(max(n_show, 2L), palette = "Dark 3")[seq_len(n_show)]
        for (k in seq_len(n_show)) {
          graphics::lines(xv, per_donor_mean[, k], col = cols[k], lwd = 1.6)
        }

        # Legend only for tractable donor counts; otherwise label most
        # extreme line at the right edge to avoid an unreadable legend.
        if (n_show <= 12L) {
          graphics::legend("topleft",
                           legend = donor_show,
                           col    = cols,
                           lty    = 1, lwd = 1.6,
                           bg     = grDevices::adjustcolor("white", alpha.f = 0.85),
                           box.col = "gray70",
                           cex    = 0.8,
                           ncol   = if (n_show > 6L) 2L else 1L)
        }
      }

      invisible(NULL)
    },

    #' @description
    #' Posterior density plots of the main scalar Bayesian parameters that the
    #' model actually estimated. Depending on the configuration, this may
    #' include:
    #' \itemize{
    #'   \item \code{rho}             -- Step-1 spatial autocorrelation (SAR/SDM).
    #'   \item \code{theta[k]}        -- Step-1 SDM spillover coefficients
    #'                                   (one panel per covariate).
    #'   \item \code{sigma_step1}     -- Step-1 residual SD (SAR or SDM).
    #'   \item \code{lambda}          -- Step-2 NASC penalty strength
    #'                                   (when \code{nasc_penalty = TRUE}).
    #'   \item \code{sigma_step2}     -- Step-2 residual SD on the standardized
    #'                                   outcome scale.
    #'   \item \code{bias_correction} -- Posterior of the multiplicative bias
    #'                                   factor \eqn{1/(1 - \langle w, s\rangle)}
    #'                                   (when \code{bias_correction = TRUE}).
    #'                                   Skipped when bias correction is off,
    #'                                   since draws are constant at 1.
    #' }
    #' Each parameter is shown in its own panel as a color-filled density curve.
    posteriorPlot = function() {
      if (is.null(private$fitted) && is.null(private$y_synth_draws)) {
        stop("Run $fit() before calling posteriorPlot().")
      }

      # ----------------------------------------------------------------
      # Collect (label, draws) pairs from wherever each parameter lives.
      # ----------------------------------------------------------------
      panels <- list()

      add_panel <- function(label, draws) {
        if (is.null(draws)) return(invisible(NULL))
        x <- as.numeric(draws)
        x <- x[is.finite(x)]
        if (length(x) < 2L || stats::sd(x) == 0) return(invisible(NULL))
        panels[[length(panels) + 1L]] <<- list(label = label, draws = x)
      }

      # Step-1 scalars (rho, sigma_step1) and theta vector live in private$fitted
      # whenever Step 1 ran. private$fitted may also be a Step-2 fit when Step 1
      # was skipped, so we extract defensively rather than by branch.
      step1_draws <- if (!is.null(private$fitted)) {
        tryCatch(
          rstan::extract(private$fitted, permuted = TRUE),
          error = function(e) list()
        )
      } else {
        list()
      }

      if (!is.null(step1_draws$rho)) {
        add_panel("rho", step1_draws$rho)
      }

      # theta is a matrix (n_draws x K_pred) in the SDM model. Plot one density
      # per component so the visualization is meaningful.
      if (!is.null(step1_draws$theta)) {
        th <- step1_draws$theta
        if (is.matrix(th)) {
          for (k in seq_len(ncol(th))) {
            add_panel(sprintf("theta[%d]", k), th[, k])
          }
        } else {
          add_panel("theta", th)
        }
      }

      # Step-1 sigma is named differently per model: sigma_sar vs sigma_sdm.
      # Surface either under a unified label.
      if (!is.null(step1_draws$sigma_sar)) {
        add_panel("sigma_step1", step1_draws$sigma_sar)
      } else if (!is.null(step1_draws$sigma_sdm)) {
        add_panel("sigma_step1", step1_draws$sigma_sdm)
      }

      # Step-2 scalars (lambda, sigma_step2, bias_correction) live in
      # private$y_synth_draws. lambda is only present when nasc_penalty = TRUE
      # (stan_2_NASC). bias_correction is a generated quantity equal to
      # 1 / (1 - <w, s>) when bias correction is on, and the constant 1.0 when
      # it is off; the constant case is auto-skipped by add_panel's sd > 0 guard.
      if (!is.null(private$y_synth_draws)) {
        if (!is.null(private$y_synth_draws$lambda)) {
          add_panel("lambda", private$y_synth_draws$lambda)
        }
        # sigma_sc (penalty model) or sigma (model1); whichever is present.
        if (!is.null(private$y_synth_draws$sigma_sc)) {
          add_panel("sigma_step2", private$y_synth_draws$sigma_sc)
        } else if (!is.null(private$y_synth_draws$sigma)) {
          add_panel("sigma_step2", private$y_synth_draws$sigma)
        }
        if (!is.null(private$y_synth_draws$bias_correction)) {
          add_panel("bias_correction", private$y_synth_draws$bias_correction)
        }
      }

      if (length(panels) == 0L) {
        stop("No estimated scalar parameters available to plot.")
      }

      # ----------------------------------------------------------------
      # Layout: stack vertically; for many panels, switch to a 2-column
      # grid so densities stay readable.
      # ----------------------------------------------------------------
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

    #' @description
    #' Posterior density plot of the average treatment effect on the treated
    #' (ATT), pooled across all post-treatment periods. When the model uses
    #' a network, a second panel is added showing the density of the
    #' average total indirect (spillover) effect, defined as the per-draw
    #' average across post-periods of \eqn{\sum_i \delta_{i,t}^{NASC}}
    #' (the spillover analog of ATT, Proposition 6.2).
    #'
    #' @param indirect Logical. If \code{TRUE} (default when the model uses
    #'   a network), append a second panel for the average indirect effect.
    attPlot = function(indirect = NULL) {
      if (is.null(private$y_synth_draws)) stop("Run $fit() before calling attPlot.")

      # Reconstruct tau draws
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

      avg_total_draws <- NULL
      if (indirect) {
        ind <- tryCatch(.nasc_indirect_draws(self), error = function(e) NULL)
        if (!is.null(ind)) {
          avg_total_draws <- ind$avg_total
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
        d_ind <- stats::density(avg_total_draws, na.rm = TRUE)
        plot(d_ind, main = "Average indirect effect (total)",
             xlab = expression(bar(delta)^{total}), ylab = "density")
        graphics::polygon(d_ind,
                          col = grDevices::adjustcolor("indianred", alpha.f = 0.3),
                          border = "indianred")
        graphics::abline(v = 0, lty = 3, col = "gray40")
      }

      invisible(NULL)
    },

    #' @description
    #' Posterior density plots of the period-by-period treatment effects
    #' for all post-treatment periods. When the model uses a network (i.e.
    #' \code{uses_rho = TRUE}), two columns are drawn: the left column is
    #' the direct effect \eqn{\tau_{1t}} and the right column is the
    #' total indirect (spillover) effect
    #' \eqn{\delta_t^{total} = \sum_{i=2}^{J+1}\delta_{i,t}^{NASC}}, with
    #' one row per post period. When no network is in use, only the
    #' direct-effect column is drawn (single-column layout, preserving
    #' previous behavior). Panels within a column share a common x-axis;
    #' the two columns use independent x-axes since direct and indirect
    #' effects are typically on very different scales.
    #'
    #' @param indirect Logical. If \code{TRUE} (default when the model
    #'   uses a network), also draw the per-period total indirect-effect
    #'   density column. Set to \code{FALSE} to suppress the indirect
    #'   panels even when a network is in use.
    tauPlot = function(indirect = NULL) {
      if (is.null(private$y_synth_draws)) stop("Run $fit() before calling tauPlot.")

      # Reconstruct tau draws (same logic as before).
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

      # ----------------------------------------------------------------
      # Decide whether to draw the indirect-effect column.
      #
      # Default: ON when the model uses a network and we can compute
      # spillover draws. Off otherwise. The user can force-disable it
      # via indirect = FALSE.
      # ----------------------------------------------------------------
      indirect_default <- isTRUE(private$uses_rho)
      if (is.null(indirect)) indirect <- indirect_default
      stopifnot(is.logical(indirect), length(indirect) == 1L)

      delta_total_draws <- NULL
      if (indirect) {
        ind <- tryCatch(.nasc_indirect_draws(self), error = function(e) NULL)
        if (!is.null(ind)) {
          delta_total_draws <- ind$delta_total      # [n_draws x T_post]
        } else if (indirect_default) {
          # Only warn when the user (or the default) actually expected
          # indirect draws. If the model has no rho, falling back to the
          # single-column layout silently is the desired behavior.
          warning("Indirect-effect draws unavailable; drawing direct effect only.")
        }
        if (is.null(delta_total_draws)) indirect <- FALSE
      }

      # Direct-effect densities (always drawn)
      dens_dir <- lapply(seq_len(n_t),
                         function(i) stats::density(tau_draws[, i], na.rm = TRUE))
      xr_dir <- range(sapply(dens_dir, function(d) d$x))

      # Indirect-effect densities (drawn only if indirect = TRUE)
      dens_ind <- NULL
      xr_ind   <- NULL
      if (indirect) {
        dens_ind <- lapply(seq_len(n_t),
                           function(i) stats::density(delta_total_draws[, i],
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

      # When mfrow = c(n_t, 2), R fills row-by-row, so for each period i
      # we draw the direct panel first, then the indirect panel. Column
      # titles are placed via mtext() in the outer top margin.
      for (i in seq_len(n_t)) {
        is_last <- i == n_t

        # Direct
        d_i <- dens_dir[[i]]
        plot(d_i,
             main = paste0("Period ", time_post[i]),
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
               main = paste0("Period ", time_post[i]),
               xlab = if (is_last) expression(delta ~ "(indirect, total)") else "",
               ylab = "density",
               xlim = xr_ind,
               xaxt = if (is_last) "s" else "n")
          graphics::polygon(d_j,
                            col = grDevices::adjustcolor("indianred", alpha.f = 0.3),
                            border = "indianred")
          graphics::abline(v = 0, lty = 3, col = "gray40")
        }
      }

      # Column headings in the outer top margin
      if (indirect) {
        graphics::mtext("Direct effect",   side = 3, outer = TRUE,
                        at = 0.25, line = 0.6, font = 2)
        graphics::mtext("Indirect effect (total)", side = 3, outer = TRUE,
                        at = 0.75, line = 0.6, font = 2)
      } else {
        graphics::mtext("",
                        side = 3, outer = TRUE, line = 0.5, font = 2)
      }

      invisible(NULL)
    },

    #' @description
    #' Ridgeline plot of the posterior weight density per donor. Densities
    #' are stacked vertically, each labelled by donor name, and adjacent
    #' ridges overlap by a configurable fraction so that each individual
    #' density can be drawn taller without compressing the layout.
    #'
    #' @param overlap Numeric in \code{[0, 1)}. Fraction of vertical
    #'   overlap between adjacent ridges. \code{0} reproduces the old
    #'   non-overlapping layout; \code{0.5} (default) makes each ridge
    #'   reach halfway into the row above; values close to \code{1}
    #'   approach a fully-stacked look. Default \code{0.5}.
    #' @param scale Numeric > 0. Multiplicative height of every ridge
    #'   relative to the baseline step. Combined with \code{overlap} this
    #'   controls how prominent each density is. Default \code{1.4}.
    #' @param fill_alpha Numeric in \code{(0, 1]}. Multiplier on the
    #'   default fill transparency (\eqn{\approx 0.4} alpha at
    #'   \code{fill_alpha = 1}). Lower values make overlapping ridges
    #'   more transparent. Default \code{0.85}, matching the look of
    #'   \code{posteriorPlot()}.
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

      # Use the canonical donor ordering stored at fit time. Falling back
      # to setdiff(levels(id), treated_id) silently scrambles labels when
      # the long data isn't sorted by id.
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

      # Per-donor density. We compute on a common x-grid so the shapes
      # are directly comparable; densities are NOT renormalized per ridge,
      # so a peaked donor is visibly taller than a flat one.
      dens <- lapply(seq_len(n), function(i) {
        stats::density(w_mat[, i], na.rm = TRUE)
      })

      x_range <- range(unlist(lapply(dens, function(d) d$x)))
      max_y   <- max(vapply(dens, function(d) max(d$y), numeric(1)))
      if (!is.finite(max_y) || max_y <= 0) max_y <- 1  # degenerate safety

      # ---------------------------------------------------------------
      # Layout geometry.
      #
      # `step` is the vertical distance between adjacent donor baselines.
      # `ridge_h` is the maximum drawn height of any single density.
      #
      # The relationship overlap = 1 - step/ridge_h means:
      #   overlap = 0   -> step = ridge_h     (no overlap; old behaviour)
      #   overlap = 0.5 -> step = 0.5*ridge_h (half-overlap; ridges
      #                                        reach halfway into row above)
      #   overlap -> 1  -> step -> 0          (fully stacked; not useful)
      #
      # `scale` then expands ridge_h relative to the implied baseline,
      # letting the user make ridges taller without changing the spacing
      # logic. We anchor on max_y so the tallest density gets exactly
      # `scale * max_y` of vertical room.
      # ---------------------------------------------------------------
      ridge_h <- scale * max_y
      step    <- ridge_h * (1 - overlap)
      # Top of the topmost ridge sits at step*n + ridge_h; pad a little.
      ylim_top <- step * n + ridge_h * 1.05
      # Bottom: a small negative cushion so the lowest baseline isn't
      # flush with the x-axis frame.
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

      # Draw from TOP donor down to BOTTOM donor. In a stacked plot the
      # last polygon drawn sits in front of earlier ones, so to get the
      # canonical ridgeline look (lower rows appearing in front), the
      # bottom ridge must be drawn last. seq(n, 1) does exactly that.
      for (i in seq(n, 1L, by = -1L)) {
        d <- dens[[i]]
        baseline <- step * i
        # Each density is scaled to ridge_h at its own peak * (this peak /
        # global max), preserving relative heights across donors.
        y <- baseline + d$y * (ridge_h / max_y)
        graphics::polygon(
          x = c(d$x, rev(d$x)),
          y = c(y, rep(baseline, length(d$x))),
          col    = cols[i],
          border = border
        )
        # Thin baseline rule for visual anchoring.
        graphics::segments(x_range[1], baseline, x_range[2], baseline,
                           col = "gray60", lwd = 0.5)
      }

      invisible(NULL)
    },

    #' @description
    #' Plot correlations between weights across draws.
    weightCorr = function() {
      if (is.null(private$fitted)) stop("Run $fit() before calling weightCorr().")

      w_mat <- private$y_synth_draws$w

      # Use canonical donor ordering (see weightDraws() for rationale).
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

      # Diverging palette red -> white -> green, midpoint 0
      pal <- grDevices::colorRampPalette(c("red", "white", "green"))(101)
      # Plot rows top-to-bottom to mirror typical heatmap orientation
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

      # Simple color legend on the right
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
