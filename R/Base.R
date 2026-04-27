# Bayesian NASC estimator
# BASE code

#' Bayesian Network-Aware Synthetic Control
#'
#' @description
#' An R6 class implementing Bayesian Synthetic Control estimators with optional
#' network/spatial awareness. Supports predictor matching, covariate adjustment,
#' and SAR/SDM spatial models via Stan.
#'
#' @export
nascSynth <- R6::R6Class(
  classname = "nascSynth",
  private = list(
    data = NULL,
    covariates = NULL,
    vs = NULL,
    W = NULL,
    spatial_model = NULL,
    time = NULL,
    id = NULL,
    treated = NULL,
    outcome = NULL,
    ci_width = NULL,
    intervention = NULL,
    fitted = NULL,
    plot_data = NULL,
    time_tiles_cache = NULL,
    stan_model = NULL,
    y_synth_draws = NULL,
    lift_draws = NULL,
    treated_ids = NULL,
    mcmc_checks = NULL
  ),
  active = list(
    #' @field timeTiles ggplot2 tile chart of treatment status over time.
    timeTiles = function() { return(private$time_tiles_cache) },

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
    },

    #' @field checks MCMC diagnostics from the fitted Stan model.
    checks = function() { private$mcmc_checks },

    #' @field lift Lift draws (if computed).
    lift = function() { private$lift_draws }
  ),
  public = list(

    #' @description
    #' Create a new nascSynth object.
    initialize = function(data, time, id, treated, outcome, ci_width = 0.75,
                          covariates = NULL, vs = NULL,
                          W = NULL, spatial_model = "none") {

      stopifnot(ci_width > 0 & ci_width < 1)

      if (!spatial_model %in% c("none", "SAR", "SDM")) {
        stop("spatial_model must be 'none', 'SAR', or 'SDM'.")
      }
      if (spatial_model %in% c("SAR", "SDM") && is.null(W)) {
        stop("A spatial weights matrix 'W' is required for SAR/SDM models.")
      }

      private$time    <- rlang::enquo(time)
      private$id      <- rlang::enquo(id)
      private$treated <- rlang::enquo(treated)
      private$outcome <- rlang::enquo(outcome)
      private$ci_width     <- ci_width
      private$covariates   <- covariates
      private$vs           <- vs
      private$W            <- W
      private$spatial_model <- spatial_model

      message(sprintf(
        "Transforming data for %s model",
        ifelse(
          spatial_model == "none",
          "standard SC predictor-match",
          paste("NASC", spatial_model)
        )
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

      private$time_tiles_cache <- .time_tiles(
        data   = private$data,
        time   = private$time,
        id     = private$id,
        status = rlang::quo(status)
      )

      # Build model lists for the two-step pipeline
      if (spatial_model == "SAR") {
        private$stan_model <- list(
          step1 = stanmodels$stan_1_SAR_NASC,
          step2 = stanmodels$stan_2_NASC
        )
      } else if (spatial_model == "SDM") {
        private$stan_model <- list(
          step1 = stanmodels$stan_1_SDM_NASC,
          step2 = stanmodels$stan_2_NASC
        )
      } else {
        private$stan_model <- stanmodels$model1
      }
    },

    #' @description
    #' Fit the Stan model via MCMC.
    #' @param n_samples Number of spatial draws to pass to NASC (default 100).
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

      if (!is.null(private$covariates)) {
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

      if (is.null(private$vs)) {
        private$vs <- rep(1, nrow(X))
      }

      if (private$spatial_model %in% c("SAR", "SDM")) {

        # ---- Build covariate matrices from private$covariates ONLY ----
        # X0: K_pred x J (donor predictors), X1: K_pred-vector (treated predictors)
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

        donor_ids   <- colnames(X)                      # X is the wide outcome panel
        treated_id  <- as.character(private$treated_ids)

        # Predictor matrices: rows = predictors, cols = units
        X0_mat <- as.matrix(cov_wide[, donor_ids,  drop = FALSE])
        X1_vec <- as.numeric(cov_wide[[treated_id]])
        K_pred <- length(pred_names)

        if (length(X1_vec) != K_pred) {
          stop("Treated unit predictor vector length does not match K_pred.")
        }
        if (anyNA(X0_mat) || anyNA(X1_vec)) {
          stop("Missing values found in predictor matrix; please impute or drop.")
        }

        # ---- Predictor importance weights ----
        if (is.null(private$vs)) {
          private$vs <- rep(1, K_pred)
        } else if (length(private$vs) != K_pred) {
          stop(sprintf("Length of 'vs' (%d) must equal number of predictors (%d).",
                       length(private$vs), K_pred))
        }
        if (any(private$vs < 0)) stop("'vs' must be non-negative.")

        # ---- Spatial weights matrix ----
        if (is.null(rownames(private$W)) || is.null(colnames(private$W))) {
          all_ids <- levels(private$data[[rlang::as_name(private$id)]])
          rownames(private$W) <- colnames(private$W) <- all_ids
        }
        W_ordered <- private$W[c(donor_ids, treated_id), c(donor_ids, treated_id)]
        J <- length(donor_ids)

        # ---- W eigenvalue check ----
        ev <- eigen(W_ordered, only.values = TRUE)$values
        if (max(abs(Im(ev))) > 1e-8) {
          warning("W has non-trivial complex eigenvalues; SAR/SDM identification ",
                  "may be unreliable.")
        }
        lambda_W <- Re(ev)

        # ---- Stan data ----
        stan_data <- list(
          K_pred   = K_pred,
          X1       = X1_vec,
          J        = J,
          X0       = X0_mat,
          vs       = as.numeric(private$vs),
          T0       = nrow(pre_data),
          Y_panel  = t(as.matrix(
            pre_data |> dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
          )),
          W        = W_ordered,
          lambda_W = lambda_W,
          W_J      = W_ordered[1:J, 1:J],
          w_J1     = as.vector(W_ordered[1:J, J + 1]),
          T_post   = nrow(post_data),
          Y0_post  = as.matrix(X_pred[, donor_ids, drop = FALSE]),
          Y1_post  = post_data |> dplyr::pull(!!private$outcome)
        )

        # ---------------------------------------------------------------------
        # TWO-STEP PARALLEL EXECUTION
        # ---------------------------------------------------------------------
        message(sprintf("Running Step 1: Estimating %s parameters...",
                        private$spatial_model))
        private$fitted <- rstan::sampling(
          private$stan_model$step1,
          data  = stan_data,
          cores = cores,
          ...
        )

        rho_draws    <- rstan::extract(private$fitted, pars = "rho")$rho
        sampled_rhos <- sample(rho_draws, size = min(n_samples, length(rho_draws)),
                               replace = FALSE)

        message(sprintf("Running Step 2: Parallel NASC across %d cores...", cores))
        cl <- parallel::makeCluster(
          cores,
          type = if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
        )
        on.exit(parallel::stopCluster(cl), add = TRUE)

        # Worker: now extracts ALL relevant quantities (w, lambda, bias_correction)
        run_nasc_worker <- function(single_rho, base_data, step2_mod) {
          require(rstan)
          worker_data        <- base_data
          worker_data$rho    <- single_rho

          fit_step2 <- rstan::sampling(
            object  = step2_mod,
            data    = worker_data,
            chains  = 1,
            iter    = 1000,
            warmup  = 500,
            refresh = 0
          )
          draws <- rstan::extract(
            fit_step2,
            pars = c("tau_nasc", "y_sim_pre", "w", "lambda",
                     "sigma_sc", "bias_correction")
          )
          list(
            tau_nasc        = draws$tau_nasc,
            y_sim_pre       = draws$y_sim_pre,
            w               = draws$w,
            lambda          = draws$lambda,
            sigma_sc        = draws$sigma_sc,
            bias_correction = draws$bias_correction,
            rho_used        = single_rho
          )
        }

        results_list <- parallel::parLapply(
          cl, sampled_rhos, run_nasc_worker,
          base_data = stan_data,
          step2_mod = private$stan_model$step2
        )

        # Aggregate
        final_tau_nasc        <- do.call(rbind, lapply(results_list, \(x) x$tau_nasc))
        final_y_sim_pre       <- do.call(rbind, lapply(results_list, \(x) x$y_sim_pre))
        final_w               <- do.call(rbind, lapply(results_list, \(x) x$w))
        final_lambda          <- do.call(c,     lapply(results_list, \(x) x$lambda))
        final_sigma_sc        <- do.call(c,     lapply(results_list, \(x) x$sigma_sc))
        final_bias_correction <- do.call(c,     lapply(results_list, \(x) x$bias_correction))
        rhos_used             <- vapply(results_list, \(x) x$rho_used, numeric(1))

        private$y_synth_draws <- list(
          tau_nasc        = final_tau_nasc,
          y_sim_pre       = final_y_sim_pre,
          w               = final_w,
          lambda          = final_lambda,
          sigma_sc        = final_sigma_sc,
          bias_correction = final_bias_correction,
          rhos_used       = rhos_used
        )

        private$plot_data <- .get_nasc_results(
          tau_draws       = final_tau_nasc,
          y_sim_pre_draws = final_y_sim_pre,
          pre_data        = pre_data,
          post_data       = post_data,
          time            = private$time,
          outcome         = private$outcome,
          ci              = private$ci_width
        )

      } else {
        # Standard Model (non-spatial) FIX
        # Map the R variables to the exact names expected by model1.stan
        stan_data <- list(
          N       = nrow(X),       # pre-intervention periods (+ covariates)
          y       = X1,            # pre-intervention outcome for treated
          K       = ncol(X),       # number of donors
          X       = as.matrix(X),  # pre-intervention outcome for donors
          N_pred  = nrow(X_pred),  # post-intervention periods
          X_pred  = as.matrix(X_pred) # post-intervention outcome for donors
        )

        private$fitted <- rstan::sampling(
          private$stan_model,
          data = stan_data,
          cores = cores,
          ...
        )

        # FIX: Use .get_synth_draws instead of .get_synth_draws_predictor_match
        # because model1.stan outputs `y_sim` and `y_pred`
        private$y_synth_draws <- .get_synth_draws(
          fit       = private$fitted,
          pre_data  = pre_data,
          post_data = post_data,
          time      = private$time,
          outcome   = private$outcome
        )

        private$plot_data <- .get_plot_df(
          y_synth_draws = private$y_synth_draws,
          pre_data      = pre_data,
          post_data     = post_data,
          ci            = private$ci_width,
          time          = private$time,
          outcome       = private$outcome
        )
      }
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

      # Branch depending on whether y_synth_draws is a list (spatial) or df (standard)
      if (private$spatial_model %in% c("SAR", "SDM")) {
        private$plot_data <- .get_nasc_results(
          tau_draws       = private$y_synth_draws$tau_nasc,
          y_sim_pre_draws = private$y_synth_draws$y_sim_pre,
          pre_data        = pre_data,
          post_data       = post_data,
          time            = private$time,
          outcome         = private$outcome,
          ci              = private$ci_width
        )
      } else {
        private$plot_data <- .get_plot_df(
          y_synth_draws = private$y_synth_draws,
          pre_data      = pre_data,
          post_data     = post_data,
          ci            = private$ci_width,
          time          = private$time,
          outcome       = private$outcome
        )
      }
    },

    #' @description
    #' Summarise the posterior lift (requires lift draws to have been computed).
    summarizeLift = function() {
      if (is.null(private$lift_draws)) {
        stop("Run liftDraws() first.")
      }
      c(
        point       = mean(private$lift_draws$lift),
        lower_bound = stats::quantile(
          private$lift_draws$lift, (1 - private$ci_width) / 2
        ),
        upper_bound = stats::quantile(
          private$lift_draws$lift, 1 - (1 - private$ci_width) / 2
        )
      )
    },

    #' @description
    #' Plot the estimated treatment effect (tau) over time.
    effectPlot = function(facet = TRUE, subset = NULL) {
      .plot_tau(
        data       = private$plot_data,
        x          = private$time,
        y          = tau,
        ymin       = tau_LB,
        ymax       = tau_UB,
        xintercept = private$intervention
      )
    }
  )
)
