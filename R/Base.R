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
                          covariates = NULL,
                          W = NULL, spatial_model = "none") {

      stopifnot(ci_width > 0 & ci_width < 1)

      if (!spatial_model %in% c("none", "SAR", "SDM")) {
        stop("spatial_model must be 'none', 'SAR', or 'SDM'.")
      }
      if (spatial_model %in% c("SAR", "SDM") && is.null(W)) {
        stop("A spatial weights matrix 'W' is required for SAR/SDM models.")
      }
      if (!is.null(W)) {
        W_mat <- as.matrix(W)
        rsums <- rowSums(W_mat, na.rm = TRUE)
        is_row_standardized <- all(abs(rsums - 1) < 1e-6 | abs(rsums) < 1e-6)
        if (!is_row_standardized) {
          stop("The spatial weights matrix 'W' must be row-standardized (row sums must equal 1, or 0 for isolated units).")
        }
      }

      private$time    <- rlang::enquo(time)
      private$id      <- rlang::enquo(id)
      private$treated <- rlang::enquo(treated)
      private$outcome <- rlang::enquo(outcome)
      private$ci_width     <- ci_width
      private$covariates   <- covariates
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

        # --- STRATIFIED SAMPLING ACROSS CHAINS ---
        # permuted = FALSE returns an array of [iterations, chains, parameters]
        rho_array <- rstan::extract(private$fitted, pars = "rho", permuted = FALSE)

        # Get dimensions
        n_iters  <- dim(rho_array)[1]
        n_chains <- dim(rho_array)[2]

        # Calculate exactly how many samples to draw per chain
        samples_per_chain <- ceiling(n_samples / n_chains)
        sampled_rhos <- numeric(0)

        # Draw evenly across all chains
        for (c_idx in seq_len(n_chains)) {
          # Extract the rho vector for this specific chain
          chain_draws <- rho_array[, c_idx, 1]

          # Sample without replacement
          take_n <- min(samples_per_chain, n_iters)
          sampled_rhos <- c(sampled_rhos, sample(chain_draws, size = take_n, replace = FALSE))
        }

        # If ceiling() caused slight over-sampling (e.g., 100/3 = 34 * 3 = 102),
        # trim randomly to exactly n_samples
        if (length(sampled_rhos) > n_samples) {
          sampled_rhos <- sample(sampled_rhos, size = n_samples, replace = FALSE)
        }
        # -----------------------------------------

        message(sprintf("Running Step 2: Parallel NASC across %d cores using furrr...", cores))

        # Save old plan and ensure it gets restored on exit
        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)

        # Set up future plan
        future::plan(future::multisession, workers = cores)

        # Worker: now extracts ALL relevant quantities (w, lambda, bias_correction)
        run_nasc_worker <- function(single_rho, base_data, step2_mod) {
          # Silently load rstan on the worker
          suppressPackageStartupMessages(require(rstan, quietly = TRUE))

          worker_data        <- base_data
          worker_data$rho    <- single_rho

          # Silently run the MCMC sampling
          fit_step2 <- rstan::sampling(
              object        = step2_mod,
              data          = worker_data,
              chains        = 1,
              iter          = 1000,
              warmup        = 500,
              refresh       = 0,
              show_messages = FALSE
            )


          draws <- rstan::extract(
            fit_step2,
            pars = c("y_counterfactual", "y_sim_pre", "w", "lambda",
                     "sigma_sc", "bias_correction")
          )

          list(
            y_counterfactual = draws$y_counterfactual,
            y_sim_pre        = draws$y_sim_pre,
            w                = draws$w,
            lambda           = draws$lambda,
            sigma_sc         = draws$sigma_sc,
            bias_correction  = draws$bias_correction,
            rho_used         = single_rho
          )
        }

        # Activate global handlers so the progress bar shows in the console
        progressr::handlers(global = TRUE)

        # Wrap the future_map call in with_progress
        results_list <- progressr::with_progress({
          # Initialize the progressor with the total number of steps
          p <- progressr::progressor(steps = length(sampled_rhos))

          furrr::future_map(
            sampled_rhos,
            function(rho) {
              # 1. Run the worker
              res <- run_nasc_worker(
                single_rho = rho,
                base_data  = stan_data,
                step2_mod  = private$stan_model$step2
              )
              # 2. Signal the progress bar that this step is done
              p()
              # 3. Return the result
              return(res)
            },
            .options = furrr::furrr_options(seed = TRUE, packages = "rstan")
          )
        })

        # Aggregate
        final_y_cf            <- do.call(rbind, lapply(results_list, \(x) x$y_counterfactual))
        final_y_sim_pre       <- do.call(rbind, lapply(results_list, \(x) x$y_sim_pre))
        final_w               <- do.call(rbind, lapply(results_list, \(x) x$w))
        final_lambda          <- do.call(c,     lapply(results_list, \(x) x$lambda))
        final_sigma_sc        <- do.call(c,     lapply(results_list, \(x) x$sigma_sc))
        final_bias_correction <- do.call(c,     lapply(results_list, \(x) x$bias_correction))
        rhos_used             <- vapply(results_list, \(x) x$rho_used, numeric(1))

        private$y_synth_draws <- list(
          y_counterfactual= final_y_cf, # <--- Changed key and value
          y_sim_pre       = final_y_sim_pre,
          w               = final_w,
          lambda          = final_lambda,
          sigma_sc        = final_sigma_sc,
          bias_correction = final_bias_correction,
          rhos_used       = rhos_used
        )

        private$plot_data <- .get_nasc_results(
          y_counterfactual_draws = final_y_cf,            # <--- Update parameter names
          bias_correction_draws  = final_bias_correction, # <--- Add bias draw parameter
          y_sim_pre_draws        = final_y_sim_pre,
          pre_data               = pre_data,
          post_data              = post_data,
          time                   = private$time,
          outcome                = private$outcome,
          ci                     = private$ci_width
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
          y_counterfactual_draws = private$y_synth_draws$y_counterfactual,
          bias_correction_draws  = private$y_synth_draws$bias_correction,
          y_sim_pre_draws        = private$y_synth_draws$y_sim_pre,
          pre_data               = pre_data,
          post_data              = post_data,
          time                   = private$time,
          outcome                = private$outcome,
          ci                     = private$ci_width
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
    },

    #' @description
    #' Plot implicit weight distribution across draws.
    #' @return ggplot object with weight distribution per unit.
    weightDraws = function() {
      if (is.null(private$fitted)) stop("Run fit() first.")

      # Extract weights depending on the model type
      if (private$spatial_model %in% c("SAR", "SDM")) {
        w_mat <- private$y_synth_draws$w
      } else {
        w_mat <- rstan::extract(private$fitted, pars = "w")$w
      }

      # Get donor names (all unique IDs minus the treated ID)
      treated_id <- as.character(private$treated_ids)
      all_ids <- levels(private$data[[rlang::as_name(private$id)]])
      donor_names <- setdiff(all_ids, treated_id)

      if (ncol(w_mat) != length(donor_names)) {
        stop("Mismatch between the number of weights and donor names.")
      }

      colnames(w_mat) <- donor_names

      # Reshape using modern tidyr::pivot_longer
      melt_w <- as.data.frame(w_mat) |>
        tidyr::pivot_longer(
          cols      = dplyr::everything(),
          names_to  = "ID",
          values_to = "weight"
        )

      # Generate Ridge Plot
      melt_w |>
        ggplot2::ggplot(ggplot2::aes(x = weight, y = ID, fill = ID)) +
        ggridges::geom_density_ridges(alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = "Donor Weight", y = "Donor Unit")
    },

    #' @description
    #' Plots correlations between weights across draws.
    #' @return ggplot heatmap object with correlations.
    weightCorr = function() {
      if (is.null(private$fitted)) stop("Run fit() first.")

      # Extract weights depending on the model type
      if (private$spatial_model %in% c("SAR", "SDM")) {
        w_mat <- private$y_synth_draws$w
      } else {
        w_mat <- rstan::extract(private$fitted, pars = "w")$w
      }

      # Get donor names
      treated_id <- as.character(private$treated_ids)
      all_ids <- levels(private$data[[rlang::as_name(private$id)]])
      donor_names <- setdiff(all_ids, treated_id)

      colnames(w_mat) <- donor_names

      # Compute correlation matrix
      cormat <- round(cor(w_mat), 3)
      diag(cormat) <- NA # Set diagonal to NA so it doesn't skew the color scale

      # Reshape without depending on reshape2::melt
      melted_cormat <- as.data.frame(cormat) |>
        dplyr::mutate(X1 = factor(rownames(cormat), levels = donor_names)) |>
        tidyr::pivot_longer(
          cols      = -X1,
          names_to  = "X2",
          values_to = "value"
        ) |>
        dplyr::mutate(X2 = factor(X2, levels = donor_names))

      # Generate Heatmap
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
