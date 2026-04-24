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
    covariates = NULL,     # long-format data frame: time, id, covariate cols
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
    #'
    #' @param data A long-format data frame with one row per unit per time
    #'   period.
    #' @param time Unquoted name of the time column.
    #' @param id Unquoted name of the unit identifier column.
    #' @param treated Unquoted name of the binary treatment indicator (1/0).
    #' @param outcome Unquoted name of the outcome column.
    #' @param ci_width Credible interval width, between 0 and 1. Default 0.75.
    #' @param covariates Optional long-format data frame of time-varying
    #'   covariates (must contain the same `time` and `id` columns as `data`).
    #'   Pre-treatment means are used as additional matching predictors.
    #' @param vs Optional numeric vector of predictor importance weights.
    #'   Length must equal the number of pre-treatment periods plus the number
    #'   of covariate columns if `covariates` is supplied.
    #' @param W Optional spatial weights matrix (required for SAR/SDM models).
    #'   Must have row/column names matching the unit IDs.
    #' @param spatial_model One of `"none"` (default), `"SAR"`, or `"SDM"`.
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

      # Build treatment tile chart (pass quosures directly)
      private$time_tiles_cache <- .time_tiles(
        data   = private$data,
        time   = private$time,
        id     = private$id,
        status = rlang::quo(status)
      )

      # Select Stan model
      if (spatial_model == "SAR") {
        private$stan_model <- stanmodels$nasc_SAR
      } else if (spatial_model == "SDM") {
        private$stan_model <- stanmodels$nasc_SDM
      } else {
        private$stan_model <- stanmodels$model1_gammaOmega
      }
    },

    #' @description
    #' Fit the Stan model via MCMC.
    #' @param ... Additional arguments forwarded to [rstan::sampling()].
    fit = function(...) {
      wide_df <- .makeWide(
        data      = private$data,
        id        = private$id,
        time      = private$time,
        outcome   = private$outcome,
        treatment = private$treated
      )

      pre_data  <- wide_df |> dplyr::filter(!!private$time <  private$intervention)
      post_data <- wide_df |> dplyr::filter(!!private$time >= private$intervention)

      # Base predictor matrices (outcome time series)
      X      <- pre_data  |> dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)
      X1     <- pre_data  |> dplyr::pull(!!private$outcome)
      X_pred <- post_data |> dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)

      # Append pre-treatment covariate means as extra predictor rows
      if (!is.null(private$covariates)) {
        cov_names <- setdiff(
          names(private$covariates),
          c(rlang::as_name(private$time), rlang::as_name(private$id))
        )

        # Compute per-unit pre-treatment means, then pivot to wide (rows = covariates)
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

      # Build Stan data list
      if (private$spatial_model %in% c("SAR", "SDM")) {
        donor_ids  <- colnames(X)
        treated_id <- as.character(private$treated_ids)

        # Auto-set dimnames from factor levels if not provided
        if (is.null(rownames(private$W)) || is.null(colnames(private$W))) {
          all_ids <- levels(private$data[[rlang::as_name(private$id)]])
          if (nrow(private$W) != length(all_ids)) {
            stop("W dimensions (", nrow(private$W), ") do not match the ",
                 "number of units (", length(all_ids), ").")
          }
          rownames(private$W) <- colnames(private$W) <- all_ids
        }
        W_ordered <- private$W[c(donor_ids, treated_id), c(donor_ids, treated_id)]
        J <- length(donor_ids)

        W_J  <- W_ordered[1:J,     1:J    ]
        w_J1 <- W_ordered[1:J,     J + 1  ]

        Y_panel_df <- pre_data |>
          dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
        Y_panel <- t(as.matrix(Y_panel_df))

        stan_data <- list(
          K_pred  = nrow(X),
          X1      = X1,
          J       = J,
          X0      = X,
          vs      = private$vs,
          T0      = nrow(pre_data),
          Y_panel = Y_panel,
          W       = W_ordered,
          W_J     = W_J,
          w_J1    = as.vector(w_J1),
          T_post  = nrow(post_data),
          Y0_post = as.matrix(X_pred),
          Y1_post = post_data |> dplyr::pull(!!private$outcome)
        )
      } else {
        stan_data <- list(
          K      = nrow(X),
          X1     = X1,
          J      = ncol(X),
          X0     = X,
          T_post = nrow(X_pred),
          X0_pred = X_pred,
          vs     = private$vs
        )
      }

      private$fitted <- rstan::sampling(
        private$stan_model,
        data = stan_data,
        ...
      )

      if (private$spatial_model %in% c("SAR", "SDM")) {
        # SAR/SDM models output tau_nasc directly — no synthetic draws
        private$plot_data <- .get_nasc_results(
          fit      = private$fitted,
          pre_data = pre_data,
          post_data = post_data,
          time     = private$time,
          outcome  = private$outcome,
          ci       = private$ci_width
        )
      } else {
        private$y_synth_draws <- .get_synth_draws_predictor_match(
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
    #' @param ci_width New width, between 0 and 1.
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
      private$plot_data <- .get_plot_df(
        y_synth_draws = private$y_synth_draws,
        pre_data      = pre_data,
        post_data     = post_data,
        ci            = private$ci_width,
        time          = private$time,
        outcome       = private$outcome
      )
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
    #' @param facet Logical; facet by unit (ignored for single treated unit).
    #' @param subset Optional character vector of unit IDs to highlight.
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
