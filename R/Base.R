# Bayesian NASC estimator
# BASE code

#' Bayesian Network-Aware Synthetic Control
#'
#' @description
#' An R6 class implementing Bayesian Synthetic Control estimators. Supports
#' single and multiple treated units, optional covariate adjustment, and
#' predictor matching. Stan models are selected automatically based on the
#' structure of the supplied data.
#'
#' @export
bayesianSynth <- R6::R6Class(
  classname = "bayesianSynth",
  private = list(
    data = NULL,
    covariates = NULL,
    predictor_match_covariates0 = NULL,
    predictor_match_covariates1 = NULL,
    predictor_match = NULL,
    vs = NULL,
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
    mcmc_checks = NULL
  ),
  active = list(

    #' @field plotData returns tibble with the observed outcome and the
    #'    counterfactual data.
    plotData = function() {
      return(private$plot_data)
    },

    #' @field interventionTime returns intervention time period (e.g., year)
    #'    in which the treatment occurred.
    interventionTime = function() {
      return(private$intervention)
    },

    #' @field synthetic returns ggplot2 object that shows the
    #'    observed and counterfactual outcomes over time.
    synthetic = function() {
      df_plot <- private$plot_data |>
        dplyr::rename(
          Observed = !!private$outcome,
          Synthetic = y_synth
        )

      if (length(private$treated_ids) == 1) {
        df_plot <- df_plot |>
          dplyr::select(!!private$time, Observed, Synthetic, LB, UB)
      } else {
        df_plot <- df_plot |>
          dplyr::select(
            !!private$id, !!private$time, Observed,
            Synthetic, LB, UB
          ) |>
          dplyr::filter(!!private$id != "Average")
      }
      df_plot <- df_plot |>
        tidyr::pivot_longer(cols = c(Observed, Synthetic))

      synthetic_plot <- ggplot2::ggplot(
        data = df_plot,
        ggplot2::aes(x = !!private$time)
      ) +
        ggplot2::geom_line(
          ggplot2::aes(y = value, linetype = name)
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = LB, ymax = UB),
          color = "gray",
          alpha = 0.2
        ) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          legend.position = c(0.9, .1),
          legend.background = ggplot2::element_rect(
            fill = ggplot2::alpha("white", 0)
          ),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line()
        ) +
        ggplot2::geom_vline(
          xintercept = private$intervention,
          linetype = "dashed"
        )

      if (length(private$treated_ids) > 1) {
        synthetic_plot <- synthetic_plot +
          ggplot2::facet_wrap(dplyr::vars(!!private$id))
      }

      return(synthetic_plot)
    },

    #' @field checks returns MCMC checks.
    checks = function() {
      private$mcmc_checks
    }
  ),
  public = list(
    #' @description
    #' Create a new bayesianSynth object.
    #'
    #' @param data A data frame in long format.
    #' @param time Unquoted name of the time column.
    #' @param id Unquoted name of the unit identifier column.
    #' @param treated Unquoted name of the binary treatment indicator column
    #'   (1 = treated, 0 = untreated).
    #' @param outcome Unquoted name of the outcome column.
    #' @param ci_width Credible interval width, between 0 and 1. Default 0.75.
    #' @param covariates Optional data frame of covariates.
    #' @param predictor_match Logical; whether to use predictor matching.
    #'   Default FALSE.
    #' @param predictor_match_covariates0 Optional donor predictor matrix.
    #' @param predictor_match_covariates1 Optional treated unit predictor vector.
    #' @param vs Optional vector of predictor importance weights.
    initialize = function(data, time, id, treated, outcome, ci_width = 0.75,
                          covariates = NULL,
                          predictor_match = FALSE,
                          predictor_match_covariates0 = NULL,
                          predictor_match_covariates1 = NULL, vs = NULL) {

      stopifnot((ci_width > 0 & ci_width < 1))
      private$time <- rlang::enquo(time)
      private$id <- rlang::enquo(id)
      private$treated <- rlang::enquo(treated)
      private$predictor_match <- predictor_match
      private$predictor_match_covariates0 <- predictor_match_covariates0
      private$predictor_match_covariates1 <- predictor_match_covariates1
      private$vs <- vs

      message("Transforming data")

      if (!setequal(
        data |>
        dplyr::select(!!private$treated) |>
        dplyr::distinct() |>
        dplyr::pull({{ treated }}),
        c(1, 0)
      )) {
        stop("Treated identifier not binary 1 - 0.")
      }

      private$data <- data |>
        dplyr::mutate(
          status = dplyr::case_when(
            {{ treated }} == 1 ~ "Treated",
            {{ treated }} == 0 ~ "Untreated",
            is.na({{ treated }}) ~ "N/A"
          ),
          status = factor(status, levels = c(
            "Treated",
            "Untreated",
            "N/A"
          ))
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

      private$intervention <- private$data |>
        dplyr::filter(status == "Treated") |>
        dplyr::summarise(!!rlang::as_label(private$time) := min(!!private$time)) |>
        dplyr::pull(!!private$time)

      private$outcome <- rlang::enquo(outcome)
      private$ci_width <- ci_width

      if (length(private$treated_ids) == 1) {
        if (is.null(covariates)) {
          if (predictor_match) {
            private$stan_model <- stanmodels$model1_gammaOmega
          } else {
            private$stan_model <- stanmodels$model1
          }
        } else {
          private$stan_model <- stanmodels$model3
          private$covariates <- covariates |>
            dplyr::arrange(!!private$time)
        }
      } else {
        if (is.null(covariates)) {
          private$stan_model <- stanmodels$model5
        } else {
          private$stan_model <- stanmodels$model6
          private$covariates <- covariates |>
            dplyr::arrange(!!private$time)
        }
      }
    },

    #' @description
    #' Fit Stan model.
    #' @param ... Additional arguments passed to [rstan::sampling()].
    fit = function(...) {
      if (length(private$treated_ids) == 1) {
        wide_df <- .makeWide(
          data = private$data,
          id = private$id,
          time = private$time,
          outcome = private$outcome,
          treatment = private$treated
        )
        pre_data <- wide_df |>
          dplyr::filter(!!private$time < private$intervention)

        post_data <- wide_df |>
          dplyr::filter(!!private$time >= private$intervention)

        X <- pre_data |>
          dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)
        X_pred <- post_data |>
          dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)

        if (!is.null(private$covariates)) {
          M <- private$covariates |>
            dplyr::filter(!!private$time < private$intervention)

          M_pred <- private$covariates |>
            dplyr::filter(!!private$time >= private$intervention)

          stan_data <- list(
            N = nrow(pre_data),
            y = pre_data |> dplyr::pull(!!private$outcome),
            K = ncol(X),
            X = X,
            M_K = ncol(M),
            M = M,
            N_pred = nrow(X_pred),
            X_pred = X_pred,
            M_pred = M_pred
          )
        } else {
          if (private$predictor_match) {
            X1 <- pre_data |> dplyr::pull(!!private$outcome)
            if (!is.null(private$predictor_match_covariates0)) {
              X <- rbind(X, setNames(
                private$predictor_match_covariates0,
                names(X)
              ))
              X1 <- c(
                pre_data |> dplyr::pull(!!private$outcome),
                private$predictor_match_covariates1
              )
            }
            if (is.null(private$vs)) {
              private$vs <- rep(1, nrow(X))
            }

            stan_data <- list(
              K = nrow(X),
              X1 = X1,
              J = ncol(X),
              X0 = X,
              T_post = nrow(X_pred),
              X0_pred = X_pred,
              vs = private$vs
            )
          } else {
            stan_data <- list(
              N = nrow(pre_data),
              y = pre_data |> dplyr::pull(!!private$outcome),
              K = ncol(X),
              X = X,
              N_pred = nrow(X_pred),
              X_pred = X_pred
            )
          }
        }
      } else {
        Y <- private$data |>
          dplyr::filter(!!private$id %in% private$treated_ids) |>
          dplyr::select(!!private$id, !!private$time, !!private$outcome) |>
          dplyr::filter(!!private$time < private$intervention) |>
          tidyr::pivot_wider(
            names_from = !!private$id,
            values_from = !!private$outcome
          ) |>
          dplyr::arrange(!!private$time) |>
          dplyr::select(-!!private$time) |>
          as.matrix() |>
          t()

        donors <- private$data |>
          dplyr::filter(!(!!private$id %in% private$treated_ids)) |>
          dplyr::select(!!private$time, !!private$outcome, !!private$id)

        X <- donors |>
          dplyr::filter(!!private$time < private$intervention) |>
          dplyr::arrange(!!private$time) |>
          tidyr::pivot_wider(
            names_from = !!private$id,
            values_from = !!private$outcome
          ) |>
          dplyr::select(-!!private$time)

        X_pred <- donors |>
          dplyr::filter(!!private$time >= private$intervention) |>
          dplyr::arrange(!!private$time) |>
          tidyr::pivot_wider(
            names_from = !!private$id,
            values_from = !!private$outcome
          ) |>
          dplyr::select(-!!private$time)

        if (!is.null(private$covariates)) {
          M <- private$covariates |>
            dplyr::filter(
              !!private$id %in% private$treated_ids,
              !!private$time < private$intervention
            ) |>
            dplyr::arrange(!!private$id, !!private$time) |>
            dplyr::select(-!!private$time) |>
            dplyr::group_by(!!private$id) |>
            tidyr::nest() |>
            dplyr::pull(data) |>
            purrr::map(as.matrix)

          M_pred <- private$covariates |>
            dplyr::filter(
              !!private$id %in% private$treated_ids,
              !!private$time >= private$intervention
            ) |>
            dplyr::arrange(!!private$id, !!private$time) |>
            dplyr::select(-!!private$time) |>
            dplyr::group_by(!!private$id) |>
            tidyr::nest() |>
            dplyr::pull(data) |>
            purrr::map(as.matrix)

          stan_data <- list(
            N = ncol(Y),
            I = nrow(Y),
            y = Y,
            K = ncol(X),
            X = X,
            M_K = ncol(M[[1]]),
            M = M,
            N_pred = nrow(X_pred),
            X_pred = X_pred,
            M_pred = M_pred
          )
        } else {
          stan_data <- list(
            N = ncol(Y),
            I = nrow(Y),
            y = Y,
            K = ncol(X),
            X = X,
            N_pred = nrow(X_pred),
            X_pred = X_pred
          )
        }
      }

      private$fitted <-
        rstan::sampling(
          private$stan_model,
          data = stan_data,
          ...
        )

      if (length(private$treated_ids) == 1) {
        if (private$predictor_match) {
          private$y_synth_draws <- .get_synth_draws_predictor_match(
            fit = private$fitted,
            pre_data = pre_data,
            post_data = post_data,
            time = private$time,
            outcome = private$outcome
          )
        } else {
          private$y_synth_draws <- .get_synth_draws(
            fit = private$fitted,
            pre_data = pre_data,
            post_data = post_data,
            time = private$time,
            outcome = private$outcome
          )
        }

        private$plot_data <-
          .get_plot_df(
            y_synth_draws = private$y_synth_draws,
            pre_data = pre_data,
            post_data = post_data,
            ci = private$ci_width,
            time = private$time,
            outcome = private$outcome
          )
      } else {
        private$y_synth_draws <- .get_synth_draws3d(
          fit = private$fitted,
          data = private$data,
          id = private$id,
          treated_ids = private$treated_ids,
          time = private$time,
          outcome = private$outcome,
          intervention = private$intervention
        )

        private$plot_data <- .get_plot_df2(
          y_synth_draws = private$y_synth_draws,
          data = private$data,
          treated_ids = private$treated_ids,
          id = private$id,
          time = private$time,
          outcome = private$outcome,
          ci = private$ci_width
        )
      }
    },

    #' @description
    #' Update the width of the credible interval.
    #' @param ci_width New credible interval width, between 0 and 1.
    updateWidth = function(ci_width = 0.75) {
      stopifnot(exprs = {
        ci_width > 0
        ci_width < 1
      })
      private$ci_width <- ci_width

      wide_df <- .makeWide(
        data = private$data,
        id = private$id,
        time = private$time,
        outcome = private$outcome,
        treatment = private$treated
      )

      pre_data <- wide_df |>
        dplyr::filter(!!private$time < private$intervention)

      post_data <- wide_df |>
        dplyr::filter(!!private$time >= private$intervention)

      private$plot_data <-
        .get_plot_df(
          y_synth_draws = private$y_synth_draws,
          pre_data = pre_data,
          post_data = post_data,
          ci = private$ci_width,
          time = private$time,
          outcome = private$outcome
        )
    },

    #' @description
    #' Plot the estimated treatment effect over time.
    #' @param facet Logical; whether to facet by unit when multiple treated
    #'   units are present. Default TRUE.
    #' @param subset Optional character vector of unit IDs to include.
    effectPlot = function(facet = TRUE, subset = NULL) {
      if (length(private$treated_ids) == 1) {
        tau_plot <- .plot_tau(
          data = private$plot_data,
          x = private$time,
          y = tau,
          ymin = tau_LB,
          ymax = tau_UB,
          xintercept = private$intervention
        )
      } else {
        if (is.null(subset)) {
          tau_plot <- .plot_tau(
            data = private$plot_data,
            x = private$time,
            y = tau,
            ymin = tau_LB,
            ymax = tau_UB,
            xintercept = private$intervention,
            facet = private$id
          )
        } else {
          if (facet) {
            tau_plot <- .plot_tau(
              data = private$plot_data,
              x = private$time,
              y = tau,
              ymin = tau_LB,
              ymax = tau_UB,
              xintercept = private$intervention,
              facet = private$id,
              id = private$id,
              subset = subset
            )
          } else {
            tau_plot <- .plot_tau(
              data = private$plot_data,
              x = private$time,
              y = tau,
              ymin = tau_LB,
              ymax = tau_UB,
              xintercept = private$intervention,
              id = private$id,
              subset = subset
            )
          }
        }
      }
      return(tau_plot)
    },

    #' @description
    #' Run a placebo test by re-fitting the model with a fake intervention time.
    #' @param periods Number of pre-intervention periods to treat as the
    #'   placebo post-period.
    #' @param ... Additional arguments passed to [rstan::sampling()].
    placeboPlot = function(periods, ...) {
      stopifnot(periods > 0)

      keep <- private$data |>
        dplyr::filter(!!private$time < private$intervention) |>
        dplyr::select(!!private$time) |>
        dplyr::distinct() |>
        dplyr::slice_max(n = periods, order_by = !!private$time) |>
        dplyr::pull(!!private$time)

      wide_df <- .makeWide(
        data = private$data,
        id = private$id,
        time = private$time,
        outcome = private$outcome,
        treatment = private$treated
      ) |>
        dplyr::filter(!!private$time < private$intervention) |>
        dplyr::mutate(
          !!rlang::as_label(private$treated) :=
            dplyr::case_when(!!private$time %in% keep ~ 1, TRUE ~ 0)
        )

      pre_data <- wide_df |>
        dplyr::filter(!!private$treated == 0)
      post_data <- wide_df |>
        dplyr::filter(!!private$treated == 1)

      X <- pre_data |>
        dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)
      X_pred <- post_data |>
        dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)

      if (!is.null(private$covariates)) {
        max_t_pre <- pre_data |>
          dplyr::pull(!!private$time) |>
          max()
        max_t_post <- post_data |>
          dplyr::pull(!!private$time) |>
          max()

        M <- private$covariates |>
          dplyr::filter(!!private$time <= max_t_pre)
        M_pred <- private$covariates |>
          dplyr::filter(
            !!private$time > max_t_pre,
            !!private$time <= max_t_post
          )

        stan_data <- list(
          N = nrow(pre_data),
          y = pre_data |> dplyr::pull(!!private$outcome),
          K = ncol(X),
          X = X,
          M_K = ncol(M),
          M = M,
          N_pred = nrow(X_pred),
          X_pred = X_pred,
          M_pred = M_pred
        )
      } else {
        stan_data <- list(
          N = nrow(pre_data),
          y = pre_data |> dplyr::pull(!!private$outcome),
          K = ncol(X),
          X = X,
          N_pred = nrow(X_pred),
          X_pred = X_pred
        )
      }

      placebo_fitted <-
        rstan::sampling(
          private$stan_model,
          data = stan_data,
          ...
        )

      placebo_y_synth_draws <- .get_synth_draws(
        fit = placebo_fitted,
        pre_data = pre_data,
        post_data = post_data,
        time = private$time,
        outcome = private$outcome
      )

      plot_placebo_data <-
        .get_plot_df(
          y_synth_draws = placebo_y_synth_draws,
          pre_data = pre_data,
          post_data = post_data,
          ci = private$ci_width,
          time = private$time,
          outcome = private$outcome
        )

      tau_plot <- .plot_tau(
        data = plot_placebo_data,
        x = private$time,
        y = tau,
        ymin = tau_LB,
        ymax = tau_UB,
        xintercept = min(keep)
      )
      return(tau_plot)
    },

    #' @description
    #' Plot the posterior distribution of donor weights.
    weightDraws = function() {
      betas <- private$fitted |>
        as.data.frame() |>
        dplyr::select(dplyr::contains("beta"))

      beta_names <- private$data |>
        dplyr::filter(!!private$treated == 0) |>
        dplyr::select(!!private$id) |>
        unique() |>
        dplyr::pull()

      treated_name <- private$data |>
        dplyr::filter(!!private$treated == 1) |>
        dplyr::select(!!private$id) |>
        unique() |>
        dplyr::pull()

      donor_names <- beta_names[beta_names %in% treated_name == FALSE]

      names(betas) <- donor_names
      melt_betas <- tidyr::gather(betas, ID, weight)

      melt_betas |>
        ggplot2::ggplot(ggplot2::aes(x = weight, y = ID, fill = ID)) +
        ggridges::geom_density_ridges() +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
    }
  )
)
