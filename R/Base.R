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
nascSynth <- R6::R6Class(
  classname = "nascSynth",
  private = list(
    data = NULL,
    covariates0 = NULL,
    covariates1 = NULL,
    vs = NULL,
    W = NULL,
    spatial_model = NULL, # NEW: Store the chosen model type
    time = NULL,
    id = NULL,
    treated = NULL,
    outcome = NULL,
    ci_width = NULL,
    intervention = NULL,
    fitted = NULL,
    plot_data = NULL,
    time_tiles = NULL,
    stan_model = NULL,
    y_synth_draws = NULL,
    lift_draws = NULL,
    treated_ids = NULL,
    mcmc_checks = NULL
  ),
  active = list(
    timeTiles = function() { return(private$time_tiles) },
    plotData = function() { return(private$plot_data) },
    interventionTime = function() { return(private$intervention) },
    synthetic = function() {
      df_plot <- private$plot_data %>%
        dplyr::rename(Observed = !!private$outcome, Synthetic = y_synth) %>%
        dplyr::select(!!private$time, Observed, Synthetic, LB, UB) %>%
        tidyr::pivot_longer(cols = c(Observed, Synthetic))

      ggplot2::ggplot(data = df_plot, ggplot2::aes(x = !!private$time)) +
        ggplot2::geom_line(ggplot2::aes(y = value, linetype = name)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = LB, ymax = UB), color = "gray", alpha = 0.2) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          legend.position = c(0.9, .1),
          legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0)),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line()
        ) +
        ggplot2::geom_vline(xintercept = private$intervention, linetype = "dashed")
    },
    checks = function() { private$mcmc_checks },
    lift = function() { private$lift_draws }
  ),
  public = list(
    initialize = function(data, time, id, treated, outcome, ci_width = 0.75,
                          covariates0 = NULL, covariates1 = NULL, vs = NULL,
                          W = NULL, spatial_model = "none") {

      stopifnot((ci_width > 0 & ci_width < 1))

      # Validate Spatial Inputs
      if (!spatial_model %in% c("none", "SAR", "SDM")) {
        stop("spatial_model must be 'none', 'SAR', or 'SDM'.")
      }
      if (spatial_model %in% c("SAR", "SDM") && is.null(W)) {
        stop("A spatial weights matrix 'W' is required for NASC models.")
      }

      private$time <- rlang::enquo(time)
      private$id <- rlang::enquo(id)
      private$treated <- rlang::enquo(treated)
      private$outcome <- rlang::enquo(outcome)
      private$ci_width <- ci_width

      private$covariates0 <- covariates0
      private$covariates1 <- covariates1
      private$vs <- vs
      private$W <- W
      private$spatial_model <- spatial_model

      message(sprintf("Transforming data for %s Model",
                      ifelse(spatial_model == "none", "Standard SC Predictor Match", paste("NASC", spatial_model))))

      if (!setequal(data %>% dplyr::select(!!private$treated) %>% dplyr::distinct() %>% dplyr::pull({{ treated }}), c(1, 0))) {
        stop("Treated identifier not binary 1 - 0.")
      }

      private$data <- data %>%
        dplyr::mutate(
          status = dplyr::case_when(
            {{ treated }} == 1 ~ "Treated",
            {{ treated }} == 0 ~ "Untreated",
            is.na({{ treated }}) ~ "N/A"
          ),
          status = factor(status, levels = c("Treated", "Untreated", "N/A"))
        )

      if (!(data %>% dplyr::pull({{ id }}) %>% is.factor())) {
        private$data <- private$data %>% dplyr::mutate(!!rlang::quo_name(private$id) := factor({{ id }}))
      }

      private$treated_ids <- private$data %>%
        dplyr::filter(!!private$treated == 1) %>%
        dplyr::select({{ id }}) %>%
        dplyr::distinct() %>%
        dplyr::pull({{ id }})

      if (length(private$treated_ids) != 1) {
        stop("Multiple treated units detected.")
      }

      private$time_tiles <- time_tiles(data = private$data, time = !!private$time, id = !!private$id, status = status)

      private$intervention <- private$data %>%
        dplyr::filter(status == "Treated") %>%
        dplyr::summarise(!!rlang::as_label(private$time) := min(!!private$time)) %>%
        dplyr::pull(!!private$time)

      # Assign Stan Model
      if (spatial_model == "SAR") {
        private$stan_model <- stanmodels$nasc_SAR
      } else if (spatial_model == "SDM") {
        private$stan_model <- stanmodels$nasc_SDM
      } else {
        private$stan_model <- stanmodels$model1_gammaOmega
      }
    },

    fit = function(...) {
      wide_df <- .makeWide(
        data = private$data,
        id = private$id,
        time = private$time,
        outcome = private$outcome,
        treatment = private$treated
      )

      pre_data <- wide_df %>% dplyr::filter(!!private$time < private$intervention)
      post_data <- wide_df %>% dplyr::filter(!!private$time >= private$intervention)

      # Extract base matrices
      X <- pre_data %>% dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)
      X1 <- pre_data %>% dplyr::pull(!!private$outcome)
      X_pred <- post_data %>% dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)

      # Append Time-Independent Covariates
      if (!is.null(private$covariates0)) {
        X <- rbind(X, setNames(private$covariates0, names(X)))
        X1 <- c(X1, private$covariates1)
      }

      if (is.null(private$vs)) {
        private$vs <- rep(1, nrow(X))
      }

      # --- Branch Data Preparation based on Model Type ---
      if (private$spatial_model %in% c("SAR", "SDM")) {
        # 1. Identify Donors and Treated Unit
        donor_ids <- colnames(pre_data %>% dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome))
        treated_id <- as.character(private$treated_ids)

        # 2. Reorder W
        if (is.null(rownames(private$W)) || is.null(colnames(private$W))) {
          stop("The spatial weights matrix 'W' must have row and column names matching the unit IDs.")
        }
        W_ordered <- private$W[c(donor_ids, treated_id), c(donor_ids, treated_id)]
        J <- length(donor_ids)

        # 3. Extract Submatrices
        W_J <- W_ordered[1:J, 1:J]
        w_J1 <- W_ordered[1:J, J + 1]

        # 4. Prepare Y_panel
        Y_panel_df <- pre_data %>% dplyr::select(dplyr::all_of(donor_ids), !!private$outcome)
        Y_panel <- t(as.matrix(Y_panel_df))

        # Build NASC Data List
        stan_data <- list(
          K_pred = nrow(X),
          X1 = X1,
          J = J,
          X0 = X,
          vs = private$vs,
          T0 = nrow(pre_data),
          Y_panel = Y_panel,
          W = W_ordered,
          W_J = W_J,
          w_J1 = as.vector(w_J1),
          T_post = nrow(post_data),
          Y0_post = as.matrix(X_pred),
          Y1_post = post_data %>% dplyr::pull(!!private$outcome)
        )
      } else {
        # Build Standard SC Predictor Match Data List
        stan_data <- list(
          K = nrow(X),
          X1 = X1,
          J = ncol(X),
          X0 = X,
          T_post = nrow(X_pred),
          X0_pred = X_pred,
          vs = private$vs
        )
      }

      # Fit the model
      private$fitted <- rstan::sampling(private$stan_model, data = stan_data, ...)

      # Process Draws
      private$y_synth_draws <- .get_synth_draws_predictor_match(
        fit = private$fitted,
        pre_data = pre_data,
        post_data = post_data,
        time = private$time,
        outcome = private$outcome
      )

      private$plot_data <- .get_plot_df(
        y_synth_draws = private$y_synth_draws,
        pre_data = pre_data,
        post_data = post_data,
        ci = private$ci_width,
        time = private$time,
        outcome = private$outcome
      )
    },

    updateWidth = function(ci_width = 0.75) {
      stopifnot(exprs = { ci_width > 0; ci_width < 1 })
      private$ci_width <- ci_width
      wide_df <- .makeWide(data = private$data, id = private$id, time = private$time, outcome = private$outcome, treatment = private$treated)
      pre_data <- wide_df %>% dplyr::filter(!!private$time < private$intervention)
      post_data <- wide_df %>% dplyr::filter(!!private$time >= private$intervention)
      private$plot_data <- .get_plot_df(y_synth_draws = private$y_synth_draws, pre_data = pre_data, post_data = post_data, ci = private$ci_width, time = private$time, outcome = private$outcome)
    },

    summarizeLift = function() {
      if (is.null(private$lift_draws)) stop("You first need to run the `liftDraws()` method.")
      return(c(
        point = mean(private$lift_draws$lift),
        lower_bound = stats::quantile(private$lift_draws$lift, (1 - private$ci_width) / 2),
        upper_bound = stats::quantile(private$lift_draws$lift, 1 - (1 - private$ci_width) / 2)
      ))
    },

    effectPlot = function(facet = TRUE, subset = NULL) {
      return(.plot_tau(
        data = private$plot_data,
        x = private$time, y = tau, ymin = tau_LB, ymax = tau_UB,
        xintercept = private$intervention
      ))
    }
  )
)

