# Helper Functions

# ----------------------------------------------------------------------------
# Internal helper: run Step 2 across one or many rho draws.
#
# Used by both engines (stan_2_NASC sets rho_field = "rho", model1 sets
# rho_field = "rho_bc"). When length(rhos) == 1, runs a single Stan fit with
# default chains/iter. When > 1, runs the parallel furrr loop with reduced
# chains (1) per worker.
# ----------------------------------------------------------------------------
.run_step2_loop <- function(rhos, base_data, step2_mod, rho_field, cores,
                            extra_args, extract_pars) {

  if (length(rhos) == 1L) {

    worker_data <- base_data
    worker_data[[rho_field]] <- rhos

    fit_args <- c(
      list(object = step2_mod, data = worker_data, cores = cores),
      extra_args
    )
    fit <- do.call(rstan::sampling, fit_args)

    draws <- rstan::extract(fit, pars = extract_pars)

    # Replicate rho across draws for bookkeeping
    n_draws <- length(draws[[extract_pars[length(extract_pars)]]])
    if (is.matrix(draws[[1]])) n_draws <- nrow(draws[[1]])

    out <- draws
    out$rhos_used <- rep(rhos, n_draws)
    out$last_fit  <- fit
    return(out)

  } else {

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)

    run_worker <- function(single_rho, base_data, step2_mod, rho_field,
                           extract_pars) {
      suppressPackageStartupMessages(require(rstan, quietly = TRUE))
      worker_data <- base_data
      worker_data[[rho_field]] <- single_rho

      fit_step2 <- rstan::sampling(
        object        = step2_mod,
        data          = worker_data,
        chains        = 1,
        iter          = 1000,
        warmup        = 500,
        refresh       = 0,
        show_messages = FALSE
      )

      draws <- rstan::extract(fit_step2, pars = extract_pars)
      draws$rho_used <- single_rho
      draws
    }

    progressr::handlers(global = TRUE)
    results_list <- progressr::with_progress({
      p <- progressr::progressor(steps = length(rhos))
      furrr::future_map(
        rhos,
        function(rho) {
          res <- run_worker(
            single_rho   = rho,
            base_data    = base_data,
            step2_mod    = step2_mod,
            rho_field    = rho_field,
            extract_pars = extract_pars
          )
          p()
          res
        },
        .options = furrr::furrr_options(seed = TRUE, packages = "rstan")
      )
    })

    # Aggregate. Matrix-valued params (multi-dim draws) are rbound; scalar
    # vector params are concatenated.
    out <- list()
    for (par in extract_pars) {
      first <- results_list[[1]][[par]]
      if (is.null(dim(first)) || length(dim(first)) == 1L) {
        out[[par]] <- do.call(c, lapply(results_list, \(x) x[[par]]))
      } else {
        out[[par]] <- do.call(rbind, lapply(results_list, \(x) x[[par]]))
      }
    }
    out$rhos_used <- vapply(results_list, \(x) x$rho_used, numeric(1))
    out$last_fit  <- NULL  # workers run in subprocesses; no last_fit retained
    return(out)
  }
}


# Creates a tile chart of treatment status over time and unit.
# `time`, `id`, and `status` must be quosures.
.time_tiles <- function(data, time, id, status) {
  ggplot2::ggplot(data, ggplot2::aes(x = !!time, y = !!id, fill = !!status)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(
      values = c("Treated" = "#E74C3C", "Untreated" = "#BDC3C7", "N/A" = "gray90"),
      drop   = FALSE
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::labs(
      x = rlang::as_name(time),
      y = rlang::as_name(id)
    )
}

.get_par_long <- function(fit, par) {
  par <- rlang::enquo(par)
  long_tlb <- fit |>
    as.data.frame(pars = rlang::quo_text(par)) |>
    dplyr::mutate(draw = 1:dplyr::n()) |>
    tidyr::pivot_longer(
      names_to = "idx",
      values_to = rlang::quo_text(par), -draw
    ) |>
    dplyr::mutate(idx = as.numeric(gsub(
      glue::glue("{rlang::quo_text(par)}\\[(\\d+)\\]"),
      "\\1",
      idx
    ))) |>
    dplyr::as_tibble()
  return(long_tlb)
}



# Unified post-processing: extracts y_counterfactual and y_sim_pre draws,
# multiplies tau by bias_correction (which is 1.0 when bias correction is off),
# and builds plot_data. Used for ALL model paths after the harmonization of
# generated-quantity names across model1.stan and stan_2_NASC.stan.
.get_nasc_results <- function(y_counterfactual_draws, bias_correction_draws, y_sim_pre_draws, pre_data, post_data, time, outcome, ci = 0.75) {
  # --- 1. POST-TREATMENT PERIOD ---
  post_times <- post_data |> dplyr::pull(!!time)
  Y1_post    <- post_data |> dplyr::pull(!!outcome)

  # Replicate the post-treatment outcome to match the dimensions of the draws matrix
  # (Rows = draws, Columns = time periods)
  Y1_mat <- matrix(Y1_post, nrow = nrow(y_counterfactual_draws), ncol = length(Y1_post), byrow = TRUE)

  # Calculate tau draws dynamically: tau = (Y1_post - y_counterfactual) * bias_correction
  # bias_correction is 1.0 when use_bias_correction = 0, so this is the naive
  # difference in that case.
  # Note: Due to R's column-major matrix memory layout, multiplying an N x T matrix
  # by an N-length vector correctly multiplies each column element-wise by the vector.
  tau_draws <- (Y1_mat - y_counterfactual_draws) * bias_correction_draws

  tau_summary <- tibble::tibble(
    !!rlang::as_name(time) := post_times,
    tau    = apply(tau_draws, 2, mean),
    tau_LB = apply(tau_draws, 2, \(x) stats::quantile(x, (1 - ci) / 2)),
    tau_UB = apply(tau_draws, 2, \(x) stats::quantile(x, 1 - (1 - ci) / 2))
  )

  # Compute y_synth = observed - tau (invert for ribbon: LB/UB swap)
  post_outcome <- post_data |> dplyr::select(!!time, !!outcome)
  post_plot <- dplyr::inner_join(tau_summary, post_outcome,
                                 by = rlang::as_name(time)
  ) |>
    dplyr::mutate(
      y_synth = !!outcome - tau,
      LB      = !!outcome - tau_UB,
      UB      = !!outcome - tau_LB
    ) |>
    dplyr::select(!!time, !!outcome, y_synth, LB, UB, tau, tau_LB, tau_UB)

  # --- 2. PRE-TREATMENT PERIOD ---
  pre_times <- pre_data |> dplyr::pull(!!time)

  pre_summary <- tibble::tibble(
    !!rlang::as_name(time) := pre_times,
    y_synth = apply(y_sim_pre_draws, 2, mean),
    LB      = apply(y_sim_pre_draws, 2, \(x) stats::quantile(x, (1 - ci) / 2)),
    UB      = apply(y_sim_pre_draws, 2, \(x) stats::quantile(x, 1 - (1 - ci) / 2))
  )

  pre_outcome <- pre_data |> dplyr::select(!!time, !!outcome)
  pre_plot <- dplyr::inner_join(pre_summary, pre_outcome,
                                by = rlang::as_name(time)) |>
    dplyr::mutate(
      tau    = !!outcome - y_synth,
      tau_LB = !!outcome - UB,
      tau_UB = !!outcome - LB
    ) |>
    dplyr::select(!!time, !!outcome, y_synth, LB, UB, tau, tau_LB, tau_UB)

  # Combine pre and post data
  dplyr::bind_rows(pre_plot, post_plot)
}




.makeWide <- function(data, id, time, outcome, treatment) {
  data <- data |>
    dplyr::mutate(.tmp_id = as.integer(as.factor(!!id)))

  treated <- data |>
    dplyr::filter(status == "Treated") |>
    dplyr::select(.tmp_id) |>
    dplyr::distinct() |>
    dplyr::pull(.tmp_id)

  wide_df_treated <- data |>
    dplyr::filter(.tmp_id %in% treated) |>
    dplyr::select(!!time, !!treatment, !!outcome)

  wide_df_untreated <- data |>
    dplyr::filter(!(.tmp_id %in% treated)) |>
    dplyr::select(!!time, !!outcome, !!id) |>
    tidyr::pivot_wider(
      names_from = !!id,
      values_from = !!outcome
    )

  dplyr::inner_join(wide_df_treated, wide_df_untreated,
    by = rlang::as_name(time)
  )
}




.plot_tau <- function(data, x, y, ymin, ymax, xintercept,
                      facet, id, subset = NULL) {
  if (!is.null(subset)) {
    data <- data |>
      dplyr::filter(!!id %in% subset)
  }
  tau_plot <- ggplot2::ggplot(data = data, ggplot2::aes(x = !!x)) +
    ggplot2::geom_line(ggplot2::aes(y = {{ y }})) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = {{ ymin }}, ymax = {{ ymax }}),
      color = "gray",
      alpha = 0.2
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line()
    ) +
    ggplot2::geom_vline(xintercept = xintercept, linetype = "dashed")

  if (!missing(facet)) {
    tau_plot <- tau_plot + ggplot2::facet_grid(cols = dplyr::vars(!!facet))
  }
  return(tau_plot)
}
