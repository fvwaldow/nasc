# Helper Functions

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



# Post-processing for SAR/SDM models: extracts y_counterfactual draws,
# calculates tau dynamically, and builds plot_data.
.get_nasc_results <- function(y_counterfactual_draws, bias_correction_draws, y_sim_pre_draws, pre_data, post_data, time, outcome, ci = 0.75) {
  # --- 1. POST-TREATMENT PERIOD ---
  post_times <- post_data |> dplyr::pull(!!time)
  Y1_post    <- post_data |> dplyr::pull(!!outcome)

  # Replicate the post-treatment outcome to match the dimensions of the draws matrix
  # (Rows = draws, Columns = time periods)
  Y1_mat <- matrix(Y1_post, nrow = nrow(y_counterfactual_draws), ncol = length(Y1_post), byrow = TRUE)

  # Calculate tau draws dynamically: tau_nasc = (Y1_post - y_counterfactual) * bias_correction
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




.get_synth_draws <- function(fit, pre_data, post_data, time, outcome) {
  y_sim_draws <- .get_par_long(fit = fit, par = y_sim)
  dateXwalk <- pre_data |>
    dplyr::mutate(idx = 1:dplyr::n()) |>
    dplyr::select(idx, !!time)
  y_hat <- dplyr::inner_join(y_sim_draws, dateXwalk, by = "idx") |>
    dplyr::rename(y_synth = y_sim)

  y_pred_draws <- .get_par_long(fit = fit, par = y_pred)

  dateXwalk <- post_data |>
    dplyr::mutate(idx = 1:dplyr::n()) |>
    dplyr::select(idx, !!time)
  y_pred_hat <- dplyr::inner_join(y_pred_draws, dateXwalk, by = "idx") |>
    dplyr::rename(y_synth = y_pred)
  y_synth <- dplyr::bind_rows(y_hat, y_pred_hat) |>
    dplyr::select(-idx)

  pre_outcome <- pre_data |>
    dplyr::select(!!outcome, !!time)
  post_outcome <- post_data |>
    dplyr::select(!!outcome, !!time)

  y <- dplyr::bind_rows(pre_outcome, post_outcome)
  y_synth <- dplyr::full_join(y_synth, y, by = rlang::as_name(time))
  return(y_synth)
}




.get_plot_df <- function(y_synth_draws, pre_data,
                         post_data, time, outcome, ci = 0.75) {
  y_synth <- y_synth_draws |>
    dplyr::group_by(!!time) |>
    dplyr::summarise(
      LB = stats::quantile(y_synth, (1 - ci) / 2),
      UB = stats::quantile(y_synth, 1 - (1 - ci) / 2),
      y_synth = mean(y_synth)
    )

  all_data <- dplyr::bind_rows(pre_data, post_data)

  df_plot_all <- dplyr::inner_join(y_synth, all_data,
    by = rlang::as_name(time)
  ) |>
    dplyr::mutate(
      tau = !!outcome - y_synth,
      tau_LB = !!outcome - UB,
      tau_UB = !!outcome - LB
    ) |>
    dplyr::select(!!time, !!outcome, y_synth, LB, UB, tau, tau_LB, tau_UB)
  return(df_plot_all)
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
