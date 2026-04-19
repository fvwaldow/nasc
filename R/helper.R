# Helper Functions

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




.get_synth_draws_predictor_match <- function(fit, pre_data,
                                             post_data, time, outcome) {
  X1_sim_draws <- .get_par_long(fit = fit, par = X1_sim)
  dateXwalk <- pre_data |>
    dplyr::mutate(idx = 1:dplyr::n()) |>
    dplyr::select(idx, !!time)
  y_hat <- dplyr::inner_join(X1_sim_draws, dateXwalk, by = "idx") |>
    dplyr::rename(y_synth = X1_sim)

  X1_pred_draws <- .get_par_long(fit = fit, par = X1_pred)
  dateXwalk <- post_data |>
    dplyr::mutate(idx = 1:dplyr::n()) |>
    dplyr::select(idx, !!time)
  X1_pred_hat <- dplyr::inner_join(X1_pred_draws, dateXwalk, by = "idx") |>
    dplyr::rename(y_synth = X1_pred)

  y_synth <- dplyr::bind_rows(y_hat, X1_pred_hat) |>
    dplyr::select(-idx)

  pre_outcome <- pre_data |>
    dplyr::select(!!outcome, !!time)
  post_outcome <- post_data |>
    dplyr::select(!!outcome, !!time)

  y <- dplyr::bind_rows(pre_outcome, post_outcome)
  y_synth <- dplyr::full_join(y_synth, y, by = rlang::as_name(time))
  return(y_synth)
}



.get_synth_draws3d <- function(fit, data, id, treated_ids, time, outcome,
                               intervention) {
  y_sim_draws <-
    .get_draws3d(
      fit = fit,
      data = data,
      id = id,
      treated_ids = treated_ids,
      time = time,
      outcome = outcome,
      intervention = intervention,
      period = "pre"
    )

  y_pred_draws <-
    .get_draws3d(
      fit = fit,
      data = data,
      id = id,
      treated_ids = treated_ids,
      time = time,
      outcome = outcome,
      intervention = intervention,
      period = "post"
    )

  y_draws <- dplyr::bind_rows(y_sim_draws, y_pred_draws)

  return(y_draws)
}



.get_draws3d <- function(fit, data, id, treated_ids, time, outcome,
                         intervention, period = c("pre", "post")) {
  period <- match.arg(period)
  if (period == "pre") {
    y_sim_draws <- rstan::extract(fit, pars = "y_sim")[[1]]
    data <- data |>
      dplyr::filter(!!time < intervention)
  } else {
    y_sim_draws <- rstan::extract(fit, pars = "y_pred")[[1]]
    data <- data |>
      dplyr::filter(!!time >= intervention)
  }

  dimnames(y_sim_draws) <- list(
    "draw" = seq(1, dim(y_sim_draws)[[1]]),
    "i_idx" = seq(1, dim(y_sim_draws)[[2]]),
    "t_idx" = seq(1, dim(y_sim_draws)[[3]])
  )

  y_sim_draws <- y_sim_draws |>
    cubelyr::as.tbl_cube() |>
    tibble::as_tibble() |>
    dplyr::rename(y_hat = 4)

  wide_df_treated <- data |>
    dplyr::filter(!!id %in% treated_ids) |>
    dplyr::select(
      !!id,
      !!time,
      !!outcome
    ) |>
    tidyr::pivot_wider(names_from = !!id, values_from = !!outcome) |>
    dplyr::arrange(!!time)

  iXwalk <- wide_df_treated |>
    dplyr::select(-!!time) |>
    colnames() |>
    dplyr::tibble(id = .) |>
    dplyr::mutate(i_idx = 1:dplyr::n())

  tXwalk <- wide_df_treated |>
    dplyr::select(!!time) |>
    dplyr::mutate(t_idx = 1:dplyr::n())

  y_sim_draws <- y_sim_draws |>
    dplyr::inner_join(iXwalk, by = "i_idx") |>
    dplyr::inner_join(tXwalk, by = "t_idx") |>
    dplyr::select(-i_idx, -t_idx)

  return(y_sim_draws)
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




.get_plot_df2 <- function(y_synth_draws, data, treated_ids,
                          id, time, outcome, ci = 0.75) {
  data_treated <- data |>
    dplyr::filter(!!id %in% treated_ids) |>
    dplyr::select(!!id, !!time, !!outcome) |>
    dplyr::mutate(!!rlang::as_label(id) := as.character(!!id))

  ate <- data_treated |>
    dplyr::inner_join(y_synth_draws, by = c(
      rlang::as_name(id),
      rlang::as_name(time)
    )) |>
    dplyr::mutate(diff = !!outcome - y_hat) |>
    dplyr::group_by(!!time, draw) |>
    dplyr::summarise(diff_draw = mean(diff)) |>
    dplyr::group_by(!!time) |>
    dplyr::summarise(
      tau = mean(diff_draw),
      tau_LB = stats::quantile(diff_draw, (1 - ci) / 2),
      tau_UB = stats::quantile(diff_draw, 1 - (1 - ci) / 2)
    ) |>
    dplyr::mutate(id = "Average")

  y_synth_i <- y_synth_draws |>
    dplyr::group_by(!!time, !!id) |>
    dplyr::summarise(
      LB = stats::quantile(y_hat, (1 - ci) / 2),
      UB = stats::quantile(y_hat, 1 - (1 - ci) / 2),
      y_synth = mean(y_hat)
    )

  df_plot_i <- y_synth_i |>
    dplyr::inner_join(data_treated, by = c(
      rlang::as_name(id),
      rlang::as_name(time)
    )) |>
    dplyr::mutate(
      tau = !!outcome - y_synth,
      tau_LB = !!outcome - UB,
      tau_UB = !!outcome - LB
    )
  df_plot <- dplyr::bind_rows(df_plot_i, ate)
  return(df_plot)
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
