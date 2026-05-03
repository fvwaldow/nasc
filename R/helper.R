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
                            extra_args, extract_pars,
                            worker_iter = 2000L, worker_warmup = 1000L) {

  stopifnot(worker_iter > worker_warmup, worker_warmup > 0)
  worker_iter   <- as.integer(worker_iter)
  worker_warmup <- as.integer(worker_warmup)

  if (length(rhos) == 1L) {

    worker_data <- base_data
    worker_data[[rho_field]] <- rhos

    fit_args <- c(
      list(object = step2_mod, data = worker_data, cores = cores),
      extra_args
    )
    fit <- do.call(rstan::sampling, fit_args)

    draws <- rstan::extract(fit, pars = extract_pars)

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
                           extract_pars, worker_iter, worker_warmup) {
      suppressPackageStartupMessages(require(rstan, quietly = TRUE))
      worker_data <- base_data
      worker_data[[rho_field]] <- single_rho

      fit_step2 <- rstan::sampling(
        object        = step2_mod,
        data          = worker_data,
        chains        = 1,
        iter          = worker_iter,
        warmup        = worker_warmup,
        refresh       = 0,
        show_messages = FALSE
      )

      draws <- rstan::extract(fit_step2, pars = extract_pars)
      draws$rho_used <- single_rho

      diag_one <- tryCatch({
        s <- rstan::summary(fit_step2)$summary
        s <- s[!rownames(s) %in% "lp__", , drop = FALSE]
        sp <- rstan::get_sampler_params(fit_step2, inc_warmup = FALSE)
        sp1 <- if (length(sp)) sp[[1]] else NULL
        max_td <- 10L
        if (length(fit_step2@stan_args)) {
          ctrl <- fit_step2@stan_args[[1]]$control
          if (!is.null(ctrl$max_treedepth)) max_td <- ctrl$max_treedepth
        }
        list(
          max_rhat    = if (nrow(s)) suppressWarnings(max(s[, "Rhat"],  na.rm = TRUE)) else NA_real_,
          min_n_eff   = if (nrow(s)) suppressWarnings(min(s[, "n_eff"], na.rm = TRUE)) else NA_real_,
          n_divergent = if (!is.null(sp1)) as.integer(sum(sp1[, "divergent__"])) else NA_integer_,
          n_max_td    = if (!is.null(sp1) && "treedepth__" %in% colnames(sp1)) {
            as.integer(sum(sp1[, "treedepth__"] >= max_td))
          } else NA_integer_,
          rho_used    = single_rho
        )
      }, error = function(e) {
        list(max_rhat = NA_real_, min_n_eff = NA_real_,
             n_divergent = NA_integer_, n_max_td  = NA_integer_,
             rho_used = single_rho)
      })

      draws$diagnostics_one <- diag_one
      draws
    }

    progressr::handlers(global = TRUE)
    results_list <- progressr::with_progress({
      p <- progressr::progressor(steps = length(rhos))
      furrr::future_map(
        rhos,
        function(rho) {
          res <- run_worker(
            single_rho    = rho,
            base_data     = base_data,
            step2_mod     = step2_mod,
            rho_field     = rho_field,
            extract_pars  = extract_pars,
            worker_iter   = worker_iter,
            worker_warmup = worker_warmup
          )
          p()
          res
        },
        .options = furrr::furrr_options(seed = TRUE, packages = "rstan")
      )
    })

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

    diags <- lapply(results_list, \(x) x$diagnostics_one)
    rhat_vec <- vapply(diags, \(d) d$max_rhat, numeric(1))
    rhat_threshold <- 1.05
    n_finite_rhat <- sum(is.finite(rhat_vec))
    n_converged   <- sum(is.finite(rhat_vec) & rhat_vec <= rhat_threshold)
    out$worker_diagnostics <- list(
      converged_share   = if (n_finite_rhat > 0) n_converged / n_finite_rhat else NA_real_,
      n_converged       = n_converged,
      n_finite_rhat     = n_finite_rhat,
      rhat_threshold    = rhat_threshold,
      max_rhat_observed = suppressWarnings(max(rhat_vec, na.rm = TRUE)),
      min_n_eff         = suppressWarnings(min(vapply(diags, \(d) d$min_n_eff,   numeric(1)),  na.rm = TRUE)),
      total_divergent   = sum(vapply(diags, \(d) d$n_divergent,                  integer(1)),  na.rm = TRUE),
      total_max_td      = sum(vapply(diags, \(d) d$n_max_td,                     integer(1)),  na.rm = TRUE),
      n_workers         = length(results_list),
      iter_per_worker   = worker_iter,
      warmup_per_worker = worker_warmup
    )

    out$last_fit  <- NULL
    return(out)
  }
}


.get_nasc_results <- function(y_counterfactual_draws, bias_correction_draws, y_sim_pre_draws, pre_data, post_data, time, outcome, ci = 0.75) {
  # Defensive: when every Step-2 worker fails to create the sampler, the
  # parallel loop returns empty arrays that rbind down to NULL or a 0-row
  # matrix. Without this check the caller eventually hits an opaque
  # "non-numeric matrix extent" error inside matrix() far below; with it,
  # the user gets an actionable message.
  if (is.null(y_counterfactual_draws) ||
      !is.matrix(y_counterfactual_draws) ||
      nrow(y_counterfactual_draws) == 0L) {
    stop("Step 2 produced no posterior draws")
  }

  post_times <- post_data |> dplyr::pull(!!time)
  Y1_post    <- post_data |> dplyr::pull(!!outcome)

  Y1_mat <- matrix(Y1_post, nrow = nrow(y_counterfactual_draws), ncol = length(Y1_post), byrow = TRUE)

  bc_vec <- as.numeric(bias_correction_draws)
  bc_mat <- matrix(bc_vec, nrow = nrow(y_counterfactual_draws),
                   ncol = ncol(y_counterfactual_draws), byrow = FALSE)
  tau_draws <- (Y1_mat - y_counterfactual_draws) * bc_mat

  tau_summary <- tibble::tibble(
    !!rlang::as_name(time) := post_times,
    tau    = apply(tau_draws, 2, mean),
    tau_LB = apply(tau_draws, 2, \(x) stats::quantile(x, (1 - ci) / 2)),
    tau_UB = apply(tau_draws, 2, \(x) stats::quantile(x, 1 - (1 - ci) / 2))
  )

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

  pre_times <- pre_data |> dplyr::pull(!!time)

  if (ncol(y_sim_pre_draws) != length(pre_times)) {
    stop(sprintf(
      "y_sim_pre_draws has %d columns but pre_data has %d rows. ",
      ncol(y_sim_pre_draws), length(pre_times)
    ),
    "")
  }

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


.plot_tau <- function(data, x, y, ymin, ymax, xintercept) {
  data <- as.data.frame(data)

  # Resolve column names whether passed as symbol/quosure/string
  resolve <- function(val, quo) {
    if (rlang::is_quosure(val))            return(rlang::as_name(val))
    if (is.character(val) && length(val) == 1L) return(val)
    if (rlang::quo_is_symbol(quo))         return(rlang::as_name(quo))
    rlang::as_name(quo)
  }
  xn    <- resolve(x,    rlang::enquo(x))
  yn    <- resolve(y,    rlang::enquo(y))
  yminn <- resolve(ymin, rlang::enquo(ymin))
  ymaxn <- resolve(ymax, rlang::enquo(ymax))

  xv  <- data[[xn]]
  yv  <- data[[yn]]
  lbv <- data[[yminn]]
  ubv <- data[[ymaxn]]

  ord <- order(xv)
  xv <- xv[ord]; yv <- yv[ord]; lbv <- lbv[ord]; ubv <- ubv[ord]

  yrng <- range(c(yv, lbv, ubv), na.rm = TRUE)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(bty = "l")

  plot(xv, yv, type = "n", ylim = yrng, xlab = xn, ylab = yn)
  graphics::grid(lty = "dotted", col = "gray80")
  graphics::polygon(c(xv, rev(xv)), c(lbv, rev(ubv)),
                    col = grDevices::adjustcolor("gray", alpha.f = 0.2),
                    border = NA)
  graphics::lines(xv, yv, lwd = 2)
  graphics::abline(v = xintercept, lty = 2)
  graphics::abline(h = 0, lty = 1, col = "black")
  invisible(NULL)
}



`%||%` <- function(a, b) if (is.null(a)) b else a

.ci_probs <- function(ci) c((1 - ci) / 2, 1 - (1 - ci) / 2)

.posterior_summary <- function(x, ci) {
  q <- stats::quantile(x, .ci_probs(ci), names = FALSE, na.rm = TRUE)
  c(
    mean   = mean(x, na.rm = TRUE),
    sd     = stats::sd(x, na.rm = TRUE),
    lower  = q[1],
    upper  = q[2],
    p_pos  = mean(x > 0, na.rm = TRUE)
  )
}

.mcmc_diagnostics <- function(fit) {
  if (is.null(fit)) return(NULL)
  diag <- tryCatch({
    s <- rstan::summary(fit)$summary
    keep <- !rownames(s) %in% c("lp__")
    s <- s[keep, , drop = FALSE]

    sp <- tryCatch(
      rstan::get_sampler_params(fit, inc_warmup = FALSE),
      error = function(e) NULL
    )
    n_div <- if (!is.null(sp)) {
      sum(vapply(sp, \(x) sum(x[, "divergent__"]), numeric(1)))
    } else NA_integer_

    list(
      source     = "single_fit",
      max_rhat   = max(s[, "Rhat"],  na.rm = TRUE),
      min_n_eff  = min(s[, "n_eff"], na.rm = TRUE),
      n_divergent = n_div,
      n_max_td   = NA_integer_,
      n_chains   = length(fit@stan_args),
      n_iter     = if (length(fit@stan_args)) fit@stan_args[[1]]$iter else NA_integer_,
      n_warmup   = if (length(fit@stan_args)) fit@stan_args[[1]]$warmup else NA_integer_,
      n_workers  = NA_integer_
    )
  }, error = function(e) NULL)
  diag
}

.worker_to_diagnostics <- function(wd) {
  if (is.null(wd)) return(NULL)
  list(
    source            = "worker_loop",
    max_rhat          = wd$max_rhat_observed,
    converged_share   = wd$converged_share,
    n_converged       = wd$n_converged,
    n_finite_rhat     = wd$n_finite_rhat,
    rhat_threshold    = wd$rhat_threshold,
    min_n_eff         = wd$min_n_eff,
    n_divergent       = wd$total_divergent,
    n_max_td          = wd$total_max_td,
    n_chains          = wd$n_workers,
    n_iter            = wd$iter_per_worker,
    n_warmup          = wd$warmup_per_worker,
    n_workers         = wd$n_workers
  )
}

.nasc_summary_stats <- function(parts) {

  ci         <- parts$ci_width
  draws      <- parts$y_synth_draws
  plot_data  <- parts$plot_data
  intervention <- parts$intervention
  outcome_nm <- rlang::as_name(parts$outcome)
  time_nm    <- rlang::as_name(parts$time)

  if (is.null(draws) || is.null(plot_data)) {
    stop("Run $fit() before calling summary().")
  }

  post_data <- plot_data |>
    dplyr::filter(!!parts$time >= intervention)
  Y1_post   <- post_data[[outcome_nm]]
  ycf       <- draws$y_counterfactual
  bc_vec    <- as.numeric(draws$bias_correction)
  Y1_mat    <- matrix(Y1_post, nrow = nrow(ycf), ncol = length(Y1_post),
                      byrow = TRUE)
  bc_mat    <- matrix(bc_vec,  nrow = nrow(ycf), ncol = ncol(ycf),
                      byrow = FALSE)
  tau_draws <- (Y1_mat - ycf) * bc_mat

  att_draws <- rowMeans(tau_draws)
  att <- .posterior_summary(att_draws, ci)

  per_period <- tibble::tibble(
    !!time_nm := post_data[[time_nm]],
    mean   = apply(tau_draws, 2, mean),
    sd     = apply(tau_draws, 2, stats::sd),
    lower  = apply(tau_draws, 2, \(x) stats::quantile(x, .ci_probs(ci)[1], names = FALSE)),
    upper  = apply(tau_draws, 2, \(x) stats::quantile(x, .ci_probs(ci)[2], names = FALSE)),
    p_pos  = apply(tau_draws, 2, \(x) mean(x > 0))
  )

  pre_data <- plot_data |>
    dplyr::filter(!!parts$time < intervention)
  if (nrow(pre_data) > 0L) {
    pre_resid <- pre_data[[outcome_nm]] - pre_data$y_synth
    pre_rmse  <- sqrt(mean(pre_resid^2))
    ss_res <- sum(pre_resid^2)
    ss_tot <- sum((pre_data[[outcome_nm]] - mean(pre_data[[outcome_nm]]))^2)
    pre_r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  } else {
    pre_rmse <- NA_real_
    pre_r2   <- NA_real_
  }

  post_resid <- post_data[[outcome_nm]] - post_data$y_synth
  post_rmse  <- if (length(post_resid)) sqrt(mean(post_resid^2)) else NA_real_
  rmspe_ratio <- if (!is.na(pre_rmse) && !is.na(post_rmse) && pre_rmse > 0) {
    post_rmse / pre_rmse
  } else NA_real_

  w_mat <- draws$w
  treated_id  <- as.character(parts$treated_ids)
  # Canonical donor ordering = colnames(X_pred) at fit time. Falling back
  # to setdiff(levels(id), treated_id) silently scrambles labels whenever
  # the long data isn't sorted by id.
  donor_ids <- parts$donor_ids
  if (is.null(donor_ids)) {
    warning("parts$donor_ids missing")
    all_ids   <- levels(parts$id_levels)
    donor_ids <- setdiff(all_ids, treated_id)
  }
  donor_names <- donor_ids
  if (!is.null(w_mat) && ncol(w_mat) == length(donor_names)) {
    colnames(w_mat) <- donor_names
  }
  weight_tbl <- tibble::tibble(
    donor = colnames(w_mat) %||% paste0("donor_", seq_len(ncol(w_mat))),
    mean  = apply(w_mat, 2, mean),
    sd    = apply(w_mat, 2, stats::sd),
    lower = apply(w_mat, 2, \(x) stats::quantile(x, .ci_probs(ci)[1], names = FALSE)),
    upper = apply(w_mat, 2, \(x) stats::quantile(x, .ci_probs(ci)[2], names = FALSE))
  ) |>
    dplyr::arrange(dplyr::desc(mean))

  w_mean <- weight_tbl$mean
  w_mean <- w_mean / sum(w_mean)
  entropy <- -sum(ifelse(w_mean > 0, w_mean * log(w_mean), 0))
  eff_donors <- exp(entropy)

  param_rows <- list()
  sigma_draws <- draws$sigma %||% draws$sigma_sc
  if (!is.null(sigma_draws)) {
    param_rows[["sigma"]] <- .posterior_summary(as.numeric(sigma_draws), ci)
  }
  if (!is.null(draws$lambda)) {
    param_rows[["lambda"]] <- .posterior_summary(as.numeric(draws$lambda), ci)
  }
  if (parts$bias_correction && !is.null(draws$bias_correction)) {
    param_rows[["bias_correction"]] <- .posterior_summary(bc_vec, ci)
  }
  rhos <- draws$rhos_used
  if (parts$uses_rho && !all(is.na(rhos))) {
    param_rows[["rho"]] <- .posterior_summary(as.numeric(rhos), ci)
  }
  param_tbl <- if (length(param_rows)) {
    tibble::as_tibble(do.call(rbind, param_rows), rownames = "parameter")
  } else NULL

  diag <- if (!is.null(draws$worker_diagnostics)) {
    .worker_to_diagnostics(draws$worker_diagnostics)
  } else {
    .mcmc_diagnostics(parts$fitted)
  }

  # ----------------------------------------------------------------
  # Indirect (spillover) effects -- Proposition 6.2.
  #
  # Computed only when (i) a fitted model object is available on parts
  # (so we can rebuild s from rho and W) AND (ii) the configuration uses
  # rho. Without rho, spillover is identically zero by construction and
  # we omit the indirect blocks rather than print zeros.
  # ----------------------------------------------------------------
  indirect_per_period_total <- NULL  # period-total spillover (sum over donors)
  indirect_per_donor        <- NULL  # average across periods, per donor
  indirect_avg_total        <- NULL  # scalar: average total spillover (analog ATT)
  if (isTRUE(parts$uses_rho) && !is.null(parts$model)) {
    ind <- tryCatch(
      .nasc_indirect_draws(parts$model),
      error = function(e) NULL
    )
    if (!is.null(ind)) {
      delta_total_draws <- ind$delta_total                  # [n_draws x T_post]
      avg_per_donor_dr  <- ind$avg_per_donor                # [n_draws x J]
      avg_total_draws   <- ind$avg_total                    # [n_draws]

      indirect_per_period_total <- tibble::tibble(
        !!time_nm := ind$time_post,
        mean   = apply(delta_total_draws, 2, mean),
        sd     = apply(delta_total_draws, 2, stats::sd),
        lower  = apply(delta_total_draws, 2, \(x) stats::quantile(x, .ci_probs(ci)[1], names = FALSE)),
        upper  = apply(delta_total_draws, 2, \(x) stats::quantile(x, .ci_probs(ci)[2], names = FALSE)),
        p_pos  = apply(delta_total_draws, 2, \(x) mean(x > 0))
      )

      indirect_per_donor <- tibble::tibble(
        donor = ind$donor_names,
        mean  = apply(avg_per_donor_dr, 2, mean),
        sd    = apply(avg_per_donor_dr, 2, stats::sd),
        lower = apply(avg_per_donor_dr, 2, \(x) stats::quantile(x, .ci_probs(ci)[1], names = FALSE)),
        upper = apply(avg_per_donor_dr, 2, \(x) stats::quantile(x, .ci_probs(ci)[2], names = FALSE)),
        p_pos = apply(avg_per_donor_dr, 2, \(x) mean(x > 0))
      ) |>
        dplyr::arrange(dplyr::desc(abs(mean)))

      indirect_avg_total <- .posterior_summary(avg_total_draws, ci)
    }
  }

  out <- list(
    header = list(
      spatial_model   = parts$spatial_model,
      bias_correction = parts$bias_correction,
      nasc_penalty    = parts$nasc_penalty,
      rho_source      = parts$rho_source,
      treated_unit    = treated_id,
      intervention    = intervention,
      n_pre           = nrow(pre_data),
      n_post          = nrow(post_data),
      n_donors        = length(donor_names),
      n_draws         = nrow(ycf),
      ci_width        = ci,
      outcome         = outcome_nm,
      time            = time_nm
    ),
    att         = att,
    per_period  = per_period,
    # Indirect effect blocks: NULL when no rho is in use.
    indirect_per_period = indirect_per_period_total,
    indirect_per_donor  = indirect_per_donor,
    indirect_avg        = indirect_avg_total,
    pre_rmse    = pre_rmse,
    pre_r2      = pre_r2,
    post_rmse   = post_rmse,
    rmspe_ratio = rmspe_ratio,
    weights     = weight_tbl,
    eff_donors  = eff_donors,
    parameters  = param_tbl,
    mcmc        = diag
  )
  class(out) <- c("summary.nascSynth", "list")
  out
}

#' Print method for nascSynth summary objects
#'
#' @param x A \code{summary.nascSynth} object.
#' @param digits Number of significant digits to display. Default \code{3}.
#' @param max_donors Maximum number of donor weights and per-donor indirect
#'   effects to display. Defaults to \code{Inf}, which shows ALL donors. Pass
#'   a finite integer (e.g. \code{10}) to restrict the output to the top
#'   donors by posterior mean (for the weights table) or by absolute
#'   posterior mean (for the per-donor indirect effects table). Also accepts
#'   \code{NULL} as an alias for "show all".
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.summary.nascSynth <- function(x, digits = 3, max_donors = Inf, ...) {
  # Allow NULL as a synonym for "show all"; coerce to Inf so the same
  # min(max_donors, nrow(.)) logic works for both bounded and unbounded
  # requests without needing branching everywhere below.
  if (is.null(max_donors)) max_donors <- Inf

  h <- x$header
  cat("Network-Aware Synthetic Control\n")
  cat("\n")
  cat(sprintf("  Outcome         : %s\n", h$outcome))
  cat(sprintf("  Treated unit    : %s\n", h$treated_unit))
  cat(sprintf("  J               : %d\n", h$n_donors))
  cat(sprintf("  T_0             : %d\n", h$n_pre))
  cat(sprintf("  T               : %d\n", h$n_post))
  cat(sprintf("  Spatial model   : %s\n", h$spatial_model))
  cat(sprintf("  Bias correction : %s\n", h$bias_correction))
  cat(sprintf("  NASC penalty    : %s\n", h$nasc_penalty))
  cat(sprintf("  Rho source      : %s\n", h$rho_source))
  cat(sprintf("  Posterior draws : %d (post-warmup, all chains pooled)\n", h$n_draws))
  cat(sprintf("  Credible level  : %.0f%%\n", 100 * h$ci_width))
  cat("\n")

  ci_pct <- 100 * h$ci_width
  ci_lo <- sprintf("l-%g%% CrI", ci_pct)
  ci_hi <- sprintf("u-%g%% CrI", ci_pct)

  # Per-period TE (direct effect)
  has_indirect <- !is.null(x$indirect_per_period)
  cat(if (has_indirect) "Per-period Direct Treatment Effect\n" else "Per-period Treatment Effect\n")
  pp <- x$per_period
  pp_print <- data.frame(
    period = format(pp[[h$time]]),
    mean   = formatC(pp$mean,  digits = digits, format = "f"),
    sd     = formatC(pp$sd,    digits = digits, format = "f"),
    lower  = formatC(pp$lower, digits = digits, format = "f"),
    upper  = formatC(pp$upper, digits = digits, format = "f"),
    `Pr>0` = formatC(pp$p_pos, digits = 3, format = "f"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(pp_print)[1] <- "Period"
  names(pp_print)[2] <- "Estimate"
  names(pp_print)[3] <- "Est.Error"
  names(pp_print)[4] <- ci_lo
  names(pp_print)[5] <- ci_hi
  print(pp_print, row.names = FALSE, right = TRUE)
  cat("\n")

  # ATT (direct)
  att <- x$att
  att_p_dir <- if (att["mean"] >= 0) att["p_pos"] else 1 - att["p_pos"]
  cat(if (has_indirect) "Average direct Treatment Effect \n" else "Average Treatment Effect\n")
  att_print <- data.frame(
    blank   = "",
    mean    = formatC(att["mean"],  digits = digits, format = "f"),
    sd      = formatC(att["sd"],    digits = digits, format = "f"),
    lower   = formatC(att["lower"], digits = digits, format = "f"),
    upper   = formatC(att["upper"], digits = digits, format = "f"),
    `Pr>0`  = formatC(att_p_dir,    digits = max(digits, 3), format = "f"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(att_print)[1] <- ""
  names(att_print)[2] <- "Estimate"
  names(att_print)[3] <- "Est.Error"
  names(att_print)[4] <- ci_lo
  names(att_print)[5] <- ci_hi
  print(att_print, row.names = FALSE, right = TRUE)
  cat("\n")

  # ----------------------------------------------------------------
  # Indirect (spillover) effects (Proposition 6.2).
  #
  # Three blocks, all gated on x$indirect_per_period being non-NULL
  # (which itself is gated on uses_rho = TRUE in .nasc_summary_stats):
  #   1. Per-period total spillover  -- sum over donors per t
  #   2. Average total spillover     -- analog of ATT for spillovers
  #   3. Per-donor average spillover -- ALL donors by default (or top-N
  #      if the user passes a finite max_donors)
  # ----------------------------------------------------------------
  if (!is.null(x$indirect_per_period)) {
    cat("Per-period total indirect Treatment Effect (summed across donors)\n")
    ip <- x$indirect_per_period
    ip_print <- data.frame(
      period = format(ip[[h$time]]),
      mean   = formatC(ip$mean,  digits = digits, format = "f"),
      sd     = formatC(ip$sd,    digits = digits, format = "f"),
      lower  = formatC(ip$lower, digits = digits, format = "f"),
      upper  = formatC(ip$upper, digits = digits, format = "f"),
      `Pr>0` = formatC(ip$p_pos, digits = 3, format = "f"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    names(ip_print)[1] <- "Period"
    names(ip_print)[2] <- "Estimate"
    names(ip_print)[3] <- "Est.Error"
    names(ip_print)[4] <- ci_lo
    names(ip_print)[5] <- ci_hi
    print(ip_print, row.names = FALSE, right = TRUE)
    cat("\n")
  }

  if (!is.null(x$indirect_avg)) {
    ia <- x$indirect_avg
    ia_p_dir <- if (ia["mean"] >= 0) ia["p_pos"] else 1 - ia["p_pos"]
    cat("Average total indirect Tratment Effect (averaged across post-periods)\n")
    ia_print <- data.frame(
      blank  = "",
      mean   = formatC(ia["mean"],  digits = digits, format = "f"),
      sd     = formatC(ia["sd"],    digits = digits, format = "f"),
      lower  = formatC(ia["lower"], digits = digits, format = "f"),
      upper  = formatC(ia["upper"], digits = digits, format = "f"),
      `Pr>0` = formatC(ia_p_dir,    digits = max(digits, 3), format = "f"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    names(ia_print)[1] <- ""
    names(ia_print)[2] <- "Estimate"
    names(ia_print)[3] <- "Est.Error"
    names(ia_print)[4] <- ci_lo
    names(ia_print)[5] <- ci_hi
    print(ia_print, row.names = FALSE, right = TRUE)
    cat("\n")
  }

  if (!is.null(x$indirect_per_donor) && nrow(x$indirect_per_donor) > 0) {
    n_total <- nrow(x$indirect_per_donor)
    n_show  <- min(max_donors, n_total)
    if (is.infinite(max_donors) || n_show >= n_total) {
      cat(sprintf(
        "Per-donor indirect effect (all %d donors, sorted by |posterior mean|, averaged across post-periods)\n",
        n_total
      ))
    } else {
      cat(sprintf(
        "Per-donor indirect effect (top %d of %d by |posterior mean|, averaged across post-periods)\n",
        n_show, n_total
      ))
    }
    pd <- utils::head(x$indirect_per_donor, n_show)
    pd_print <- data.frame(
      donor = pd$donor,
      mean  = formatC(pd$mean,  digits = digits, format = "f"),
      sd    = formatC(pd$sd,    digits = digits, format = "f"),
      lower = formatC(pd$lower, digits = digits, format = "f"),
      upper = formatC(pd$upper, digits = digits, format = "f"),
      `Pr>0`= formatC(pd$p_pos, digits = 3, format = "f"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    names(pd_print)[1] <- "Donor"
    names(pd_print)[2] <- "Estimate"
    names(pd_print)[3] <- "Est.Error"
    names(pd_print)[4] <- ci_lo
    names(pd_print)[5] <- ci_hi
    # max = Inf disables the row-count truncation print() applies to long
    # data.frames, so all donors are actually rendered.
    print(pd_print, row.names = FALSE, right = TRUE, max = .Machine$integer.max)
    if (n_total > n_show) {
      cat(sprintf("  ... %d more donors not shown\n", n_total - n_show))
    }
    cat("\n")
  }

  # Model parameters
  if (!is.null(x$parameters) && nrow(x$parameters) > 0) {
    cat("Estimated model parameters\n")
    p <- x$parameters
    p_print <- data.frame(
      parameter = p$parameter,
      mean      = formatC(p$mean,  digits = digits, format = "f"),
      sd        = formatC(p$sd,    digits = digits, format = "f"),
      lower     = formatC(p$lower, digits = digits, format = "f"),
      upper     = formatC(p$upper, digits = digits, format = "f"),
      stringsAsFactors = FALSE
    )
    names(p_print)[1] <- "Parameter"
    names(p_print)[2] <- "Estimate"
    names(p_print)[3] <- "Est.Error"
    names(p_print)[4] <- ci_lo
    names(p_print)[5] <- ci_hi
    print(p_print, row.names = FALSE, right = TRUE)
    cat("\n")
  }

  # Pre-treatment fit
  cat("Pre-treatment fit\n")
  cat(sprintf("  Pre-period R^2    : %s\n",
              if (is.na(x$pre_r2)) "NA"
              else formatC(x$pre_r2, digits = digits, format = "f")))
  cat(sprintf("  Pre-period RMSE   : %s\n",
              formatC(x$pre_rmse, digits = digits, format = "f")))
  cat(sprintf("  Post-period RMSPE : %s\n",
              if (is.na(x$post_rmse)) "NA"
              else formatC(x$post_rmse, digits = digits, format = "f")))
  cat(sprintf("  RMSPE ratio       : %s\n",
              if (is.na(x$rmspe_ratio)) "NA"
              else formatC(x$rmspe_ratio, digits = digits, format = "f")))
  cat("\n")

  # Donor weights -- show all by default, or top-N when max_donors finite.
  n_w_total <- nrow(x$weights)
  n_w_show  <- min(max_donors, n_w_total)
  if (is.infinite(max_donors) || n_w_show >= n_w_total) {
    cat(sprintf("Donor weights (all %d donors, sorted by posterior mean)\n",
                n_w_total))
  } else {
    cat(sprintf("Donor weights (top %d of %d by posterior mean)\n",
                n_w_show, n_w_total))
  }
  w <- utils::head(x$weights, n_w_show)
  w_print <- data.frame(
    donor = w$donor,
    mean  = formatC(w$mean,  digits = digits, format = "f"),
    sd    = formatC(w$sd,    digits = digits, format = "f"),
    lower = formatC(w$lower, digits = digits, format = "f"),
    upper = formatC(w$upper, digits = digits, format = "f"),
    stringsAsFactors = FALSE
  )
  names(w_print)[1] <- "Donor"
  names(w_print)[2] <- "Estimate"
  names(w_print)[3] <- "Est.Error"
  names(w_print)[4] <- ci_lo
  names(w_print)[5] <- ci_hi
  # See note above on max = .Machine$integer.max.
  print(w_print, row.names = FALSE, right = TRUE, max = .Machine$integer.max)
  if (n_w_total > n_w_show) {
    cat(sprintf("  ... %d more donors not shown\n", n_w_total - n_w_show))
  }
  cat("\n")

  # MCMC diagnostics
  cat("MCMC diagnostics\n")
  if (is.null(x$mcmc)) {
    cat("  (not available -- workers ran in subprocesses or fit not retained)\n")
  } else {
    m <- x$mcmc
    src <- m$source %||% "single_fit"
    if (identical(src, "worker_loop")) {
      cat(sprintf("  Source                : per-worker (%d workers, 1 chain each)\n",
                  m$n_workers))
      cat(sprintf("  Iter / warmup (each)  : %s / %s\n",
                  m$n_iter, m$n_warmup))
      cat(sprintf("  Converged chains      : %d / %d (%.1f%%) at split-Rhat <= %s\n",
                  m$n_converged %||% NA_integer_,
                  m$n_finite_rhat %||% NA_integer_,
                  100 * (m$converged_share %||% NA_real_),
                  formatC(m$rhat_threshold %||% 1.05, digits = 2, format = "f")))
      cat(sprintf("  Worst split-Rhat      : %s\n",
                  formatC(m$max_rhat, digits = 3, format = "f")))
      cat(sprintf("  Min n_eff (worst)     : %s\n",
                  formatC(m$min_n_eff, digits = 0, format = "f")))
      cat(sprintf("  Divergent (total)     : %s\n",
                  if (is.na(m$n_divergent)) "NA" else as.character(m$n_divergent)))
      cat(sprintf("  Max-treedepth (total) : %s\n",
                  if (is.na(m$n_max_td)) "NA" else as.character(m$n_max_td)))
    } else {
      cat(sprintf("  Chains / iter / warmup : %s / %s / %s\n",
                  m$n_chains, m$n_iter, m$n_warmup))
      cat(sprintf("  Max Rhat               : %s\n",
                  formatC(m$max_rhat, digits = 3, format = "f")))
      cat(sprintf("  Min n_eff              : %s\n",
                  formatC(m$min_n_eff, digits = 0, format = "f")))
      cat(sprintf("  Divergent transitions  : %s\n",
                  if (is.na(m$n_divergent)) "NA" else as.character(m$n_divergent)))
    }
    if (identical(src, "worker_loop") && !is.na(m$converged_share) &&
        m$converged_share < 1) {
      cat(sprintf(
        "  ! Only %.1f%% of workers converged (split-Rhat <= %s); the\n",
        100 * m$converged_share,
        formatC(m$rhat_threshold, digits = 2, format = "f")
      ))
      cat("    rho ensemble may include poorly-mixed chains.\n")
    }
    if (!is.na(m$max_rhat) && m$max_rhat > 1.05) {
      cat("  ! Rhat > 1.05 suggests convergence problems.\n")
    }
    if (!is.na(m$n_divergent) && m$n_divergent > 0) {
      cat("  ! Divergent transitions detected; consider raising adapt_delta.\n")
    }
    if (!is.na(m$n_max_td) && m$n_max_td > 0) {
      cat("  ! Max treedepth saturated; consider raising max_treedepth.\n")
    }
  }
  invisible(x)
}

#' Summary method for nascSynth objects
#'
#' Dispatches to the \code{$summary()} method of the
#' \code{\link{nascSynth}} R6 object.
#'
#' @param object A \code{nascSynth} object.
#' @param ... Additional arguments passed to \code{object$summary()}.
#' @return Invisibly returns a \code{summary.nascSynth} list.
#' @export
summary.nascSynth <- function(object, ...) {
  invisible(object$summary(...))
}



# ----------------------------------------------------------------------------
# Internal: pull what we need out of a fitted nascSynth object and rebuild
# the posterior draws of s. Returns a list with donor_names, treated_id,
# w_mat (draws x J), s_mat (draws x J), W_full, rhos_used.
# ----------------------------------------------------------------------------
.nasc_contamination_draws <- function(model) {

  if (!inherits(model, "nascSynth")) {
    stop("'model' must be a fitted nascSynth object.")
  }
  priv <- model$.__enclos_env__$private
  if (is.null(priv$y_synth_draws)) {
    stop("Run $fit() before requesting contamination plots.")
  }
  if (!isTRUE(priv$uses_rho) || is.null(priv$W)) {
    stop(
      "Contamination plots require a spatial weights matrix and a non-trivial ",
      "rho. Re-fit with bias_correction = TRUE or nasc_penalty = TRUE."
    )
  }

  w_mat <- priv$y_synth_draws$w
  if (is.null(w_mat)) {
    stop("Posterior draws of donor weights 'w' are missing from the fit.")
  }

  # Canonical donor ordering: exactly what Stan saw at fit time. Stored on
  # private state by $fit() so post-hoc helpers don't have to re-derive it
  # (and get it wrong by using factor-level order instead of column order
  # from pivot_wider).
  treated_id <- as.character(priv$treated_ids)
  donor_ids  <- priv$donor_ids
  if (is.null(donor_ids)) {
    warning("priv$donor_ids not stored on this fit; falling back to ",
            "factor-level ordering. Donor labels and contamination values ",
            "may be misaligned for older fitted objects -- re-fit to be safe.")
    all_ids   <- levels(priv$data[[rlang::as_name(priv$id)]])
    donor_ids <- setdiff(all_ids, treated_id)
  }

  if (ncol(w_mat) != length(donor_ids)) {
    stop("Mismatch: w_mat has ", ncol(w_mat), " columns but ",
         length(donor_ids), " donors were used at fit time.")
  }
  colnames(w_mat) <- donor_ids

  # Reorder W to match Stan's ordering exactly: donors first (in donor_ids
  # order), treated last. Same logic as in $fit().
  W_full <- as.matrix(priv$W)
  if (is.null(rownames(W_full)) || is.null(colnames(W_full))) {
    all_ids <- levels(priv$data[[rlang::as_name(priv$id)]])
    rownames(W_full) <- colnames(W_full) <- all_ids
  }
  W_full <- W_full[c(donor_ids, treated_id), c(donor_ids, treated_id)]
  J     <- length(donor_ids)
  W_J   <- W_full[seq_len(J), seq_len(J), drop = FALSE]
  # IMPORTANT: w_J1 is the donor-to-treated COLUMN of W (a fixed property
  # of the network), not the SC simplex weights from Stan. The previous
  # implementation used w_mat[d, ] here, which produced a meaningless
  # quantity that happened to have the right shape.
  w_J1  <- as.numeric(W_full[seq_len(J), J + 1])

  # Posterior draws of rho. May come in three shapes depending on the fit
  # path -- see comments inline.
  rhos    <- priv$y_synth_draws$rhos_used
  n_draws <- nrow(w_mat)

  rhos_per_draw <-
    if (length(rhos) == 1L) {
      # Exogenous rho or single-fit Step 2 with one fixed rho.
      rep(rhos, n_draws)
    } else if (length(rhos) == n_draws) {
      # One rho per posterior draw.
      rhos
    } else if (n_draws %% length(rhos) == 0L) {
      # Multi-rho parallel loop: each worker ran with one fixed rho and
      # contributed n_draws/length(rhos) rows to w_mat. Expand each
      # worker's rho across its draws.
      draws_per_worker <- n_draws %/% length(rhos)
      rep(rhos, each = draws_per_worker)
    } else {
      warning("rhos_used has length ", length(rhos),
              " but w has ", n_draws, " draws; using the posterior mean of ",
              "rho for all draws.")
      rep(mean(rhos, na.rm = TRUE), n_draws)
    }

  # The contamination vector
  #   s = rho * (I_J - rho * W_J)^{-1} %*% w_J1
  # depends ONLY on rho (and on the fixed network terms W_J, w_J1), not on
  # the SC weights. Solve once per UNIQUE rho rather than once per
  # posterior draw -- in a multi-rho parallel run with, say, 100 workers
  # this is two orders of magnitude cheaper than the previous loop.
  I_J         <- diag(J)
  unique_rhos <- unique(rhos_per_draw)
  s_by_rho    <- vapply(unique_rhos, function(r) {
    as.numeric(r * solve(I_J - r * W_J, w_J1))
  }, numeric(J))
  if (J == 1L) {
    # vapply returns a vector when FUN.VALUE is length-1; promote to
    # matrix so the lookup below works uniformly.
    s_by_rho <- matrix(s_by_rho, nrow = 1L)
  }

  s_mat <- matrix(NA_real_, nrow = n_draws, ncol = J,
                  dimnames = list(NULL, donor_ids))
  rho_idx <- match(rhos_per_draw, unique_rhos)
  for (d in seq_len(n_draws)) {
    s_mat[d, ] <- s_by_rho[, rho_idx[d]]
  }

  list(
    donor_names = donor_ids,
    treated_id  = treated_id,
    w_mat       = w_mat,
    s_mat       = s_mat,
    W_full      = W_full,
    rhos_used   = rhos_per_draw
  )
}


# ----------------------------------------------------------------------------
# Internal: posterior draws of the indirect (spillover) treatment effect.
#
# By Proposition 6.2 of the proposal, the per-donor spillover at post-period
# t is
#     delta_{i,t}^NASC = s_i / (1 - gamma' s) * tau^SC_{1t}
#                      = s_i * tau^NASC_{1t},
# i.e. the same multiplicative factor s applied to the (bias-corrected) direct
# effect. We reconstruct tau draws here exactly as `.nasc_summary_stats()` and
# `$tauPlot()` do, then multiply by the per-draw contamination vector s to
# obtain the [n_draws x T_post x J] tensor of donor-by-period spillover
# draws.
#
# We also return the period-totals delta_t^total = sum_i delta_{i,t} and the
# per-donor average across post-periods, which are the natural scalars to
# summarize and plot.
#
# Returns NULL when the model is configured without a network (uses_rho =
# FALSE) -- callers should treat this as "indirect effect not defined".
# ----------------------------------------------------------------------------
.nasc_indirect_draws <- function(model) {

  if (!inherits(model, "nascSynth")) {
    stop("'model' must be a fitted nascSynth object.")
  }
  priv <- model$.__enclos_env__$private
  if (is.null(priv$y_synth_draws)) {
    stop("Run $fit() before requesting indirect-effect draws.")
  }
  if (!isTRUE(priv$uses_rho) || is.null(priv$W)) {
    # No network in use -> spillover is identically zero. Returning NULL
    # lets callers cleanly degrade (skip indirect panels in plots, omit
    # indirect blocks from summary output) rather than emit zero-valued
    # noise.
    return(NULL)
  }

  bits <- .nasc_contamination_draws(model)
  s_mat       <- bits$s_mat            # [n_draws x J]
  donor_names <- bits$donor_names

  # Reconstruct tau draws exactly as in .nasc_summary_stats() / tauPlot().
  ycf <- priv$y_synth_draws$y_counterfactual
  bc  <- priv$y_synth_draws$bias_correction
  if (is.null(bc)) bc <- rep(1, ncol(ycf))

  wide_df <- .makeWide(
    data      = priv$data,
    id        = priv$id,
    time      = priv$time,
    outcome   = priv$outcome,
    treatment = priv$treated
  )
  post_data <- wide_df |>
    dplyr::filter(!!priv$time >= priv$intervention)
  Y1_post   <- post_data[[rlang::as_name(priv$outcome)]]
  time_post <- post_data[[rlang::as_name(priv$time)]]

  Y1_mat <- matrix(Y1_post, nrow = nrow(ycf), ncol = length(Y1_post),
                   byrow = TRUE)
  bc_mat <- matrix(as.numeric(bc), nrow = nrow(ycf), ncol = ncol(ycf),
                   byrow = FALSE)
  tau_draws <- (Y1_mat - ycf) * bc_mat   # [n_draws x T_post]

  if (nrow(tau_draws) != nrow(s_mat)) {
    # Should never happen -- both are aligned on Step-2 posterior draws --
    # but if it does, surfacing a clear error beats producing nonsense.
    stop("Internal: tau_draws (", nrow(tau_draws),
         " draws) and s_mat (", nrow(s_mat),
         " draws) are misaligned; cannot compute indirect effect.")
  }

  n_draws <- nrow(tau_draws)
  T_post  <- ncol(tau_draws)
  J       <- ncol(s_mat)

  # 3D tensor: delta[d, t, i] = s_mat[d, i] * tau_draws[d, t].
  # Implementation: outer product per draw via rep + multiplication; an
  # explicit loop over draws is simpler and just as fast at the sizes we
  # see in practice (n_draws on the order of 10^3 - 10^4, T_post and J
  # both small).
  delta_arr <- array(NA_real_, dim = c(n_draws, T_post, J),
                     dimnames = list(NULL, NULL, donor_names))
  for (d in seq_len(n_draws)) {
    delta_arr[d, , ] <- tcrossprod(tau_draws[d, ], s_mat[d, ])
  }

  # Period totals: sum across donors for each (draw, period). Equivalent
  # to tau_draws[d, t] * sum(s_mat[d, ]).
  delta_total <- tau_draws * matrix(rowSums(s_mat),
                                    nrow = n_draws, ncol = T_post,
                                    byrow = FALSE)

  # Average across post-periods, draw by draw -> [n_draws x J] matrix of
  # per-donor "average indirect effect" (analog of ATT for spillovers).
  if (T_post == 1L) {
    avg_per_donor <- matrix(delta_arr[, 1, ], nrow = n_draws, ncol = J,
                            dimnames = list(NULL, donor_names))
  } else {
    # apply(., c(1, 3), mean) collapses dimension 2 (time).
    avg_per_donor <- apply(delta_arr, c(1L, 3L), mean)
    if (!is.matrix(avg_per_donor)) {
      avg_per_donor <- matrix(avg_per_donor, nrow = n_draws, ncol = J)
    }
    colnames(avg_per_donor) <- donor_names
  }

  # Average total indirect effect across post-periods (a single scalar
  # per draw -- the spillover analog of ATT).
  avg_total <- if (T_post == 1L) {
    as.numeric(delta_total[, 1])
  } else {
    rowMeans(delta_total)
  }

  list(
    donor_names    = donor_names,
    time_post      = time_post,
    delta_arr      = delta_arr,       # [n_draws x T_post x J]
    delta_total    = delta_total,     # [n_draws x T_post] (sum over donors)
    avg_per_donor  = avg_per_donor,   # [n_draws x J]
    avg_total      = avg_total        # [n_draws]
  )
}
