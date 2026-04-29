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
  #
  # bias_correction is declared as a scalar `real` in both Stan files, so
  # rstan::extract() returns it as a 1-D array of length n_draws (with a
  # `dim` attribute). A 1-D array does not conform to a 2-D matrix under
  # `*`, so we explicitly broadcast it to an n_draws x n_post matrix and
  # multiply elementwise. This works regardless of whether the input is a
  # plain numeric vector or a 1-D array.
  bc_vec <- as.numeric(bias_correction_draws)  # strip dim attribute
  bc_mat <- matrix(bc_vec, nrow = nrow(y_counterfactual_draws),
                   ncol = ncol(y_counterfactual_draws), byrow = FALSE)
  tau_draws <- (Y1_mat - y_counterfactual_draws) * bc_mat

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

  # Defensive sanity check: y_sim_pre_draws must have one column per real
  # pre-treatment time period. The model1 covariate-stacking feature in
  # fit() inflates y_sim_pre with extra rows for predictor matching, which
  # are stripped before reaching this helper. If the trim wasn't applied
  # (or the caller passed something inconsistent) we surface a clear error
  # here rather than letting tibble::tibble() complain about column sizes.
  if (ncol(y_sim_pre_draws) != length(pre_times)) {
    stop(sprintf(
      "y_sim_pre_draws has %d columns but pre_data has %d rows. ",
      ncol(y_sim_pre_draws), length(pre_times)
    ),
    "If covariates were supplied, y_sim_pre may include rows from the ",
    "Abadie-style predictor matching stack and needs to be trimmed to ",
    "real pre-treatment periods before post-processing.")
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




.plot_tau <- function(data, x, y, ymin, ymax, xintercept) {
  ggplot2::ggplot(data = data, ggplot2::aes(x = !!x)) +
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
}




# ----------------------------------------------------------------------------
# Summary infrastructure
#
# .nasc_summary_stats() builds the structured summary list from a fitted
# nascSynth object (called via the public $summary() method, which forwards
# private state through `parts`). The result is a list with class
# "summary.nascSynth"; print.summary.nascSynth formats it for the console.
# ----------------------------------------------------------------------------

# Small null-coalescing operator (avoids a hard dep on rlang's %||%).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Quantile probabilities from a central credible interval width (e.g. 0.95
# -> c(0.025, 0.975)).
.ci_probs <- function(ci) c((1 - ci) / 2, 1 - (1 - ci) / 2)

# Posterior-summary helper: mean, sd, central CrI, and Pr(x > 0).
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

# MCMC diagnostics from the rstan fit. Returns NULL if anything goes wrong
# (e.g. workers were spawned in subprocesses and last_fit is NULL).
.mcmc_diagnostics <- function(fit) {
  if (is.null(fit)) return(NULL)
  diag <- tryCatch({
    s <- rstan::summary(fit)$summary
    keep <- !rownames(s) %in% c("lp__")
    s <- s[keep, , drop = FALSE]

    # Divergence count from sampler params (post-warmup).
    sp <- tryCatch(
      rstan::get_sampler_params(fit, inc_warmup = FALSE),
      error = function(e) NULL
    )
    n_div <- if (!is.null(sp)) {
      sum(vapply(sp, \(x) sum(x[, "divergent__"]), numeric(1)))
    } else NA_integer_

    list(
      max_rhat   = max(s[, "Rhat"],  na.rm = TRUE),
      min_n_eff  = min(s[, "n_eff"], na.rm = TRUE),
      n_divergent = n_div,
      n_chains   = length(fit@stan_args),
      n_iter     = if (length(fit@stan_args)) fit@stan_args[[1]]$iter else NA_integer_,
      n_warmup   = if (length(fit@stan_args)) fit@stan_args[[1]]$warmup else NA_integer_
    )
  }, error = function(e) NULL)
  diag
}

# Build the summary structure. `parts` is a named list assembled in the
# public method (private fields are not visible from outside the R6 class).
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

  # --- 1. Post-treatment tau draws (n_draws x n_post) ---
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

  # --- 2. ATT (average over post-treatment periods, per draw) ---
  att_draws <- rowMeans(tau_draws)
  att <- .posterior_summary(att_draws, ci)

  # --- 3. Per-period tau table ---
  per_period <- tibble::tibble(
    !!time_nm := post_data[[time_nm]],
    mean   = apply(tau_draws, 2, mean),
    sd     = apply(tau_draws, 2, stats::sd),
    lower  = apply(tau_draws, 2, \(x) stats::quantile(x, .ci_probs(ci)[1], names = FALSE)),
    upper  = apply(tau_draws, 2, \(x) stats::quantile(x, .ci_probs(ci)[2], names = FALSE)),
    p_pos  = apply(tau_draws, 2, \(x) mean(x > 0))
  )

  # --- 4. Pre-treatment fit (RMSE between observed and posterior-mean y_synth) ---
  pre_data <- plot_data |>
    dplyr::filter(!!parts$time < intervention)
  pre_rmse <- if (nrow(pre_data) > 0L) {
    sqrt(mean((pre_data[[outcome_nm]] - pre_data$y_synth)^2))
  } else NA_real_

  # --- 5. Donor weights (n_draws x n_donors -> per-donor summaries) ---
  w_mat <- draws$w
  treated_id  <- as.character(parts$treated_ids)
  all_ids     <- levels(parts$id_levels)
  donor_names <- setdiff(all_ids, treated_id)
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

  # Effective number of donors via exp(entropy) of the posterior-mean weight
  # vector (a soft sparsity measure: equals K under uniform weights, 1 under
  # a one-hot weight).
  w_mean <- weight_tbl$mean
  w_mean <- w_mean / sum(w_mean)  # guard against tiny numerical drift
  entropy <- -sum(ifelse(w_mean > 0, w_mean * log(w_mean), 0))
  eff_donors <- exp(entropy)

  # --- 6. Other model parameters ---
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

  # --- 7. MCMC diagnostics ---
  diag <- .mcmc_diagnostics(parts$fitted)

  out <- list(
    header = list(
      engine          = if (parts$nasc_penalty) "stan_2_NASC" else "model1",
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
    att        = att,
    per_period = per_period,
    pre_rmse   = pre_rmse,
    weights    = weight_tbl,
    eff_donors = eff_donors,
    parameters = param_tbl,
    mcmc       = diag
  )
  class(out) <- c("summary.nascSynth", "list")
  out
}

# S3 print method for the summary object. Plain base-R formatting; no extra
# package dependencies.
#' @export
print.summary.nascSynth <- function(x, digits = 3, max_donors = 10, ...) {
  h <- x$header
  bar <- strrep("=", 72)
  rule <- strrep("-", 72)

  cat(bar, "\n", sep = "")
  cat("Bayesian Network-Aware Synthetic Control\n")
  cat(bar, "\n", sep = "")
  cat(sprintf("  Engine          : %s\n", h$engine))
  cat(sprintf("  Spatial model   : %s\n", h$spatial_model))
  cat(sprintf("  NASC penalty    : %s\n", h$nasc_penalty))
  cat(sprintf("  Bias correction : %s\n", h$bias_correction))
  cat(sprintf("  Rho source      : %s\n", h$rho_source))
  cat(sprintf("  Outcome         : %s\n", h$outcome))
  cat(sprintf("  Treated unit    : %s\n", h$treated_unit))
  cat(sprintf("  Intervention at : %s = %s\n", h$time, format(h$intervention)))
  cat(sprintf("  Periods         : %d pre / %d post\n", h$n_pre, h$n_post))
  cat(sprintf("  Donors          : %d\n", h$n_donors))
  cat(sprintf("  Posterior draws : %d\n", h$n_draws))
  cat(sprintf("  Credible level  : %.0f%%\n", 100 * h$ci_width))
  cat(rule, "\n", sep = "")

  # ---- ATT ----
  ci_lo <- 100 * (1 - h$ci_width) / 2
  ci_hi <- 100 - ci_lo
  att_dir   <- if (x$att["mean"] >= 0) ">" else "<"
  att_p_dir <- if (x$att["mean"] >= 0) x$att["p_pos"] else 1 - x$att["p_pos"]
  cat("Average treatment effect on the treated (ATT)\n")
  cat(sprintf("  Posterior mean   : %s\n",  formatC(x$att["mean"], digits = digits, format = "f")))
  cat(sprintf("  Posterior SD     : %s\n",  formatC(x$att["sd"],   digits = digits, format = "f")))
  cat(sprintf("  %.1f%% CrI         : [%s, %s]\n",
              100 * h$ci_width,
              formatC(x$att["lower"], digits = digits, format = "f"),
              formatC(x$att["upper"], digits = digits, format = "f")))
  cat(sprintf("  Pr(ATT %s 0)      : %s\n",
              att_dir,
              formatC(att_p_dir, digits = max(digits, 3), format = "f")))
  cat(rule, "\n", sep = "")

  # ---- Per-period tau ----
  cat("Per-period treatment effect (tau)\n")
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
  names(pp_print)[1] <- h$time
  names(pp_print)[4] <- sprintf("%.1f%%", ci_lo)
  names(pp_print)[5] <- sprintf("%.1f%%", ci_hi)
  print(pp_print, row.names = FALSE, right = TRUE)
  cat(rule, "\n", sep = "")

  # ---- Pre-treatment fit ----
  cat("Pre-treatment fit\n")
  cat(sprintf("  RMSE (observed vs synthetic) : %s\n",
              formatC(x$pre_rmse, digits = digits, format = "f")))
  cat(rule, "\n", sep = "")

  # ---- Donor weights ----
  cat(sprintf("Donor weights (top %d by posterior mean)\n",
              min(max_donors, nrow(x$weights))))
  w <- utils::head(x$weights, max_donors)
  w_print <- data.frame(
    donor = w$donor,
    mean  = formatC(w$mean,  digits = digits, format = "f"),
    sd    = formatC(w$sd,    digits = digits, format = "f"),
    lower = formatC(w$lower, digits = digits, format = "f"),
    upper = formatC(w$upper, digits = digits, format = "f"),
    stringsAsFactors = FALSE
  )
  names(w_print)[4] <- sprintf("%.1f%%", ci_lo)
  names(w_print)[5] <- sprintf("%.1f%%", ci_hi)
  print(w_print, row.names = FALSE, right = TRUE)
  cat(sprintf("  Effective # of donors (exp-entropy) : %s\n",
              formatC(x$eff_donors, digits = digits, format = "f")))
  if (nrow(x$weights) > max_donors) {
    cat(sprintf("  ... %d more donors not shown\n",
                nrow(x$weights) - max_donors))
  }
  cat(rule, "\n", sep = "")

  # ---- Other parameters ----
  if (!is.null(x$parameters) && nrow(x$parameters) > 0) {
    cat("Other model parameters\n")
    p <- x$parameters
    p_print <- data.frame(
      parameter = p$parameter,
      mean      = formatC(p$mean,  digits = digits, format = "f"),
      sd        = formatC(p$sd,    digits = digits, format = "f"),
      lower     = formatC(p$lower, digits = digits, format = "f"),
      upper     = formatC(p$upper, digits = digits, format = "f"),
      stringsAsFactors = FALSE
    )
    names(p_print)[4] <- sprintf("%.1f%%", ci_lo)
    names(p_print)[5] <- sprintf("%.1f%%", ci_hi)
    print(p_print, row.names = FALSE, right = TRUE)
    cat(rule, "\n", sep = "")
  }

  # ---- MCMC diagnostics ----
  cat("MCMC diagnostics\n")
  if (is.null(x$mcmc)) {
    cat("  (not available -- workers ran in subprocesses or fit not retained)\n")
  } else {
    m <- x$mcmc
    cat(sprintf("  Chains / iter / warmup : %s / %s / %s\n",
                m$n_chains, m$n_iter, m$n_warmup))
    cat(sprintf("  Max Rhat               : %s\n",
                formatC(m$max_rhat, digits = 3, format = "f")))
    cat(sprintf("  Min n_eff              : %s\n",
                formatC(m$min_n_eff, digits = 0, format = "f")))
    cat(sprintf("  Divergent transitions  : %s\n",
                if (is.na(m$n_divergent)) "NA" else as.character(m$n_divergent)))
    # Quick warnings
    if (!is.na(m$max_rhat) && m$max_rhat > 1.05) {
      cat("  ! Rhat > 1.05 suggests convergence problems.\n")
    }
    if (!is.na(m$n_divergent) && m$n_divergent > 0) {
      cat("  ! Divergent transitions detected; consider raising adapt_delta.\n")
    }
  }
  cat(bar, "\n", sep = "")
  invisible(x)
}

# S3 method dispatch for summary() on nascSynth R6 objects. This calls back
# into the R6 method so the implementation stays in one place. Returns
# invisibly because the R6 method already prints by default; otherwise the
# REPL would auto-print the returned list a second time.
#' @export
summary.nascSynth <- function(object, ...) {
  invisible(object$summary(...))
}
