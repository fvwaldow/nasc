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
#' @param max_donors Maximum number of donor weights to display.
#'   Default \code{10}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.summary.nascSynth <- function(x, digits = 3, max_donors = 10, ...) {
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

  # Per-period TE
  cat("Per-period TE\n")
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

  # ATT
  att <- x$att
  att_p_dir <- if (att["mean"] >= 0) att["p_pos"] else 1 - att["p_pos"]
  cat("ATT\n")
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

  # Donor weights
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
  names(w_print)[1] <- "Donor"
  names(w_print)[2] <- "Estimate"
  names(w_print)[3] <- "Est.Error"
  names(w_print)[4] <- ci_lo
  names(w_print)[5] <- ci_hi
  print(w_print, row.names = FALSE, right = TRUE)
  if (nrow(x$weights) > max_donors) {
    cat(sprintf("  ... %d more donors not shown\n",
                nrow(x$weights) - max_donors))
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
