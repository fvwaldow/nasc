# Helper Functions


# run module 2 across one or many rho draws
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
      suppressWarnings(
        suppressPackageStartupMessages(require(rstan, quietly = TRUE))
      )
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

    # Filter only the "package built under R version XXX"
    .is_build_version_warning <- function(w) {
      msg <- conditionMessage(w)
      grepl("built under R version", msg, fixed = TRUE)
    }

    results_list <- withCallingHandlers(
      furrr::future_map(
        rhos,
        function(rho) {
          run_worker(
            single_rho    = rho,
            base_data     = base_data,
            step2_mod     = step2_mod,
            rho_field     = rho_field,
            extract_pars  = extract_pars,
            worker_iter   = worker_iter,
            worker_warmup = worker_warmup
          )
        },
        .options = furrr::furrr_options(seed = TRUE, packages = "rstan")
      ),
      warning = function(w) {
        if (.is_build_version_warning(w)) invokeRestart("muffleWarning")
      }
    )

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

# Construct tidy plot data
.get_nasc_results <- function(y_counterfactual_draws, bias_correction_draws, y_sim_pre_draws, pre_data, post_data, time, outcome, ci = 0.75) {
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

# Reshape the long panel into wide
.makeWide <- function(data, id, time, outcome, treatment) {
  id_name <- rlang::as_name(id)
  treated_lab <- data |>
    dplyr::filter(status == "Treated") |>
    dplyr::select(!!id) |>
    dplyr::distinct() |>
    dplyr::pull(!!id)
  treated_lab <- as.character(treated_lab)

  wide_df_treated <- data |>
    dplyr::filter(as.character(!!id) %in% treated_lab) |>
    dplyr::select(!!time, !!treatment, !!outcome)

  wide_df_untreated <- data |>
    dplyr::filter(!(as.character(!!id) %in% treated_lab)) |>
    dplyr::select(!!time, !!outcome, !!id) |>
    tidyr::pivot_wider(
      names_from  = !!id,
      values_from = !!outcome
    )

  dplyr::inner_join(wide_df_treated, wide_df_untreated,
                    by = rlang::as_name(time))
}

# construct module 2 predictor-matching matrix
.build_predictor_matrix <- function(data,
                                    covariates,
                                    id,
                                    time,
                                    outcome,
                                    treated_id,
                                    donor_ids,
                                    intervention,
                                    predictors_op = "mean",
                                    special_predictors = NULL,
                                    time_pred_prior = NULL) {

  id_name      <- rlang::as_name(id)
  time_name    <- rlang::as_name(time)
  outcome_name <- rlang::as_name(outcome)

  if (!is.character(predictors_op) || length(predictors_op) != 1L ||
      is.na(predictors_op) || !nzchar(predictors_op)) {
    stop("'predictors.op' must be a single non-empty character string ",
         "(e.g. \"mean\", \"median\").")
  }
  op_fun <- tryCatch(match.fun(predictors_op),
                     error = function(e) {
                       stop("predictors.op = '", predictors_op,
                            "' could not be resolved to a function: ",
                            conditionMessage(e))
                     })
  agg_fun <- function(x) op_fun(x[!is.na(x)])

  unit_order <- c(donor_ids, treated_id)

  reg_names <- if (is.null(covariates)) {
    character(0)
  } else {
    setdiff(names(covariates), c(time_name, id_name))
  }

  X0_reg <- matrix(numeric(0), nrow = 0L, ncol = length(donor_ids),
                   dimnames = list(NULL, donor_ids))
  X1_reg <- numeric(0)
  reg_labels <- character(0)

  if (length(reg_names) > 0L) {
    cov_pre <- covariates |>
      dplyr::filter(!!time < intervention)
    if (!is.null(time_pred_prior)) {
      cov_pre <- cov_pre |> dplyr::filter(!!time %in% time_pred_prior)
      if (nrow(cov_pre) == 0L) {
        stop("'time.predictors.prior' selects no rows in the covariate panel.")
      }
    }

    cov_agg <- cov_pre |>
      dplyr::group_by(!!id) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(reg_names), agg_fun),
        .groups = "drop"
      )

    cov_wide <- cov_agg |>
      tidyr::pivot_longer(
        cols      = dplyr::all_of(reg_names),
        names_to  = ".predictor",
        values_to = ".value"
      ) |>
      tidyr::pivot_wider(
        names_from  = !!id,
        values_from = .value
      )

    if (!all(unit_order %in% names(cov_wide))) {
      missing <- setdiff(unit_order, names(cov_wide))
      stop("Covariates are missing values for unit(s): ",
           paste(missing, collapse = ", "))
    }

    cov_wide <- cov_wide |> dplyr::arrange(match(.predictor, reg_names))
    reg_labels <- as.character(cov_wide$.predictor)

    X0_reg <- as.matrix(cov_wide[, donor_ids,  drop = FALSE])
    X1_reg <- as.numeric(unlist(cov_wide[, treated_id, drop = TRUE]))

    if (anyNA(X0_reg) || anyNA(X1_reg)) {
      stop("NA values produced when aggregating regular predictors. ",
           "Check the covariate panel and 'time.predictors.prior'.")
    }
  }

  X0_sp <- matrix(numeric(0), nrow = 0L, ncol = length(donor_ids),
                  dimnames = list(NULL, donor_ids))
  X1_sp <- numeric(0)
  sp_labels <- character(0)

  if (!is.null(special_predictors)) {
    if (!is.list(special_predictors)) {
      stop("'special.predictors' must be a list of length-3 lists.")
    }

    .colname_from_ref <- function(ref, src_df, what) {
      if (is.character(ref) && length(ref) == 1L) {
        if (!ref %in% names(src_df)) {
          stop(what, ": column '", ref, "' not found in `data` ",
               "(or `covariates`).")
        }
        return(ref)
      }
      if (is.numeric(ref) && length(ref) == 1L) {
        idx <- as.integer(ref)
        if (idx < 1L || idx > ncol(src_df)) {
          stop(what, ": column index ", idx, " out of range.")
        }
        return(names(src_df)[idx])
      }
      stop(what, ": predictor reference must be a column name (character) ",
           "or column number (numeric scalar).")
    }

    sp_X0_list <- vector("list", length(special_predictors))
    sp_X1_list <- numeric(length(special_predictors))

    for (i in seq_along(special_predictors)) {
      entry <- special_predictors[[i]]
      what  <- sprintf("special.predictors[[%d]]", i)

      if (!is.list(entry) || length(entry) != 3L) {
        stop(what, ": each entry must be a list of length 3 ",
             "(predictor, time periods, operator).")
      }
      ref       <- entry[[1]]
      sp_times  <- entry[[2]]
      sp_op_str <- entry[[3]]

      if (!is.character(sp_op_str) || length(sp_op_str) != 1L ||
          is.na(sp_op_str) || !nzchar(sp_op_str)) {
        stop(what, ": operator must be a single non-empty character string.")
      }
      sp_op_fun <- tryCatch(match.fun(sp_op_str),
                            error = function(e) {
                              stop(what, ": operator '", sp_op_str,
                                   "' could not be resolved to a function.")
                            })
      sp_agg <- function(x) sp_op_fun(x[!is.na(x)])

      if (!is.numeric(sp_times) || length(sp_times) < 1L ||
          anyNA(sp_times)) {
        stop(what, ": time periods must be a non-empty numeric vector ",
             "without NAs.")
      }
      if (any(sp_times >= intervention)) {
        warning(what, ": time periods include post-intervention periods; ",
                "Synth conventionally uses only pre-intervention periods.")
      }

      colname <- NULL
      src_df  <- NULL
      if (!is.null(covariates) && is.character(ref) && length(ref) == 1L &&
          ref %in% names(covariates)) {
        colname <- ref
        src_df  <- covariates
      } else if (is.character(ref) && length(ref) == 1L &&
                 ref %in% names(data)) {
        colname <- ref
        src_df  <- data
      } else {
        if (is.numeric(ref)) {
          colname <- .colname_from_ref(ref, data, what)
          src_df  <- data
        } else {
          stop(what, ": predictor '", ref, "' not found in `data` or `covariates`.")
        }
      }

      sp_pre <- src_df |>
        dplyr::filter(!!time %in% sp_times) |>
        dplyr::select(!!id, !!time, dplyr::all_of(colname))

      if (nrow(sp_pre) == 0L) {
        stop(what, ": no rows match the requested time periods (",
             paste(range(sp_times), collapse = "-"), ").")
      }

      sp_agg_df <- sp_pre |>
        dplyr::group_by(!!id) |>
        dplyr::summarise(.value = sp_agg(.data[[colname]]),
                         .groups = "drop")

      vals <- setNames(sp_agg_df$.value,
                       as.character(sp_agg_df[[id_name]]))
      if (!all(unit_order %in% names(vals))) {
        missing <- setdiff(unit_order, names(vals))
        stop(what, ": missing values for unit(s): ",
             paste(missing, collapse = ", "))
      }

      donor_vec   <- as.numeric(vals[donor_ids])
      treated_val <- as.numeric(vals[treated_id])
      if (anyNA(donor_vec) || is.na(treated_val)) {
        stop(what, ": NA aggregate produced for some units.")
      }

      sp_X0_list[[i]] <- donor_vec
      sp_X1_list[i]   <- treated_val

      tt <- sort(unique(sp_times))
      period_lab <- if (length(tt) == 1L) {
        as.character(tt)
      } else if (all(diff(tt) == 1)) {
        sprintf("%s_%s", tt[1], tt[length(tt)])
      } else {
        paste(tt, collapse = "_")
      }
      sp_labels[i] <- sprintf("%s_%s_%s", colname, sp_op_str, period_lab)
    }

    if (length(sp_X0_list) > 0L) {
      X0_sp <- do.call(rbind, sp_X0_list)
      colnames(X0_sp) <- donor_ids
      X1_sp <- sp_X1_list
    }
  }

  X0 <- rbind(X0_reg, X0_sp)
  X1 <- c(X1_reg, X1_sp)
  nm <- c(reg_labels, sp_labels)

  if (length(nm) > 0L) {
    rownames(X0) <- nm
    names(X1)    <- nm
  }

  list(X1 = X1, X0 = X0, names = nm)
}


# transform predictor weights to numeric vector
.resolve_predictor_weights <- function(predictor_weights, pred_names) {
  K <- length(pred_names)
  if (is.null(predictor_weights)) {
    return(rep(1.0, K))
  }
  if (!is.numeric(predictor_weights) ||
      any(!is.finite(predictor_weights)) ||
      any(predictor_weights < 0)) {
    stop("predictor_weights must be a finite, non-negative numeric vector.")
  }
  pw_names <- names(predictor_weights)
  if (!is.null(pw_names) && all(nzchar(pw_names))) {
    if (!all(pred_names %in% pw_names)) {
      missing <- setdiff(pred_names, pw_names)
      stop("predictor_weights is missing entries for: ",
           paste(missing, collapse = ", "))
    }
    return(as.numeric(predictor_weights[pred_names]))
  }
  if (length(predictor_weights) != K) {
    stop(sprintf(
      "predictor_weights has length %d but %d predictor row(s) were built.",
      length(predictor_weights), K
    ))
  }
  as.numeric(predictor_weights)
}

# Plot TE trajectory with CrI
#
# manage_par = TRUE snapshots and restores the graphics parameters. That
# snapshot includes `mfg`, so restoring it inside a multi-panel layout
# (par(mfrow = ...)) would rewind the panel cursor and the next plot would
# overwrite this one. Callers that set up their own layout pass FALSE.
.plot_tau <- function(data, x, y, ymin, ymax, xintercept, manage_par = TRUE) {
  data <- as.data.frame(data)

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

  if (isTRUE(manage_par)) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))
    graphics::par(bty = "l")
  }

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

# Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# Lower/upper probabilities for CrI
.ci_probs <- function(ci) c((1 - ci) / 2, 1 - (1 - ci) / 2)

# Summarize posterior draw vector
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

# Extract MCMC diagnostics
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

# Repack parallel worker diagnostics
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

# Compute summary statistics for fitted model
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
  )
  weight_tbl <- weight_tbl[rev(seq_len(nrow(weight_tbl))), , drop = FALSE]

  w_mean <- weight_tbl$mean
  w_mean <- w_mean / sum(w_mean)
  entropy <- -sum(ifelse(w_mean > 0, w_mean * log(w_mean), 0))
  eff_donors <- exp(entropy)

  param_rows <- list()
  sigma_draws <- draws$sigma %||% draws$sigma_sc
  if (!is.null(sigma_draws)) {
    param_rows[["sigma"]] <- .posterior_summary(as.numeric(sigma_draws), ci)
  }
  if (!is.null(draws$lambda_tilde)) {
    param_rows[["lambda_tilde"]] <- .posterior_summary(as.numeric(draws$lambda_tilde), ci)
  }
  if (parts$bias_correction && !is.null(draws$bias_correction)) {
    param_rows[["bias_correction"]] <- .posterior_summary(bc_vec, ci)
  }
  rhos <- draws$rhos_used
  if (parts$uses_rho && !all(is.na(rhos))) {
    param_rows[["rho"]] <- .posterior_summary(as.numeric(rhos), ci)
  }

  # Re-scale module 1 coefficients
  .add_step1_coef_rows <- function(par_name, label_prefix) {
    if (is.null(parts$fitted)) return(invisible(NULL))
    arr <- tryCatch(
      rstan::extract(parts$fitted, pars = par_name, permuted = TRUE)[[par_name]],
      error = function(e) NULL
    )
    if (is.null(arr)) return(invisible(NULL))
    if (!is.matrix(arr)) arr <- matrix(arr, ncol = 1L)
    K_b <- ncol(arr)
    if (K_b == 0L) return(invisible(NULL))
    cn <- parts$cov_names
    labs <- if (!is.null(cn) && length(cn) == K_b) {
      sprintf("%s[%s]", label_prefix, cn)
    } else if (K_b > 1L) {
      sprintf("%s[%d]", label_prefix, seq_len(K_b))
    } else {
      label_prefix
    }
    keep <- if (!is.null(parts$beta_identified) &&
                length(parts$beta_identified) == K_b) {
      parts$beta_identified
    } else {
      rep(TRUE, K_b)
    }
    for (k in seq_len(K_b)) {
      if (!isTRUE(keep[k])) next
      col_k <- arr[, k]
      col_k <- col_k[is.finite(col_k)]
      if (length(col_k) < 2L) next
      param_rows[[labs[k]]] <<- .posterior_summary(col_k, ci)
    }
    invisible(NULL)
  }
  .add_step1_coef_rows("beta_orig",  "beta")
  .add_step1_coef_rows("theta_orig", "theta")

  if (!is.null(parts$fitted)) {
    .extract_one <- function(par_name) {
      tryCatch(
        rstan::extract(parts$fitted, pars = par_name, permuted = TRUE)[[par_name]],
        error = function(e) NULL
      )
    }
    sig_draws <- .extract_one("sigma_sar") %||% .extract_one("sigma_sdm")
    if (!is.null(sig_draws)) {
      v <- as.numeric(sig_draws)
      v <- v[is.finite(v)]
      if (length(v) >= 2L) {
        param_rows[["sigma_step1"]] <- .posterior_summary(v, ci)
      }
    }
  }

  param_tbl <- if (length(param_rows)) {
    tibble::as_tibble(do.call(rbind, param_rows), rownames = "parameter")
  } else NULL

  diag <- if (!is.null(draws$worker_diagnostics)) {
    .worker_to_diagnostics(draws$worker_diagnostics)
  } else {
    .mcmc_diagnostics(parts$fitted)
  }


  # Indirect effects
  indirect_per_period_avg   <- NULL  # per-period donor-average spillover
  indirect_per_donor        <- NULL  # average across periods per donor
  indirect_avg              <- NULL  # average across periods of donor-average spillover
  if (isTRUE(parts$uses_rho) && !is.null(parts$model)) {
    ind <- tryCatch(
      .nasc_indirect_draws(parts$model),
      error = function(e) NULL
    )
    if (!is.null(ind)) {
      delta_total_draws <- ind$delta_total                  # n_draws x T_post
      avg_per_donor_dr  <- ind$avg_per_donor                # n_draws x J
      J                 <- ncol(avg_per_donor_dr)

      # Per-period donor-average spillover
      delta_avg_draws <- delta_total_draws / J

      # average indirect effect
      avg_indirect_draws <- rowMeans(delta_avg_draws)         # n_draws

      indirect_per_period_avg <- tibble::tibble(
        !!time_nm := ind$time_post,
        mean   = apply(delta_avg_draws, 2, mean),
        sd     = apply(delta_avg_draws, 2, stats::sd),
        lower  = apply(delta_avg_draws, 2, \(x) stats::quantile(x, .ci_probs(ci)[1], names = FALSE)),
        upper  = apply(delta_avg_draws, 2, \(x) stats::quantile(x, .ci_probs(ci)[2], names = FALSE)),
        p_pos  = apply(delta_avg_draws, 2, \(x) mean(x > 0))
      )

      indirect_per_donor <- tibble::tibble(
        donor = ind$donor_names,
        mean  = apply(avg_per_donor_dr, 2, mean),
        sd    = apply(avg_per_donor_dr, 2, stats::sd),
        lower = apply(avg_per_donor_dr, 2, \(x) stats::quantile(x, .ci_probs(ci)[1], names = FALSE)),
        upper = apply(avg_per_donor_dr, 2, \(x) stats::quantile(x, .ci_probs(ci)[2], names = FALSE)),
        p_pos = apply(avg_per_donor_dr, 2, \(x) mean(x > 0))
      )
      # Sort by descending unit index
      indirect_per_donor <- indirect_per_donor[
        rev(seq_len(nrow(indirect_per_donor))), , drop = FALSE
      ]

      indirect_avg <- .posterior_summary(avg_indirect_draws, ci)
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
    indirect_per_period = indirect_per_period_avg,
    indirect_per_donor  = indirect_per_donor,
    indirect_avg        = indirect_avg,
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

# Print nascSynth summary objects
print.summary.nascSynth <- function(x, digits = 3, ...) {

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

  # Per-period direct TE
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

  # direct ATT
  att <- x$att
  att_p_dir <- if (att["mean"] >= 0) att["p_pos"] else 1 - att["p_pos"]
  cat(if (has_indirect) "Average direct Treatment Effect \n" else "Average Treatment Effect\n")
  att_print <- data.frame(
    blank   = " ",
    mean    = formatC(att["mean"],  digits = digits, format = "f"),
    sd      = formatC(att["sd"],    digits = digits, format = "f"),
    lower   = formatC(att["lower"], digits = digits, format = "f"),
    upper   = formatC(att["upper"], digits = digits, format = "f"),
    `Pr>0`  = formatC(att_p_dir,    digits = max(digits, 3), format = "f"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(att_print)[1] <- " "
  names(att_print)[2] <- "Estimate"
  names(att_print)[3] <- "Est.Error"
  names(att_print)[4] <- ci_lo
  names(att_print)[5] <- ci_hi
  print(att_print, row.names = FALSE, right = TRUE)
  cat("\n")

  if (!is.null(x$indirect_per_period)) {
    ip <- x$indirect_per_period

    cat("Per-period average indirect Treatment Effect (averaged across donors)\n")
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
    cat("Average indirect Treatment Effect (averaged across donors and post-periods)\n")
    ia_print <- data.frame(
      blank  = " ",
      mean   = formatC(ia["mean"],  digits = digits, format = "f"),
      sd     = formatC(ia["sd"],    digits = digits, format = "f"),
      lower  = formatC(ia["lower"], digits = digits, format = "f"),
      upper  = formatC(ia["upper"], digits = digits, format = "f"),
      `Pr>0` = formatC(ia_p_dir,    digits = max(digits, 3), format = "f"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    names(ia_print)[1] <- " "
    names(ia_print)[2] <- "Estimate"
    names(ia_print)[3] <- "Est.Error"
    names(ia_print)[4] <- ci_lo
    names(ia_print)[5] <- ci_hi
    print(ia_print, row.names = FALSE, right = TRUE)
    cat("\n")
  }

  if (!is.null(x$indirect_per_donor) && nrow(x$indirect_per_donor) > 0) {
    cat("Per-donor indirect effect (averaged across post-periods)\n")
    pd <- x$indirect_per_donor
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
    print(pd_print, row.names = FALSE, right = TRUE, max = .Machine$integer.max)
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

  # Donor weights -- always show every donor.
  cat("Donor weights \n")
  w <- x$weights
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

  print(w_print, row.names = FALSE, right = TRUE, max = .Machine$integer.max)
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

# Summary method for nascSynth objects
summary.nascSynth <- function(object, ...) {
  invisible(object$summary(...))
}



# Rebuild posterior draws of contamination vector
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

  W_full <- as.matrix(priv$W)
  if (is.null(rownames(W_full)) || is.null(colnames(W_full))) {
    all_ids <- levels(priv$data[[rlang::as_name(priv$id)]])
    rownames(W_full) <- colnames(W_full) <- all_ids
  }
  W_full <- W_full[c(donor_ids, treated_id), c(donor_ids, treated_id)]
  J     <- length(donor_ids)
  W_J   <- W_full[seq_len(J), seq_len(J), drop = FALSE]
  w_J1  <- as.numeric(W_full[seq_len(J), J + 1])

  rhos    <- priv$y_synth_draws$rhos_used
  n_draws <- nrow(w_mat)

  rhos_per_draw <-
    if (length(rhos) == 1L) {
      # Exogenous rho or single-fit module 2 with one fixed rho.
      rep(rhos, n_draws)
    } else if (length(rhos) == n_draws) {
      # One rho per posterior draw.
      rhos
    } else if (n_draws %% length(rhos) == 0L) {
      # Multi-rho parallel loop
      draws_per_worker <- n_draws %/% length(rhos)
      rep(rhos, each = draws_per_worker)
    } else {
      warning("rhos_used has length ", length(rhos),
              " but w has ", n_draws, " draws; using the posterior mean of ",
              "rho for all draws.")
      rep(mean(rhos, na.rm = TRUE), n_draws)
    }

  I_J         <- diag(J)
  unique_rhos <- unique(rhos_per_draw)
  s_by_rho    <- vapply(unique_rhos, function(r) {
    as.numeric(r * solve(I_J - r * W_J, w_J1))
  }, numeric(J))
  if (J == 1L) {
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


# Rebuild direct-effect draws from a fitted model
.indirect_get_tau_draws <- function(model) {
  priv <- model$.__enclos_env__$private
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
  Y1_post <- post_data[[rlang::as_name(priv$outcome)]]

  Y1_mat <- matrix(Y1_post, nrow = nrow(ycf), ncol = length(Y1_post),
                   byrow = TRUE)
  bc_mat <- matrix(as.numeric(bc), nrow = nrow(ycf), ncol = ncol(ycf),
                   byrow = FALSE)
  (Y1_mat - ycf) * bc_mat
}

# Posterior draws of the indirect (spillover) effect
.nasc_indirect_draws <- function(model) {

  if (!inherits(model, "nascSynth")) {
    stop("'model' must be a fitted nascSynth object.")
  }
  priv <- model$.__enclos_env__$private
  if (is.null(priv$y_synth_draws)) {
    stop("Run $fit() before requesting indirect-effect draws.")
  }
  if (!isTRUE(priv$uses_rho) || is.null(priv$W)) {
    return(NULL)
  }

  bits <- .nasc_contamination_draws(model)
  s_mat       <- bits$s_mat            # n_draws x J
  donor_names <- bits$donor_names

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
  tau_draws <- (Y1_mat - ycf) * bc_mat   # n_draws x T_post

  if (nrow(tau_draws) != nrow(s_mat)) {
    stop("Internal: tau_draws (", nrow(tau_draws),
         " draws) and s_mat (", nrow(s_mat),
         " draws) are misaligned; cannot compute indirect effect.")
  }

  n_draws <- nrow(tau_draws)
  T_post  <- ncol(tau_draws)
  J       <- ncol(s_mat)

  delta_arr <- array(NA_real_, dim = c(n_draws, T_post, J),
                     dimnames = list(NULL, NULL, donor_names))
  for (d in seq_len(n_draws)) {
    delta_arr[d, , ] <- tcrossprod(tau_draws[d, ], s_mat[d, ])
  }

  delta_total <- tau_draws * matrix(rowSums(s_mat),
                                    nrow = n_draws, ncol = T_post,
                                    byrow = FALSE)

  if (T_post == 1L) {
    avg_per_donor <- matrix(delta_arr[, 1, ], nrow = n_draws, ncol = J,
                            dimnames = list(NULL, donor_names))
  } else {
    avg_per_donor <- apply(delta_arr, c(1L, 3L), mean)
    if (!is.matrix(avg_per_donor)) {
      avg_per_donor <- matrix(avg_per_donor, nrow = n_draws, ncol = J)
    }
    colnames(avg_per_donor) <- donor_names
  }

  avg_total <- if (T_post == 1L) {
    as.numeric(delta_total[, 1])
  } else {
    rowMeans(delta_total)
  }

  list(
    donor_names    = donor_names,
    time_post      = time_post,
    delta_arr      = delta_arr,       # n_draws x T_post x J
    delta_total    = delta_total,     # n_draws x T_post
    avg_per_donor  = avg_per_donor,   # n_draws x J
    avg_total      = avg_total        # n_draws
  )
}


# ------------------------------------------------------------------------------
# .nasc_indirect_matrix: posterior-mean indirect (spillover) effect per donor
# and period -- the per-donor reduction of .nasc_indirect_draws()$delta_arr
# that tauPlot()/attPlot() collapse to a donor average.
#
# `pre` controls the pre-treatment periods:
#   "zero"    (default) delta_j(t) = tau(t) * s_j with tau(t) = 0 before
#             treatment, so the spillover is EXACTLY zero. This is imposed by
#             the structure, not estimated -- there is no uncertainty to show.
#   "drop"    post-treatment periods only.
#   "placebo" the pre-treatment fit residual carried through the same
#             contamination vector, delta_j(t) = tau_pre(t) * s_j. This is the
#             exact analogue of what the DIRECT panel shows before treatment
#             (a posterior-predictive residual from y_sim_pre, i.e. a balance
#             diagnostic rather than an effect), but since s_j is a constant
#             per donor it is that same residual rescaled -- it carries no
#             information the direct panel does not already show.
#
# Returns a list with
#   time   : periods (length T)
#   donors : donor labels (length J)
#   mean   : T x J matrix of posterior means (rows = periods)
#   avg    : length-T donor average
#   pre    : the `pre` mode used
# plus `lower`/`upper` matrices when `ci` is supplied. Returns NULL when the fit
# carries no spatial information (uses_rho = FALSE or no W).
# ------------------------------------------------------------------------------
.nasc_indirect_matrix <- function(model, ci = NULL,
                                  pre = c("zero", "drop", "placebo")) {
  pre <- match.arg(pre)

  ind <- tryCatch(.nasc_indirect_draws(model), error = function(e) NULL)
  if (is.null(ind)) return(NULL)

  delta  <- ind$delta_arr                       # n_draws x T_post x J
  T_post <- dim(delta)[2L]
  J      <- dim(delta)[3L]
  donors <- ind$donor_names

  times    <- ind$time_post
  mean_mat <- matrix(apply(delta, c(2L, 3L), mean),
                     nrow = T_post, ncol = J,
                     dimnames = list(NULL, donors))

  lo_mat <- up_mat <- NULL
  if (!is.null(ci)) {
    p <- .ci_probs(ci)
    lo_mat <- matrix(apply(delta, c(2L, 3L), stats::quantile,
                           probs = p[1], names = FALSE),
                     nrow = T_post, ncol = J, dimnames = list(NULL, donors))
    up_mat <- matrix(apply(delta, c(2L, 3L), stats::quantile,
                           probs = p[2], names = FALSE),
                     nrow = T_post, ncol = J, dimnames = list(NULL, donors))
  }

  n_pre <- 0L
  if (!identical(pre, "drop")) {
    priv <- model$.__enclos_env__$private
    wide_df <- .makeWide(
      data      = priv$data,
      id        = priv$id,
      time      = priv$time,
      outcome   = priv$outcome,
      treatment = priv$treated
    )
    pre_data  <- wide_df |>
      dplyr::filter(!!priv$time < priv$intervention)
    pre_times <- pre_data[[rlang::as_name(priv$time)]]
    n_pre     <- length(pre_times)

    if (n_pre > 0L) {
      if (identical(pre, "zero")) {

        pre_mean <- matrix(0, nrow = n_pre, ncol = J,
                           dimnames = list(NULL, donors))
        pre_lo <- pre_up <- pre_mean

      } else {

        y_sim_pre <- priv$y_synth_draws$y_sim_pre       # n_draws x T_pre
        if (is.null(y_sim_pre) || ncol(y_sim_pre) != n_pre) {
          stop("Pre-treatment draws 'y_sim_pre' are missing or do not match ",
               "the pre-period; use pre = 'zero' or pre = 'drop'.")
        }
        s_mat  <- .nasc_contamination_draws(model)$s_mat  # n_draws x J
        Y1_pre <- pre_data[[rlang::as_name(priv$outcome)]]

        # same sign convention as .get_nasc_results(): tau = Y - y_synth
        tau_pre <- matrix(Y1_pre, nrow = nrow(y_sim_pre), ncol = n_pre,
                          byrow = TRUE) - y_sim_pre
        if (nrow(tau_pre) != nrow(s_mat)) {
          stop("Internal: pre-period draws (", nrow(tau_pre),
               ") and contamination draws (", nrow(s_mat),
               ") are misaligned.")
        }

        # E[tau_pre(t) * s_j] over draws, without materialising the array
        pre_mean <- crossprod(tau_pre, s_mat) / nrow(tau_pre)
        dimnames(pre_mean) <- list(NULL, donors)

        if (!is.null(ci)) {
          p <- .ci_probs(ci)
          pre_lo <- pre_up <- matrix(NA_real_, nrow = n_pre, ncol = J,
                                     dimnames = list(NULL, donors))
          for (t in seq_len(n_pre)) {
            slice <- tau_pre[, t] * s_mat                 # n_draws x J
            pre_lo[t, ] <- apply(slice, 2, stats::quantile,
                                 probs = p[1], names = FALSE)
            pre_up[t, ] <- apply(slice, 2, stats::quantile,
                                 probs = p[2], names = FALSE)
          }
        } else {
          pre_lo <- pre_up <- NULL
        }
      }

      times    <- c(pre_times, times)
      mean_mat <- rbind(pre_mean, mean_mat)
      if (!is.null(ci)) {
        lo_mat <- rbind(pre_lo, lo_mat)
        up_mat <- rbind(pre_up, up_mat)
      }
    }
  }

  out <- list(
    time   = times,
    donors = donors,
    mean   = mean_mat,
    avg    = rowMeans(mean_mat),
    pre    = pre,
    n_pre  = n_pre
  )
  if (!is.null(ci)) {
    out$lower <- lo_mat
    out$upper <- up_mat
  }
  out
}


# Empty panel for indirect-effect trajectories; the caller adds the lines.
.plot_indirect_panel <- function(xlim, ylim, xintercept,
                                 xlab = "time",
                                 ylab = "indirect effect",
                                 main = NULL) {
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  graphics::grid(lty = "dotted", col = "gray80")
  graphics::abline(h = 0, lty = 1, col = "black")
  graphics::abline(v = xintercept, lty = 3, col = "gray40")
  if (!is.null(main)) graphics::title(main = main, font.main = 1)
  invisible(NULL)
}


# ------------------------------------------------------------------------------
# .simplex_project: exact Euclidean projection onto the probability simplex
# (Duchi, Shalev-Shwartz, Singer & Chandra, 2008). Always returns a valid
# simplex point.
# ------------------------------------------------------------------------------
.simplex_project <- function(v) {
  n   <- length(v)
  u   <- sort(v, decreasing = TRUE)
  css <- cumsum(u) - 1
  idx <- seq_len(n)
  ok  <- (u - css / idx) > 0
  if (!any(ok)) return(rep(1 / n, n))
  k <- max(idx[ok])
  pmax(v - css[k] / k, 0)
}

# ------------------------------------------------------------------------------
# .penalized_sc_weights: MAP donor weights of the Step-2 model, solved exactly.
#
# The w-dependent part of the Stan log-posterior is
#   -(1/sigma^2) * [ 0.5*||Xs w - ys||^2
#                    + 0.5*sum_k v_k (C1k - C0k'w)^2
#                    + lambda*n_eff*(w's_abs) ]
# so sigma is a positive scale factor and CANCELS from argmin_w: the MAP weights
# are the minimizer of the bracket, a convex QP over the simplex. Solved here by
# FISTA with exact simplex projection. This is deterministic, needs no sigma and
# no Stan optimizer, and ALWAYS returns a valid simplex point -- which is why the
# CV table can no longer contain NA rows. (rstan::optimizing could fail its line
# search whenever the design has an exactly flat direction, e.g. duplicated donor
# columns or lambda = 0 with T0 < J, leaving NA behind.)
# ------------------------------------------------------------------------------
.penalized_sc_weights <- function(Xs, ys, C0, C1, v_cov, s_abs, lam, n_eff,
                                  max_iter = 5000L, tol = 1e-11) {
  J <- ncol(Xs)
  A <- crossprod(Xs)
  b <- as.numeric(crossprod(Xs, ys))
  if (length(v_cov)) {
    for (k in seq_along(v_cov)) {
      if (isTRUE(v_cov[k] > 0)) {
        A <- A + v_cov[k] * tcrossprod(C0[k, ])
        b <- b + v_cov[k] * C1[k] * C0[k, ]
      }
    }
  }
  lin  <- lam * n_eff * s_abs
  L    <- max(abs(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
  if (!is.finite(L) || L <= 0) L <- 1
  step <- 1 / L

  w <- rep(1 / J, J); z <- w; tk <- 1
  for (i in seq_len(max_iter)) {
    g  <- as.numeric(A %*% z) - b + lin
    w1 <- .simplex_project(z - step * g)
    t1 <- (1 + sqrt(1 + 4 * tk^2)) / 2
    z  <- w1 + ((tk - 1) / t1) * (w1 - w)
    if (max(abs(w1 - w)) < tol) { w <- w1; break }
    w <- w1; tk <- t1
  }
  w
}

# ------------------------------------------------------------------------------
# .compute_sigma_ref: reference noise scale used to calibrate the penalty.
#
# The residual scale of the UNPENALIZED (lambda = 0) simplex fit over the full
# pre-period, on the same standardized scale as the Stan model's sigma_sc, and
# including the covariate-matching terms so it matches sigma's MLE there:
#   sigma_ref^2 = ( SSR + sum_k v_k (C1k - C0k'w0)^2 ) / n_eff
# Floored at the model's sigma_floor. Passed to Stan as data, so the penalty is
# calibrated at the scale the data actually exhibit while leaving sigma's own
# conditional untouched.
# ------------------------------------------------------------------------------
.compute_sigma_ref <- function(base_data, sigma_floor = 0.05) {
  J  <- base_data$J
  T0 <- base_data$T0

  Y_panel <- base_data$Y_panel
  y_pre   <- as.numeric(Y_panel[J + 1L, ])
  X_don   <- t(Y_panel[seq_len(J), , drop = FALSE])     # T0 x J

  mean_y <- mean(y_pre)
  sd_y   <- stats::sd(y_pre)
  if (!is.finite(sd_y) || sd_y <= 0) sd_y <- 1
  ys <- (y_pre - mean_y) / sd_y
  Xs <- (X_don - mean_y) / sd_y

  K_cov <- base_data$K_cov
  v_cov <- as.numeric(base_data$v_cov)
  C0 <- matrix(0, nrow = max(K_cov, 1L), ncol = J)
  C1 <- numeric(max(K_cov, 1L))
  if (K_cov > 0) {
    for (k in seq_len(K_cov)) {
      m_k <- mean(base_data$X_cov0[k, ])
      s_k <- stats::sd(base_data$X_cov0[k, ])
      if (is.finite(s_k) && s_k > 1e-12) {
        C0[k, ] <- (base_data$X_cov0[k, ] - m_k) / s_k
        C1[k]   <- (base_data$X_cov1[k]   - m_k) / s_k
      }
    }
  } else {
    v_cov <- numeric(0)
  }

  n_eff <- T0 + sum(v_cov)
  w0 <- tryCatch(
    .penalized_sc_weights(Xs, ys, C0, C1, v_cov, rep(0, J), 0, n_eff),
    error = function(e) rep(1 / J, J))

  ss <- sum((ys - as.numeric(Xs %*% w0))^2)
  if (length(v_cov)) {
    for (k in seq_along(v_cov)) {
      if (isTRUE(v_cov[k] > 0))
        ss <- ss + v_cov[k] * (C1[k] - sum(C0[k, ] * w0))^2
    }
  }
  out <- sqrt(ss / n_eff)
  if (!is.finite(out)) out <- sigma_floor
  max(sigma_floor, out)
}

# ------------------------------------------------------------------------------
# .cv_lambda: single hold-out validation for the NASC penalty strength lambda.
#
# Partition: ONE split of the pre-period. The first `train_frac` of the
# pre-periods are the training window; the remaining periods -- the block
# IMMEDIATELY BEFORE the intervention -- are the validation block. With
# T0 = 20 and train_frac = 0.8 that is train 1-16, validate 17-20.
#
# Rationale: the estimator's real task is to extrapolate the treated unit from
# the pre-period into the post-period, so the validation block should be the
# periods adjacent to treatment. Earlier folds of a rolling-origin scheme
# validate on periods far from the intervention and dilute that signal.
#
# Prediction uses the raw weighted donor combination (no bias correction): there
# is no treatment in the pre-period, so contamination-free forecasting accuracy
# is the right criterion.
#
# base_data: assembled Step-2 data list (penalty branch), BEFORE `lambda` is set.
# rho:       single rho for the contamination vector (exogenous, or the Step-1
#            posterior mean).
# grid:      candidate lambda values.
# train_frac: share of the pre-period used for training.
# ------------------------------------------------------------------------------
.cv_lambda <- function(base_data, rho,
                       grid       = seq(from=0,to=100,by=1),
                       train_frac = 0.8,
                       one_se     = FALSE) {

  J  <- base_data$J
  T0 <- base_data$T0

  # ---- single hold-out split, validation adjacent to the intervention --------
  tr_end <- max(2L, floor(train_frac * T0))
  if (tr_end >= T0) tr_end <- T0 - 1L
  ho_idx <- (tr_end + 1L):T0
  if (tr_end < 2L || !length(ho_idx))
    stop("Pre-period too short to form a train/validation split.")

  Y_panel <- base_data$Y_panel                        # (J+1) x T0
  y_tr    <- as.numeric(Y_panel[J + 1L, ])            # treated pre path
  X_don   <- t(Y_panel[seq_len(J), , drop = FALSE])   # T0 x J

  # contamination vector, same construction as the Stan model
  s_abs <- abs(as.numeric(rho * solve(diag(J) - rho * base_data$W_J,
                                      base_data$w_J1)))

  # standardization exactly as in the Stan transformed-data block, computed on
  # the TRAINING window only
  y_train <- y_tr[seq_len(tr_end)]
  mean_y  <- mean(y_train)
  sd_y    <- stats::sd(y_train)
  if (!is.finite(sd_y) || sd_y <= 0) sd_y <- 1
  ys <- (y_train - mean_y) / sd_y
  Xs <- (X_don[seq_len(tr_end), , drop = FALSE] - mean_y) / sd_y

  K_cov <- base_data$K_cov
  v_cov <- as.numeric(base_data$v_cov)
  C0 <- matrix(0, nrow = max(K_cov, 1L), ncol = J)
  C1 <- numeric(max(K_cov, 1L))
  if (K_cov > 0) {
    for (k in seq_len(K_cov)) {
      m_k <- mean(base_data$X_cov0[k, ])
      s_k <- stats::sd(base_data$X_cov0[k, ])
      if (is.finite(s_k) && s_k > 1e-12) {
        C0[k, ] <- (base_data$X_cov0[k, ] - m_k) / s_k
        C1[k]   <- (base_data$X_cov1[k]   - m_k) / s_k
      }
    }
  } else {
    v_cov <- numeric(0)
  }

  n_eff <- tr_end + sum(v_cov)          # n_eff on the TRAINING window

  # ---- evaluate the grid on the hold-out block ------------------------------
  X_ho <- X_don[ho_idx, , drop = FALSE]
  y_ho <- y_tr[ho_idx]

  n_g   <- length(grid)
  rmse  <- rep(NA_real_, n_g)
  se    <- rep(NA_real_, n_g)
  contam <- rep(NA_real_, n_g)

  for (g in seq_len(n_g)) {
    w <- tryCatch(
      .penalized_sc_weights(Xs, ys, C0, C1, v_cov, s_abs, grid[g], n_eff),
      error = function(e) NULL)
    if (is.null(w) || anyNA(w)) next
    resid     <- y_ho - as.numeric(X_ho %*% w)
    rmse[g]   <- sqrt(mean(resid^2))
    contam[g] <- sum(w * s_abs)
    # hold-out SE of the RMSE, delta method from the squared errors. With a
    # single fold there is no across-fold SD to use, so the one-SE rule draws
    # its scale from within the validation block instead.
    e2 <- resid^2
    se[g] <- if (length(e2) > 1L && mean(e2) > 0)
      (stats::sd(e2) / sqrt(length(e2))) / (2 * sqrt(mean(e2))) else 0
  }

  tbl <- data.frame(lambda = grid, rmse = rmse, se = se,
                    contamination = contam, picked = FALSE)

  if (all(!is.finite(rmse))) {
    warning("Hold-out CV failed for every lambda candidate; using lambda = 0.")
    return(list(lambda = 0, rmse = NA_real_, train = seq_len(tr_end),
                validation = ho_idx, table = tbl))
  }

  # One-standard-error rule, breaking ties toward MORE regularization: among
  # candidates within one hold-out SE of the best RMSE, take the largest lambda.
  # The pre-fit curve is often nearly flat in lambda, so a plain argmin would
  # select lambda ~ 0 on noise alone.
  g_min  <- which.min(rmse)
  thresh <- if (isTRUE(one_se)) rmse[g_min] + se[g_min] else rmse[g_min]
  g_pick <- max(which(is.finite(rmse) & rmse <= thresh))
  tbl$picked[g_pick] <- TRUE

  list(lambda     = grid[g_pick],
       rmse       = rmse[g_pick],
       train      = seq_len(tr_end),
       validation = ho_idx,
       table      = tbl)
}
