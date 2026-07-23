# Comparison Plots of Multiple NASC Models

nascPlot <- function(models, show_ci = FALSE, indirect = TRUE, show_avg = TRUE) {

  if (!is.list(models) || is.null(names(models))) {
    stop("'models' must be a named list of fitted nascSynth objects.")
  }

  mod1 <- models[[1]]
  if (is.null(mod1$plotData)) {
    stop("The models must be fitted (run $fit()) before plotting.")
  }

  intervention_time <- mod1$interventionTime

  combined_list <- lapply(names(models), function(m_name) {
    df <- models[[m_name]]$plotData
    if (is.null(df)) {
      stop(sprintf("Model '%s' is missing plotData. Has it been fitted?", m_name))
    }
    df <- as.data.frame(df)
    colnames(df)[1:2] <- c("time_var", "outcome_var")
    df$Model <- m_name
    df <- df[order(df$time_var), , drop = FALSE]
    df
  })
  names(combined_list) <- names(models)

  n_models <- length(models)
  cols <- grDevices::hcl.colors(max(n_models, 2), palette = "Dark 3")[seq_len(n_models)]


  # Plot Synthetic vs Observed
  obs_df <- combined_list[[1]][, c("time_var", "outcome_var")]

  y_vals <- c(obs_df$outcome_var,
              unlist(lapply(combined_list, function(d) d$y_synth)))
  if (isTRUE(show_ci)) {
    y_vals <- c(y_vals,
                unlist(lapply(combined_list, function(d) c(d$LB, d$UB))))
  }
  yrng <- range(y_vals, na.rm = TRUE)
  xrng <- range(obs_df$time_var, na.rm = TRUE)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(bty = "l")

  plot(obs_df$time_var, obs_df$outcome_var, type = "n",
       xlim = xrng, ylim = yrng,
       xlab = "time", ylab = "Y")
  graphics::grid(lty = "dotted", col = "gray80")

  if (isTRUE(show_ci)) {
    for (i in seq_len(n_models)) {
      d <- combined_list[[i]]
      graphics::polygon(c(d$time_var, rev(d$time_var)),
                        c(d$LB, rev(d$UB)),
                        col = grDevices::adjustcolor(cols[i], alpha.f = 0.15),
                        border = NA)
    }
  }

  graphics::lines(obs_df$time_var, obs_df$outcome_var,
                  lwd = 2, col = "black", lty = 1)
  for (i in seq_len(n_models)) {
    d <- combined_list[[i]]
    graphics::lines(d$time_var, d$y_synth,
                    col = cols[i], lwd = 2, lty = 2)
  }
  graphics::abline(v = intervention_time, lty = 3, col = "gray40")

  graphics::legend("topleft",
                   legend = c("Observed", names(models)),
                   col    = c("black", cols),
                   lty    = c(1, rep(2, n_models)),
                   lwd    = 2,
                   bg     = grDevices::adjustcolor("white", alpha.f = 0.85),
                   box.col = "gray70")


  # Plot Direct treatment effect (tau)
  y_vals2 <- unlist(lapply(combined_list, function(d) d$tau))
  if (isTRUE(show_ci)) {
    y_vals2 <- c(y_vals2,
                 unlist(lapply(combined_list, function(d) c(d$tau_LB, d$tau_UB))))
  }
  yrng2 <- range(c(y_vals2, 0), na.rm = TRUE)

  plot(combined_list[[1]]$time_var, combined_list[[1]]$tau, type = "n",
       xlim = xrng, ylim = yrng2,
       xlab = "time", ylab = "effect")
  graphics::grid(lty = "dotted", col = "gray80")

  if (isTRUE(show_ci)) {
    for (i in seq_len(n_models)) {
      d <- combined_list[[i]]
      graphics::polygon(c(d$time_var, rev(d$time_var)),
                        c(d$tau_LB, rev(d$tau_UB)),
                        col = grDevices::adjustcolor(cols[i], alpha.f = 0.15),
                        border = NA)
    }
  }

  for (i in seq_len(n_models)) {
    d <- combined_list[[i]]
    graphics::lines(d$time_var, d$tau, col = cols[i], lwd = 2)
  }
  graphics::abline(h = 0, lty = 1, col = "black")
  graphics::abline(v = intervention_time, lty = 3, col = "gray40")

  graphics::legend("topleft",
                   legend = names(models),
                   col    = cols,
                   lty    = 1, lwd = 2,
                   bg     = grDevices::adjustcolor("white", alpha.f = 0.85),
                   box.col = "gray70")


  # Plot indirect (spillover) effects: one thin line per untreated unit,
  # coloured by model; the model's donor average is drawn bold in the same
  # colour. Models fitted without a rho/W carry no spillovers and are skipped.
  if (isTRUE(indirect)) {
    ind_list <- lapply(models, function(m) {
      if (!inherits(m, "nascSynth")) return(NULL)
      .nasc_indirect_matrix(m)
    })
    names(ind_list) <- names(models)
    has_ind <- !vapply(ind_list, is.null, logical(1))

    if (!any(has_ind)) {
      message("nascPlot: no model carries indirect effects ",
              "(bias_correction / nasc_penalty = FALSE); ",
              "skipping the spillover panel.")
    } else {
      if (any(!has_ind)) {
        message("nascPlot: no indirect effects for model(s): ",
                paste(names(models)[!has_ind], collapse = ", "), ".")
      }

      yrng3 <- range(
        c(0, unlist(lapply(ind_list[has_ind], function(d) as.numeric(d$mean)))),
        na.rm = TRUE
      )

      plot(NA, xlim = xrng, ylim = yrng3,
           xlab = "time", ylab = "indirect effect")
      graphics::grid(lty = "dotted", col = "gray80")

      for (i in seq_len(n_models)) {
        d <- ind_list[[i]]
        if (is.null(d)) next
        ltype <- if (length(d$time_post) == 1L) "p" else "l"
        faded <- grDevices::adjustcolor(cols[i], alpha.f = 0.55)
        for (j in seq_along(d$donors)) {
          graphics::lines(d$time_post, d$mean[, j],
                          col = faded, lwd = 1, type = ltype, pch = 16)
        }
        if (isTRUE(show_avg)) {
          graphics::lines(d$time_post, d$avg,
                          col = cols[i], lwd = 2.5, type = ltype, pch = 16)
        }
      }

      graphics::abline(h = 0, lty = 1, col = "black")
      graphics::abline(v = intervention_time, lty = 3, col = "gray40")

      graphics::legend("topleft",
                       legend = names(models)[has_ind],
                       col    = cols[has_ind],
                       lty    = 1, lwd = 2,
                       bg     = grDevices::adjustcolor("white", alpha.f = 0.85),
                       box.col = "gray70")
    }
  }

  invisible(NULL)
}


# Comparison Plot of Posterior Donor-Weight Distributions Across Models

nascWeight <- function(models,
                       overlap     = 0.5,
                       scale       = 1.4,
                       fill_alpha  = 0.45,
                       max_donors  = NULL) {

  if (!is.list(models) || is.null(names(models)) || any(names(models) == "")) {
    stop("'models' must be a named list of fitted nascSynth objects.")
  }
  if (any(duplicated(names(models)))) {
    stop("Names of 'models' must be unique; got duplicates: ",
         paste(unique(names(models)[duplicated(names(models))]),
               collapse = ", "), ".")
  }
  if (!all(vapply(models, inherits, logical(1), what = "nascSynth"))) {
    stop("Every element of 'models' must be a fitted nascSynth object.")
  }
  stopifnot(
    is.numeric(overlap), length(overlap) == 1L, overlap >= 0, overlap < 1,
    is.numeric(scale),   length(scale)   == 1L, scale > 0,
    is.numeric(fill_alpha), length(fill_alpha) == 1L,
    fill_alpha > 0, fill_alpha <= 1
  )
  if (!is.null(max_donors)) {
    stopifnot(is.numeric(max_donors), length(max_donors) == 1L,
              max_donors >= 1)
    max_donors <- as.integer(max_donors)
  }

  per_model <- lapply(names(models), function(m_name) {
    mod  <- models[[m_name]]
    priv <- mod$.__enclos_env__$private
    if (is.null(priv$y_synth_draws)) {
      stop(sprintf("Model '%s' has not been fitted (run $fit() first).", m_name))
    }
    w_mat <- priv$y_synth_draws$w
    if (is.null(w_mat)) {
      stop(sprintf("Model '%s' has no posterior draws of donor weights.",
                   m_name))
    }

    donor_ids <- priv$donor_ids
    if (is.null(donor_ids)) {
      treated_id  <- as.character(priv$treated_ids)
      all_ids     <- levels(priv$data[[rlang::as_name(priv$id)]])
      donor_ids   <- setdiff(all_ids, treated_id)
      warning(sprintf(
        "Model '%s' lacks priv$donor_ids (older fit?); falling back to factor-level ordering. Donor labels may be misaligned.",
        m_name))
    }
    if (ncol(w_mat) != length(donor_ids)) {
      stop(sprintf("Model '%s': %d weight columns vs %d donors.",
                   m_name, ncol(w_mat), length(donor_ids)))
    }
    colnames(w_mat) <- donor_ids
    list(name = m_name, w_mat = w_mat, donors = donor_ids)
  })

  all_donors <- unique(unlist(lapply(per_model, `[[`, "donors")))

  num_attempt <- suppressWarnings(as.numeric(all_donors))
  if (!any(is.na(num_attempt))) {
    donor_order <- all_donors[order(num_attempt, decreasing = FALSE)]
  } else {
    donor_order <- all_donors[order(as.character(all_donors), decreasing = FALSE)]
  }
  if (!is.null(max_donors) && length(donor_order) > max_donors) {
    keep_n <- as.integer(max_donors)
    donor_order <- donor_order[(length(donor_order) - keep_n + 1L):length(donor_order)]
  }

  donor_plot <- donor_order
  n_d <- length(donor_plot)

  if (n_d == 0L) {
    stop("No donors to plot.")
  }

  n_m <- length(per_model)
  dens_grid <- vector("list", n_d)
  for (i in seq_len(n_d)) {
    donor <- donor_plot[i]
    dens_grid[[i]] <- lapply(per_model, function(pm) {
      if (donor %in% pm$donors) {
        stats::density(pm$w_mat[, donor], na.rm = TRUE)
      } else {
        NULL
      }
    })
  }

  all_x <- unlist(lapply(dens_grid, function(row) {
    unlist(lapply(row, function(d) if (!is.null(d)) d$x else NULL))
  }))
  x_range <- if (length(all_x)) range(all_x) else c(0, 1)
  max_y <- max(vapply(dens_grid, function(row) {
    ys <- vapply(row, function(d) if (!is.null(d)) max(d$y) else 0,
                 numeric(1))
    if (length(ys)) max(ys) else 0
  }, numeric(1)))
  if (!is.finite(max_y) || max_y <= 0) max_y <- 1

  ridge_h  <- scale * max_y
  step     <- ridge_h * (1 - overlap)
  ylim_top <- step * n_d + ridge_h * 1.05
  ylim_bot <- -ridge_h * 0.05

  cols <- grDevices::hcl.colors(max(n_m, 2), palette = "Dark 3")[seq_len(n_m)]
  fills <- grDevices::adjustcolor(cols, alpha.f = fill_alpha)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mar = c(4, 6, 2, 1), bty = "l")

  plot(NA,
       xlim = x_range, ylim = c(ylim_bot, ylim_top),
       xlab = "Donor Weight", ylab = "Donor Unit",
       yaxt = "n")
  graphics::axis(2, at = step * seq_len(n_d), labels = donor_plot, las = 1)

  for (i in seq(n_d, 1L, by = -1L)) {
    baseline <- step * i
    graphics::segments(x_range[1], baseline, x_range[2], baseline,
                       col = "gray60", lwd = 0.5)
    for (m in seq_len(n_m)) {
      d <- dens_grid[[i]][[m]]
      if (is.null(d)) next
      y <- baseline + d$y * (ridge_h / max_y)
      graphics::polygon(
        x = c(d$x, rev(d$x)),
        y = c(y, rep(baseline, length(d$x))),
        col    = fills[m],
        border = cols[m],
        lwd    = 1
      )
    }
  }

  graphics::legend("topright",
                   legend = names(models),
                   fill   = fills,
                   border = cols,
                   bg     = grDevices::adjustcolor("white", alpha.f = 0.85),
                   box.col = "gray70",
                   cex     = 0.9)

  summary_rows <- list()
  for (donor in rev(donor_order)) {
    for (pm in per_model) {
      if (donor %in% pm$donors) {
        col <- pm$w_mat[, donor]
        summary_rows[[length(summary_rows) + 1L]] <- tibble::tibble(
          donor      = donor,
          model      = pm$name,
          mean       = mean(col, na.rm = TRUE),
          sd         = stats::sd(col, na.rm = TRUE),
          share_gt0  = mean(col > 1e-3, na.rm = TRUE)
        )
      }
    }
  }
  out <- if (length(summary_rows)) {
    do.call(rbind, summary_rows)
  } else {
    tibble::tibble(donor = character(), model = character(),
                   mean = numeric(), sd = numeric(),
                   share_gt0 = numeric())
  }

  invisible(out)
}
