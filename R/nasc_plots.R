# compare_plots.R

#' Compare Multiple NASC Models
#'
#' @description
#' Produces comparative base R plots overlaying the results of multiple fitted
#' \code{nascSynth} objects. Useful for comparing configurations such as
#' Conventional SC, Bias-Corrected SC, NASC, and Fully NASC.
#'
#' @param models A named list of fitted \code{nascSynth} objects. The names of
#'   the list are used as the legend labels.
#' @param show_ci Logical. If \code{TRUE}, adds shaded credible-interval bands
#'   (\code{LB}, \code{UB}, \code{tau_LB}, \code{tau_UB}). Defaults to
#'   \code{FALSE} because overlapping bands can be hard to read.
#'
#' @details
#' Draws plots directly to the active graphics device using base R. Two plots
#' are produced in sequence (synthetic, then effect). Use
#' \code{par(mfrow = c(1, 2))} or \code{dev.new()} between calls if you want
#' them side by side or in separate windows.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   drawing plots.
#'
#' @export
nascPlot <- function(models, show_ci = FALSE) {

  if (!is.list(models) || is.null(names(models))) {
    stop("'models' must be a named list of fitted nascSynth objects.")
  }

  mod1 <- models[[1]]
  if (is.null(mod1$plotData)) {
    stop("The models must be fitted (run $fit()) before plotting.")
  }

  intervention_time <- mod1$interventionTime

  # Standardize each model's data
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

  # ---------------------------------------------------------------
  # Plot 1: Synthetic vs Observed
  # ---------------------------------------------------------------
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
       xlab = "time", ylab = "value")
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

  # ---------------------------------------------------------------
  # Plot 2: Treatment effect (tau)
  # ---------------------------------------------------------------
  y_vals2 <- unlist(lapply(combined_list, function(d) d$tau))
  if (isTRUE(show_ci)) {
    y_vals2 <- c(y_vals2,
                 unlist(lapply(combined_list, function(d) c(d$tau_LB, d$tau_UB))))
  }
  yrng2 <- range(c(y_vals2, 0), na.rm = TRUE)

  plot(combined_list[[1]]$time_var, combined_list[[1]]$tau, type = "n",
       xlim = xrng, ylim = yrng2,
       xlab = "time", ylab = "tau")
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

  invisible(NULL)
}


#' Compare Posterior Donor-Weight Distributions Across Models
#'
#' @description
#' Ridgeline plot of posterior donor-weight densities, with one row per
#' donor and one overlaid translucent density per model. Designed as the
#' multi-model counterpart of \code{$weightDraws()} and as a companion
#' to \code{nascPlot()}: model colours are drawn from the same
#' \code{hcl.colors(palette = "Dark 3")} palette so that lines in
#' \code{nascPlot()} and densities here can be cross-referenced at a
#' glance.
#'
#' @details
#' Donors are listed on the y-axis sorted by their mean posterior weight
#' across the supplied models (heaviest at the top). When models use
#' different donor pools, the union of donors is shown; a model that
#' does not include a particular donor simply has no density drawn in
#' that row.
#'
#' @param models A named list of fitted \code{nascSynth} objects. The
#'   list names are used as legend labels and must be unique.
#' @param overlap Numeric in \code{[0, 1)}. Fraction of vertical overlap
#'   between adjacent donor rows. \code{0} reproduces a non-overlapping
#'   layout; \code{0.5} (default) makes each ridge reach halfway into
#'   the row above. See \code{$weightDraws()} for the same parameter.
#' @param scale Numeric > 0. Multiplicative height of every ridge
#'   relative to the baseline step. Default \code{1.4}.
#' @param fill_alpha Numeric in \code{(0, 1]}. Alpha applied to the
#'   fill colours; lower values make overlapping model ridges more
#'   transparent. Default \code{0.45}, calibrated for readable overlap
#'   at 2-4 models.
#' @param max_donors Integer or \code{NULL}. If a positive integer, only
#'   the \code{max_donors} donors with the largest aggregate mean weight
#'   are plotted. Useful when the donor pool is large and only the
#'   leading contributors matter. Default \code{NULL} (plot all donors).
#'
#' @return Invisibly returns a tibble with one row per (donor, model)
#'   combination giving the posterior mean weight and the share of
#'   draws above \code{1e-3}. Sorted by aggregate descending mean.
#'
#' @export
nascWeight <- function(models,
                          overlap     = 0.5,
                          scale       = 1.4,
                          fill_alpha  = 0.45,
                          max_donors  = NULL) {

  # ---- input validation ---------------------------------------------
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

  # ---- pull each model's draws + donor labels -----------------------
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
    # Canonical donor ordering = colnames(X_pred) at fit time.
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

  # ---- aggregate donor ordering: union sorted by mean weight --------
  all_donors <- unique(unlist(lapply(per_model, `[[`, "donors")))

  # Per-donor aggregate score = mean of per-model means (donors absent
  # in a model get 0 there, which correctly downweights donors that
  # only appear in one configuration).
  donor_score <- vapply(all_donors, function(donor) {
    vals <- vapply(per_model, function(pm) {
      if (donor %in% pm$donors) mean(pm$w_mat[, donor]) else 0
    }, numeric(1))
    mean(vals)
  }, numeric(1))

  donor_order <- all_donors[order(donor_score, decreasing = TRUE)]
  if (!is.null(max_donors) && length(donor_order) > max_donors) {
    donor_order <- donor_order[seq_len(max_donors)]
  }
  # Reverse for plotting: largest weight at TOP of plot. We assign
  # baseline = step * i with i = 1..n, so i = n must be the heaviest.
  donor_plot <- rev(donor_order)
  n_d <- length(donor_plot)

  if (n_d == 0L) {
    stop("No donors to plot.")
  }

  # ---- compute densities --------------------------------------------
  # Index by [donor, model]; entries are NULL where the donor is absent
  # from a given model.
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

  # Global x-range and per-row max y for a coherent layout.
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

  # ---- layout geometry (same idea as $weightDraws()) ----------------
  ridge_h  <- scale * max_y
  step     <- ridge_h * (1 - overlap)
  ylim_top <- step * n_d + ridge_h * 1.05
  ylim_bot <- -ridge_h * 0.05

  # ---- colours: same palette as nascPlot() --------------------------
  cols <- grDevices::hcl.colors(max(n_m, 2), palette = "Dark 3")[seq_len(n_m)]
  fills <- grDevices::adjustcolor(cols, alpha.f = fill_alpha)

  # ---- draw -----------------------------------------------------------
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mar = c(4, 6, 2, 1), bty = "l")

  plot(NA,
       xlim = x_range, ylim = c(ylim_bot, ylim_top),
       xlab = "Donor Weight", ylab = "Donor Unit",
       yaxt = "n")
  graphics::axis(2, at = step * seq_len(n_d), labels = donor_plot, las = 1)

  # Top-down draw order so the bottommost (heaviest) donor's ridge is
  # in front. Within a row, draw all models in the order given so the
  # legend's color order is preserved visually.
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

  # ---- summary tibble for programmatic use ---------------------------
  summary_rows <- list()
  for (donor in donor_order) {  # in descending order, not reversed
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
