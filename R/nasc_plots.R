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
