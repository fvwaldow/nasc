# Network graph coloured by donor contamination

contaminationGraph <- function(model,
                               signed         = FALSE,
                               edge_threshold = 1e-3,
                               layout         = NULL,
                               vertex_size    = 12,
                               label_cex      = 0.8,
                               directed       = FALSE,
                               seed           = 1L) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for contaminationPlot(). ",
         "Install it via install.packages(\"igraph\").")
  }

  bits        <- .nasc_contamination_draws(model)
  donor_names <- bits$donor_names
  treated_id  <- bits$treated_id
  s_mat       <- bits$s_mat
  W_full      <- bits$W_full

  s_mean     <- colMeans(s_mat, na.rm = TRUE)
  s_abs_mean <- colMeans(abs(s_mat), na.rm = TRUE)

  # Build graph from the donors+treated block of W.
  A <- W_full
  A[abs(A) < edge_threshold] <- 0

  mode <- if (directed) "directed" else "undirected"
  if (!directed) {
    # Symmetrise via mean of W and W' so an edge survives if it appears
    # in either direction (row-standardised matrices are rarely symmetric).
    A <- (A + t(A)) / 2
    A[abs(A) < edge_threshold] <- 0
  }

  g <- igraph::graph_from_adjacency_matrix(
    A, mode = mode, weighted = TRUE, diag = FALSE
  )

  v_names <- igraph::V(g)$name
  is_treated <- v_names == treated_id

  # Build per-donor contamination value used for colouring.
  v_score <- numeric(length(v_names))
  names(v_score) <- v_names

  donor_match <- match(v_names, donor_names)
  if (signed) {
    v_score[!is.na(donor_match)] <- s_mean[donor_match[!is.na(donor_match)]]
  } else {
    v_score[!is.na(donor_match)] <- s_abs_mean[donor_match[!is.na(donor_match)]]
  }

  # ----------------------------------------------------------------
  # Map donor scores to colours. The range is set to the actual
  # [min, max] of the donor scores (treated unit is excluded -- it is
  # rendered in a flat grey marker), so the legend interval shrinks
  # to what is really observed instead of being padded out
  # symmetrically. For the signed palette we still anchor white at
  # zero so the diverging colours retain their meaning; the bar
  # end-points are the true min and max.
  # ----------------------------------------------------------------
  donor_score  <- v_score[!is_treated]
  finite_donor <- donor_score[is.finite(donor_score)]
  if (signed) {
    v_min <- if (length(finite_donor)) min(finite_donor) else -1
    v_max <- if (length(finite_donor)) max(finite_donor) else  1
    # Guard against a degenerate constant score (all equal): give the
    # bar a tiny symmetric spread so colorRampPalette stays well-defined.
    if (!is.finite(v_min) || !is.finite(v_max) || v_min == v_max) {
      eps   <- if (is.finite(v_min) && v_min != 0) abs(v_min) else 1
      v_min <- v_min - eps
      v_max <- v_max + eps
    }
    pal  <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
    # Map score to [-1, 1] with white pinned at 0 on the score scale.
    # If the observed range is one-sided (e.g. all positive), the bar
    # simply starts in the white-to-red half -- still honest.
    span <- max(abs(v_min), abs(v_max))
    z    <- pmax(-1, pmin(1, v_score / span))             # in [-1, 1]
    idx  <- round((z + 1) / 2 * 100) + 1                  # 1..101
  } else {
    v_min <- if (length(finite_donor)) min(finite_donor) else 0
    v_max <- if (length(finite_donor)) max(finite_donor) else 1
    if (!is.finite(v_min) || !is.finite(v_max) || v_min == v_max) {
      v_min <- 0
      v_max <- if (is.finite(v_max) && v_max > 0) v_max else 1
    }
    pal  <- grDevices::hcl.colors(101, palette = "YlOrRd", rev = TRUE)
    span <- v_max - v_min
    z    <- pmax(0, pmin(1, (v_score - v_min) / span))
    idx  <- round(z * 100) + 1
  }
  v_col <- pal[idx]
  v_col[is_treated] <- "#444444"                         # treated marker

  v_shape <- ifelse(is_treated, "square", "circle")
  v_size  <- ifelse(is_treated, vertex_size * 1.2, vertex_size)

  e_width <- 1

  # ----------------------------------------------------------------
  # Layout. FR (and most force-directed layouts) use a random initial
  # configuration, so successive calls yield different pictures. We
  # materialise the layout matrix here under a fixed seed and then
  # pass the matrix -- not the function -- to plot.igraph, which
  # makes the result reproducible across calls. Pass `seed = NULL`
  # to opt out and get the original stochastic behaviour.
  # ----------------------------------------------------------------
  if (is.null(layout)) layout <- igraph::layout_with_fr
  if (is.function(layout)) {
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        if (is.null(old_seed) &&
            exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        } else if (!is.null(old_seed)) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed)
    }
    layout <- layout(g)
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  # Right margin reserved for the colour bar (was 6; the new bar is compact).
  graphics::par(mar = c(1, 1, 2, 4))

  igraph::plot.igraph(
    g,
    layout            = layout,
    vertex.color      = v_col,
    vertex.size       = v_size,
    vertex.shape      = v_shape,
    vertex.frame.color = "gray30",
    vertex.label      = v_names,
    vertex.label.cex  = label_cex,
    vertex.label.color = "black",
    vertex.label.family = "sans",
    edge.width        = e_width,
    edge.color        = grDevices::adjustcolor("gray40", alpha.f = 0.5),
    edge.arrow.size   = if (directed) 0.4 else 0,
    main = ""
  )

  # ----------------------------------------------------------------
  # Compact colour-bar legend, anchored in the upper right margin.
  #
  # Geometry choices:
  #   * The bar sits *outside* the plot area (in the reserved right
  #     margin) so donor nodes and labels never collide with it.
  #   * Bar height = 35% of the plot height, top-aligned. This is
  #     deliberately short -- a tall bar dominates whitespace and
  #     competes visually with the network.
  #   * Bar width is set in inches (independent of x-coordinate
  #     scale, which igraph rescales unpredictably).
  # ----------------------------------------------------------------
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")          # plot region in inches
  cxy <- graphics::par("cxy")          # one character cell in usr units

  # Convert "0.18 inch" bar width into usr x-units.
  bar_width_in <- 0.18
  ux_per_in    <- (usr[2] - usr[1]) / pin[1]
  bar_w        <- bar_width_in * ux_per_in

  # Place bar just outside the right edge, with a small gap.
  gap   <- 0.5 * cxy[1]
  xl    <- usr[2] + gap
  xr    <- xl + bar_w

  # Vertical: 35% of plot height, anchored 5% below the top.
  bar_h_frac <- 0.35
  top_pad    <- 0.05
  yt <- usr[4] - top_pad * (usr[4] - usr[3])
  yb_bot <- yt - bar_h_frac * (usr[4] - usr[3])
  # Slice the palette to the [v_min, v_max] sub-range. For the signed
  # case the palette is parameterised on [-span, span] with white at
  # zero, so we map [v_min, v_max] into palette indices accordingly.
  # For the unsigned case the palette already spans [v_min, v_max].
  if (signed) {
    lo_frac <- (v_min / span + 1) / 2
    hi_frac <- (v_max / span + 1) / 2
  } else {
    lo_frac <- 0
    hi_frac <- 1
  }
  lo_idx  <- max(1L, round(lo_frac * 100) + 1L)
  hi_idx  <- min(length(pal), round(hi_frac * 100) + 1L)
  bar_pal <- pal[lo_idx:hi_idx]

  yb <- seq(yb_bot, yt, length.out = length(bar_pal) + 1)
  graphics::rect(xl, yb[-length(yb)], xr, yb[-1],
                 col = bar_pal, border = NA, xpd = TRUE)
  graphics::rect(xl, yb_bot, xr, yt,
                 col = NA, border = "gray40", xpd = TRUE)

  # Tick labels: bottom = observed min, top = observed max, middle =
  # midpoint of the observed range (or zero for the signed case when
  # zero falls inside the range, which keeps the diverging anchor
  # legible).
  if (signed && v_min < 0 && v_max > 0) {
    mid_val <- 0
    mid_at  <- yb_bot + (0 - v_min) / (v_max - v_min) * (yt - yb_bot)
  } else {
    mid_val <- (v_min + v_max) / 2
    mid_at  <- (yb_bot + yt) / 2
  }
  lab_txt <- formatC(c(v_min, mid_val, v_max), digits = 2, format = "g")
  lab_at  <- c(yb_bot, mid_at, yt)
  graphics::text(xr, lab_at, labels = lab_txt,
                 pos = 4, cex = 0.7, xpd = TRUE)

  # Title above the bar.
  graphics::text(mean(c(xl, xr)), yt,
                 labels = if (signed) "s" else "|s|",
                 pos = 3, cex = 0.85, xpd = TRUE)

  out <- tibble::tibble(
    donor       = donor_names,
    s_mean      = s_mean,
    abs_s_mean  = s_abs_mean
  )
  invisible(out)
}


# Network graph coloured by direct and indirect treatment effects

effectGraph <- function(model,
                        signed         = TRUE,
                        edge_threshold = 1e-3,
                        layout         = NULL,
                        vertex_size    = 12,
                        label_cex      = 0.8,
                        directed       = FALSE,
                        show_values    = FALSE,
                        digits         = 2,
                        seed           = 1L) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for effectGraph(). ",
         "Install it via install.packages(\"igraph\").")
  }
  stopifnot(
    is.logical(signed),      length(signed)      == 1L,
    is.logical(directed),    length(directed)    == 1L,
    is.logical(show_values), length(show_values) == 1L,
    is.numeric(edge_threshold), length(edge_threshold) == 1L,
    edge_threshold >= 0,
    is.numeric(vertex_size), length(vertex_size) == 1L, vertex_size > 0,
    is.numeric(label_cex),   length(label_cex)   == 1L, label_cex   > 0,
    is.numeric(digits),      length(digits)      == 1L, digits      >= 1,
    is.null(seed) || (is.numeric(seed) && length(seed) == 1L)
  )

  # ----------------------------------------------------------------
  # The network and donor labels still come from the contamination
  # helper. We allow uses_rho = FALSE here (unlike the previous
  # spillover-only graph): if there is no network the treated node
  # still carries the ATT, which is informative on its own. We just
  # warn so the caller knows donors will all be zero.
  # ----------------------------------------------------------------
  bits        <- tryCatch(.nasc_contamination_draws(model),
                          error = function(e) NULL)
  if (is.null(bits)) {
    stop(
      "effectGraph() requires a model with a spatial weights matrix W. ",
      "Re-fit with bias_correction = TRUE or nasc_penalty = TRUE so a ",
      "network is in use, or supply W explicitly."
    )
  }
  donor_names <- bits$donor_names
  treated_id  <- bits$treated_id
  W_full      <- bits$W_full

  # Per-donor average indirect effect: zero everywhere when uses_rho is
  # FALSE (spillover is identically zero by construction). Otherwise
  # pull from the indirect-effect helper.
  ind <- .nasc_indirect_draws(model)
  if (is.null(ind)) {
    delta_mean     <- stats::setNames(numeric(length(donor_names)), donor_names)
    delta_abs_mean <- delta_mean
    warning("Model has no rho in use; donor spillovers are identically ",
            "zero. Only the treated unit's ATT carries a non-trivial colour.")
  } else {
    delta_mean     <- colMeans(ind$avg_per_donor, na.rm = TRUE)
    delta_abs_mean <- colMeans(abs(ind$avg_per_donor), na.rm = TRUE)
    names(delta_mean)     <- donor_names
    names(delta_abs_mean) <- donor_names
  }

  # Treated-unit effect: posterior-mean ATT (direct effect averaged
  # over post-periods).
  att_mean <- mean(rowMeans(.indirect_get_tau_draws(model)), na.rm = TRUE)

  # ----------------------------------------------------------------
  # Build graph (same logic as contaminationPlot for consistency).
  # ----------------------------------------------------------------
  A <- W_full
  A[abs(A) < edge_threshold] <- 0

  mode <- if (directed) "directed" else "undirected"
  if (!directed) {
    A <- (A + t(A)) / 2
    A[abs(A) < edge_threshold] <- 0
  }

  g <- igraph::graph_from_adjacency_matrix(
    A, mode = mode, weighted = TRUE, diag = FALSE
  )

  v_names    <- igraph::V(g)$name
  is_treated <- v_names == treated_id

  # ----------------------------------------------------------------
  # Per-node colouring score: signed effect or |effect|. Treated node
  # gets the ATT; donors get delta_mean (or its absolute value).
  # ----------------------------------------------------------------
  v_score <- numeric(length(v_names))
  names(v_score) <- v_names
  donor_match <- match(v_names, donor_names)
  has_donor   <- !is.na(donor_match)
  if (signed) {
    v_score[has_donor]  <- delta_mean[donor_match[has_donor]]
    v_score[is_treated] <- att_mean
  } else {
    v_score[has_donor]  <- delta_abs_mean[donor_match[has_donor]]
    v_score[is_treated] <- abs(att_mean)
  }

  # ----------------------------------------------------------------
  # Single shared palette covering ALL nodes (donors + treated). The
  # range is set to the actual [min, max] of the plotted scores, so
  # the legend interval shrinks to what is really observed instead
  # of being padded out symmetrically. For the signed palette we
  # still anchor white at zero (otherwise the diverging colours lose
  # their meaning), but the bar end-points are the true min and max.
  # ----------------------------------------------------------------
  finite_score <- v_score[is.finite(v_score)]
  if (signed) {
    v_min <- if (length(finite_score)) min(finite_score) else -1
    v_max <- if (length(finite_score)) max(finite_score) else  1
    # Guard against a degenerate constant score (all equal): give the
    # bar a tiny symmetric spread so colorRampPalette stays well-defined.
    if (!is.finite(v_min) || !is.finite(v_max) || v_min == v_max) {
      eps   <- if (is.finite(v_min) && v_min != 0) abs(v_min) else 1
      v_min <- v_min - eps
      v_max <- v_max + eps
    }
    pal <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
    # Map score to [0, 1] with white pinned at 0 on the score scale.
    # If the observed range is one-sided (e.g. all positive), the bar
    # simply starts in the white-to-red half -- still honest.
    span <- max(abs(v_min), abs(v_max))
    z    <- pmax(-1, pmin(1, v_score / span))
    idx  <- round((z + 1) / 2 * 100) + 1
  } else {
    v_min <- if (length(finite_score)) min(finite_score) else 0
    v_max <- if (length(finite_score)) max(finite_score) else 1
    if (!is.finite(v_min) || !is.finite(v_max) || v_min == v_max) {
      v_min <- 0
      v_max <- if (is.finite(v_max) && v_max > 0) v_max else 1
    }
    pal  <- grDevices::hcl.colors(101, palette = "YlOrRd", rev = TRUE)
    span <- v_max - v_min
    z    <- pmax(0, pmin(1, (v_score - v_min) / span))
    idx  <- round(z * 100) + 1
  }
  v_col <- pal[idx]

  # Treated node still distinguished by shape (square) and size (1.2x);
  # the colour now reflects the ATT rather than the grey-marker hack.
  v_shape <- ifelse(is_treated, "square", "circle")
  v_size  <- ifelse(is_treated, vertex_size * 1.2, vertex_size)

  # Optional numeric annotations on node labels.
  v_label <- v_names
  if (isTRUE(show_values)) {
    fmt <- function(x) formatC(x, digits = digits, format = "g")
    for (k in which(has_donor)) {
      v_label[k] <- sprintf("%s\n%s", v_names[k],
                            fmt(delta_mean[donor_match[k]]))
    }
    v_label[is_treated] <- sprintf("%s\nATT=%s", treated_id, fmt(att_mean))
  }

  e_width <- 1

  # ----------------------------------------------------------------
  # Layout. FR (and most force-directed layouts) use a random initial
  # configuration, so successive calls yield different pictures. We
  # materialise the layout matrix here under a fixed seed and then
  # pass the matrix -- not the function -- to plot.igraph, which
  # makes the result reproducible across calls. Pass `seed = NULL`
  # to opt out and get the original stochastic behaviour.
  # ----------------------------------------------------------------
  if (is.null(layout)) layout <- igraph::layout_with_fr
  if (is.function(layout)) {
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        if (is.null(old_seed) &&
            exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        } else if (!is.null(old_seed)) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed)
    }
    layout <- layout(g)
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(1, 1, 2, 4))

  igraph::plot.igraph(
    g,
    layout             = layout,
    vertex.color       = v_col,
    vertex.size        = v_size,
    vertex.shape       = v_shape,
    vertex.frame.color = "gray30",
    vertex.label       = v_label,
    vertex.label.cex   = label_cex,
    vertex.label.color = "black",
    vertex.label.family = "sans",
    edge.width         = e_width,
    edge.color         = grDevices::adjustcolor("gray40", alpha.f = 0.5),
    edge.arrow.size    = if (directed) 0.4 else 0,
    main = ""
  )

  # ----------------------------------------------------------------
  # Compact colour-bar legend (geometry matches contaminationPlot()
  # so the two figures align side by side). The bar now covers only
  # the slice of the palette that lies between the observed min and
  # max scores -- i.e. it shrinks to the data, instead of running
  # the full -rng..+rng span.
  # ----------------------------------------------------------------
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  cxy <- graphics::par("cxy")

  bar_width_in <- 0.18
  ux_per_in    <- (usr[2] - usr[1]) / pin[1]
  bar_w        <- bar_width_in * ux_per_in

  gap   <- 0.5 * cxy[1]
  xl    <- usr[2] + gap
  xr    <- xl + bar_w

  bar_h_frac <- 0.35
  top_pad    <- 0.05
  yt <- usr[4] - top_pad * (usr[4] - usr[3])
  yb_bot <- yt - bar_h_frac * (usr[4] - usr[3])

  # Slice the palette to the [v_min, v_max] sub-range. For the signed
  # case the palette is parameterised on [-span, span] with white at
  # zero, so we map [v_min, v_max] into palette indices accordingly.
  # For the unsigned case the palette already spans [v_min, v_max].
  if (signed) {
    lo_frac <- (v_min / span + 1) / 2
    hi_frac <- (v_max / span + 1) / 2
  } else {
    lo_frac <- 0
    hi_frac <- 1
  }
  lo_idx  <- max(1L, round(lo_frac * 100) + 1L)
  hi_idx  <- min(length(pal), round(hi_frac * 100) + 1L)
  bar_pal <- pal[lo_idx:hi_idx]

  yb <- seq(yb_bot, yt, length.out = length(bar_pal) + 1)
  graphics::rect(xl, yb[-length(yb)], xr, yb[-1],
                 col = bar_pal, border = NA, xpd = TRUE)
  graphics::rect(xl, yb_bot, xr, yt,
                 col = NA, border = "gray40", xpd = TRUE)

  # Tick labels: bottom = observed min, top = observed max, middle =
  # midpoint of the observed range (or zero for the signed case when
  # zero falls inside the range, which keeps the diverging anchor
  # legible).
  if (signed && v_min < 0 && v_max > 0) {
    mid_val <- 0
    mid_at  <- yb_bot + (0 - v_min) / (v_max - v_min) * (yt - yb_bot)
  } else {
    mid_val <- (v_min + v_max) / 2
    mid_at  <- (yb_bot + yt) / 2
  }
  lab_txt <- formatC(c(v_min, mid_val, v_max), digits = 2, format = "g")
  lab_at  <- c(yb_bot, mid_at, yt)
  graphics::text(xr, lab_at, labels = lab_txt,
                 pos = 4, cex = 0.7, xpd = TRUE)

  # Title above the bar: tau-bar marks "average treatment effect"
  # (ATT for the treated unit, delta-bar for donors -- both share
  # the scale so a single label covers them).
  graphics::text(mean(c(xl, xr)), yt,
                 labels = if (signed)
                   expression(bar(tau))
                 else
                   expression(bar("|") * tau * bar("|")),
                 pos = 3, cex = 0.95, xpd = TRUE)

  # Return tibble: one row per node, donors first then treated. The
  # `effect_mean` column carries delta-bar for donors and ATT for the
  # treated unit -- exactly the quantity plotted.
  out <- tibble::tibble(
    node            = c(donor_names, treated_id),
    effect_mean     = c(delta_mean,     att_mean),
    abs_effect_mean = c(delta_abs_mean, abs(att_mean)),
    is_treated      = c(rep(FALSE, length(donor_names)), TRUE)
  )
  invisible(out)
}


# ----------------------------------------------------------------------------
# Internal: rebuild tau (direct-effect) draws from a fitted nascSynth.
#
# Used by effectGraph() to surface the ATT as the treated-unit colour
# and annotation. The same reconstruction is implemented inline in
# $tauPlot(), $attPlot(), $effectPlot() and .nasc_indirect_draws();
# centralizing it here would be a good follow-up refactor but is out
# of scope for this patch.
# ----------------------------------------------------------------------------
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


# Mean contamination vs mean weight scatter plot

contaminationScatter <- function(model,
                                 signed = FALSE,
                                 label  = TRUE,
                                 top_n  = 10L) {

  bits        <- .nasc_contamination_draws(model)
  donor_names <- bits$donor_names
  w_mat       <- bits$w_mat
  s_mat       <- bits$s_mat

  w_mean     <- colMeans(w_mat, na.rm = TRUE)
  s_mean     <- colMeans(s_mat, na.rm = TRUE)
  s_abs_mean <- colMeans(abs(s_mat), na.rm = TRUE)

  contrib <- w_mean * abs(s_mean)
  ord <- order(contrib, decreasing = TRUE)

  out <- tibble::tibble(
    donor      = donor_names,
    w_mean     = w_mean,
    s_mean     = s_mean,
    abs_s_mean = s_abs_mean,
    contrib    = contrib
  )[ord, ]

  x_vals <- if (signed) s_mean else s_abs_mean
  y_vals <- w_mean

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mar = c(4, 5, 2, 1), bty = "l")

  xlab <- if (signed) expression(bar(s)[j]) else expression(bar("|s|")[j])
  ylab <- expression(bar(w)[j])

  plot(x_vals, y_vals, type = "n",
       xlab = xlab, ylab = ylab,
       main = "")
  graphics::grid(lty = "dotted", col = "gray80")
  if (signed) graphics::abline(v = 0, lty = 2, col = "gray50")
  graphics::abline(h = 0, lty = 2, col = "gray50")

  graphics::points(x_vals, y_vals,
                   pch = 1, col = "gray30", cex = 1.4, lwd = 1.2)

  # Label every marker uniformly in normal weight.
  # To reduce overlap without adding a dependency on
  # ggrepel, points in the upper half of the y-range are labelled below
  # the marker (pos = 1) and points in the lower half above (pos = 3).
  if (isTRUE(label)) {
    n_donor <- length(donor_names)
    if (n_donor > 0L) {

      # Stagger above/below by y-position to spread labels vertically.
      y_mid <- mean(range(y_vals, na.rm = TRUE))
      pos   <- ifelse(y_vals >= y_mid, 1L, 3L)

      graphics::text(x_vals, y_vals,
                     labels = donor_names,
                     pos    = pos,
                     cex    = 0.75, # Uniform size for all labels
                     font   = 1L,   # Uniform normal weight for all labels
                     offset = 0.4,
                     xpd    = TRUE)
    }
  }

  invisible(out)
}
