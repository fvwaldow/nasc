# nasc_contamination_plots.R
#
#' Network graph coloured by donor contamination
#'
#' @description
#' Plots the spatial-weights network underlying a fitted \code{nascSynth}
#' object as a graph, with each donor node coloured by its posterior mean
#' contamination \eqn{|s_j|}. The treated unit is drawn as a distinct
#' marker and is not coloured on the contamination scale (its row of
#' \eqn{W} feeds donors but it has no weight of its own).
#'
#' Edges follow the row-standardised weights matrix \eqn{W}. By default the
#' graph is treated as undirected and edges weaker than \code{edge_threshold}
#' are dropped to keep the plot readable.
#'
#' @param model A fitted \code{nascSynth} object.
#' @param signed Logical. If \code{TRUE}, nodes are coloured by signed
#'   \eqn{s_j} on a diverging palette; if \code{FALSE} (default), by
#'   \eqn{|s_j|} on a sequential palette.
#' @param edge_threshold Numeric. Edges with absolute weight below this
#'   value are not drawn. Default \code{1e-3}.
#' @param layout An \code{igraph} layout function or a two-column matrix
#'   of node coordinates. Default \code{igraph::layout_with_fr}.
#' @param vertex_size Numeric scalar for node size. Default \code{12}.
#' @param label_cex Numeric scalar for label size. Default \code{0.8}.
#' @param directed Logical. If \code{TRUE}, treats \eqn{W} as a directed
#'   graph and draws arrows. Default \code{FALSE}.
#'
#' @return Invisibly returns a tibble with one row per donor giving the
#'   posterior mean of \eqn{s_j} and \eqn{|s_j|}.
#'
#' @export
contaminationPlot <- function(model,
                              signed         = FALSE,
                              edge_threshold = 1e-3,
                              layout         = NULL,
                              vertex_size    = 12,
                              label_cex      = 0.8,
                              directed       = FALSE) {

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
  v_names <- igraph::V(g)$name
  is_treated <- v_names == treated_id

  # Build per-donor contamination value used for colouring.
  v_score <- numeric(length(v_names))
  names(v_score) <- v_names

  # Identify which vertices are donors
  donor_v <- v_names[!is_treated]

  # Safely assign by exact name matching
  if (signed) {
    v_score[donor_v] <- s_mean[donor_v]
  } else {
    v_score[donor_v] <- s_abs_mean[donor_v]
  }

  donor_match <- match(v_names, donor_names)
  if (signed) {
    v_score[!is.na(donor_match)] <- s_mean[donor_match[!is.na(donor_match)]]
  } else {
    v_score[!is.na(donor_match)] <- s_abs_mean[donor_match[!is.na(donor_match)]]
  }

  # Map scores to colours.
  if (signed) {
    rng <- max(abs(v_score[!is_treated]), na.rm = TRUE)
    rng <- if (is.finite(rng) && rng > 0) rng else 1
    pal <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
    z   <- pmax(-1, pmin(1, v_score / rng))             # in [-1, 1]
    idx <- round((z + 1) / 2 * 100) + 1                  # 1..101
  } else {
    rng <- max(v_score[!is_treated], na.rm = TRUE)
    rng <- if (is.finite(rng) && rng > 0) rng else 1
    pal <- grDevices::hcl.colors(101, palette = "YlOrRd", rev = TRUE)
    z   <- pmax(0, pmin(1, v_score / rng))
    idx <- round(z * 100) + 1
  }
  v_col <- pal[idx]
  v_col[is_treated] <- "#444444"                         # treated marker

  v_shape <- ifelse(is_treated, "square", "circle")
  v_size  <- ifelse(is_treated, vertex_size * 1.2, vertex_size)

  e_width <- 1

  if (is.null(layout)) layout <- igraph::layout_with_fr

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
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
  yb <- seq(yb_bot, yt, length.out = length(pal) + 1)

  graphics::rect(xl, yb[-length(yb)], xr, yb[-1],
                 col = pal, border = NA, xpd = TRUE)
  graphics::rect(xl, yb_bot, xr, yt,
                 col = NA, border = "gray40", xpd = TRUE)

  # Three tick labels (min, mid, max), placed to the right of the bar.
  if (signed) {
    lab_txt <- formatC(c(-rng, 0, rng), digits = 2, format = "g")
  } else {
    lab_txt <- formatC(c(0, rng / 2, rng), digits = 2, format = "g")
  }
  lab_at <- c(yb_bot, (yb_bot + yt) / 2, yt)
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


#' Network graph coloured by direct and indirect treatment effects
#'
#' @description
#' Plots the spatial-weights network underlying a fitted \code{nascSynth}
#' object as a graph, with each donor node coloured by its posterior
#' mean indirect (spillover) effect \eqn{\bar\delta_j} averaged over
#' the post-treatment periods, and the treated node coloured by the
#' posterior-mean ATT (the direct effect averaged over post-treatment
#' periods). By Proposition 6.2,
#' \eqn{\delta_{j,t}^{NASC} = s_j \cdot \tau_{1t}^{NASC}}, so this plot
#' is the effect-side analog of \code{contaminationPlot()}: where
#' \code{contaminationPlot} shows how \emph{exposed} each donor is to
#' the treated unit through the network, this plot shows the
#' \emph{realized} treatment effect at every node -- direct at the
#' treated unit, indirect at every donor.
#'
#' Direct and indirect effects share a single colour scale so the
#' relative magnitudes are directly comparable on the figure. Because
#' the ATT is typically larger than any single donor's spillover, the
#' treated node usually appears strongly coloured while donors near
#' the periphery of the network appear pale -- the contrast itself
#' communicates the share of the total effect that is direct vs
#' indirect. The treated unit is also drawn as a slightly larger
#' square so it remains visually identifiable.
#'
#' Edges follow the row-standardised weights matrix \eqn{W}. By default
#' the graph is treated as undirected and edges weaker than
#' \code{edge_threshold} are dropped to keep the plot readable.
#'
#' @param model A fitted \code{nascSynth} object. The model must use a
#'   network (\code{uses_rho = TRUE}); otherwise spillover is identically
#'   zero by construction and only the treated unit would carry a
#'   non-trivial colour.
#' @param signed Logical. If \code{TRUE} (default), nodes are coloured
#'   by signed effects on a diverging blue-white-red palette centered at
#'   zero -- the natural choice for treatment effects, where the sign
#'   distinguishes positive from negative effects. If \code{FALSE}, by
#'   absolute effects on a sequential palette.
#' @param edge_threshold Numeric. Edges with absolute weight below this
#'   value are not drawn. Default \code{1e-3}.
#' @param layout An \code{igraph} layout function or a two-column matrix
#'   of node coordinates. Default \code{igraph::layout_with_fr}.
#' @param vertex_size Numeric scalar for node size. The treated unit is
#'   drawn 1.2x this size for visual identification. Default \code{12}.
#' @param label_cex Numeric scalar for label size. Default \code{0.8}.
#' @param directed Logical. If \code{TRUE}, treats \eqn{W} as a directed
#'   graph and draws arrows. Default \code{FALSE}.
#' @param show_values Logical. If \code{TRUE}, append the numeric
#'   effect to each node label (\eqn{\bar\delta_j} for donors, ATT for
#'   the treated unit). Useful for small networks; can clutter larger
#'   ones. Default \code{FALSE}.
#' @param digits Integer. Number of significant digits used when
#'   \code{show_values = TRUE}. Default \code{2}.
#'
#' @return Invisibly returns a tibble with one row per node giving the
#'   posterior mean effect (\code{effect_mean} -- \eqn{\bar\delta_j}
#'   for donors, ATT for the treated unit), its absolute value
#'   (\code{abs_effect_mean}), and a logical flag \code{is_treated}.
#'
#' @export
effectGraph <- function(model,
                        signed         = TRUE,
                        edge_threshold = 1e-3,
                        layout         = NULL,
                        vertex_size    = 12,
                        label_cex      = 0.8,
                        directed       = FALSE,
                        show_values    = FALSE,
                        digits         = 2) {

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
    is.numeric(digits),      length(digits)      == 1L, digits      >= 1
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
    delta_mean     <- setNames(numeric(length(donor_names)), donor_names)
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
  # range is computed jointly so direct and indirect effects sit on
  # the same scale -- this is the whole point of plotting them
  # together. With the ATT typically dominating, the treated node is
  # usually saturated and donors trail off; the contrast itself is
  # informative.
  # ----------------------------------------------------------------
  if (signed) {
    rng <- max(abs(v_score), na.rm = TRUE)
    rng <- if (is.finite(rng) && rng > 0) rng else 1
    pal <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
    z   <- pmax(-1, pmin(1, v_score / rng))
    idx <- round((z + 1) / 2 * 100) + 1
  } else {
    rng <- max(v_score, na.rm = TRUE)
    rng <- if (is.finite(rng) && rng > 0) rng else 1
    pal <- grDevices::hcl.colors(101, palette = "YlOrRd", rev = TRUE)
    z   <- pmax(0, pmin(1, v_score / rng))
    idx <- round(z * 100) + 1
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

  if (is.null(layout)) layout <- igraph::layout_with_fr

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
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
  # so the two figures align side by side). Tick labels are picked
  # from the joint donor+treated range so the bar covers everything
  # actually plotted.
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
  yb <- seq(yb_bot, yt, length.out = length(pal) + 1)

  graphics::rect(xl, yb[-length(yb)], xr, yb[-1],
                 col = pal, border = NA, xpd = TRUE)
  graphics::rect(xl, yb_bot, xr, yt,
                 col = NA, border = "gray40", xpd = TRUE)

  if (signed) {
    lab_txt <- formatC(c(-rng, 0, rng), digits = 2, format = "g")
  } else {
    lab_txt <- formatC(c(0, rng / 2, rng), digits = 2, format = "g")
  }
  lab_at <- c(yb_bot, (yb_bot + yt) / 2, yt)
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


#' Mean contamination vs mean weight scatter plot
#'
#' @description
#' Per-donor scatter plot of posterior-mean synthetic-control weight
#' (\eqn{\bar w_j}) against posterior-mean contamination
#' (\eqn{\bar s_j} or \eqn{\overline{|s_j|}}). Donors that combine large
#' weight with large contamination contribute the most to the NASC penalty
#' \eqn{\langle w, |s|\rangle} and to the bias-correction term \eqn{w's}
#' and are the diagnostic targets of this plot.
#'
#' @param model A fitted \code{nascSynth} object.
#' @param signed Logical. If \code{TRUE}, plot signed \eqn{\bar s_j} on the
#'   x-axis with a vertical reference line at zero; if \code{FALSE}
#'   (default), plot \eqn{\overline{|s_j|}}.
#' @param label Logical. If \code{TRUE} (default), every marker is labelled
#'   with its donor name. The top \code{top_n} contributors are emphasized
#'   in bold. Set to \code{FALSE} to suppress all labels (useful for very
#'   dense donor pools where labels would be unreadable).
#' @param top_n Integer. When \code{label = TRUE}, the \code{top_n} donors
#'   with largest \eqn{\bar w_j \cdot |\bar s_j|} (their contribution to
#'   the penalty) are drawn in bold; remaining donors are still labelled,
#'   in normal weight. Default \code{10}.
#'
#' @return Invisibly returns a tibble with columns \code{donor},
#'   \code{w_mean}, \code{s_mean}, \code{abs_s_mean}, sorted by descending
#'   contribution \eqn{\bar w_j |\bar s_j|}.
#'
#' @export
contaminationScatter <- function(model,
                                 signed = FALSE,
                                 label  = TRUE,
                                 top_n  = 10L) { # Kept to maintain backward compatibility with existing calls

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
