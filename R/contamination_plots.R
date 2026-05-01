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
  v_score <- numeric(length(v_names))
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

  # Edge widths proportional to |W_ij|.
  e_w <- abs(igraph::E(g)$weight)
  if (length(e_w) > 0) {
    e_width <- 0.3 + 3 * (e_w / max(e_w, na.rm = TRUE))
  } else {
    e_width <- 1
  }

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

  # Label every marker. Top contributors are drawn in bold; the rest in
  # normal weight. To reduce overlap without adding a dependency on
  # ggrepel, points in the upper half of the y-range are labelled below
  # the marker (pos = 1) and points in the lower half above (pos = 3).
  if (isTRUE(label)) {
    n_donor <- length(donor_names)
    if (n_donor > 0L) {
      tn <- suppressWarnings(as.integer(top_n))
      if (length(tn) != 1L || is.na(tn)) tn <- 0L
      n_top  <- min(max(0L, tn), n_donor)
      is_top <- logical(n_donor)
      if (n_top > 0L) is_top[ord[seq_len(n_top)]] <- TRUE

      # Stagger above/below by y-position to spread labels vertically.
      y_mid <- mean(range(y_vals, na.rm = TRUE))
      pos   <- ifelse(y_vals >= y_mid, 1L, 3L)

      graphics::text(x_vals, y_vals,
                     labels = donor_names,
                     pos    = pos,
                     cex    = ifelse(is_top, 0.85, 0.7),
                     font   = ifelse(is_top, 2L, 1L),
                     offset = 0.4,
                     xpd    = TRUE)
    }
  }

  invisible(out)
}
