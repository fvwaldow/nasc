# Network graph by donor contamination

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

  v_names <- igraph::V(g)$name
  is_treated <- v_names == treated_id

  v_score <- numeric(length(v_names))
  names(v_score) <- v_names

  donor_match <- match(v_names, donor_names)
  if (signed) {
    v_score[!is.na(donor_match)] <- s_mean[donor_match[!is.na(donor_match)]]
  } else {
    v_score[!is.na(donor_match)] <- s_abs_mean[donor_match[!is.na(donor_match)]]
  }

  donor_score  <- v_score[!is_treated]
  finite_donor <- donor_score[is.finite(donor_score)]
  if (signed) {
    v_min <- if (length(finite_donor)) min(finite_donor) else -1
    v_max <- if (length(finite_donor)) max(finite_donor) else  1
    if (!is.finite(v_min) || !is.finite(v_max) || v_min == v_max) {
      eps   <- if (is.finite(v_min) && v_min != 0) abs(v_min) else 1
      v_min <- v_min - eps
      v_max <- v_max + eps
    }
    pal  <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
    span <- max(abs(v_min), abs(v_max))
    z    <- pmax(-1, pmin(1, v_score / span))
    idx  <- round((z + 1) / 2 * 100) + 1
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
  v_col[is_treated] <- "#444444"                         # treated unit

  v_shape <- ifelse(is_treated, "square", "circle")
  v_size  <- ifelse(is_treated, vertex_size * 1.2, vertex_size)

  e_width <- 1

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


# Network graph by direct and indirect treatment effects

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

  # Per-donor average indirect effect
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

  # Treated-unit effect:
  att_mean <- mean(rowMeans(.indirect_get_tau_draws(model)), na.rm = TRUE)

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

  finite_score <- v_score[is.finite(v_score)]
  if (signed) {
    v_min <- if (length(finite_score)) min(finite_score) else -1
    v_max <- if (length(finite_score)) max(finite_score) else  1
    if (!is.finite(v_min) || !is.finite(v_max) || v_min == v_max) {
      eps   <- if (is.finite(v_min) && v_min != 0) abs(v_min) else 1
      v_min <- v_min - eps
      v_max <- v_max + eps
    }
    pal <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)

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

  v_shape <- ifelse(is_treated, "square", "circle")
  v_size  <- ifelse(is_treated, vertex_size * 1.2, vertex_size)

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

  graphics::text(mean(c(xl, xr)), yt,
                 labels = if (signed)
                   expression(bar(tau))
                 else
                   expression(bar("|") * tau * bar("|")),
                 pos = 3, cex = 0.95, xpd = TRUE)

  out <- tibble::tibble(
    node            = c(donor_names, treated_id),
    effect_mean     = c(delta_mean,     att_mean),
    abs_effect_mean = c(delta_abs_mean, abs(att_mean)),
    is_treated      = c(rep(FALSE, length(donor_names)), TRUE)
  )
  invisible(out)
}
