

rm(list = ls())
library(nasc)
library(dplyr)


source("C:/Users/frede/Desktop/Master Thesis/Simulation_R/NASC Estimator/DGP functions 2.R")



# -----------------------------------------------------------------------------
# 1. Generate one dataset
# -----------------------------------------------------------------------------
sim <- generate_data_v2(
  type             = "SAR",
  N                = 20,
  T                = 30,
  T_0              = 20,
  treated_idx      = 1,
  # Effects
  beta             = c(0.5, 1.0, 0.8),
  tau_global       = 0.05,
  delta            = 5,
  # Spatial
  rho              = -0.8,
  # Noise & factor structure (tuned for high pre-period R²)
  sigma_e          = 0.2,
  n_factors        = 3,
  factor_sd        = 1.0,
  loading_sd       = 1.0,
  # Best-donors design
  n_best_donors            = 6,
  include_isolated_in_best = TRUE,
  # Network
  connectivity     = "low",
  treated_position = "peripheral",
  allow_isolated   = TRUE,
  n_isolated       = 3,
  seed             = 123
)


# -----------------------------------------------------------------------------
# 2. Fit both estimators
# -----------------------------------------------------------------------------
df <- sim$df
W  <- sim$weights$W
covariates <- df %>%
  dplyr::select(time, id, X1, X2)


# Standard SC
fit_plain <- nascSynth$new(
  data            = df,
  time            = time,
  id              = id,
  treated         = D,
  outcome         = Y,
  covariates = covariates,
  W               = W,
  spatial_model   = "SAR",     # Step 1 still estimates rho, but Step 2 ignores it
  bias_correction = FALSE,
  nasc_penalty    = FALSE,
  ci_width   = 0.95,
)
fit_plain$fit(n_samples = 100, cores = 4)
sum_plain <- fit_plain$summary(print = FALSE)
fit_plain$syntheticPlot()

# Full NASC
fit_nasc <- nascSynth$new(
  data            = df,
  time            = time,
  id              = id,
  treated         = D,
  outcome         = Y,
  covariates = covariates,
  W               = W,
  spatial_model   = "SAR",
  bias_correction = TRUE,
  nasc_penalty    = TRUE,
  ci_width   = 0.95
)
fit_nasc$fit(n_samples = 100, cores = 4)
sum_nasc <- fit_nasc$summary(print = FALSE)
fit_nasc$syntheticPlot()

models_list <- list(
  "Conventional SC" = fit_plain,
  "Fully NASC"      = fit_nasc
)
nascPlot(models_list, show_ci = TRUE)
contaminationPlot(fit_nasc)
fit_nasc$weightDraws()
fit_nasc$posteriorPlot()

fit_plain$summary()
fit_nasc$summary()
cat(sprintf("True ATT: %.4f\n\n", mean(sim$true_att)))






  plot_network <- function(weights,
                           treated_idx     = 1,
                           directed        = FALSE,
                           edge_by_weight  = TRUE,
                           show_isolated   = TRUE,
                           main            = "Spatial weight matrix") {

    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Package 'igraph' is required. Install with install.packages('igraph').")
    }

    W      <- weights$W
    coords <- weights$coords
    N      <- nrow(W)

    # Identify isolated units. Use the constructor's diagnostic if present;
    # otherwise compute on the fly from row/col sums of the (unweighted) graph.
    iso <- weights$isolated_units
    if (is.null(iso)) {
      A <- (W > 0) * 1
      iso <- which(rowSums(A) == 0 & colSums(A) == 0)
    }

    # Build the graph
    mode <- if (directed) "directed" else "undirected"
    g <- igraph::graph_from_adjacency_matrix(W, mode = mode, weighted = TRUE,
                                             diag = FALSE)

    # Vertex styling -----------------------------------------------------------
    igraph::V(g)$shape       <- "circle"
    igraph::V(g)$label       <- seq_len(N)
    igraph::V(g)$label.cex   <- 0.85
    igraph::V(g)$label.color <- "black"
    igraph::V(g)$frame.color <- "darkgray"
    igraph::V(g)$size        <- 12
    igraph::V(g)$color       <- "lightblue"

    # Treated unit: highlight prominently
    igraph::V(g)$color[treated_idx]       <- "royalblue"
    igraph::V(g)$size[treated_idx]        <- 16
    igraph::V(g)$frame.color[treated_idx] <- "black"

    # Isolated donors: distinct color so the design is visible
    if (show_isolated && length(iso) > 0) {
      iso_donors <- setdiff(iso, treated_idx)  # treated should never be isolated
      igraph::V(g)$color[iso_donors]       <- "salmon"
      igraph::V(g)$frame.color[iso_donors] <- "darkred"
      igraph::V(g)$shape[iso_donors]       <- "square"
    }

    # Edge styling -------------------------------------------------------------
    igraph::E(g)$width <- 1.2
    if (edge_by_weight && igraph::ecount(g) > 0) {
      w_edge <- igraph::E(g)$weight
      # Map edge weight to gray intensity: small weight = light gray, large = black
      w_norm <- (w_edge - min(w_edge)) / (max(w_edge) - min(w_edge) + 1e-12)
      gray_levels <- 0.7 - 0.65 * w_norm   # 0.05 (dark) .. 0.7 (light)
      igraph::E(g)$color <- gray(gray_levels)
      igraph::E(g)$width <- 0.8 + 2.0 * w_norm
    } else {
      igraph::E(g)$color <- "gray60"
    }
    if (directed) igraph::E(g)$arrow.size <- 0.4

    # Layout: use the actual spatial coordinates --------------------------------
    # Pad the plot window slightly so labels don't get clipped.
    x_range <- range(coords[, 1]); y_range <- range(coords[, 2])
    pad_x   <- 0.08 * diff(x_range); pad_y <- 0.08 * diff(y_range)

    plot(g,
         layout    = coords,
         rescale   = FALSE,
         xlim      = c(x_range[1] - pad_x, x_range[2] + pad_x),
         ylim      = c(y_range[1] - pad_y, y_range[2] + pad_y),
         main      = main)

    # Legend -------------------------------------------------------------------
    legend_labels <- c(sprintf("Treated unit (%d)", treated_idx),
                       "Connected donors")
    legend_cols   <- c("royalblue", "lightblue")
    legend_pch    <- c(21, 21)
    if (show_isolated && length(iso) > 0) {
      legend_labels <- c(legend_labels,
                         sprintf("Isolated donors (n=%d)", length(iso)))
      legend_cols   <- c(legend_cols, "salmon")
      legend_pch    <- c(legend_pch, 22)
    }
    legend("bottom",
           legend = legend_labels,
           pch    = legend_pch,
           pt.bg  = legend_cols,
           col    = "darkgray",
           pt.cex = 2,
           bty    = "n",
           horiz  = TRUE,
           inset  = c(0, -0.05),
           xpd    = TRUE)

    invisible(g)
  }


w <- construct_weight_matrix_v2(
N = 20, treated_idx = 1, donor_idx = 2:20,
connectivity = "low", allow_isolated = TRUE, n_isolated = 3
)
 plot_network(w, main = "NASC test network: low connectivity, 3 isolated")

