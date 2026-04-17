# ==============================================================================
# 05_threshold_sensitivity.R
# Purpose: Sweep the gap-threshold parameter to find the optimal value
#          across different tissue complexities. This directly addresses
#          the observation that threshold=0.15 works for simple tissues
#          but is too aggressive for complex ones.
# Output:  analysis/results/threshold_sensitivity.csv
# ==============================================================================

devtools::load_all(".")

library(oneiric)
library(mclust)

set.seed(42)

dir.create("analysis/results", recursive = TRUE, showWarnings = FALSE)

# --- Configuration -----------------------------------------------------------

thresholds <- c(0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30, 0.50)

architectures <- list(
    list(name = "3_territories", n_terr = 3, pattern = "circle"),
    list(name = "5_territories", n_terr = 5, pattern = "circle"),
    list(name = "8_territories", n_terr = 8, pattern = "circle")
)

n_cells      <- 1000
n_replicates <- 5
signal       <- 5   # mid-range signal

# --- Helper ------------------------------------------------------------------

make_data <- function(arch, seed) {
    set.seed(seed)
    layout <- simulate_spatial(n_cells = n_cells,
                               n_territories = arch$n_terr,
                               pattern = arch$pattern)[[1]]
    n_actual    <- nrow(layout)
    n_sig_genes <- arch$n_terr * 10
    n_noise     <- 50
    total       <- n_sig_genes + n_noise

    counts <- matrix(rpois(n_actual * total, lambda = 2),
                     nrow = total, ncol = n_actual)
    colnames(counts) <- layout$barcodes
    rownames(counts) <- paste0("Gene", seq_len(total))

    for (t in seq(0, arch$n_terr - 1)) {
        g_start <- t * 10 + 1
        g_end   <- (t + 1) * 10
        idx     <- which(layout$Territory == t)
        counts[g_start:g_end, idx] <- counts[g_start:g_end, idx] + signal
    }

    coords <- layout[, c("barcodes", "x", "y")]
    list(coords = coords, counts = counts, layout = layout,
         total_genes = total)
}

# --- Main sweep --------------------------------------------------------------

results <- data.frame()

for (arch in architectures) {
    message(sprintf("\n=== %s ===", arch$name))

    for (rep_i in seq_len(n_replicates)) {
        dat <- make_data(arch, seed = rep_i * 100 + arch$n_terr)

        # Build and process once per replicate
        ves <- tryCatch({
            v <- build_vesalius_assay(dat$coords, dat$counts, verbose = FALSE)
            v <- generate_embeddings(v, dim_reduction = "PCA",
                                     nfeatures = dat$total_genes,
                                     verbose = FALSE)
            smooth_image(v, sigma = 2, iter = 5, verbose = FALSE)
        }, error = function(e) {
            message(sprintf("  SKIP rep %d: %s", rep_i, e$message))
            NULL
        })
        if (is.null(ves)) next

        dims_to_test <- seq_len(min(arch$n_terr + 5, 15))
        embed <- vesalius:::check_embedding_selection(ves, "last",
                                                       max(dims_to_test))

        # Baseline: all dims, no selection
        km_all <- kmeans(embed[, dims_to_test, drop = FALSE],
                         centers = arch$n_terr, nstart = 10, iter.max = 100)
        ari_all <- adjustedRandIndex(km_all$cluster, dat$layout$Territory)

        # Sweep thresholds
        for (thr in thresholds) {
            ves_scored <- tryCatch(
                select_dimension(ves, dimensions = dims_to_test,
                                 metric = "morans_I",
                                 threshold = thr,
                                 verbose = FALSE),
                error = function(e) NULL
            )
            if (is.null(ves_scored)) next

            sel <- ves_scored@meta$selected_dims

            if (length(sel) >= 1 && max(sel) <= ncol(embed)) {
                km_sel <- kmeans(embed[, sel, drop = FALSE],
                                 centers = arch$n_terr,
                                 nstart = 10, iter.max = 100)
                ari_sel <- adjustedRandIndex(km_sel$cluster,
                                              dat$layout$Territory)
            } else {
                ari_sel <- NA
            }

            row <- data.frame(
                architecture  = arch$name,
                n_territories = arch$n_terr,
                replicate     = rep_i,
                threshold     = thr,
                n_selected    = length(sel),
                ari_baseline  = round(ari_all, 4),
                ari_selected  = round(ari_sel, 4),
                improvement   = round(ari_sel - ari_all, 4),
                stringsAsFactors = FALSE
            )
            results <- rbind(results, row)
        }
    }
}

# --- Save results -------------------------------------------------------------
write.csv(results, "analysis/results/threshold_sensitivity.csv",
          row.names = FALSE)

message(sprintf("\nThreshold sweep complete: %d rows saved.", nrow(results)))

# --- Summary ------------------------------------------------------------------
message("\n=== MEAN ARI BY THRESHOLD & ARCHITECTURE ===")
agg <- aggregate(cbind(ari_selected, n_selected, improvement) ~
                     threshold + architecture,
                 data = results, FUN = mean, na.rm = TRUE)
agg <- agg[order(agg$architecture, agg$threshold), ]
print(agg, row.names = FALSE)

message("\n=== BEST THRESHOLD PER ARCHITECTURE ===")
for (a in unique(agg$architecture)) {
    sub <- agg[agg$architecture == a, ]
    best <- sub[which.max(sub$ari_selected), ]
    message(sprintf("  %s: threshold=%.2f -> ARI=%.3f (%d dims selected)",
                    a, best$threshold, best$ari_selected, round(best$n_selected)))
}

message("\n=== OVERALL BEST THRESHOLD (mean across architectures) ===")
overall <- aggregate(ari_selected ~ threshold, data = results,
                     FUN = mean, na.rm = TRUE)
best_overall <- overall[which.max(overall$ari_selected), ]
message(sprintf("  Best overall: threshold=%.2f -> mean ARI=%.3f",
                best_overall$threshold, best_overall$ari_selected))