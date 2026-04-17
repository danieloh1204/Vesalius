# ==============================================================================
# 02_preprocessing_impact.R
# Purpose: Compare how different preprocessing pipelines affect the
#          spatial information scores across ALL dimensions (not just max).
# Output:  analysis/results/preprocessing_impact.csv
# ==============================================================================

devtools::load_all(".")

library(oneiric)

set.seed(42)

dir.create("analysis/results", recursive = TRUE, showWarnings = FALSE)

# --- Configuration -----------------------------------------------------------

pipelines <- list(
    list(name = "raw",       smooth = FALSE, equalize = FALSE),
    list(name = "smoothed",  smooth = TRUE,  equalize = FALSE),
    list(name = "smooth_eq", smooth = TRUE,  equalize = TRUE)
)

n_territories <- 4
n_cells       <- 1000
dims_to_test  <- 1:10
n_replicates  <- 5
metrics       <- c("morans_I", "gearys_C", "spatial_entropy")

# --- Helper: build one synthetic dataset -------------------------------------

make_data <- function(seed) {
    set.seed(seed)
    layout <- simulate_spatial(n_cells = n_cells,
                               n_territories = n_territories,
                               pattern = "circle")[[1]]
    n_actual    <- nrow(layout)
    n_sig_genes <- n_territories * 10
    n_noise     <- 50
    total_genes <- n_sig_genes + n_noise

    counts <- matrix(rpois(n_actual * total_genes, lambda = 2),
                     nrow = total_genes, ncol = n_actual)
    colnames(counts) <- layout$barcodes
    rownames(counts) <- paste0("Gene", seq_len(total_genes))

    for (t in seq(0, n_territories - 1)) {
        g_start <- t * 10 + 1
        g_end   <- (t + 1) * 10
        idx     <- which(layout$Territory == t)
        counts[g_start:g_end, idx] <- counts[g_start:g_end, idx] + 5
    }

    coords <- layout[, c("barcodes", "x", "y")]
    list(coords = coords, counts = counts, total_genes = total_genes)
}


# --- Main loop ----------------------------------------------------------------

results <- data.frame()

for (rep_i in seq_len(n_replicates)) {
    dat <- make_data(seed = rep_i * 42)

    for (pipe in pipelines) {
        ves <- build_vesalius_assay(dat$coords, dat$counts, verbose = FALSE)
        ves <- generate_embeddings(ves, dim_reduction = "PCA",
                                   nfeatures = dat$total_genes,
                                   verbose = FALSE)

        if (pipe$equalize) {
            ves <- tryCatch(
                equalize_image(ves, verbose = FALSE),
                error = function(e) ves  # skip if not available
            )
        }
        if (pipe$smooth) {
            ves <- smooth_image(ves, sigma = 2, iter = 5, verbose = FALSE)
        }

        for (met in metrics) {
            ves_scored <- tryCatch(
                select_dimension(ves, dimensions = dims_to_test,
                                 metric = met, verbose = FALSE),
                error = function(e) {
                    message(sprintf("ERROR %s/%s: %s", pipe$name, met,
                                    e$message))
                    NULL
                }
            )
            if (is.null(ves_scored)) next

            scores_df <- ves_scored@meta$dim_scores

            for (r in seq_len(nrow(scores_df))) {
                row <- data.frame(
                    replicate  = rep_i,
                    pipeline   = pipe$name,
                    metric     = met,
                    dimension  = scores_df$dimension[r],
                    score      = scores_df$score[r],
                    rank       = r,
                    stringsAsFactors = FALSE
                )
                results <- rbind(results, row)
            }
        }
    }
}

# --- Save results -------------------------------------------------------------
write.csv(results, "analysis/results/preprocessing_impact.csv",
          row.names = FALSE)

message(sprintf("\nPreprocessing impact analysis complete: %d rows saved.",
                nrow(results)))

# --- Summary ------------------------------------------------------------------
message("\n=== MEAN SCORE BY PIPELINE & METRIC ===")
agg <- aggregate(score ~ pipeline + metric, data = results,
                 FUN = mean, na.rm = TRUE)
agg <- agg[order(agg$metric, -agg$score), ]
print(agg)

message("\n=== TOP DIMENSION SCORE (RANK 1) BY PIPELINE ===")
top <- results[results$rank == 1, ]
agg_top <- aggregate(score ~ pipeline + metric, data = top,
                     FUN = function(x) round(mean(x), 4))
print(agg_top)