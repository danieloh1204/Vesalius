# ==============================================================================
# 07_dimreduction_comparison.R
# Purpose: Compare how different dimensionality reduction methods
#          respond to spatial dimension selection.
#          NMF excluded due to Seurat v5 bug.
# Output:  analysis/results/dimreduction_comparison.csv
# ==============================================================================

devtools::load_all(".")

library(oneiric)
library(mclust)

set.seed(42)

dir.create("analysis/results", recursive = TRUE, showWarnings = FALSE)

# --- Configuration -----------------------------------------------------------

# NMF excluded: Seurat v5 Assay5 incompatibility
# UMAP excluded: requires PCA first, test separately if needed
dim_methods <- c("PCA", "PCA_L")

architectures <- list(
    list(name = "3_territories", n_terr = 3, pattern = "circle"),
    list(name = "5_territories", n_terr = 5, pattern = "circle")
)

n_cells      <- 1000
n_replicates <- 5
signal       <- 5

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

# --- Main loop ----------------------------------------------------------------

results <- data.frame()

for (arch in architectures) {
    message(sprintf("\n=== %s ===", arch$name))

    for (rep_i in seq_len(n_replicates)) {
        dat <- make_data(arch, seed = rep_i * 100 + arch$n_terr)

        for (dr in dim_methods) {
            message(sprintf("  %s rep %d: %s", arch$name, rep_i, dr))

            ves <- tryCatch({
                v <- build_vesalius_assay(dat$coords, dat$counts,
                                         verbose = FALSE)
                v <- generate_embeddings(v, dim_reduction = dr,
                                         nfeatures = dat$total_genes,
                                         verbose = FALSE)
                smooth_image(v, sigma = 2, iter = 5, verbose = FALSE)
            }, error = function(e) {
                message(sprintf("    SKIP: %s", e$message))
                NULL
            })
            if (is.null(ves)) next

            dims_to_test <- seq_len(min(arch$n_terr + 5, 10))
            embed <- tryCatch(
                vesalius:::check_embedding_selection(ves, "last",
                                                      max(dims_to_test)),
                error = function(e) NULL
            )
            if (is.null(embed)) next

            truth <- dat$layout$Territory

            # Baseline: all dims
            km_all <- kmeans(embed[, dims_to_test, drop = FALSE],
                             centers = arch$n_terr,
                             nstart = 10, iter.max = 100)
            ari_all <- adjustedRandIndex(km_all$cluster, truth)

            # Selected dims
            ves_scored <- tryCatch(
                select_dimension(ves, dimensions = dims_to_test,
                                 metric = "morans_I", verbose = FALSE),
                error = function(e) NULL
            )
            if (is.null(ves_scored)) next

            sel <- ves_scored@meta$selected_dims
            if (length(sel) >= 1 && max(sel) <= ncol(embed)) {
                km_sel <- kmeans(embed[, sel, drop = FALSE],
                                 centers = arch$n_terr,
                                 nstart = 10, iter.max = 100)
                ari_sel <- adjustedRandIndex(km_sel$cluster, truth)
            } else {
                ari_sel <- NA
            }

            scores_df <- ves_scored@meta$dim_scores

            row <- data.frame(
                architecture  = arch$name,
                n_territories = arch$n_terr,
                replicate     = rep_i,
                dim_reduction = dr,
                n_selected    = length(sel),
                selected_dims = paste(sel, collapse = ";"),
                ari_baseline  = round(ari_all, 4),
                ari_selected  = round(ari_sel, 4),
                improvement   = round(ari_sel - ari_all, 4),
                top_score     = round(scores_df$score[1], 4),
                bottom_score  = round(scores_df$score[nrow(scores_df)], 4),
                stringsAsFactors = FALSE
            )
            results <- rbind(results, row)
        }
    }
}

# --- Save results -------------------------------------------------------------
write.csv(results, "analysis/results/dimreduction_comparison.csv",
          row.names = FALSE)

message(sprintf("\nDim reduction comparison complete: %d rows saved.",
                nrow(results)))

# --- Summary ------------------------------------------------------------------
message("\n=== MEAN ARI BY METHOD & ARCHITECTURE ===")
agg <- aggregate(cbind(ari_baseline, ari_selected, improvement) ~
                     dim_reduction + architecture,
                 data = results, FUN = mean, na.rm = TRUE)
print(agg, row.names = FALSE)

message("\n=== OVERALL BY DIM REDUCTION METHOD ===")
agg2 <- aggregate(cbind(ari_baseline, ari_selected, improvement) ~
                      dim_reduction,
                  data = results, FUN = mean, na.rm = TRUE)
print(agg2, row.names = FALSE)