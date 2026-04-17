# ==============================================================================
# 08_real_data_all.R
# Purpose: Run select_dimension on Vesalius built-in real data
#          (Slide-seq V2 hippocampus subset, 1052 beads).
#          Full-resolution runs go on HPC.
# Output:  analysis/results/real_data_all.csv
# ==============================================================================

devtools::load_all(".")

set.seed(42)

dir.create("analysis/results", recursive = TRUE, showWarnings = FALSE)

metrics <- c("morans_I", "gearys_C", "spatial_entropy")
all_results <- data.frame()

# ============================================================================
# Vesalius built-in dataset (1052 beads, Slide-seq V2 dentate gyrus)
# ============================================================================
message("\n========================================")
message("Vesalius built-in dataset (1052 beads)")
message("========================================")

data(vesalius, package = "vesalius")

ves <- build_vesalius_assay(coordinates, counts, verbose = FALSE)
ves <- generate_embeddings(ves, dim_reduction = "PCA",
                           nfeatures = 2000, verbose = FALSE)
ves <- smooth_image(ves, sigma = 5, iter = 10, verbose = FALSE)
ves <- equalize_image(ves, sleft = 5, sright = 5, verbose = FALSE)

for (met in metrics) {
    message(sprintf("  Scoring with %s...", met))
    ves_scored <- tryCatch(
        select_dimension(ves, dimensions = 1:15,
                         metric = met, verbose = FALSE),
        error = function(e) {
            message(sprintf("  ERROR: %s", e$message))
            NULL
        }
    )
    if (is.null(ves_scored)) next

    scores <- ves_scored@meta$dim_scores
    scores$dataset    <- "builtin_1052"
    scores$metric     <- met
    scores$selected   <- scores$dimension %in% ves_scored@meta$selected_dims
    all_results <- rbind(all_results, scores)
}

# --- Territory recovery comparison ---
message("\n  Running territory recovery comparison...")
ves_scored <- select_dimension(ves, dimensions = 1:15,
                                metric = "morans_I", verbose = FALSE)
sel_dims <- ves_scored@meta$selected_dims

# Full pipeline with ALL dims (default Vesalius behavior)
ves_all <- segment_image(ves, method = "kmeans",
                          col_resolution = 5, verbose = FALSE)
ves_all <- isolate_territories(ves_all, verbose = FALSE)

# Full pipeline with SELECTED dims only
ves_sel <- smooth_image(ves, dimensions = sel_dims,
                         sigma = 5, iter = 10, verbose = FALSE)
ves_sel <- equalize_image(ves_sel, dimensions = sel_dims,
                           sleft = 5, sright = 5, verbose = FALSE)
ves_sel <- segment_image(ves_sel, dimensions = sel_dims,
                          method = "kmeans", col_resolution = 5,
                          verbose = FALSE)
ves_sel <- isolate_territories(ves_sel, verbose = FALSE)

message(sprintf("\n  Selected dims: [%s]", paste(sel_dims, collapse = ", ")))

# Count territories recovered
ter_all <- get_territories(ves_all)
ter_sel <- get_territories(ves_sel)
n_ter_all <- length(unique(ter_all[, ncol(ter_all)])) - 1
n_ter_sel <- length(unique(ter_sel[, ncol(ter_sel)])) - 1

message(sprintf("  Territories with ALL 15 dims: %d", n_ter_all))
message(sprintf("  Territories with SELECTED dims [%s]: %d",
                paste(sel_dims, collapse = ","), n_ter_sel))

# ============================================================================
# SAVE & SUMMARIZE
# ============================================================================
write.csv(all_results, "analysis/results/real_data_all.csv",
          row.names = FALSE)

message(sprintf("\nReal data analysis complete: %d rows saved.", nrow(all_results)))

message("\n=== SELECTED DIMS PER METRIC ===")
sel <- all_results[all_results$selected == TRUE, ]
if (nrow(sel) > 0) {
    agg <- aggregate(dimension ~ metric, data = sel,
                     FUN = function(x) paste(sort(x), collapse = ", "))
    colnames(agg)[2] <- "selected_dims"
    print(agg, row.names = FALSE)
}

message("\n=== FULL SCORE TABLE (Moran's I) ===")
moran_scores <- all_results[all_results$metric == "morans_I",
                             c("dimension", "score", "selected")]
moran_scores <- moran_scores[order(-moran_scores$score), ]
print(moran_scores, row.names = FALSE)